#!/usr/bin/env python3
"""Established AMR determinant baseline for the frozen E. coli research cohort.

Commands:
  call   resumable AMRFinderPlus calls on the PR #6 genome cache, through the
         PR #4 caller wrapper (bin/amrscan.py) and its provenance checks
  build  determinant-to-antimicrobial interpretation with a versioned mapping,
         join to PR #5 phenotype endpoints, metrics, stratified counts and
         discordance tables

The mapping is fixed before calls are read. PR #5 phenotypes are used as
frozen, except where a versioned correction table named in the config
swaps R and S for a whole submitting project (see decision 0004).
Nothing here trains a model.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import types

sys.path.insert(0, str(Path(__file__).resolve().parent))
import amrscan  # noqa: E402
import population_structure as pop  # noqa: E402

VERSION = "1.1.0"
NA = "NA"

CALL_POLICY = {"ident_min": -1, "coverage_min": 0.5, "translation_table": 11,
               "report_all_equal": True, "organism": "Escherichia", "plus": False}
DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "benchmark/cohorts/ecoli_research_v1/baseline_config.json"


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------- caller layer

def database_fingerprint(database: Path) -> str:
    """Same fingerprint as amrscan.call: SHA-256 of the sorted file-hash map."""
    database = Path(database)
    files = {p.relative_to(database).as_posix(): amrscan.sha256(p)
             for p in sorted(database.rglob("*")) if p.is_file()}
    if not files:
        raise ValueError("AMRFinderPlus database directory is empty")
    return "sha256:" + hashlib.sha256(amrscan.json_text(files).encode()).hexdigest()


def call_status(call_dir: Path, isolate: str, genome_sha: str, database_id: str,
                caller_version: str | None = None) -> tuple[str, str]:
    """Classify an existing call directory without changing it.

    Returns (status, detail); status is one of missing, success, caller_failed,
    provenance_mismatch, harmonization_failed.
    """
    call_dir = Path(call_dir)
    prov_path = call_dir / "provenance.json"
    if not prov_path.exists():
        return "missing", "no provenance"
    try:
        prov = json.loads(prov_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return "provenance_mismatch", f"unreadable provenance: {error}"
    if prov.get("sample_id") != isolate or prov.get("tool") != "AMRFinderPlus":
        return "provenance_mismatch", "sample or tool differs"
    if prov.get("status") != "success":
        return "caller_failed", f"caller status {prov.get('status')}"
    if prov.get("input_sha256") != genome_sha:
        return "provenance_mismatch", "caller input hash differs from verified genome"
    if prov.get("database_id") != database_id:
        return "provenance_mismatch", "database fingerprint differs"
    if caller_version and prov.get("caller_version") != caller_version:
        return "provenance_mismatch", f"caller version {prov.get('caller_version')}"
    if prov.get("thresholds") != CALL_POLICY:
        return "provenance_mismatch", "caller options differ from the call policy"
    raw = call_dir / "amrfinder.tsv"
    if not raw.is_file() or prov.get("raw_sha256") != amrscan.sha256(raw):
        return "provenance_mismatch", "raw report hash differs"
    harmonized = call_dir / "harmonized.tsv"
    if not harmonized.is_file():
        return "harmonization_failed", "harmonized report missing"
    with tempfile.TemporaryDirectory(prefix="amrscan-harmonize-check-") as tmp:
        rebuilt = Path(tmp) / "expected.tsv"
        try:
            amrscan.harmonize(raw, prov_path, isolate, rebuilt)
        except (ValueError, KeyError, OSError) as error:
            return "harmonization_failed", str(error)
        if amrscan.sha256(rebuilt) != amrscan.sha256(harmonized):
            return "harmonization_failed", "harmonized report differs from raw report"
    return "success", ""


def run_one(row: dict, cache: Path, calls: Path, database: Path, database_id: str,
            caller_version: str, executable: str, threads: int, database_source: str,
            caller=None, harmonizer=None) -> dict:
    """Call one isolate unless a valid call exists. Invalid earlier output is
    moved to calls/.quarantine rather than overwritten."""
    caller = caller or amrscan.call
    harmonizer = harmonizer or amrscan.harmonize
    isolate, assembly, genome_sha = row["isolate_id"], row["assembly_accession"], row["genome_sha256"]
    out = Path(calls) / isolate
    status, detail = call_status(out, isolate, genome_sha, database_id, caller_version)
    if status == "success":
        return {"isolate_id": isolate, "status": "success", "action": "reused", "detail": ""}
    genome = Path(cache) / assembly / "genome.fna"
    if not genome.is_file() or amrscan.sha256(genome) != genome_sha:
        return {"isolate_id": isolate, "status": "genome_unavailable", "action": "skipped",
                "detail": "cached genome missing or differs from PR #6 hash"}
    if out.exists():
        quarantine = Path(calls) / ".quarantine" / f"{isolate}.{datetime.now(timezone.utc):%Y%m%dT%H%M%S%f}"
        quarantine.parent.mkdir(parents=True, exist_ok=True)
        out.rename(quarantine)
    out.mkdir(parents=True)
    args = types.SimpleNamespace(
        sample=isolate, input=genome, database=Path(database), database_source=database_source,
        executable=executable, threads=threads, identity=CALL_POLICY["ident_min"],
        coverage=CALL_POLICY["coverage_min"], organism=CALL_POLICY["organism"], plus=CALL_POLICY["plus"],
        log=out / "amrfinder.log", raw=out / "amrfinder.tsv", provenance=out / "provenance.json")
    try:
        caller(args)
    except Exception as error:  # the wrapper has already written failed provenance
        return {"isolate_id": isolate, "status": "caller_failed", "action": "called", "detail": str(error)[:300]}
    try:
        harmonizer(out / "amrfinder.tsv", out / "provenance.json", isolate, out / "harmonized.tsv")
    except (ValueError, KeyError, OSError) as error:
        return {"isolate_id": isolate, "status": "harmonization_failed", "action": "called", "detail": str(error)}
    status, detail = call_status(out, isolate, genome_sha, database_id, caller_version)
    return {"isolate_id": isolate, "status": status, "action": "called", "detail": detail}


def call_cohort(population: Path, cache: Path, calls: Path, database: Path, executable: str = "amrfinder",
                workers: int = 4, threads: int = 1, limit: int | None = None, database_source: str | None = None,
                caller=None, harmonizer=None, version_probe=None) -> dict:
    rows = pop.read_tsv(Path(population) / "population_metadata.tsv")
    rows = sorted(rows, key=lambda r: r["isolate_id"])
    if limit is not None:
        rows = rows[:limit]
    database = Path(database)
    database_id = database_fingerprint(database)
    if version_probe is None:
        import subprocess
        resolved = shutil.which(executable)
        if not resolved:
            raise ValueError(f"AMRFinderPlus executable not found: {executable}")
        caller_version = subprocess.run([resolved, "--version"], check=True, text=True,
                                        capture_output=True).stdout.strip()
    else:
        caller_version = version_probe()
    source = database_source or f"NCBI AMRFinderPlus database {database.resolve().name}, software {caller_version}"
    calls = Path(calls)
    calls.mkdir(parents=True, exist_ok=True)
    log_path = calls / "call_log.jsonl"
    counts = {}
    with ThreadPoolExecutor(max_workers=workers) as pool, log_path.open("a", encoding="utf-8") as log:
        futures = [pool.submit(run_one, r, cache, calls, database, database_id, caller_version, executable,
                               threads, source, caller, harmonizer) for r in rows]
        for future in as_completed(futures):
            result = future.result()
            result["logged_at"] = now()
            log.write(json.dumps(result, sort_keys=True) + "\n")
            log.flush()
            key = f"{result['status']}:{result['action']}"
            counts[key] = counts.get(key, 0) + 1
    return {"isolates": len(rows), "database_id": database_id, "caller_version": caller_version,
            "results": dict(sorted(counts.items()))}




# ---------------------------------------------------------------- strict tables

def read_table(path: Path, required=()) -> list[dict]:
    """TSV with a unique header; every row must have the header's width."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or not any(h.strip() for h in header):
            raise ValueError(f"{path}: missing header")
        if len(header) != len(set(header)):
            raise ValueError(f"{path}: duplicate column in header")
        missing = sorted(set(required) - set(header))
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        rows = []
        for line, record in enumerate(reader, start=2):
            if len(record) != len(header):
                raise ValueError(f"{path}: malformed row at line {line}")
            rows.append(dict(zip(header, record)))
    return rows


def write_table(path: Path, fields, rows) -> None:
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------- mapping

MAPPING_FIELDS = ("mapping_version organism rule_id determinant determinant_type gene_family allele_or_mutation "
                  "antimicrobial interpretation requires evidence_source evidence_terms evidence_reference "
                  "evidence_level review_status notes").split()
INTERPRETATIONS = {"predicts_resistance", "does_not_predict", "uncertain", "HUMAN_REVIEW_REQUIRED", "component"}
PARTIAL_METHODS = {"INTERNAL_STOP", "PARTIALP", "PARTIALX"}
FLUOROQUINOLONES = ("ciprofloxacin", "levofloxacin")
VARIANTS = ("primary", "uncertain_as_positive", "uncertain_as_negative", "fq_two_mechanisms")


def load_mapping(rules_path: Path, meta_path: Path, expected_version: str | None = None,
                 organism: str = "Escherichia coli") -> dict:
    meta = json.loads(Path(meta_path).read_text(encoding="utf-8"))
    version = meta.get("mapping_version")
    if not version:
        raise ValueError("mapping metadata lacks mapping_version")
    if expected_version and version != expected_version:
        raise ValueError(f"mapping version {version} differs from expected {expected_version}")
    if meta.get("rules_sha256") != amrscan.sha256(rules_path):
        raise ValueError("mapping rules checksum differs from mapping metadata")
    study = meta.get("study_antimicrobials") or {}
    combos = meta.get("combination_rules", {})
    rules, seen = {}, {}
    for row in read_table(rules_path, MAPPING_FIELDS):
        key = (row["determinant"], row["antimicrobial"])
        if row["mapping_version"] != version:
            raise ValueError(f"rule {row['rule_id']} has mapping version {row['mapping_version']}")
        if row["organism"] != organism:
            raise ValueError(f"rule {row['rule_id']} is for organism {row['organism']}")
        if row["interpretation"] not in INTERPRETATIONS:
            raise ValueError(f"rule {row['rule_id']} has unknown interpretation {row['interpretation']}")
        if row["antimicrobial"] not in study:
            raise ValueError(f"rule {row['rule_id']} names an antimicrobial outside the study set")
        if row["rule_id"] != f"{row['determinant']}|{row['antimicrobial']}":
            raise ValueError(f"rule id {row['rule_id']} does not match determinant and antimicrobial")
        if row["interpretation"] == "component":
            parts = combos.get(row["antimicrobial"], {}).get("components", {})
            if row["requires"] not in parts:
                raise ValueError(f"rule {row['rule_id']} names an undefined component {row['requires']}")
        if key in seen:
            if seen[key] == row["interpretation"]:
                raise ValueError(f"duplicate mapping rule {row['rule_id']}")
            raise ValueError(f"conflicting mapping rules for {row['rule_id']}")
        seen[key] = row["interpretation"]
        rules.setdefault(row["antimicrobial"], {})[row["determinant"]] = row
    return {"meta": meta, "version": version, "rules": rules, "study": study, "combinations": combos,
            "rules_sha256": meta["rules_sha256"], "metadata_sha256": amrscan.sha256(meta_path)}


def amr_determinants(records: list[dict]) -> list[dict]:
    """AMR elements of one isolate, one entry per symbol; full evidence beats partial."""
    out = {}
    for r in records:
        if r["element_type"] != "AMR":
            continue
        sym = r["gene"]
        partial = r["evidence_type"] in PARTIAL_METHODS
        prev = out.get(sym)
        if prev is None or (prev["partial"] and not partial):
            out[sym] = {"symbol": sym, "partial": partial, "method": r["evidence_type"],
                        "subtype": r["element_subtype"], "class": r["subclass"]}
    return [out[k] for k in sorted(out)]


def interpret(dets: list[dict], drug_rules: dict, combination: dict | None = None) -> dict:
    """Determinant status of one isolate for one antimicrobial."""
    full = {d["symbol"] for d in dets if not d["partial"]}
    positive, uncertain, components = set(), set(), {}
    for d in dets:
        rule = drug_rules.get(d["symbol"])
        if rule is None or rule["interpretation"] == "does_not_predict":
            continue
        if rule["interpretation"] in {"uncertain", "HUMAN_REVIEW_REQUIRED"} or d["partial"]:
            uncertain.add(d["symbol"])
        elif rule["interpretation"] == "component":
            components.setdefault(rule["requires"], set()).add(d["symbol"])
        elif rule["requires"] not in {"", NA} and not set(rule["requires"].split("|")) & full:
            uncertain.add(d["symbol"])
        else:
            positive.add(d["symbol"])
    if components:
        members = set().union(*components.values())
        if combination and set(combination["components"]) <= set(components):
            positive |= members
        else:
            uncertain |= members
    status = "positive" if positive else "uncertain" if uncertain else "negative"
    return {"status": status, "positive": sorted(positive), "uncertain": sorted(uncertain - positive)}


def variant_status(result: dict, variant: str, drug: str) -> str | None:
    status = result["status"]
    if variant == "uncertain_as_positive":
        return "positive" if status == "uncertain" else status
    if variant == "uncertain_as_negative":
        return "negative" if status == "uncertain" else status
    if variant == "fq_two_mechanisms":
        if drug not in FLUOROQUINOLONES:
            return None
        n = len(result["positive"])
        return "positive" if n >= 2 else "uncertain" if n == 1 else status
    return status


# ---------------------------------------------------------------- metrics

def ratio(a, b):
    return a / b if b else None


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple:
    if n == 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def metrics(tp: int, fn: int, fp: int, tn: int) -> dict:
    """Phenotype R is the positive reference state; determinant-positive is the
    baseline positive call. Undefined values are None, never a division by zero."""
    n = tp + fn + fp + tn
    sens, spec = ratio(tp, tp + fn), ratio(tn, tn + fp)
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out = {"N": n, "N_R": tp + fn, "N_S": fp + tn, "TP": tp, "FN": fn, "FP": fp, "TN": tn,
           "sensitivity": sens, "specificity": spec, "PPV": ratio(tp, tp + fp), "NPV": ratio(tn, tn + fn),
           "accuracy": ratio(tp + tn, n),
           "balanced_accuracy": (sens + spec) / 2 if sens is not None and spec is not None else None,
           "F1": ratio(2 * tp, 2 * tp + fp + fn), "MCC": (tp * tn - fp * fn) / den if den else None}
    for name, k, total in (("sensitivity", tp, tp + fn), ("specificity", tn, tn + fp),
                           ("PPV", tp, tp + fp), ("NPV", tn, tn + fn)):
        out[f"{name}_ci_low"], out[f"{name}_ci_high"] = wilson(k, total)
    return out


def fmt(value) -> str:
    if value is None:
        return NA
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def confusion(pairs) -> tuple:
    """pairs of (phenotype, determinant status) with status positive/negative."""
    tp = fn = fp = tn = 0
    for phenotype, status in pairs:
        if phenotype == "R":
            tp += status == "positive"
            fn += status == "negative"
        elif phenotype == "S":
            fp += status == "positive"
            tn += status == "negative"
    return tp, fn, fp, tn


METRIC_FIELDS = ("N N_R N_S TP FN FP TN sensitivity sensitivity_ci_low sensitivity_ci_high specificity "
                 "specificity_ci_low specificity_ci_high PPV PPV_ci_low PPV_ci_high NPV NPV_ci_low NPV_ci_high "
                 "accuracy balanced_accuracy F1 MCC").split()


# ---------------------------------------------------------------- upstream inputs

def load_config(path: Path) -> dict:
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("mapping", "expected_upstream", "caller", "stratum_support", "temporal", "mapping_audit"):
        if key not in config:
            raise ValueError(f"baseline config lacks {key}")
    corrections = config.get("phenotype_corrections")
    if corrections is not None and not (corrections.get("table") and corrections.get("sha256")
                                        and corrections.get("decision_record")):
        raise ValueError("phenotype_corrections needs table, sha256 and decision_record")
    support = config["stratum_support"]
    if support.get("decision") != "HUMAN_REVIEW_REQUIRED" and not support.get("decision_record"):
        raise ValueError("stratum support thresholds need a decision record or HUMAN_REVIEW_REQUIRED")
    return config


SWAP_RS = {"R": "S", "S": "R"}
CORRECTION_FIELDS = ("correction_id", "project", "action", "scope", "evidence", "decision_record")


def load_corrections(path: Path, expected_sha256: str) -> dict:
    """Project-level phenotype corrections: {project: correction_id}. Only swap_R_S is defined."""
    if amrscan.sha256(path) != expected_sha256:
        raise ValueError("phenotype correction table hash differs from the config")
    out = {}
    for row in read_table(path, CORRECTION_FIELDS):
        if row["action"] != "swap_R_S" or row["scope"] != "all_binary_endpoints":
            raise ValueError(f"unsupported phenotype correction: {row['correction_id']}")
        if row["project"] in out:
            raise ValueError(f"project corrected twice: {row['project']}")
        out[row["project"]] = row["correction_id"]
    return out


def load_population(population: Path, expected_manifest_sha256: str | None, pr5_manifest_sha256: str) -> dict:
    population = Path(population)
    sums = {}
    for line in (population / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            sums[name.strip()] = digest
    for name in ("population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json"):
        if name not in sums or amrscan.sha256(population / name) != sums[name]:
            raise ValueError(f"PR #6 population file missing or changed: {name}")
    manifest_sha = amrscan.sha256(population / "population_manifest.json")
    if expected_manifest_sha256 and manifest_sha != expected_manifest_sha256:
        raise ValueError("PR #6 population manifest hash differs from the expected value")
    manifest = json.loads((population / "population_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("pr5_cohort_manifest_sha256") != pr5_manifest_sha256:
        raise ValueError("PR #6 population metadata was built from a different PR #5 cohort")
    meta = {r["isolate_id"]: r for r in read_table(population / "population_metadata.tsv",
                                                   ("isolate_id", "assembly_accession", "genome_sha256"))}
    groups = {r["isolate_id"]: r for r in read_table(population / "evaluation_groups.tsv", ("isolate_id",))}
    if set(meta) != set(groups):
        raise ValueError("PR #6 metadata and evaluation groups list different isolates")
    return {"manifest": manifest, "manifest_sha256": manifest_sha, "meta": meta, "groups": groups,
            "population_id": manifest.get("population_id", NA)}


def load_calls(population: dict, calls: Path, database_id: str, caller_version: str, workers: int = 4) -> dict:
    """Per isolate: call status and AMR determinants. Never turns a failure into no-hit."""
    def one(item):
        isolate, row = item
        if pop.missing(row.get("genome_sha256")):
            return isolate, {"status": "genome_unavailable", "detail": "no verified genome", "dets": []}
        if row.get("genome_qc_status") == "fail":
            return isolate, {"status": "genome_qc_failure", "detail": row.get("genome_qc_flags", ""), "dets": []}
        call_dir = Path(calls) / isolate
        status, detail = call_status(call_dir, isolate, row["genome_sha256"], database_id, caller_version)
        if status == "missing":
            return isolate, {"status": "caller_missing", "detail": detail, "dets": []}
        if status != "success":
            return isolate, {"status": status, "detail": detail, "dets": []}
        prov = json.loads((call_dir / "provenance.json").read_text(encoding="utf-8"))
        records = read_table(call_dir / "harmonized.tsv", amrscan.SCHEMA)
        if any(r["sample_id"] != isolate for r in records):
            return isolate, {"status": "provenance_mismatch", "detail": "harmonized sample differs", "dets": []}
        return isolate, {"status": "success", "detail": "", "dets": amr_determinants(records),
                         "raw_sha256": prov["raw_sha256"], "records": len(records),
                         "element_types": sorted({r["element_type"] for r in records})}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(pool.map(one, sorted(population["meta"].items())))


def call_times(calls: Path) -> str | None:
    """Latest time a call in this directory was made (from call_log.jsonl)."""
    path = Path(calls) / "call_log.jsonl"
    if not path.exists():
        return None
    latest = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rec = json.loads(line)
            if rec.get("action") == "called" and rec.get("status") == "success":
                latest = max(latest or rec["logged_at"], rec["logged_at"])
    return latest


# ---------------------------------------------------------------- baseline

ENDPOINT_FIELDS = ("cohort_id isolate_id pdt_accession biosample_accession assembly_accession antimicrobial "
                   "study_tier phenotype phenotype_pr5 phenotype_correction phenotype_status eligible_binary_rs population_eligible genotype_status "
                   "determinant_status determinant_count determinant_ids determinant_types "
                   "uncertain_determinant_ids concordance_category fq_two_mechanisms_status mlst_st lineage_id "
                   "lineage_component phylogroup ncbi_snp_cluster country collection_year period geography_group "
                   "source_category project mapping_version caller_version database_id endpoint_eligible "
                   "exclusion_reason").split()
GENOTYPE_UNAVAILABLE = ("genome_unavailable", "genome_qc_failure", "caller_missing", "caller_failed",
                        "provenance_mismatch", "harmonization_failed")
DENOM_FIELDS = ("antimicrobial study_tier total_endpoints binary_rs_endpoints excluded_non_binary_or_conflict "
                "population_eligible_binary excluded_population_ineligible genome_available caller_successful "
                "mapping_interpretable determinant_positive determinant_negative determinant_uncertain "
                "excluded_not_study_antimicrobial excluded_genotype_unavailable excluded_mapping_uncertain "
                "metric_eligible reconciled").split()
STRATUM_FIELDS = ("dimension stratum antimicrobial study_tier N R S determinant_positive determinant_negative "
                  "R_DET_POS R_DET_NEG S_DET_POS S_DET_NEG r_det_neg_fraction top_lineage top_lineage_fraction "
                  "top_project top_project_fraction top_country top_country_fraction support_status flags "
                  "sensitivity specificity PPV NPV balanced_accuracy F1 MCC").split()
DIMENSIONS = {
    "lineage": ("lineage_baseline.tsv", lambda e: e["lineage_id"]),
    "phylogroup": ("lineage_baseline.tsv", lambda e: e["phylogroup"]),
    "country": ("country_baseline.tsv", lambda e: e["country"]),
    "geography_group": ("country_baseline.tsv", lambda e: e["geography_group"]),
    "period": ("temporal_baseline.tsv", lambda e: e["period"]),
    "collection_year": ("temporal_baseline.tsv", lambda e: e["collection_year"]),
    "source": ("source_baseline.tsv", lambda e: e["source_category"]),
    "project": ("error_distribution.tsv", lambda e: e["project"]),
    "snp_cluster": ("error_distribution.tsv", lambda e: e["ncbi_snp_cluster"]),
}
CANDIDATE_FIELDS = ("cohort_id isolate_id pdt_accession biosample_accession assembly_accession antimicrobial "
                    "study_tier phenotype determinant_status determinant_ids determinant_types "
                    "uncertain_determinant_ids other_amr_determinants mlst_st lineage_id lineage_component phylogroup "
                    "ncbi_snp_cluster country collection_year period source_category project genome_qc_status "
                    "caller_version database_id raw_report_sha256 mapping_version").split()
AUDIT_FIELDS = ("antimicrobial determinant determinant_type interpretation requires evidence_level review_status "
                "carriers carriers_R carriers_S r_fraction sole_support_R sole_support_S lineages_json "
                "countries_json flag").split()
UNMAPPED_FIELDS = "determinant subtype subclass isolates population_eligible_isolates".split()


def stratum_label(value: str) -> str:
    return "unknown" if pop.missing(value) or value == "unknown" else value


def share(counter: Counter) -> tuple:
    total = sum(counter.values())
    if not total:
        return NA, NA
    top, n = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[0]
    return top, fmt(n / total)


def stratum_rows(eligible: list[dict], drugs: dict, support: dict) -> dict:
    """Descriptive counts and metrics per dimension x stratum x antimicrobial."""
    out = defaultdict(list)
    min_r, min_s = support["min_R"], support["min_S"]
    for dim, (filename, getter) in DIMENSIONS.items():
        groups = defaultdict(list)
        for e in eligible:
            groups[(stratum_label(getter(e)), e["antimicrobial"])].append(e)
        for (stratum, drug), rows in sorted(groups.items()):
            cats = Counter(e["concordance_category"] for e in rows)
            r, s = cats["R_DET_POS"] + cats["R_DET_NEG"], cats["S_DET_POS"] + cats["S_DET_NEG"]
            top_lin, lin_frac = share(Counter(e["lineage_id"] for e in rows))
            top_proj, proj_frac = share(Counter(e["project"] for e in rows))
            top_ctry, ctry_frac = share(Counter(e["country"] for e in rows))
            flags = []
            if r == 0 or s == 0:
                flags.append("only_R" if s == 0 else "only_S")
            elif min(r, s) / (r + s) < 0.05:
                flags.append("extreme_R_S_imbalance")
            if dim != "lineage" and lin_frac != NA and float(lin_frac) >= 0.5 and len(rows) >= min_r:
                flags.append("one_lineage_majority")
            if dim != "project" and proj_frac != NA and float(proj_frac) >= 0.5 and len(rows) >= min_r:
                flags.append("one_project_majority")
            row = {"dimension": dim, "stratum": stratum, "antimicrobial": drug,
                   "study_tier": drugs[drug], "N": len(rows), "R": r, "S": s,
                   "determinant_positive": cats["R_DET_POS"] + cats["S_DET_POS"],
                   "determinant_negative": cats["R_DET_NEG"] + cats["S_DET_NEG"],
                   **{k: cats[k] for k in ("R_DET_POS", "R_DET_NEG", "S_DET_POS", "S_DET_NEG")},
                   "r_det_neg_fraction": fmt(ratio(cats["R_DET_NEG"], r)),
                   "top_lineage": top_lin, "top_lineage_fraction": lin_frac,
                   "top_project": top_proj, "top_project_fraction": proj_frac,
                   "top_country": top_ctry, "top_country_fraction": ctry_frac,
                   "flags": ";".join(flags) or NA}
            if r >= min_r and s >= min_s:
                m = metrics(cats["R_DET_POS"], cats["R_DET_NEG"], cats["S_DET_POS"], cats["S_DET_NEG"])
                row["support_status"] = "reported"
                row.update({k: fmt(m[k]) for k in ("sensitivity", "specificity", "PPV", "NPV",
                                                    "balanced_accuracy", "F1", "MCC")})
            else:
                row["support_status"] = f"below_support: R {r}, S {s}; need R>={min_r} and S>={min_s}"
                row.update({k: NA for k in ("sensitivity", "specificity", "PPV", "NPV",
                                            "balanced_accuracy", "F1", "MCC")})
            out[filename].append(row)
    return out


def build(processed: Path, population: Path, calls: Path, database: Path, config_path: Path, outdir: Path,
          workers: int = 4, code_commit: str | None = None) -> dict:
    config = load_config(config_path)
    repo = Path(__file__).resolve().parents[1]
    mapping = load_mapping(repo / config["mapping"]["rules"], repo / config["mapping"]["metadata"],
                           config["mapping"]["version"])
    upstream = config["expected_upstream"]
    cohort = pop.load_cohort(processed, upstream.get("pr5_cohort_manifest_sha256"))
    popn = load_population(population, upstream.get("pr6_population_manifest_sha256"), cohort["manifest_sha256"])
    database_id = database_fingerprint(database)
    if database_id != config["caller"]["database_id"]:
        raise ValueError("AMRFinderPlus database fingerprint differs from the configured database")
    calls_by_isolate = load_calls(popn, calls, database_id, config["caller"]["version"], workers)
    candidates = {r["isolate_id"]: r for r in cohort["candidates"]}
    if set(candidates) != set(popn["meta"]):
        raise ValueError("PR #5 candidates and PR #6 population metadata list different isolates")
    drugs = mapping["study"]
    corr_cfg = config.get("phenotype_corrections")
    corrections = load_corrections(repo / corr_cfg["table"], corr_cfg["sha256"]) if corr_cfg else {}
    corrected = Counter()

    endpoints, seen = [], set()
    for e in sorted(cohort["endpoints"], key=lambda r: (r["isolate_id"], r["antimicrobial"])):
        iso, drug = e["isolate_id"], e["antimicrobial"]
        if (iso, drug) in seen:
            raise ValueError(f"duplicate PR #5 endpoint {iso} {drug}")
        seen.add((iso, drug))
        p, g = popn["meta"].get(iso), popn["groups"].get(iso)
        if p is None:
            raise ValueError(f"endpoint isolate {iso} missing from PR #6 population metadata")
        for field in ("pdt_accession", "biosample_accession", "assembly_accession"):
            if e[field] != p[field]:
                raise ValueError(f"{field} mismatch for {iso}: PR #5 {e[field]}, PR #6 {p[field]}")
        call = calls_by_isolate[iso]
        phenotype, correction = e["phenotype_endpoint"], corrections.get(p["project"], NA)
        if correction != NA and phenotype in SWAP_RS:
            phenotype = SWAP_RS[phenotype]
            corrected[correction] += 1
        else:
            correction = NA
        binary = e["eligible_binary_rs"] == "true" and phenotype in {"R", "S"}
        pop_ok = p["eligible_population_analysis"] == "true"
        study = drug in drugs
        result = None
        if not study:
            gstatus = "mapping_unavailable"
        elif call["status"] != "success":
            gstatus = call["status"]
        else:
            result = interpret(call["dets"], mapping["rules"].get(drug, {}), mapping["combinations"].get(drug))
            gstatus = "determinant_" + result["status"]
        reasons = []
        if not binary:
            reasons.append("non_binary_or_conflict")
        if not pop_ok:
            reasons.append("population_ineligible")
        if not study:
            reasons.append("not_study_antimicrobial")
        elif gstatus in GENOTYPE_UNAVAILABLE:
            reasons.append("genotype_unavailable:" + gstatus)
        elif gstatus == "determinant_uncertain":
            reasons.append("mapping_uncertain")
        dstatus = result["status"] if result else NA
        category = (f"{phenotype}_DET_{'POS' if dstatus == 'positive' else 'NEG'}"
                    if binary and dstatus in {"positive", "negative"} else NA)
        rule_types = mapping["rules"].get(drug, {})
        pos = result["positive"] if result else []
        fq = variant_status(result, "fq_two_mechanisms", drug) if result else None
        endpoints.append({
            "cohort_id": e["cohort_id"], "isolate_id": iso, "pdt_accession": e["pdt_accession"],
            "biosample_accession": e["biosample_accession"], "assembly_accession": e["assembly_accession"],
            "antimicrobial": drug, "study_tier": drugs.get(drug, "not_study"),
            "phenotype": phenotype, "phenotype_pr5": e["phenotype_endpoint"], "phenotype_correction": correction,
            "phenotype_status": e["phenotype_status"],
            "eligible_binary_rs": "true" if binary else "false", "population_eligible": "true" if pop_ok else "false",
            "genotype_status": gstatus, "determinant_status": dstatus,
            "determinant_count": len(pos) if result else NA,
            "determinant_ids": ";".join(pos) or NA,
            "determinant_types": ";".join(rule_types[s]["determinant_type"] for s in pos) or NA,
            "uncertain_determinant_ids": ";".join(result["uncertain"]) if result and result["uncertain"] else NA,
            "concordance_category": category, "fq_two_mechanisms_status": fq or NA,
            "mlst_st": p["mlst_st"], "lineage_id": p["lineage_id"], "lineage_component": g["lineage_component"],
            "phylogroup": p["phylogroup"], "ncbi_snp_cluster": p["ncbi_snp_cluster"], "country": e["country"],
            "collection_year": e["collection_year"], "period": g["period"], "geography_group": g["geography_group"],
            "source_category": e["source_category"], "project": p["project"],
            "mapping_version": mapping["version"] if study else NA,
            "caller_version": config["caller"]["version"] if call["status"] == "success" else NA,
            "database_id": database_id if call["status"] == "success" else NA,
            "endpoint_eligible": "false" if reasons else "true", "exclusion_reason": reasons[0] if reasons else NA,
            "_result": result, "_call": call, "_qc": p["genome_qc_status"]})

    # Denominators: every endpoint is counted once, under its first exclusion reason.
    denominators, metric_rows = [], []
    by_drug = defaultdict(list)
    for e in endpoints:
        by_drug[e["antimicrobial"]].append(e)
    for drug in sorted(by_drug):
        rows = by_drug[drug]
        first = Counter(e["exclusion_reason"].split(":")[0] for e in rows)
        pop_binary = [e for e in rows if e["eligible_binary_rs"] == "true" and e["population_eligible"] == "true"]
        gs = Counter(e["genotype_status"] for e in pop_binary)
        d = {"antimicrobial": drug, "study_tier": drugs.get(drug, "not_study"), "total_endpoints": len(rows),
             "binary_rs_endpoints": sum(e["eligible_binary_rs"] == "true" for e in rows),
             "excluded_non_binary_or_conflict": first["non_binary_or_conflict"],
             "population_eligible_binary": len(pop_binary),
             "excluded_population_ineligible": first["population_ineligible"],
             "genome_available": sum(e["_call"]["status"] not in {"genome_unavailable"} for e in pop_binary),
             "caller_successful": sum(e["_call"]["status"] == "success" for e in pop_binary),
             "mapping_interpretable": gs["determinant_positive"] + gs["determinant_negative"],
             "determinant_positive": gs["determinant_positive"], "determinant_negative": gs["determinant_negative"],
             "determinant_uncertain": gs["determinant_uncertain"],
             "excluded_not_study_antimicrobial": first["not_study_antimicrobial"],
             "excluded_genotype_unavailable": first["genotype_unavailable"],
             "excluded_mapping_uncertain": first["mapping_uncertain"], "metric_eligible": first[NA]}
        total = sum(d[k] for k in ("excluded_non_binary_or_conflict", "excluded_population_ineligible",
                                   "excluded_not_study_antimicrobial", "excluded_genotype_unavailable",
                                   "excluded_mapping_uncertain", "metric_eligible"))
        if total != len(rows):
            raise ValueError(f"denominators do not reconcile for {drug}")
        d["reconciled"] = "true"
        denominators.append(d)
        if drug not in drugs:
            continue
        base = [e for e in pop_binary if e["_result"] is not None]
        for variant in VARIANTS:
            if variant == "fq_two_mechanisms" and drug not in FLUOROQUINOLONES:
                continue
            pairs = [(e["phenotype"], variant_status(e["_result"], variant, drug)) for e in base]
            used = [(ph, st) for ph, st in pairs if st in {"positive", "negative"}]
            m = metrics(*confusion(used))
            metric_rows.append({"antimicrobial": drug, "study_tier": drugs[drug], "rule_variant": variant,
                                "excluded_uncertain": sum(st == "uncertain" for _, st in pairs),
                                **{k: fmt(m[k]) for k in METRIC_FIELDS}})

    eligible = [e for e in endpoints if e["endpoint_eligible"] == "true"]
    strata = stratum_rows(eligible, drugs, config["stratum_support"])

    def candidate(e):
        call = e["_call"]
        others = [d["symbol"] for d in call["dets"] if d["symbol"] not in set(e["_result"]["positive"])]
        return {**e, "other_amr_determinants": ";".join(others) or NA, "genome_qc_status": e["_qc"],
                "raw_report_sha256": call.get("raw_sha256", NA)}
    unexplained = [candidate(e) for e in eligible if e["concordance_category"] == "R_DET_NEG"]
    susceptible_pos = [candidate(e) for e in eligible if e["concordance_category"] == "S_DET_POS"]

    # Mapping audit over population-eligible binary endpoints with a successful call.
    audit_cfg = config["mapping_audit"]
    audit = []
    called = [e for e in endpoints if e["_result"] is not None and e["eligible_binary_rs"] == "true"
              and e["population_eligible"] == "true"]
    carriers = defaultdict(list)
    for e in called:
        rules = mapping["rules"].get(e["antimicrobial"], {})
        for d in e["_call"]["dets"]:
            if d["symbol"] in rules:
                carriers[(e["antimicrobial"], d["symbol"])].append(e)
    for (drug, sym), rows in sorted(carriers.items()):
        rule = mapping["rules"][drug][sym]
        r = sum(e["phenotype"] == "R" for e in rows)
        s = len(rows) - r
        sole = [e for e in rows if e["_result"]["positive"] == [sym]]
        frac = ratio(r, len(rows))
        flag = NA
        if rule["interpretation"] in {"predicts_resistance", "component"} and len(rows) >= audit_cfg["min_carriers"] \
                and frac is not None and frac < audit_cfg["low_r_fraction"]:
            flag = "HUMAN_REVIEW_REQUIRED: low R fraction among carriers"
        elif rule["interpretation"] in {"uncertain", "does_not_predict"} and len(rows) >= audit_cfg["min_carriers"] \
                and frac is not None and frac >= audit_cfg["high_r_fraction"]:
            flag = "HUMAN_REVIEW_REQUIRED: high R fraction for a non-predicting rule"
        audit.append({"antimicrobial": drug, "determinant": sym, "determinant_type": rule["determinant_type"],
                      "interpretation": rule["interpretation"], "requires": rule["requires"],
                      "evidence_level": rule["evidence_level"], "review_status": rule["review_status"],
                      "carriers": len(rows), "carriers_R": r, "carriers_S": s, "r_fraction": fmt(frac),
                      "sole_support_R": sum(e["phenotype"] == "R" for e in sole),
                      "sole_support_S": sum(e["phenotype"] == "S" for e in sole),
                      "lineages_json": pop.canonical(dict(Counter(e["lineage_id"] for e in rows).most_common(5))),
                      "countries_json": pop.canonical(dict(Counter(e["country"] for e in rows).most_common(5))),
                      "flag": flag})

    mapped_symbols = {s for rules in mapping["rules"].values() for s in rules}
    unmapped = defaultdict(lambda: {"isolates": 0, "population_eligible_isolates": 0})
    for iso, call in calls_by_isolate.items():
        for d in call["dets"]:
            if d["symbol"] not in mapped_symbols:
                u = unmapped[(d["symbol"], d["subtype"], d["class"])]
                u["isolates"] += 1
                u["population_eligible_isolates"] += popn["meta"][iso]["eligible_population_analysis"] == "true"
    unmapped_rows = [{"determinant": k[0], "subtype": k[1], "subclass": k[2], **v}
                     for k, v in sorted(unmapped.items())]

    # Summary
    call_counts = Counter(c["status"] for c in calls_by_isolate.values())
    study_binary = [e for e in endpoints if e["antimicrobial"] in drugs and e["eligible_binary_rs"] == "true"
                    and e["population_eligible"] == "true"]
    qc = Counter(p["genome_qc_status"] for p in popn["meta"].values())
    summary = {
        "cohort_isolates": len(candidates), "phenotype_endpoints": len(endpoints),
        "binary_rs_endpoints": sum(e["eligible_binary_rs"] == "true" for e in endpoints),
        "unique_assemblies": len({p["assembly_accession"] for p in popn["meta"].values()}),
        "genomes_available": sum(not pop.missing(p["genome_sha256"]) for p in popn["meta"].values()),
        "genome_qc": dict(sorted(qc.items())),
        "caller_status": dict(sorted(call_counts.items())),
        "study_antimicrobials": drugs,
        "study_endpoint_genotype_status": dict(sorted(Counter(e["genotype_status"] for e in study_binary).items())),
        "metric_eligible_endpoints": len(eligible),
        "concordance_categories": dict(sorted(Counter(e["concordance_category"] for e in eligible).items())),
        "unexplained_resistance_endpoints": len(unexplained),
        "susceptible_with_determinant_endpoints": len(susceptible_pos),
        "mapping_audit_flags": sum(a["flag"] != NA for a in audit),
        "unmapped_amr_determinants": len(unmapped_rows),
        "phenotype_corrections": dict(sorted(corrected.items())),
        "stratum_support": config["stratum_support"],
        "definitions": {
            "positive_reference": "phenotype R",
            "baseline_positive_call": "determinant_status positive under the mapping version",
            "determinant_negative": "no determinant satisfying the mapping was detected in a successfully "
                                    "analysed genome; not biological susceptibility",
            "R_DET_NEG": "candidate for unexplained resistance, not evidence of a novel mechanism",
        },
    }
    generated_at = call_times(calls) or now()
    manifest = {
        "baseline_version": VERSION, "builder_sha256": amrscan.sha256(Path(__file__)),
        "code_commit": code_commit or NA, "config_sha256": amrscan.sha256(config_path),
        "generated_at": generated_at,
        "generated_at_semantics": "latest successful AMRFinderPlus call used; fixed so rebuilds are deterministic",
        "pr5_cohort_id": cohort["manifest"].get("cohort_id", "ecoli_research_v1"),
        "pr5_cohort_manifest_sha256": cohort["manifest_sha256"], "pr5_files_sha256": cohort["file_sha256"],
        "pr6_population_id": popn["population_id"], "pr6_population_manifest_sha256": popn["manifest_sha256"],
        "pr6_code_commit": popn["manifest"].get("code_commit", NA),
        "cohort_isolate_count": len(candidates), "phenotype_endpoint_count": len(endpoints),
        "amrfinderplus_version": config["caller"]["version"],
        "amrfinderplus_database_version": config["caller"]["database_version"],
        "amrfinderplus_database_fingerprint": database_id, "caller_policy": CALL_POLICY,
        "mapping_version": mapping["version"], "mapping_rules_sha256": mapping["rules_sha256"],
        "mapping_metadata_sha256": mapping["metadata_sha256"],
        "phenotype_corrections_sha256": corr_cfg["sha256"] if corr_cfg else NA,
        "phenotype_corrected_endpoints": sum(corrected.values()),
        "caller_successful_isolates": call_counts["success"],
        "mapping_interpretable_endpoints": sum(e["genotype_status"] in {"determinant_positive",
                                                                       "determinant_negative"} for e in study_binary),
        "metric_eligible_endpoints": len(eligible),
    }

    outdir = Path(outdir)
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"output already exists: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".genotype-baseline-", dir=outdir.parent))
    try:
        write_table(staging / "genotype_phenotype_endpoints.tsv", ENDPOINT_FIELDS, endpoints)
        write_table(staging / "baseline_metrics.tsv", ["antimicrobial", "study_tier", "rule_variant",
                                                       "excluded_uncertain", *METRIC_FIELDS], metric_rows)
        write_table(staging / "baseline_denominators.tsv", DENOM_FIELDS, denominators)
        for filename in sorted({f for f, _ in DIMENSIONS.values()}):
            write_table(staging / filename, STRATUM_FIELDS, strata.get(filename, []))
        write_table(staging / "unexplained_resistance.tsv", CANDIDATE_FIELDS, unexplained)
        write_table(staging / "susceptible_with_determinant.tsv", CANDIDATE_FIELDS, susceptible_pos)
        write_table(staging / "mapping_audit.tsv", AUDIT_FIELDS, audit)
        write_table(staging / "unmapped_determinants.tsv", UNMAPPED_FIELDS, unmapped_rows)
        for name, value in (("baseline_summary.json", summary), ("baseline_manifest.json", manifest)):
            (staging / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        names = sorted(p.name for p in staging.iterdir())
        (staging / "SHA256SUMS").write_text("".join(f"{amrscan.sha256(staging / n)}  {n}\n" for n in names),
                                            encoding="utf-8")
        staging.rename(outdir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return {"summary": summary, "manifest": manifest}


# ---------------------------------------------------------------- command line

def git_commit() -> str | None:
    import subprocess
    try:
        return subprocess.run(["git", "-C", str(Path(__file__).resolve().parents[1]), "rev-parse", "HEAD"],
                              check=True, text=True, capture_output=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("call", help="resumable AMRFinderPlus calls on cached genomes")
    c.add_argument("--population", type=Path, required=True, help="PR #6 population build directory")
    c.add_argument("--cache", type=Path, required=True, help="PR #6 genome cache")
    c.add_argument("--calls", type=Path, required=True)
    c.add_argument("--database", type=Path, required=True)
    c.add_argument("--database-source")
    c.add_argument("--executable", default="amrfinder")
    c.add_argument("--workers", type=int, default=4)
    c.add_argument("--threads", type=int, default=1)
    c.add_argument("--limit", type=int)
    b = sub.add_parser("build", help="genotype-phenotype baseline from validated calls")
    b.add_argument("--processed", type=Path, required=True, help="PR #5 processed cohort directory")
    b.add_argument("--population", type=Path, required=True, help="PR #6 population build directory")
    b.add_argument("--calls", type=Path, required=True)
    b.add_argument("--database", type=Path, required=True)
    b.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    b.add_argument("--outdir", type=Path, required=True)
    b.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        if args.workers < 1:
            raise ValueError("workers must be positive")
        if args.command == "call":
            if args.threads < 1:
                raise ValueError("threads must be positive")
            result = call_cohort(args.population, args.cache, args.calls, args.database, args.executable,
                                 args.workers, args.threads, args.limit, args.database_source)
        else:
            result = build(args.processed, args.population, args.calls, args.database, args.config, args.outdir,
                           args.workers, git_commit())["summary"]
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ValueError, KeyError, OSError, csv.Error, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
