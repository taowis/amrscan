#!/usr/bin/env python3
"""E. coli population structure metadata for the frozen PR #5 research cohort.

Commands:
  fetch-genomes  resumable, accession-verified NCBI assembly cache
  annotate       genome QC, MLST (mlst) or phylogroup (EzClermont) per cached assembly
  build          population metadata, lineage/cluster summaries, grouping
                 variables, holdout feasibility and review flags

Lineages and clusters come from genomes and NCBI metadata only; phenotype is
joined afterwards for description. Nothing here trains or evaluates a model.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from urllib.error import URLError

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genome_dataset  # noqa: E402

VERSION = "1.1.0"
NA = "NA"
csv.field_size_limit(sys.maxsize)
DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "benchmark/cohorts/ecoli_research_v1/population_config.json"
COUNTRY_MISSING = {"", NA, "unknown", "UNKNOWN"}
PDS = re.compile(r"^(PDS\d+)\.\d+$")
DOWNLOAD = "/download?include_annotation_type=GENOME_FASTA,SEQUENCE_REPORT&filename=assembly.zip"
COHORT_FILES = ("candidate_isolates.tsv", "isolate_drug_endpoints.tsv", "cohort_manifest.json")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    return genome_dataset.sha256(path)


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def missing(value) -> bool:
    return value is None or str(value).strip() in {"", NA}


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path}: malformed row width")
    return rows


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json_atomic(path: Path, value) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def load_config(path: Path) -> dict:
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("genome_qc", "mlst", "phylogroup", "lineage", "holdout_eligibility", "feasibility_rules", "review_flags",
                "evaluation_design"):
        if key not in config:
            raise ValueError(f"population config lacks {key}")
    if config["feasibility_rules"].get("decision") != "HUMAN_REVIEW_REQUIRED":
        raise ValueError("feasibility decision must remain HUMAN_REVIEW_REQUIRED until reviewed")
    return config


# ---------------------------------------------------------------- PR #5 cohort

def load_cohort(processed: Path, expected_manifest_sha256: str | None = None) -> dict:
    """Read the frozen PR #5 build after checking its own SHA256SUMS."""
    processed = Path(processed)
    sums = {}
    for line in (processed / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            sums[name.strip()] = digest
    for name in COHORT_FILES:
        if name not in sums or sha256(processed / name) != sums[name]:
            raise ValueError(f"PR #5 cohort file missing or changed: {name}")
    manifest_sha = sha256(processed / "cohort_manifest.json")
    if expected_manifest_sha256 and manifest_sha != expected_manifest_sha256:
        raise ValueError("PR #5 cohort manifest hash differs from the expected value")
    manifest = json.loads((processed / "cohort_manifest.json").read_text(encoding="utf-8"))
    candidates = read_tsv(processed / "candidate_isolates.tsv")
    endpoints = read_tsv(processed / "isolate_drug_endpoints.tsv")
    ids = [r["isolate_id"] for r in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate isolate in PR #5 candidates")
    assemblies = Counter(r["assembly_accession"] for r in candidates if not missing(r["assembly_accession"]))
    duplicates = sorted(a for a, n in assemblies.items() if n > 1)
    if int(manifest.get("candidate_isolate_count", -1)) != len(candidates):
        raise ValueError("PR #5 manifest isolate count differs from candidate table")
    return {"processed": processed, "manifest": manifest, "manifest_sha256": manifest_sha,
            "candidates": candidates, "endpoints": endpoints, "duplicate_assemblies": duplicates,
            "file_sha256": {name: sums[name] for name in COHORT_FILES}}


def snp_cluster(candidate: dict) -> tuple[str, str, list[str]]:
    """Unversioned NCBI SNP cluster accession, its status and the raw versions."""
    try:
        raw = json.loads(candidate.get("snp_clusters_json") or "[]")
    except json.JSONDecodeError as error:
        raise ValueError(f"malformed snp_clusters_json for {candidate.get('isolate_id')}") from error
    raw = sorted({str(x) for x in raw if not missing(x)})
    bases = sorted({PDS.match(x).group(1) if PDS.match(x) else x for x in raw})
    if not bases:
        return NA, "missing", raw
    if len(bases) > 1:
        return NA, "conflicting", raw
    return bases[0], "assigned", raw


# ---------------------------------------------------------------- genome cache

def cached_genome(cache: Path, assembly: str) -> dict | None:
    """Return verified cache provenance, or None if absent. Raise if corrupt."""
    entry = Path(cache) / assembly
    prov = entry / "provenance.json"
    if not prov.exists():
        return None
    record = json.loads(prov.read_text(encoding="utf-8"))
    genome = entry / "genome.fna"
    if record.get("assembly_accession") != assembly or not genome.is_file() or sha256(genome) != record.get("sha256"):
        raise ValueError(f"cached assembly failed checksum: {assembly}")
    return record


def fetch_one(row: dict, cache: Path, retries: int = 3, opener=None, pause: float = 5.0) -> dict:
    assembly = row["assembly_accession"]
    cache = Path(cache)
    try:
        record = cached_genome(cache, assembly)
        if record:
            return {"assembly_accession": assembly, "status": "cached"}
    except (ValueError, json.JSONDecodeError):
        shutil.rmtree(cache / assembly, ignore_errors=True)
    opener = opener or genome_dataset.request_bytes
    url = genome_dataset.API.format(assembly) + DOWNLOAD
    error = None
    for attempt in range(1, retries + 1):
        try:
            started = now()
            archive, final_url = opener(url, genome_dataset.MAX_ARCHIVE)
            fasta, member, records, bases = genome_dataset.package_identity(archive, row)
            break
        except (URLError, OSError, ValueError, zipfile.BadZipFile, UnicodeError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(pause * attempt)
    else:
        return {"assembly_accession": assembly, "status": "failed", "error": error, "attempts": retries}
    staging = Path(tempfile.mkdtemp(prefix=f".{assembly}-", dir=cache))
    try:
        (staging / "genome.fna").write_bytes(fasta)
        write_json_atomic(staging / "provenance.json", {
            "assembly_accession": assembly, "biosample_accession": row["biosample_accession"],
            "organism": row["organism"], "isolate_id": row.get("isolate_id", NA),
            "sha256": sha256(staging / "genome.fna"), "package_sha256": hashlib.sha256(archive).hexdigest(),
            "package_member": member, "sequence_records": records, "sequence_bases": bases,
            "source_url": final_url, "retrieved_at": started, "attempts": attempt,
            "retrieval": f"population_structure.py {VERSION} via genome_dataset.package_identity"})
        staging.rename(cache / assembly)
    finally:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
    return {"assembly_accession": assembly, "status": "downloaded", "attempts": attempt}


def fetch_genomes(cohort: dict, cache: Path, workers: int = 3, limit: int | None = None,
                  opener=None, pause: float = 5.0) -> dict:
    if not 1 <= workers <= 5:
        raise ValueError("workers must be between 1 and 5 to respect NCBI request limits")
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    rows = sorted((r for r in cohort["candidates"] if genome_dataset.ACCESSION.fullmatch(r["assembly_accession"])),
                  key=lambda r: r["assembly_accession"])
    seen, todo = set(), []
    for r in rows:
        if r["assembly_accession"] not in seen:
            seen.add(r["assembly_accession"])
            todo.append({"assembly_accession": r["assembly_accession"], "biosample_accession": r["biosample_accession"],
                         "organism": r["organism"], "isolate_id": r["isolate_id"]})
    if limit is not None:
        todo = todo[:limit]
    started = now()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(lambda r: fetch_one(r, cache, opener=opener, pause=pause), todo))
    summary = {"started_at": started, "finished_at": now(), "requested": len(todo),
               "status_counts": dict(sorted(Counter(r["status"] for r in results).items())),
               "failures": [r for r in results if r["status"] == "failed"]}
    with (cache / "fetch_log.jsonl").open("a", encoding="utf-8") as log:
        log.write(canonical(summary) + "\n")
    return summary


# ---------------------------------------------------------------- genome QC

def fasta_stats(path: Path) -> dict:
    lengths, counts, current, header = [], Counter(), 0, False
    with Path(path).open("rb") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if line.startswith(b">"):
                if header and current == 0:
                    raise ValueError(f"empty FASTA record before line {number}")
                if not line[1:].strip():
                    raise ValueError(f"empty FASTA header at line {number}")
                if header:
                    lengths.append(current)
                header, current = True, 0
            elif line:
                if not header or not re.fullmatch(rb"[A-Za-z*.-]+", line):
                    raise ValueError(f"invalid FASTA sequence at line {number}")
                current += len(line)
                upper = line.upper()
                for base in b"ACGT":
                    counts[base] += upper.count(base)
    if header:
        if current == 0:
            raise ValueError("empty FASTA record at end of file")
        lengths.append(current)
    if not lengths:
        raise ValueError("FASTA contains no sequence")
    total = sum(lengths)
    acgt = sum(counts[c] for c in b"ACGT")
    ordered, running, n50 = sorted(lengths, reverse=True), 0, 0
    for length in ordered:
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return {"assembly_length": total, "contig_count": len(lengths), "n50": n50, "longest_contig": ordered[0],
            "gc_fraction": round((counts[ord("G")] + counts[ord("C")]) / acgt, 6) if acgt else 0.0,
            "ambiguous_fraction": round((total - acgt) / total, 6)}


def qc_status(stats: dict | None, config: dict, problem: str | None = None) -> tuple[str, list[str]]:
    if problem:
        return "fail", [problem]
    q, flags = config["genome_qc"], []
    if not q["fail_length_min"] <= stats["assembly_length"] <= q["fail_length_max"]:
        flags.append("assembly_size_fail")
    elif not q["warn_length_min"] <= stats["assembly_length"] <= q["warn_length_max"]:
        flags.append("assembly_size_warning")
    if stats["ambiguous_fraction"] > q["fail_ambiguous_fraction_max"]:
        flags.append("excessive_ambiguous_bases_fail")
    elif stats["ambiguous_fraction"] > q["warn_ambiguous_fraction_max"]:
        flags.append("excessive_ambiguous_bases")
    if stats["contig_count"] > q["warn_contigs_max"]:
        flags.append("excessive_contigs")
    if stats["n50"] < q["warn_n50_min"]:
        flags.append("low_n50")
    if any(f.endswith("_fail") for f in flags):
        return "fail", flags
    return ("warning" if flags else "pass"), flags or ["qc_pass"]


def genome_qc(cache: Path, assembly: str, config: dict) -> dict:
    row = {"assembly_accession": assembly, "genome_sha256": NA, "assembly_length": NA, "contig_count": NA,
           "n50": NA, "longest_contig": NA, "gc_fraction": NA, "ambiguous_fraction": NA}
    try:
        record = cached_genome(cache, assembly)
    except (ValueError, json.JSONDecodeError):
        status, flags = qc_status(None, config, "assembly_checksum_failed")
        return {**row, "genome_qc_status": status, "genome_qc_flags": ";".join(flags)}
    if record is None:
        status, flags = qc_status(None, config, "assembly_missing")
        return {**row, "genome_qc_status": status, "genome_qc_flags": ";".join(flags)}
    try:
        stats = fasta_stats(Path(cache) / assembly / "genome.fna")
    except ValueError:
        status, flags = qc_status(None, config, "malformed_fasta")
        return {**row, "genome_sha256": record["sha256"], "genome_qc_status": status, "genome_qc_flags": ";".join(flags)}
    status, flags = qc_status(stats, config)
    return {**row, **{k: str(v) for k, v in stats.items()}, "genome_sha256": record["sha256"],
            "genome_qc_status": status, "genome_qc_flags": ";".join(flags)}


# ---------------------------------------------------------------- MLST / phylogroup

def tool_info(tool: str, env_bin: Path, config: dict) -> dict:
    env_bin = Path(env_bin)
    env = {**os.environ, "PATH": f"{env_bin}{os.pathsep}{os.environ.get('PATH', '')}"}
    def version(*cmd):
        out = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    if tool == "mlst":
        scheme = config["mlst"]["scheme"]
        db = env_bin.parent / "db"
        profile = db / "pubmlst" / scheme / f"{scheme}.txt"
        info_file = db / "pubmlst" / scheme / f"{scheme}_info.json"
        if not profile.is_file():
            raise ValueError(f"MLST scheme not installed: {scheme}")
        return {"tool": "mlst", "version": version(str(env_bin / "mlst"), "--version"),
                "blastn": version(str(env_bin / "blastn"), "-version"), "scheme": scheme,
                "database_version": (db / "VERSION.txt").read_text().strip() if (db / "VERSION.txt").exists() else NA,
                "scheme_info": json.loads(info_file.read_text()) if info_file.exists() else {},
                "scheme_profile_sha256": sha256(profile)}
    if tool == "phylogroup":
        return {"tool": "ezclermont", "version": version(str(env_bin / "ezclermont"), "--version"),
                "method": config["phylogroup"]["method"]}
    raise ValueError(f"unknown annotation tool: {tool}")


def parse_mlst(text: str, config: dict) -> dict:
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError("malformed mlst JSON") from error
    if not isinstance(data, list) or len(data) != 1 or not isinstance(data[0], dict):
        raise ValueError("mlst JSON must contain one result")
    item = data[0]
    scheme, loci = config["mlst"]["scheme"], config["mlst"]["loci"]
    alleles = item.get("alleles") or {}
    if item.get("scheme") not in {scheme, "-"} or (alleles and set(alleles) != set(loci)):
        raise ValueError("mlst result uses an unexpected scheme or loci")
    values = [str(alleles.get(locus, "-")) for locus in loci]
    profile = ";".join(f"{locus}({value})" for locus, value in zip(loci, values))
    st = str(item.get("sequence_type", "-"))
    exact = all(re.fullmatch(r"\d+", v) for v in values)
    if exact and re.fullmatch(r"\d+", st):
        status = "assigned"
    elif exact:
        status = "novel_profile"
    elif any(v == "-" for v in values):
        status = "missing_alleles"
    elif any("," in v for v in values):
        status = "multiple_alleles"
    else:
        status = "novel_or_partial_alleles"
    return {"mlst_scheme": scheme, "mlst_st": st if status == "assigned" else NA, "mlst_status": status,
            "mlst_profile": profile, "mlst_raw_st": st}


def parse_phylogroup(stdout: str, config: dict) -> dict:
    lines = [line for line in stdout.splitlines() if line.strip()]
    if not lines or len(lines[-1].split("\t")) != 2:
        raise ValueError("malformed EzClermont output")
    raw = lines[-1].split("\t")[1].strip()
    if raw in config["phylogroup"]["assigned"]:
        return {"phylogroup": raw, "phylogroup_status": "assigned", "phylogroup_raw": raw}
    if raw in config["phylogroup"]["cryptic"]:
        return {"phylogroup": "cryptic", "phylogroup_status": "cryptic_clade", "phylogroup_raw": raw}
    status = "control_failed" if raw == "EC_control_fail" else "unresolved"
    return {"phylogroup": NA, "phylogroup_status": status, "phylogroup_raw": raw}


def annotate_one(tool: str, assembly: str, cache: Path, env_bin: Path, info: dict, config: dict, runner=None) -> str:
    """Run one tool on one cached genome; reuse a result made with the same input and tool."""
    runner = runner or subprocess.run
    entry = Path(cache) / assembly
    record = cached_genome(cache, assembly)
    if record is None:
        return "missing_genome"
    out_path = entry / f"{tool}.json"
    fingerprint = hashlib.sha256(canonical(info).encode()).hexdigest()
    if out_path.exists():
        old = json.loads(out_path.read_text(encoding="utf-8"))
        if old.get("input_sha256") == record["sha256"] and old.get("tool_fingerprint") == fingerprint:
            return "cached"
    env = {**os.environ, "PATH": f"{Path(env_bin)}{os.pathsep}{os.environ.get('PATH', '')}"}
    genome = entry / "genome.fna"
    with tempfile.TemporaryDirectory(prefix=f".{tool}-", dir=entry) as tmp:
        if tool == "mlst":
            json_out = Path(tmp) / "mlst.json"
            cmd = [str(Path(env_bin) / "mlst"), "--scheme", config["mlst"]["scheme"], "--threads", "1",
                   "--json", str(json_out), str(genome)]
        else:
            cmd = [str(Path(env_bin) / "ezclermont"), str(genome)]
        started = now()
        proc = runner(cmd, capture_output=True, text=True, env=env, check=False)
        raw = json_out.read_text(encoding="utf-8") if tool == "mlst" and json_out.exists() else proc.stdout
    result = {"tool_fingerprint": fingerprint, "tool_info": info, "input_sha256": record["sha256"],
              "assembly_accession": assembly, "command": cmd, "returncode": proc.returncode,
              "started_at": started, "raw_output": raw, "stderr_tail": proc.stderr[-2000:]}
    try:
        if proc.returncode != 0:
            raise ValueError(f"exit code {proc.returncode}")
        result["parsed"] = parse_mlst(raw, config) if tool == "mlst" else parse_phylogroup(raw, config)
        result["status"] = "success"
    except ValueError as error:
        result["status"], result["error"] = "failed", str(error)
    write_json_atomic(out_path, result)
    return result["status"]


def annotate(tool: str, cohort: dict, cache: Path, env_bin: Path, config: dict, workers: int = 4,
             runner=None, info: dict | None = None) -> dict:
    info = info or tool_info(tool, env_bin, config)
    assemblies = sorted({r["assembly_accession"] for r in cohort["candidates"]
                         if genome_dataset.ACCESSION.fullmatch(r["assembly_accession"])})
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        statuses = list(pool.map(lambda a: annotate_one(tool, a, cache, env_bin, info, config, runner), assemblies))
    return {"tool": info, "status_counts": dict(sorted(Counter(statuses).items()))}


def read_annotation(cache: Path, assembly: str, tool: str, genome_sha: str) -> dict | None:
    path = Path(cache) / assembly / f"{tool}.json"
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("input_sha256") != genome_sha:
        raise ValueError(f"{tool} result does not match cached genome for {assembly}")
    return data


# ---------------------------------------------------------------- build

def valid_year(value: str) -> bool:
    return bool(re.fullmatch(r"\d{4}", value or ""))


def ratio(a: int, b: int) -> str:
    return f"{a / b:.6f}" if b else NA


def population_rows(cohort: dict, cache: Path, config: dict, workers: int = 4) -> tuple[list[dict], list[dict], list[dict], list[dict], dict]:
    elig = config["holdout_eligibility"]
    known_sources = set(elig["known_source_categories"])
    per_isolate = defaultdict(lambda: {"binary": 0, "drugs": set(), "eligible": False})
    for e in cohort["endpoints"]:
        p = per_isolate[e["isolate_id"]]
        p["drugs"].add(e["antimicrobial"])
        if e["eligible_overall"] == "true":
            p["eligible"] = True
            p["binary"] += 1
    rows, qc_rows, mlst_rows, phylo_rows = [], [], [], []
    tools = {"mlst": set(), "phylogroup": set()}
    candidates = sorted(cohort["candidates"], key=lambda r: r["isolate_id"])
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        qcs = list(pool.map(lambda c: genome_qc(cache, c["assembly_accession"], config), candidates))
    for c, qc in zip(candidates, qcs):
        assembly = c["assembly_accession"]
        qc_rows.append({"isolate_id": c["isolate_id"], **qc})
        annotations = {}
        for tool in ("mlst", "phylogroup"):
            data = read_annotation(cache, assembly, tool, qc["genome_sha256"]) if qc["genome_sha256"] != NA else None
            annotations[tool] = data
            if data:
                tools[tool].add(canonical(data["tool_info"]))
        m = annotations["mlst"]
        if m and m["status"] == "success":
            mp = m["parsed"]
        else:
            mp = {"mlst_scheme": config["mlst"]["scheme"], "mlst_st": NA, "mlst_profile": NA, "mlst_raw_st": NA,
                  "mlst_status": "failed" if m else "not_run"}
        mlst_rows.append({"isolate_id": c["isolate_id"], "assembly_accession": assembly, **mp,
                          "software": m["tool_info"]["tool"] if m else NA,
                          "software_version": m["tool_info"]["version"] if m else NA,
                          "database_version": m["tool_info"].get("database_version", NA) if m else NA,
                          "scheme_profile_sha256": m["tool_info"].get("scheme_profile_sha256", NA) if m else NA,
                          "input_sha256": qc["genome_sha256"]})
        ph = annotations["phylogroup"]
        if ph and ph["status"] == "success":
            pp = ph["parsed"]
        else:
            pp = {"phylogroup": NA, "phylogroup_raw": NA, "phylogroup_status": "failed" if ph else "not_run"}
        phylo_rows.append({"isolate_id": c["isolate_id"], "assembly_accession": assembly, **pp,
                           "method": config["phylogroup"]["method"],
                           "software_version": ph["tool_info"]["version"] if ph else NA,
                           "input_sha256": qc["genome_sha256"]})
        cluster, cluster_status, cluster_raw = snp_cluster(c)
        p = per_isolate[c["isolate_id"]]
        population = p["eligible"] and qc["genome_qc_status"] != "fail"
        lineage = f"ST{mp['mlst_st']}" if mp["mlst_status"] == "assigned" else "UNKNOWN"
        country_ok = c["country"] not in COUNTRY_MISSING
        rows.append({
            "cohort_id": c["cohort_id"], "isolate_id": c["isolate_id"], "pdt_accession": c["pdt_accession"],
            "biosample_accession": c["biosample_accession"], "assembly_accession": assembly,
            "organism": c["organism"], "strain": c["strain"],
            "genome_sha256": qc["genome_sha256"], "genome_qc_status": qc["genome_qc_status"],
            "genome_qc_flags": qc["genome_qc_flags"], "assembly_length": qc["assembly_length"],
            "contig_count": qc["contig_count"], "n50": qc["n50"], "gc_fraction": qc["gc_fraction"],
            "ambiguous_fraction": qc["ambiguous_fraction"],
            "mlst_scheme": mp["mlst_scheme"], "mlst_st": mp["mlst_st"], "mlst_status": mp["mlst_status"],
            "mlst_profile": mp["mlst_profile"], "phylogroup": pp["phylogroup"],
            "phylogroup_status": pp["phylogroup_status"], "ncbi_snp_cluster": cluster,
            "ncbi_snp_cluster_status": cluster_status, "ncbi_snp_cluster_versions_json": canonical(cluster_raw),
            "lineage_id": lineage, "country": c["country"], "collection_year": c["collection_year"],
            "source_category": c["source_category"], "submitter": c["submitter"], "project": c["project"],
            "binary_rs_endpoint_count": str(p["binary"]), "antimicrobial_count": str(len(p["drugs"])),
            "pr5_eligible_isolate": "true" if p["eligible"] else "false",
            "eligible_population_analysis": "true" if population else "false",
            "eligible_lineage_holdout": "true" if population and lineage != "UNKNOWN" else "false",
            "eligible_country_holdout": "true" if population and country_ok else "false",
            "eligible_temporal_holdout": "true" if population and valid_year(c["collection_year"]) else "false",
            "eligible_source_holdout": "true" if population and c["source_category"] in known_sources else "false",
        })
    return rows, qc_rows, mlst_rows, phylo_rows, {k: sorted(v) for k, v in tools.items()}


def grouped_counts(rows: list[dict], key: str, value: str) -> list[dict]:
    counts = Counter((r[key], r[value] if not missing(r[value]) else NA) for r in rows)
    return [{key: k, value: v, "isolate_count": str(n)} for (k, v), n in sorted(counts.items())]


def lineage_summary(rows: list[dict], binary: list[dict], config: dict) -> list[dict]:
    total = len(rows)
    by = defaultdict(list)
    for r in rows:
        by[r["lineage_id"]].append(r)
    endpoints = Counter(e["lineage_id"] for e in binary)
    out = []
    for lineage, group in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        countries = Counter(r["country"] for r in group if r["country"] not in COUNTRY_MISSING)
        years = sorted(int(r["collection_year"]) for r in group if valid_year(r["collection_year"]))
        sources = {r["source_category"] for r in group if not missing(r["source_category"])}
        top_country, top_n = countries.most_common(1)[0] if countries else (NA, 0)
        out.append({
            "lineage_id": lineage, "is_lineage": "false" if lineage == "UNKNOWN" else "true",
            "isolate_count": str(len(group)), "fraction_of_cohort": ratio(len(group), total),
            "population_eligible_isolates": str(sum(r["eligible_population_analysis"] == "true" for r in group)),
            "country_count": str(len(countries)), "largest_country": top_country,
            "fraction_from_largest_country": ratio(top_n, sum(countries.values())),
            "earliest_year": str(years[0]) if years else NA, "latest_year": str(years[-1]) if years else NA,
            "year_span": str(years[-1] - years[0]) if years else NA,
            "source_category_count": str(len(sources)),
            "antimicrobial_endpoint_count": str(endpoints[lineage]),
            "phylogroups_json": canonical(dict(sorted(Counter(r["phylogroup"] for r in group).items()))),
        })
    return out


def lineage_antimicrobial(binary: list[dict]) -> list[dict]:
    counts = defaultdict(Counter)
    for e in binary:
        counts[(e["antimicrobial"], e["lineage_id"])][e["phenotype_endpoint"]] += 1
    return [{"antimicrobial": drug, "lineage_id": lineage, "n_total": str(c["R"] + c["S"]), "n_R": str(c["R"]),
             "n_S": str(c["S"]), "resistance_fraction": ratio(c["R"], c["R"] + c["S"])}
            for (drug, lineage), c in sorted(counts.items())]


def snp_cluster_summary(rows: list[dict], binary: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        if r["ncbi_snp_cluster"] != NA:
            by[r["ncbi_snp_cluster"]].append(r)
    rs = defaultdict(Counter)
    for e in binary:
        if e["ncbi_snp_cluster"] != NA:
            rs[e["ncbi_snp_cluster"]][e["phenotype_endpoint"]] += 1
    out = []
    for cluster, group in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        years = sorted(int(r["collection_year"]) for r in group if valid_year(r["collection_year"]))
        out.append({"ncbi_snp_cluster": cluster, "isolate_count": str(len(group)),
                    "countries_json": canonical(dict(sorted(Counter(r["country"] for r in group).items()))),
                    "earliest_year": str(years[0]) if years else NA, "latest_year": str(years[-1]) if years else NA,
                    "lineages_json": canonical(dict(sorted(Counter(r["lineage_id"] for r in group).items()))),
                    "projects_json": canonical(dict(sorted(Counter(r["project"] for r in group).items()))),
                    "binary_endpoints_R": str(rs[cluster]["R"]), "binary_endpoints_S": str(rs[cluster]["S"])})
    return out


def project_summary(rows: list[dict], binary: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        by[(r["project"], r["submitter"])].append(r)
    rs = defaultdict(Counter)
    for e in binary:
        rs[e["project"]][e["phenotype_endpoint"]] += 1
    out = []
    for (project, submitter), group in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        lineages = Counter(r["lineage_id"] for r in group)
        out.append({"project": project, "submitter": submitter, "isolate_count": str(len(group)),
                    "country_count": str(len({r["country"] for r in group if r["country"] not in COUNTRY_MISSING})),
                    "countries_json": canonical(dict(sorted(Counter(r["country"] for r in group).items()))),
                    "lineage_count": str(len(set(lineages) - {"UNKNOWN"})),
                    "top_lineages_json": canonical(lineages.most_common(5)),
                    "binary_endpoints_R": str(rs[project]["R"]), "binary_endpoints_S": str(rs[project]["S"])})
    return out


def recent_cutoff(years: list[int], quantile: float) -> int | None:
    if not years:
        return None
    ordered = sorted(years)
    return ordered[min(len(ordered) - 1, int(math.floor(quantile * len(ordered))))]


def temporal_cutoff(years: list[int], rule: dict) -> int | None:
    """Fixed cutoff year when the config sets one, otherwise the configured quantile."""
    return rule["cutoff_year"] if "cutoff_year" in rule else recent_cutoff(years, rule["recent_quantile"])


def classify(metrics: dict, feasible: dict, marginal: dict) -> str:
    def meets(rule):
        for key, limit in rule.items():
            name = key.removeprefix("min_").removeprefix("max_")
            value = metrics[name]
            if key.startswith("min_") and value < limit or key.startswith("max_") and value > limit:
                return False
        return True
    return "feasible" if meets(feasible) else "marginal" if meets(marginal) else "not_feasible"


def holdout_feasibility(binary: list[dict], rows: list[dict], config: dict) -> list[dict]:
    rules = config["feasibility_rules"]
    both_min = rules["min_both_per_group"]
    years_all = [int(r["collection_year"]) for r in rows
                 if r["eligible_temporal_holdout"] == "true" and valid_year(r["collection_year"])]
    cutoff = temporal_cutoff(years_all, rules["temporal"])
    out = []
    by_drug = defaultdict(list)
    for e in binary:
        by_drug[e["antimicrobial"]].append(e)
    for drug, eps in sorted(by_drug.items()):
        def group_stats(field, eligible_field, skip=()):
            g = defaultdict(Counter)
            for e in eps:
                if e[eligible_field] == "true" and e[field] not in skip:
                    g[e[field]][e["phenotype_endpoint"]] += 1
            return g
        # lineage
        g = group_stats("lineage_id", "eligible_lineage_holdout", {"UNKNOWN"})
        n = sum(sum(c.values()) for c in g.values())
        nR = sum(c["R"] for c in g.values())
        largest = max((sum(c.values()) for c in g.values()), default=0)
        major = sorted(g.items(), key=lambda kv: -sum(kv[1].values()))[:5]
        metrics = {"R": nR, "S": n - nR, "sts_with_both": sum(c["R"] > 0 and c["S"] > 0 for c in g.values()),
                   "largest_st_fraction": largest / n if n else 1.0}
        out.append({"antimicrobial": drug, "holdout": "lineage", "eligible_isolates": str(n), "groups": str(len(g)),
                    "largest_group_size": str(largest), "largest_group_fraction": ratio(largest, n),
                    "groups_with_both_R_S": str(metrics["sts_with_both"]), "n_R": str(nR), "n_S": str(n - nR),
                    "detail_json": canonical({k: [c["R"], c["S"]] for k, c in major}),
                    "status": classify(metrics, rules["lineage"]["feasible"], rules["lineage"]["marginal"]),
                    "decision": rules["decision"]})
        # country and source
        for holdout, field, eligible_field in (("country", "country", "eligible_country_holdout"),
                                               ("source", "source_category", "eligible_source_holdout")):
            g = group_stats(field, eligible_field)
            n = sum(sum(c.values()) for c in g.values())
            nR = sum(c["R"] for c in g.values())
            largest = max((sum(c.values()) for c in g.values()), default=0)
            size = rules[holdout]["min_group_isolates"]
            both = sum(sum(c.values()) >= size and c["R"] >= both_min and c["S"] >= both_min for c in g.values())
            metrics = {"groups_with_both": both, "largest_fraction": largest / n if n else 1.0}
            out.append({"antimicrobial": drug, "holdout": holdout, "eligible_isolates": str(n), "groups": str(len(g)),
                        "largest_group_size": str(largest), "largest_group_fraction": ratio(largest, n),
                        "groups_with_both_R_S": str(both), "n_R": str(nR), "n_S": str(n - nR),
                        "detail_json": canonical({k: [c["R"], c["S"]] for k, c in sorted(g.items())}),
                        "status": classify(metrics, rules[holdout]["feasible"], rules[holdout]["marginal"]),
                        "decision": rules["decision"]})
        # temporal
        t = [e for e in eps if e["eligible_temporal_holdout"] == "true"]
        recent = Counter(e["phenotype_endpoint"] for e in t if cutoff is not None and int(e["collection_year"]) >= cutoff)
        earlier = Counter(e["phenotype_endpoint"] for e in t if cutoff is not None and int(e["collection_year"]) < cutoff)
        years = sorted({int(e["collection_year"]) for e in t})
        missing_year = sum(e["eligible_population_analysis"] == "true" and not valid_year(e["collection_year"]) for e in eps)
        metrics = {"recent": sum(recent.values()), "recent_R": recent["R"], "recent_S": recent["S"]}
        out.append({"antimicrobial": drug, "holdout": "temporal", "eligible_isolates": str(len(t)),
                    "groups": str(len(years)), "largest_group_size": NA, "largest_group_fraction": NA,
                    "groups_with_both_R_S": NA, "n_R": str(recent["R"] + earlier["R"]),
                    "n_S": str(recent["S"] + earlier["S"]),
                    "detail_json": canonical({"cutoff_year": cutoff, "year_range": [years[0], years[-1]] if years else None,
                                              "recent_R_S": [recent["R"], recent["S"]],
                                              "earlier_R_S": [earlier["R"], earlier["S"]],
                                              "missing_year_endpoints": missing_year}),
                    "status": classify(metrics, rules["temporal"]["feasible"], rules["temporal"]["marginal"]),
                    "decision": rules["decision"]})
    return out


def review_flags(rows: list[dict], binary: list[dict], config: dict, cutoff: int | None) -> list[dict]:
    f = config["review_flags"]
    flags = []
    def add(risk, subject, metric, value, threshold, detail):
        flags.append({"risk": risk, "subject": subject, "metric": metric, "value": f"{value:.6f}" if isinstance(value, float) else str(value),
                      "threshold": str(threshold), "detail": detail})
    population = [r for r in rows if r["eligible_population_analysis"] == "true"]
    total = len(population)
    by_lineage = defaultdict(list)
    for r in population:
        by_lineage[r["lineage_id"]].append(r)
    for lineage, group in sorted(by_lineage.items()):
        if lineage == "UNKNOWN" or len(group) < f["min_lineage_isolates"]:
            continue
        n = len(group)
        if n / total >= f["dominant_lineage_fraction"]:
            add("dominant_lineage", lineage, "fraction_of_population", n / total, f["dominant_lineage_fraction"], f"{n} isolates")
        for risk, field, limit in (("lineage_country_concentration", "country", f["lineage_single_country_fraction"]),
                                   ("lineage_source_concentration", "source_category", f["lineage_single_source_fraction"]),
                                   ("lineage_project_concentration", "project", f["lineage_single_project_fraction"])):
            known = Counter(r[field] for r in group if r[field] not in COUNTRY_MISSING)
            if known:
                top, k = known.most_common(1)[0]
                share = k / sum(known.values())
                if share >= limit:
                    add(risk, lineage, f"fraction_from_{field}", share, limit, f"{top}: {k}/{sum(known.values())}")
        years = [int(r["collection_year"]) for r in group if valid_year(r["collection_year"])]
        if cutoff is not None and years:
            late = sum(y >= cutoff for y in years) / len(years)
            if late >= f["lineage_after_cutoff_fraction"]:
                add("lineage_year_concentration", lineage, "fraction_at_or_after_cutoff", late,
                    f["lineage_after_cutoff_fraction"], f"cutoff {cutoff}")
            if 1 - late >= f["lineage_after_cutoff_fraction"] and len(years) >= f["min_lineage_isolates"]:
                add("lineage_absent_from_recent_years", lineage, "fraction_before_cutoff", 1 - late,
                    f["lineage_after_cutoff_fraction"], f"cutoff {cutoff}")
    by_drug = defaultdict(list)
    for e in binary:
        if e["phenotype_endpoint"] == "R":
            by_drug[e["antimicrobial"]].append(e)
    for drug, resistant in sorted(by_drug.items()):
        if len(resistant) < f["min_resistant_for_concentration"]:
            continue
        for risk, field, limit in (("resistance_lineage_concentration", "lineage_id", f["resistance_single_lineage_fraction"]),
                                   ("resistance_snp_cluster_concentration", "ncbi_snp_cluster", f["resistance_single_cluster_fraction"]),
                                   ("resistance_project_concentration", "project", f["resistance_single_project_fraction"])):
            known = Counter(e[field] for e in resistant if e[field] not in {NA, "UNKNOWN"})
            if known:
                top, k = known.most_common(1)[0]
                if k / len(resistant) >= limit:
                    add(risk, drug, f"fraction_of_R_from_{field}", k / len(resistant), limit, f"{top}: {k}/{len(resistant)}")
    clusters = Counter(r["ncbi_snp_cluster"] for r in population if r["ncbi_snp_cluster"] != NA)
    for cluster, n in sorted(clusters.items()):
        if n >= f["large_snp_cluster_isolates"]:
            add("large_snp_cluster", cluster, "population_isolates", n, f["large_snp_cluster_isolates"], "")
    genomes = [r for r in rows if r["genome_qc_status"] != "fail" and r["genome_sha256"] != NA]
    if genomes:
        success = sum(r["mlst_status"] == "assigned" for r in genomes) / len(genomes)
        if success < f["low_mlst_success_fraction"]:
            add("low_mlst_success", "cohort", "assigned_fraction", success, f["low_mlst_success_fraction"], "")
    return flags


def lineage_components(rows: list[dict]) -> dict:
    """Group population isolates that share an assigned ST or an NCBI SNP cluster.

    A cluster spanning two STs joins them. UNKNOWN isolates join through their
    cluster; one without a cluster has no component and is never test data.
    """
    parent = {}

    def find(key):
        parent.setdefault(key, key)
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    keys = {}
    for r in rows:
        if r["eligible_population_analysis"] != "true":
            continue
        own = ([f"ST:{r['lineage_id']}"] if r["lineage_id"] != "UNKNOWN" else []) + \
              ([f"C:{r['ncbi_snp_cluster']}"] if r["ncbi_snp_cluster"] != NA else [])
        for key in own:
            find(key)
        if len(own) == 2:
            parent[find(own[0])] = find(own[1])
        keys[r["isolate_id"]] = own
    members = defaultdict(set)
    for key in parent:
        members[find(key)].add(key)
    names = {}
    for root, group in members.items():
        sts = sorted(k[3:] for k in group if k.startswith("ST:"))
        names[root] = "+".join(sts) if sts else sorted(group)[0][2:]
    return {iid: names[find(own[0])] if own else NA for iid, own in keys.items()}


def evaluation_plan(binary: list[dict], config: dict, cutoff: int | None) -> list[dict]:
    """Per-antimicrobial tier and analyses, following the evaluation design decisions."""
    design = config["evaluation_design"]
    tiers, reference = design["drug_tiers"], design["geography"]["reference_country"]
    by_drug = defaultdict(list)
    for e in binary:
        by_drug[e["antimicrobial"]].append(e)
    out = []
    for drug, eps in sorted(by_drug.items()):
        def per(field, skip=()):
            g = defaultdict(Counter)
            for e in eps:
                if e[field] not in skip:
                    g[e[field]][e["phenotype_endpoint"]] += 1
            return g
        n = Counter(e["phenotype_endpoint"] for e in eps)
        resistant_projects = Counter(e["project"] for e in eps if e["phenotype_endpoint"] == "R")
        top_project, top_r = resistant_projects.most_common(1)[0] if resistant_projects else (NA, 0)
        projects, sts = per("project"), per("lineage_id", {"UNKNOWN"})
        metrics = {"R": n["R"], "S": n["S"], "top_project_R_fraction": top_r / n["R"] if n["R"] else 1.0,
                   "projects_with_both": sum(c["R"] >= tiers["project_both_min"] and c["S"] >= tiers["project_both_min"]
                                             for c in projects.values()),
                   "sts_with_both": sum(c["R"] >= tiers["st_both_min"] and c["S"] >= tiers["st_both_min"]
                                        for c in sts.values())}
        dated = [e for e in eps if e["eligible_temporal_holdout"] == "true"]
        recent = Counter(e["phenotype_endpoint"] for e in dated if cutoff is not None and int(e["collection_year"]) >= cutoff)
        earlier = Counter(e["phenotype_endpoint"] for e in dated if cutoff is not None and int(e["collection_year"]) < cutoff)
        countries = per("country", COUNTRY_MISSING)
        other = Counter()
        for country, c in countries.items():
            if country != reference:
                other.update(c)
        local = {k: [c["R"], c["S"]] for k, c in sorted(countries.items()) if k != reference
                 and c["R"] >= design["geography"]["min_per_country"] and c["S"] >= design["geography"]["min_per_country"]}
        status = classify(metrics, tiers["primary"], tiers["secondary"])
        if drug in tiers["excluded"]:
            tier, reason = "excluded", tiers["excluded"][drug]
        elif status == "feasible":
            tier, reason = "primary", "meets primary criteria"
        else:
            tier = "secondary" if status == "marginal" else "excluded"
            rule = tiers["primary"] if status == "marginal" else tiers["secondary"]
            failed = []
            for key, limit in rule.items():
                name = key.removeprefix("min_").removeprefix("max_")
                value = metrics[name]
                if key.startswith("min_") and value < limit or key.startswith("max_") and value > limit:
                    failed.append(f"{name} {round(value, 2)} {'<' if key.startswith('min_') else '>'} {limit}")
            reason = "; ".join(failed)
        t_min, g_min = design["temporal_min_per_period"], design["geography"]["min_pooled"]
        temporal = (tier == "primary" and min(recent["R"], earlier["R"]) >= t_min
                    and min(recent["S"], earlier["S"]) >= t_min)
        geographic = tier == "primary" and other["R"] >= g_min and other["S"] >= g_min
        yes = lambda flag: "yes" if flag else "no"  # noqa: E731
        out.append({"antimicrobial": drug, "tier": tier, "tier_reason": reason, "n_R": str(n["R"]), "n_S": str(n["S"]),
                    "top_project": top_project, "top_project_R_fraction": ratio(top_r, n["R"]),
                    "projects_with_both": str(metrics["projects_with_both"]),
                    "sts_with_both": str(metrics["sts_with_both"]),
                    "earlier_R": str(earlier["R"]), "earlier_S": str(earlier["S"]),
                    "recent_R": str(recent["R"]), "recent_S": str(recent["S"]),
                    "non_reference_R": str(other["R"]), "non_reference_S": str(other["S"]),
                    "non_reference_countries_json": canonical(local),
                    "lineage_cv": yes(tier in ("primary", "secondary")), "within_st_cluster": yes(tier == "primary"),
                    "project_holdout": yes(tier == "primary"), "geographic_transfer": yes(geographic),
                    "temporal": yes(temporal), "source_holdout": "no", "decision_record": design["decision_record"]})
    return out


POP_FIELDS = ("cohort_id isolate_id pdt_accession biosample_accession assembly_accession organism strain genome_sha256 "
              "genome_qc_status genome_qc_flags assembly_length contig_count n50 gc_fraction ambiguous_fraction "
              "mlst_scheme mlst_st mlst_status mlst_profile phylogroup phylogroup_status ncbi_snp_cluster "
              "ncbi_snp_cluster_status ncbi_snp_cluster_versions_json lineage_id country collection_year "
              "source_category submitter project binary_rs_endpoint_count antimicrobial_count pr5_eligible_isolate "
              "eligible_population_analysis eligible_lineage_holdout eligible_country_holdout "
              "eligible_temporal_holdout eligible_source_holdout").split()
QC_FIELDS = ("isolate_id assembly_accession genome_sha256 genome_qc_status genome_qc_flags assembly_length "
             "contig_count n50 longest_contig gc_fraction ambiguous_fraction").split()
MLST_FIELDS = ("isolate_id assembly_accession mlst_scheme mlst_st mlst_status mlst_profile mlst_raw_st software "
               "software_version database_version scheme_profile_sha256 input_sha256").split()
PHYLO_FIELDS = "isolate_id assembly_accession phylogroup phylogroup_status phylogroup_raw method software_version input_sha256".split()
GROUP_FIELDS = ("isolate_id lineage_group snp_cluster_group country_group year_group source_group "
                "population_eligible lineage_holdout_eligible country_holdout_eligible temporal_holdout_eligible "
                "source_holdout_eligible lineage_component lineage_cv_test_eligible split_unit project_group "
                "geography_group period").split()
PLAN_FIELDS = ("antimicrobial tier tier_reason n_R n_S top_project top_project_R_fraction projects_with_both "
               "sts_with_both earlier_R earlier_S recent_R recent_S non_reference_R non_reference_S "
               "non_reference_countries_json lineage_cv within_st_cluster project_holdout geographic_transfer "
               "temporal source_holdout decision_record").split()
LINEAGE_FIELDS = ("lineage_id is_lineage isolate_count fraction_of_cohort population_eligible_isolates country_count "
                  "largest_country fraction_from_largest_country earliest_year latest_year year_span "
                  "source_category_count antimicrobial_endpoint_count phylogroups_json").split()
FEAS_FIELDS = ("antimicrobial holdout eligible_isolates groups largest_group_size largest_group_fraction "
               "groups_with_both_R_S n_R n_S detail_json status decision").split()


def build(processed: Path, cache: Path, outdir: Path, config_path: Path = DEFAULT_CONFIG,
          expected_manifest_sha256: str | None = None, workers: int = 4) -> dict:
    config = load_config(config_path)
    cohort = load_cohort(processed, expected_manifest_sha256)
    if not cohort["candidates"]:
        raise ValueError("PR #5 cohort contains no isolates")
    if cohort["duplicate_assemblies"]:
        raise ValueError(f"assembly shared by several isolates: {cohort['duplicate_assemblies'][:5]}")
    rows, qc_rows, mlst_rows, phylo_rows, tools = population_rows(cohort, cache, config, workers)
    by_id = {r["isolate_id"]: r for r in rows}
    binary = []
    for e in cohort["endpoints"]:
        r = by_id.get(e["isolate_id"])
        if r and e["eligible_overall"] == "true" and r["eligible_population_analysis"] == "true":
            if e["phenotype_endpoint"] not in {"R", "S"}:
                raise ValueError("eligible PR #5 endpoint is not binary R/S")
            binary.append({**{k: r[k] for k in ("lineage_id", "ncbi_snp_cluster", "project", "country", "collection_year",
                                                 "source_category", "eligible_lineage_holdout", "eligible_country_holdout",
                                                 "eligible_temporal_holdout", "eligible_source_holdout",
                                                 "eligible_population_analysis")},
                           "antimicrobial": e["antimicrobial"], "phenotype_endpoint": e["phenotype_endpoint"]})
    population = [r for r in rows if r["eligible_population_analysis"] == "true"]
    years = [int(r["collection_year"]) for r in population if valid_year(r["collection_year"])]
    cutoff = temporal_cutoff(years, config["feasibility_rules"]["temporal"])
    lineages = lineage_summary(rows, binary, config)
    assigned = [r for r in rows if r["mlst_status"] == "assigned"]
    st_counts = Counter(r["lineage_id"] for r in assigned)
    genomes = [r for r in rows if r["genome_sha256"] != NA]
    summary = {
        "cohort_isolates": len(rows), "unique_assemblies": len({r["assembly_accession"] for r in rows}),
        "genomes_available": len(genomes),
        "genome_qc": dict(sorted(Counter(r["genome_qc_status"] for r in rows).items())),
        "mlst_status": dict(sorted(Counter(r["mlst_status"] for r in rows).items())),
        "mlst_assigned_fraction_of_genomes": round(len(assigned) / len(genomes), 6) if genomes else None,
        "unique_sts": len(st_counts), "singleton_sts": sum(n == 1 for n in st_counts.values()),
        "unresolved_st_isolates": len(rows) - len(assigned),
        "top_sts": st_counts.most_common(config["lineage"]["top_n"]),
        "phylogroup": dict(sorted(Counter(r["phylogroup"] if r["phylogroup"] != NA else r["phylogroup_status"]
                                          for r in rows).items())),
        "snp_cluster_status": dict(sorted(Counter(r["ncbi_snp_cluster_status"] for r in rows).items())),
        "unique_snp_clusters": len({r["ncbi_snp_cluster"] for r in rows if r["ncbi_snp_cluster"] != NA}),
        "largest_snp_clusters": Counter(r["ncbi_snp_cluster"] for r in rows if r["ncbi_snp_cluster"] != NA).most_common(10),
        "eligibility": {k: sum(r[k] == "true" for r in rows) for k in
                        ("pr5_eligible_isolate", "eligible_population_analysis", "eligible_lineage_holdout",
                         "eligible_country_holdout", "eligible_temporal_holdout", "eligible_source_holdout")},
        "binary_endpoints_population": len(binary), "temporal_recent_cutoff_year": cutoff,
        "reconciliation": {
            "rows_equal_pr5_candidates": len(rows) == int(cohort["manifest"]["candidate_isolate_count"]),
            "population_endpoints_le_pr5_eligible": len(binary) <= int(cohort["manifest"]["eligible_binary_rs_endpoint_count"]),
            "lineage_counts_sum_to_cohort": sum(int(x["isolate_count"]) for x in lineages) == len(rows),
            "qc_rows_equal_cohort": len(qc_rows) == len(rows),
        },
    }
    if not all(summary["reconciliation"].values()):
        raise ValueError(f"population counts do not reconcile: {summary['reconciliation']}")
    flags = review_flags(rows, binary, config, cutoff)
    feasibility = holdout_feasibility(binary, rows, config)
    components = lineage_components(rows)
    reference = config["evaluation_design"]["geography"]["reference_country"]
    groups = [{"isolate_id": r["isolate_id"], "lineage_group": r["lineage_id"], "snp_cluster_group": r["ncbi_snp_cluster"],
               "country_group": r["country"] if r["country"] not in COUNTRY_MISSING else NA,
               "year_group": r["collection_year"] if valid_year(r["collection_year"]) else NA,
               "source_group": r["source_category"], "population_eligible": r["eligible_population_analysis"],
               "lineage_holdout_eligible": r["eligible_lineage_holdout"],
               "country_holdout_eligible": r["eligible_country_holdout"],
               "temporal_holdout_eligible": r["eligible_temporal_holdout"],
               "source_holdout_eligible": r["eligible_source_holdout"],
               "lineage_component": components.get(r["isolate_id"], NA),
               "lineage_cv_test_eligible": "true" if components.get(r["isolate_id"], NA) != NA else "false",
               "split_unit": r["ncbi_snp_cluster"] if r["ncbi_snp_cluster"] != NA else r["isolate_id"],
               "project_group": r["project"] if not missing(r["project"]) else NA,
               "geography_group": NA if r["country"] in COUNTRY_MISSING else
               "reference" if r["country"] == reference else "non_reference",
               "period": NA if not valid_year(r["collection_year"]) or cutoff is None else
               "recent" if int(r["collection_year"]) >= cutoff else "earlier"} for r in rows]
    plan = evaluation_plan(binary, config, cutoff)
    population_components = Counter(v for v in components.values() if v != NA)
    summary["evaluation_design"] = {
        "temporal_cutoff_year": cutoff,
        "lineage_components": len(population_components),
        "lineage_components_joining_several_sts": sorted(k for k in population_components if "+" in k),
        "unknown_lineage_in_component": sum(by_id[i]["lineage_id"] == "UNKNOWN" and v != NA
                                            for i, v in components.items()),
        "population_without_component": sum(v == NA for v in components.values()),
        "drug_tiers": dict(sorted(Counter(p["tier"] for p in plan).items())),
        "drugs_per_analysis": {k: sum(p[k] == "yes" for p in plan) for k in
                               ("lineage_cv", "within_st_cluster", "project_holdout", "geographic_transfer",
                                "temporal", "source_holdout")}}
    code = "unknown"
    try:
        code = subprocess.run(["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True,
                              cwd=Path(__file__).resolve().parent).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    retrieved = sorted(json.loads((Path(cache) / r["assembly_accession"] / "provenance.json").read_text())["retrieved_at"]
                       for r in genomes)
    manifest = {
        "population_id": config["population_id"], "builder_version": VERSION, "config_sha256": sha256(Path(config_path)),
        "cohort_id": cohort["manifest"]["cohort_id"], "pr5_cohort_manifest_sha256": cohort["manifest_sha256"],
        "pr5_files_sha256": cohort["file_sha256"], "pr5_source_snapshot": cohort["manifest"].get("source_snapshot_id"),
        "pr5_source_snapshot_manifest_sha256": cohort["manifest"].get("source_snapshot_manifest_sha256"),
        "pr5_code_commit": cohort["manifest"].get("code_commit"), "code_commit": code,
        "builder_sha256": sha256(Path(__file__)),
        "evaluation_design_record": config["evaluation_design"]["decision_record"],
        "isolate_count": len(rows), "assembly_count": summary["unique_assemblies"],
        "genome_retrieval": {"source": "NCBI Datasets v2 genome package", "verified_by": "genome_dataset.package_identity",
                             "genomes": len(genomes), "retrieved_first": retrieved[0] if retrieved else None,
                             "retrieved_last": retrieved[-1] if retrieved else None,
                             "genome_sha256_list_sha256": hashlib.sha256(
                                 "\n".join(f"{r['assembly_accession']}\t{r['genome_sha256']}" for r in genomes).encode()).hexdigest()},
        "tools": {k: [json.loads(x) for x in v] for k, v in tools.items()},
        "generated_at": retrieved[-1] if retrieved else None,
        "generated_at_semantics": "latest genome retrieval time; fixed so rebuilds are deterministic",
        "annotation_counts": {"mlst": summary["mlst_status"], "phylogroup": dict(sorted(Counter(r["phylogroup_status"] for r in rows).items())),
                              "genome_qc": summary["genome_qc"]},
    }
    outdir = Path(outdir)
    if outdir.exists():
        raise ValueError(f"output already exists: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".population-", dir=outdir.parent))
    try:
        write_tsv(staging / "population_metadata.tsv", POP_FIELDS, rows)
        write_tsv(staging / "genome_qc.tsv", QC_FIELDS, qc_rows)
        write_tsv(staging / "mlst.tsv", MLST_FIELDS, mlst_rows)
        write_tsv(staging / "phylogroup.tsv", PHYLO_FIELDS, phylo_rows)
        write_tsv(staging / "evaluation_groups.tsv", GROUP_FIELDS, groups)
        write_tsv(staging / "lineage_summary.tsv", LINEAGE_FIELDS, lineages)
        write_tsv(staging / "lineage_antimicrobial_summary.tsv",
                  "antimicrobial lineage_id n_total n_R n_S resistance_fraction".split(), lineage_antimicrobial(binary))
        for name, field in (("lineage_country_summary.tsv", "country"), ("lineage_year_summary.tsv", "collection_year"),
                            ("lineage_source_summary.tsv", "source_category")):
            write_tsv(staging / name, ["lineage_id", field, "isolate_count"], grouped_counts(rows, "lineage_id", field))
        write_tsv(staging / "snp_cluster_summary.tsv",
                  ("ncbi_snp_cluster isolate_count countries_json earliest_year latest_year lineages_json projects_json "
                   "binary_endpoints_R binary_endpoints_S").split(), snp_cluster_summary(rows, binary))
        write_tsv(staging / "project_summary.tsv",
                  ("project submitter isolate_count country_count countries_json lineage_count top_lineages_json "
                   "binary_endpoints_R binary_endpoints_S").split(), project_summary(rows, binary))
        write_tsv(staging / "holdout_feasibility.tsv", FEAS_FIELDS, feasibility)
        write_tsv(staging / "evaluation_plan.tsv", PLAN_FIELDS, plan)
        write_tsv(staging / "confounding_summary.tsv", "risk subject metric value threshold detail".split(), flags)
        (staging / "population_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (staging / "population_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        sums = [f"{sha256(p)}  {p.name}" for p in sorted(staging.iterdir()) if p.is_file()]
        (staging / "SHA256SUMS").write_text("\n".join(sums) + "\n", encoding="utf-8")
        staging.rename(outdir)
    finally:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("fetch-genomes", "annotate", "build"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--cohort", type=Path, required=True, help="PR #5 ecoli_cohort.py build directory")
        cmd.add_argument("--cache", type=Path, required=True, help="genome cache directory")
        cmd.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
        cmd.add_argument("--expected-cohort-manifest-sha256")
        if name == "fetch-genomes":
            cmd.add_argument("--workers", type=int, default=3)
            cmd.add_argument("--limit", type=int)
        elif name == "annotate":
            cmd.add_argument("--tool", choices=["mlst", "phylogroup"], required=True)
            cmd.add_argument("--env-bin", type=Path, required=True, help="bin directory of the pinned tool environment")
            cmd.add_argument("--workers", type=int, default=4)
        else:
            cmd.add_argument("--outdir", type=Path, required=True)
            cmd.add_argument("--workers", type=int, default=4, help="threads for genome QC")
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            result = build(args.cohort, args.cache, args.outdir, args.config, args.expected_cohort_manifest_sha256,
                           args.workers)
        else:
            cohort = load_cohort(args.cohort, args.expected_cohort_manifest_sha256)
            if args.command == "fetch-genomes":
                result = fetch_genomes(cohort, args.cache, args.workers, args.limit)
            else:
                result = annotate(args.tool, cohort, args.cache, args.env_bin, load_config(args.config), args.workers)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ValueError, KeyError, OSError, csv.Error, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
