#!/usr/bin/env python3
"""Leakage-aware generalisation benchmark for the frozen E. coli research cohort.

Commands, run in this order:
  sketch  verify every genome against the PR #8 manifest and write a
          FracMinHash k-mer sketch per verified genome (resumable)
  matrix  build the isolate x k-mer presence matrix from the sketches
  freeze  build benchmark endpoints, cohorts and split tasks, validate them
          and write the benchmark manifest; no model is fitted
  run     fit the model for every task and write test predictions
          (resumable; one checkpoint per task)
  report  metrics, generalisation drops, leakage and shift audits,
          determinant comparison and figure data from the frozen predictions

The split rules, cohorts, model and thresholds all come from the versioned
benchmark config. Phenotypes are the PR #7 endpoint labels (with the
decision 0004 correction); the PR #7 determinant calls are only read as a
comparator.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
from scipy import sparse
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import population_structure as pop  # noqa: E402

VERSION = "1.0.0"
NA = "NA"
DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "benchmark/cohorts/ecoli_research_v1/generalisation_v1.json"
GCA = re.compile(r"^GCA_\d{9}\.\d+$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
MANIFEST_FIELDS = ("isolate_id", "assembly_accession", "biosample_accession", "genome_sha256", "sequence_records",
                   "sequence_bases")
GENOME_FAILURES = ("genome_missing", "genome_hash_mismatch", "genome_unreadable", "genome_version_mismatch")
CV_STRATEGIES = ("random", "random_cluster", "lineage")
GROUP_STRATEGIES = ("country", "project")
STRATEGIES = CV_STRATEGIES + GROUP_STRATEGIES + ("geography", "temporal")
COHORTS = ("common", "maximal")


class LeakageError(ValueError):
    """A split puts held-out information into training."""


class SingleClassError(ValueError):
    """Training data hold only one phenotype."""


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def text_sha(value) -> str:
    return sha256_bytes(pop.canonical(value).encode())


def fmt(value, digits: int = 6) -> str:
    if value is None or isinstance(value, float) and math.isnan(value):
        return NA
    if isinstance(value, float):
        return f"{round(value, digits):.{digits}f}".rstrip("0").rstrip(".") if value != int(value) else str(int(value))
    return str(value)


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with opener(tmp, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, NA)) for k in fields})
    tmp.replace(path)


# ---------------------------------------------------------------- config and upstream

def load_config(path: Path) -> dict:
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("benchmark_version", "upstream", "endpoints", "features", "model", "seeds", "splits", "cohorts",
                "support", "bootstrap", "decision_record"):
        if key not in config:
            raise ValueError(f"benchmark config lacks {key}")
    model = config["model"]
    if model["penalty"] != "l1" and config["features"].get("deduplicate_identical_columns"):
        raise ValueError("merging identical columns is only exact for an L1 penalty")
    if not model["C_grid"] or any(c <= 0 for c in model["C_grid"]):
        raise ValueError("C_grid must hold positive values")
    if not 0 < model["threshold"] < 1:
        raise ValueError("decision threshold must lie in (0, 1)")
    for key in ("split", "inner_cv", "model", "bootstrap"):
        if not isinstance(config["seeds"].get(key), int):
            raise ValueError(f"seed {key} must be a fixed integer")
    if config["splits"]["folds"] < 2:
        raise ValueError("at least two folds are needed")
    return config


def config_sha(path: Path) -> str:
    return pop.sha256(Path(path))


def check_sums(directory: Path, names: list[str]) -> None:
    """Check named files against the directory's own SHA256SUMS."""
    directory = Path(directory)
    sums = {}
    for line in (directory / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        sums[name.lstrip("*")] = digest
    for name in names:
        if name not in sums:
            raise ValueError(f"{directory}: {name} not listed in SHA256SUMS")
        if pop.sha256(directory / name) != sums[name]:
            raise ValueError(f"{directory}: {name} differs from SHA256SUMS")


def read_genome_manifest(path: Path, expected_sha256: str | None = None, expected_rows: int | None = None) -> dict:
    """PR #8 genome manifest keyed by isolate, after schema and uniqueness checks."""
    path = Path(path)
    if expected_sha256 and pop.sha256(path) != expected_sha256:
        raise ValueError("genome manifest SHA-256 differs from the benchmark config")
    rows = pop.read_tsv(path)
    if rows and any(field not in rows[0] for field in MANIFEST_FIELDS):
        raise ValueError(f"genome manifest lacks {[f for f in MANIFEST_FIELDS if f not in rows[0]]}")
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError(f"genome manifest has {len(rows)} rows, expected {expected_rows}")
    out, assemblies = {}, {}
    for row in rows:
        iid, gca, digest = row["isolate_id"], row["assembly_accession"], row["genome_sha256"]
        if pop.missing(iid) or not GCA.match(gca) or not HEX64.match(digest):
            raise ValueError(f"genome manifest row for {iid!r} has a missing ID, unversioned GCA or bad SHA-256")
        if iid in out:
            raise ValueError(f"genome manifest lists isolate {iid} twice")
        if gca in assemblies:
            raise ValueError(f"assembly {gca} is listed for {assemblies[gca]} and {iid}")
        out[iid], assemblies[gca] = row, iid
    return out


def duplicate_genomes(manifest: dict) -> list[list[str]]:
    """Isolates whose genome bytes are identical (same SHA-256 under different assemblies)."""
    by = defaultdict(list)
    for iid, row in manifest.items():
        by[row["genome_sha256"]].append(iid)
    return sorted(sorted(v) for v in by.values() if len(v) > 1)


def upstream_hashes(config: dict, population: Path, baseline: Path) -> dict:
    """Check the PR #6 and PR #7 builds the benchmark reads and return their hashes."""
    up = config["upstream"]
    check_sums(population, ["population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json"])
    check_sums(baseline, ["genotype_phenotype_endpoints.tsv", "baseline_manifest.json"])
    pop_sha = pop.sha256(Path(population) / "population_manifest.json")
    base_sha = pop.sha256(Path(baseline) / "baseline_manifest.json")
    if pop_sha != up["pr6_population_manifest_sha256"]:
        raise ValueError("PR #6 population manifest differs from the benchmark config")
    if base_sha != up["pr7_baseline_manifest_sha256"]:
        raise ValueError("PR #7 baseline manifest differs from the benchmark config")
    base = json.loads((Path(baseline) / "baseline_manifest.json").read_text(encoding="utf-8"))
    popm = json.loads((Path(population) / "population_manifest.json").read_text(encoding="utf-8"))
    if base["pr5_cohort_manifest_sha256"] != up["pr5_cohort_manifest_sha256"] or \
            popm["pr5_cohort_manifest_sha256"] != up["pr5_cohort_manifest_sha256"]:
        raise ValueError("PR #5 cohort manifest hash differs between the config and upstream builds")
    if base.get("phenotype_corrections_sha256") != up["pr7_phenotype_corrections_sha256"]:
        raise ValueError("PR #7 build used a different phenotype correction table")
    if base.get("mapping_version") != up["pr7_mapping_version"]:
        raise ValueError("PR #7 build used a different determinant mapping")
    return {"pr5_cohort_manifest_sha256": up["pr5_cohort_manifest_sha256"],
            "pr6_population_manifest_sha256": pop_sha, "pr7_baseline_manifest_sha256": base_sha,
            "pr7_endpoints_sha256": pop.sha256(Path(baseline) / "genotype_phenotype_endpoints.tsv"),
            "pr6_evaluation_groups_sha256": pop.sha256(Path(population) / "evaluation_groups.tsv"),
            "pr7_mapping_version": base["mapping_version"],
            "pr7_phenotype_corrections_sha256": base["phenotype_corrections_sha256"]}


# ---------------------------------------------------------------- genomes and sketches

def parse_fasta(data: bytes) -> list[str]:
    records, current = [], []
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            if current:
                records.append("".join(current).upper())
            current = []
        else:
            current.append(line.strip())
    if current:
        records.append("".join(current).upper())
    if not records:
        raise ValueError("no sequence records")
    return records


def verify_genome(cache: Path, row: dict, expected_assembly: str | None = None) -> tuple[str, str, bytes | None]:
    """Classify one cached genome against the PR #8 manifest row; returns its bytes when verified."""
    gca = row["assembly_accession"]
    if expected_assembly is not None and expected_assembly != gca:
        return "genome_version_mismatch", f"cohort assembly {expected_assembly}, manifest {gca}", None
    path = Path(cache) / gca / "genome.fna"
    if not path.is_file():
        return "genome_missing", str(path), None
    try:
        data = path.read_bytes()
    except OSError as error:
        return "genome_unreadable", str(error), None
    digest = sha256_bytes(data)
    if digest != row["genome_sha256"]:
        return "genome_hash_mismatch", digest, None
    return "verified", "", data


def sketch_params(config: dict) -> dict:
    f = config["features"]
    return {"method": f["method"], "tool": f["tool"], "tool_version": f["tool_version"], "ksize": f["ksize"],
            "scaled": f["scaled"], "seed": f["hash_seed"]}


def sketch_genome(data: bytes, params: dict) -> np.ndarray:
    import sourmash
    if sourmash.VERSION != params["tool_version"]:
        raise RuntimeError(f"sourmash {sourmash.VERSION} installed, config needs {params['tool_version']}")
    mh = sourmash.MinHash(n=0, ksize=params["ksize"], scaled=params["scaled"], seed=params["seed"])
    for record in parse_fasta(data):
        mh.add_sequence(record, force=True)
    return np.array(sorted(mh.hashes), dtype=np.uint64)


def sketch_one(cache: str, outdir: str, row: dict, expected_assembly: str | None, params: dict) -> dict:
    """Verify and sketch one genome; a valid earlier sketch is reused."""
    gca, iid = row["assembly_accession"], row["isolate_id"]
    meta_path, npy = Path(outdir) / f"{gca}.json", Path(outdir) / f"{gca}.npy"
    base = {"isolate_id": iid, "assembly_accession": gca, "genome_sha256": row["genome_sha256"]}
    if meta_path.exists() and npy.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("params") == params and meta.get("genome_sha256") == row["genome_sha256"] and \
                meta.get("hash_count") == int(np.load(npy, mmap_mode="r").shape[0]) and \
                (expected_assembly is None or expected_assembly == gca):
            return {**base, "status": "verified", "detail": "", "hash_count": meta["hash_count"], "action": "reused"}
    status, detail, data = verify_genome(Path(cache), row, expected_assembly)
    if status != "verified":
        return {**base, "status": status, "detail": detail, "hash_count": 0, "action": "skipped"}
    try:
        hashes = sketch_genome(data, params)
    except (ValueError, UnicodeDecodeError, OSError) as error:
        return {**base, "status": "genome_unreadable", "detail": str(error), "hash_count": 0, "action": "skipped"}
    tmp = npy.with_name(f".{npy.name}.{os.getpid()}.tmp.npy")
    np.save(tmp, hashes)
    tmp.replace(npy)
    pop.write_json_atomic(meta_path, {"isolate_id": iid, "assembly_accession": gca,
                                      "genome_sha256": row["genome_sha256"], "params": params,
                                      "hash_count": int(hashes.shape[0]), "created_at": now()})
    return {**base, "status": "verified", "detail": "", "hash_count": int(hashes.shape[0]), "action": "sketched"}


VERIFY_FIELDS = "isolate_id assembly_accession genome_sha256 status detail hash_count".split()


def sketch_all(config: dict, cache: Path, outdir: Path, cohort_assemblies: dict, workers: int = 4,
               limit: int | None = None) -> dict:
    """Verify and sketch all manifest genomes; writes genome_verification.tsv."""
    up = config["upstream"]
    root = Path(__file__).resolve().parents[1]
    manifest = read_genome_manifest(root / up["genome_manifest"], up["genome_manifest_sha256"],
                                    up["genome_manifest_rows"])
    outdir = Path(outdir)
    sketches = outdir / "sketches"
    sketches.mkdir(parents=True, exist_ok=True)
    params = sketch_params(config)
    rows = sorted(manifest.values(), key=lambda r: r["isolate_id"])[:limit]
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(sketch_one, str(cache), str(sketches), r, cohort_assemblies.get(r["isolate_id"]),
                               params) for r in rows]
        for i, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if i % 500 == 0:
                print(f"{now()} sketched {i}/{len(rows)}", flush=True)
    results.sort(key=lambda r: r["isolate_id"])
    write_tsv(outdir / "genome_verification.tsv", VERIFY_FIELDS, results)
    counts = Counter(r["status"] for r in results)
    summary = {"genomes": len(results), "status_counts": dict(sorted(counts.items())),
               "actions": dict(sorted(Counter(r["action"] for r in results).items())),
               "duplicate_genomes": duplicate_genomes(manifest), "params": params,
               "genome_manifest_sha256": up["genome_manifest_sha256"], "finished_at": now()}
    pop.write_json_atomic(outdir / "genome_verification.json", summary)
    return summary


# ---------------------------------------------------------------- presence matrix

def build_matrix(sketch_dir: Path, isolates: list[tuple[str, str]], min_present: int, min_absent: int,
                 dedup: bool, buckets: int = 16, scaled: int | None = None) -> dict:
    """Isolate x feature presence matrix from per-genome sketches.

    `isolates` is a list of (isolate_id, assembly) in row order. Features
    present in fewer than `min_present` or absent from fewer than `min_absent`
    of all rows are dropped: the per-task filter, which counts training rows
    only, would drop them in every task, so this changes nothing. With
    `dedup`, columns with identical presence over all rows are merged into one
    pattern, which is exact for an L1-penalised linear model.
    """
    arrays = []
    for _iid, gca in isolates:
        h = np.load(Path(sketch_dir) / f"{gca}.npy")
        if scaled is not None:  # downsample a finer sketch
            h = h[h <= np.uint64(max_hash(scaled))]
        arrays.append(h)
    n = len(arrays)
    edges = [np.uint64(int(i * (2 ** 64) // buckets)) for i in range(buckets)]
    signatures, pattern_rows, hash_ids, hash_pattern = {}, [], [], []
    observed = 0
    for b in range(buckets):
        lo = edges[b]
        parts, rows = [], []
        for r, h in enumerate(arrays):
            start = np.searchsorted(h, lo, side="left")
            stop = np.searchsorted(h, edges[b + 1], side="left") if b + 1 < buckets else h.shape[0]
            if stop > start:
                parts.append(h[start:stop])
                rows.append(np.full(stop - start, r, dtype=np.int32))
        if not parts:
            continue
        hv, rv = np.concatenate(parts), np.concatenate(rows)
        order = np.argsort(hv, kind="stable")
        hv, rv = hv[order], rv[order]
        starts = np.flatnonzero(np.r_[True, hv[1:] != hv[:-1]])
        counts = np.diff(np.r_[starts, hv.shape[0]])
        observed += starts.shape[0]
        keep = (counts >= min_present) & (n - counts >= min_absent)
        for s, c in zip(starts[keep], counts[keep]):
            members = rv[s:s + c]
            key = hashlib.blake2b(members.tobytes(), digest_size=16).digest() if dedup else len(pattern_rows)
            pid = signatures.get(key)
            if pid is None:
                pid = signatures[key] = len(pattern_rows)
                pattern_rows.append(members)
            hash_ids.append(hv[s])
            hash_pattern.append(pid)
    p = len(pattern_rows)
    if p == 0:
        raise ValueError("no informative features")
    sizes = np.array([len(m) for m in pattern_rows], dtype=np.int64)
    indptr = np.zeros(p + 1, dtype=index_dtype(int(sizes.sum())))
    indptr[1:] = np.cumsum(sizes)
    indices = np.concatenate(pattern_rows)
    del pattern_rows, signatures
    csc = sparse.csc_matrix((np.ones(indices.shape[0], dtype=np.int8), indices, indptr), shape=(n, p))
    del indices
    csr = csc.tocsr()
    csr.sort_indices()
    hashes = np.array(hash_ids, dtype=np.uint64)
    return {"matrix": csr, "hashes": hashes, "hash_pattern": np.array(hash_pattern, dtype=np.int64),
            "observed_features": int(observed), "kept_features": int(hashes.shape[0]), "patterns": p,
            "nnz": int(csr.nnz)}


def max_hash(scaled: int) -> int:
    """FracMinHash threshold used by sourmash for a scaled value."""
    from sourmash.minhash import _get_max_hash_for_scaled
    return _get_max_hash_for_scaled(scaled)


def index_dtype(nnz: int):
    """int32 row pointers while they fit, so scipy does not copy the indices to int64."""
    return np.int32 if nnz < 2 ** 31 else np.int64


def save_matrix(outdir: Path, built: dict, isolates: list[tuple[str, str]], meta: dict) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    m = built["matrix"]
    np.save(outdir / "indptr.npy", m.indptr.astype(index_dtype(m.nnz)))
    np.save(outdir / "indices.npy", m.indices.astype(np.int32))
    np.save(outdir / "feature_hashes.npy", built["hashes"])
    np.save(outdir / "feature_pattern.npy", built["hash_pattern"])
    (outdir / "rows.tsv").write_text("isolate_id\tassembly_accession\n" +
                                     "".join(f"{i}\t{g}\n" for i, g in isolates), encoding="utf-8")
    manifest = {**meta, "rows": len(isolates), "patterns": built["patterns"], "nnz": built["nnz"],
                "observed_features": built["observed_features"], "kept_features": built["kept_features"],
                "rows_sha256": pop.sha256(outdir / "rows.tsv"),
                "feature_hashes_sha256": sha256_bytes(built["hashes"].tobytes()),
                "feature_pattern_sha256": sha256_bytes(built["hash_pattern"].tobytes()),
                "indices_sha256": sha256_bytes(m.indices.astype(np.int32).tobytes()),
                "indptr_sha256": sha256_bytes(m.indptr.astype(index_dtype(m.nnz)).tobytes()), "created_at": now()}
    pop.write_json_atomic(outdir / "feature_manifest.json", manifest)
    return manifest


def load_matrix(directory: Path, check: bool = False) -> tuple[sparse.csr_matrix, dict, dict]:
    directory = Path(directory)
    meta = json.loads((directory / "feature_manifest.json").read_text(encoding="utf-8"))
    indptr = np.load(directory / "indptr.npy", mmap_mode="r")
    indices = np.load(directory / "indices.npy", mmap_mode="r")
    if check and (sha256_bytes(np.asarray(indices).tobytes()) != meta["indices_sha256"]
                  or sha256_bytes(np.asarray(indptr).tobytes()) != meta["indptr_sha256"]):
        raise ValueError("feature matrix differs from its manifest")
    rows = [line.split("\t")[0] for line in (directory / "rows.tsv").read_text(encoding="utf-8").splitlines()[1:]]
    m = sparse.csr_matrix((np.ones(indices.shape[0], dtype=np.int8), indices, indptr),
                          shape=(meta["rows"], meta["patterns"]))
    return m, {iid: i for i, iid in enumerate(rows)}, meta


# ---------------------------------------------------------------- benchmark endpoints

def drug_plan(binary_rows: list[dict], population_config: dict) -> dict:
    """Decision 0001 tiers recomputed on the corrected labels (decision 0005)."""
    cutoff = population_config["feasibility_rules"]["temporal"]["cutoff_year"]
    return {p["antimicrobial"]: p for p in pop.evaluation_plan(binary_rows, population_config, cutoff)}


def load_endpoints(baseline: Path, population: Path, verified: set, population_config: dict) -> tuple[list[dict], dict]:
    """Binary R/S endpoints of population-eligible isolates, with split metadata.

    Returns (endpoints of study drugs with a verified genome, plan by drug).
    """
    groups = {r["isolate_id"]: r for r in pop.read_tsv(Path(population) / "evaluation_groups.tsv")}
    meta = {r["isolate_id"]: r for r in pop.read_tsv(Path(population) / "population_metadata.tsv")}
    rows = []
    for r in pop.read_tsv(Path(baseline) / "genotype_phenotype_endpoints.tsv"):
        if r["eligible_binary_rs"] != "true" or r["population_eligible"] != "true":
            continue
        if r["phenotype"] not in ("R", "S"):
            raise ValueError(f"{r['isolate_id']} {r['antimicrobial']}: binary endpoint with phenotype {r['phenotype']}")
        g, m = groups[r["isolate_id"]], meta[r["isolate_id"]]
        if m["assembly_accession"] != r["assembly_accession"]:
            raise ValueError(f"{r['isolate_id']}: PR #6 and PR #7 assemblies differ")
        year = m["collection_year"]
        rows.append({
            "isolate_id": r["isolate_id"], "assembly_accession": r["assembly_accession"],
            "antimicrobial": r["antimicrobial"], "phenotype": r["phenotype"], "phenotype_pr5": r["phenotype_pr5"],
            "phenotype_correction": r["phenotype_correction"],
            "determinant_status": r["determinant_status"] if r["endpoint_eligible"] == "true" else NA,
            "concordance_category": r["concordance_category"] if r["endpoint_eligible"] == "true" else NA,
            "lineage_id": m["lineage_id"], "lineage_component": g["lineage_component"],
            "lineage_test_eligible": g["lineage_cv_test_eligible"] == "true",
            "snp_cluster": m["ncbi_snp_cluster"], "split_unit": g["split_unit"],
            "country": NA if m["country"] in pop.COUNTRY_MISSING else m["country"],
            "year": int(year) if pop.valid_year(year) else None, "period": g["period"],
            "source": m["source_category"], "project": m["project"], "geography_group": g["geography_group"],
            # fields evaluation_plan reads
            "phenotype_endpoint": r["phenotype"], "eligible_temporal_holdout": m["eligible_temporal_holdout"],
            "collection_year": year})
    plan = drug_plan(rows, population_config)
    study = {d for d, p in plan.items() if p["tier"] in ("primary", "secondary")}
    for e in rows:
        e["tier"] = plan[e["antimicrobial"]]["tier"]
        e["genome_verified"] = e["isolate_id"] in verified
    return [e for e in rows if e["antimicrobial"] in study], plan


# ---------------------------------------------------------------- split engine

def order_key(seed: int, *parts) -> str:
    return sha256_bytes(":".join([str(seed), *map(str, parts)]).encode())


def assign_folds(units: dict, k: int, seed: int, salt: str) -> dict:
    """Deterministic stratified assignment of units to k folds.

    `units` maps unit -> (n_R, n_S). Units are taken largest first (ties in
    hash order) and each goes to the fold where its addition leaves the
    smallest maximum of R and S load relative to the totals.
    """
    tot_r = sum(r for r, _ in units.values()) or 1
    tot_s = sum(s for _, s in units.values()) or 1
    load = [[0, 0] for _ in range(k)]
    out = {}
    for unit in sorted(units, key=lambda u: (-sum(units[u]), order_key(seed, salt, u))):
        r, s = units[unit]
        best = min(range(k), key=lambda f: (max((load[f][0] + r) / tot_r, (load[f][1] + s) / tot_s),
                                            load[f][0] + load[f][1], f))
        load[best][0] += r
        load[best][1] += s
        out[unit] = best
    return out


def counts_by(endpoints: list[dict], key: str) -> dict:
    c = defaultdict(lambda: [0, 0])
    for e in endpoints:
        c[e[key]][0 if e["phenotype"] == "R" else 1] += 1
    return {u: tuple(v) for u, v in c.items()}


def rs(endpoints) -> tuple[int, int]:
    n_r = sum(e["phenotype"] == "R" for e in endpoints)
    return n_r, len(endpoints) - n_r


def cohort_members(endpoints: list[dict], cohort: str, strategy: str) -> tuple[list[dict], dict]:
    """Endpoints a strategy may use in a cohort, and exclusion reasons for the rest."""
    out, reasons = [], Counter()
    for e in endpoints:
        reason = None
        if not e["genome_verified"]:
            reason = "genome_not_verified"
        elif cohort == "common" and not e["lineage_test_eligible"]:
            reason = "no_lineage_component"
        elif (cohort == "common" or strategy in ("country", "geography")) and e["country"] == NA:
            reason = "missing_country"
        elif (cohort == "common" or strategy == "temporal") and e["year"] is None:
            reason = "missing_year"
        if reason:
            reasons[reason] += 1
        else:
            out.append(e)
    return out, dict(reasons)


def task_id(cohort: str, drug: str, strategy: str, fold: str) -> str:
    return "__".join([cohort, drug.replace(" ", "_").replace("/", "_"), strategy, fold])


def make_task(cohort, drug, strategy, fold, train, test, note=""):
    train = sorted(train, key=lambda e: e["isolate_id"])
    test = sorted(test, key=lambda e: e["isolate_id"])
    return {"task_id": task_id(cohort, drug, strategy, fold), "cohort": cohort, "antimicrobial": drug,
            "strategy": strategy, "fold": fold, "train": train, "test": test, "note": note}


def feasibility(train, test, rule: dict) -> str | None:
    tr, ts = rs(train), rs(test)
    if tr[0] < rule["min_train_R"] or tr[1] < rule["min_train_S"]:
        return f"training {tr[0]} R / {tr[1]} S below {rule['min_train_R']} / {rule['min_train_S']}"
    if ts[0] < rule["min_test_R"] or ts[1] < rule["min_test_S"]:
        return f"test {ts[0]} R / {ts[1]} S below {rule['min_test_R']} / {rule['min_test_S']}"
    return None


def strategy_tasks(endpoints: list[dict], drug: str, cohort: str, strategy: str, config: dict,
                   plan_row: dict) -> tuple[list[dict], list[dict]]:
    """Tasks for one drug, cohort and strategy, and NOT_FEASIBLE records."""
    split = config["splits"]
    k, seed = split["folds"], config["seeds"]["split"]
    rule = split["strategies"][strategy]
    tier = plan_row["tier"]
    pool, _ = cohort_members(endpoints, cohort, strategy)
    skipped = []

    def nf(fold, reason):
        skipped.append({"cohort": cohort, "antimicrobial": drug, "strategy": strategy, "fold": fold,
                        "status": "NOT_FEASIBLE", "reason": reason})

    if tier not in rule["tiers"]:
        nf("all", f"{tier} drugs are not evaluated with {strategy} (decision 0001)")
        return [], skipped
    if rule.get("plan_flag") and plan_row.get(rule["plan_flag"]) != "yes":
        nf("all", f"evaluation plan {rule['plan_flag']} = no (decision 0001 support rule)")
        return [], skipped
    tasks = []
    if strategy in CV_STRATEGIES:
        test_ok = [e for e in pool if strategy != "lineage" or e["lineage_test_eligible"]]
        unit_key = {"random": "isolate_id", "random_cluster": "split_unit", "lineage": "lineage_component"}[strategy]
        folds = assign_folds(counts_by(test_ok, unit_key), k, seed, f"{drug}:{strategy}")
        for f in range(k):
            test = [e for e in test_ok if folds[e[unit_key]] == f]
            held = {e[unit_key] for e in test}
            if strategy == "lineage":
                # isolates without a component are training-only
                train = [e for e in pool if not (e["lineage_test_eligible"] and e["lineage_component"] in held)]
            else:
                train = [e for e in pool if e[unit_key] not in held]
            tasks.append(make_task(cohort, drug, strategy, f"fold{f + 1}", train, test))
        return tasks, skipped
    if strategy in GROUP_STRATEGIES:
        key = strategy
        groups = counts_by(pool, key)
        for group in sorted(groups):
            if group == NA:
                continue
            r, s = groups[group]
            if r < rule["min_test_R"] or s < rule["min_test_S"]:
                continue
            test = [e for e in pool if e[key] == group]
            train = [e for e in pool if e[key] != group and e[key] != NA]
            reason = feasibility(train, test, rule)
            if reason:
                nf(group, reason)
            else:
                tasks.append(make_task(cohort, drug, strategy, group, train, test))
        if not tasks and not skipped:
            nf("all", f"no {strategy} with at least {rule['min_test_R']} R and {rule['min_test_S']} S")
        return tasks, skipped
    if strategy == "geography":
        ref = rule["reference_country"]
        train = [e for e in pool if e["country"] == ref]
        test = [e for e in pool if e["country"] not in (ref, NA)]
        fold = f"{ref}_to_other"
    else:
        cutoff = rule["cutoff_year"]
        train = [e for e in pool if e["year"] < cutoff]
        test = [e for e in pool if e["year"] >= cutoff]
        fold = f"before_{cutoff}_to_{cutoff}_on"
    reason = feasibility(train, test, rule)
    if reason:
        nf(fold, reason)
    else:
        tasks.append(make_task(cohort, drug, strategy, fold, train, test))
    return tasks, skipped


def build_tasks(endpoints: list[dict], plan: dict, config: dict) -> tuple[list[dict], list[dict]]:
    by_drug = defaultdict(list)
    for e in endpoints:
        by_drug[e["antimicrobial"]].append(e)
    tasks, skipped = [], []
    for drug in sorted(by_drug):
        for cohort in COHORTS:
            for strategy in STRATEGIES:
                if cohort not in config["splits"]["strategies"][strategy]["cohorts"]:
                    continue
                t, s = strategy_tasks(by_drug[drug], drug, cohort, strategy, config, plan[drug])
                tasks += t
                skipped += s
    return tasks, skipped


def task_digest(task: dict) -> str:
    return text_sha({"train": [e["isolate_id"] for e in task["train"]],
                     "test": [e["isolate_id"] for e in task["test"]]})


# ---------------------------------------------------------------- validation and leakage audit

def known(values):
    return {v for v in values if v not in (NA, "UNKNOWN", None)}


def validate_task(task: dict, config: dict) -> None:
    """Raise LeakageError when a task breaks its strategy's separation rule."""
    train, test, strategy = task["train"], task["test"], task["strategy"]
    tr_ids, ts_ids = [e["isolate_id"] for e in train], [e["isolate_id"] for e in test]
    if len(set(tr_ids)) != len(tr_ids) or len(set(ts_ids)) != len(ts_ids):
        raise LeakageError(f"{task['task_id']}: duplicate isolate within a partition")
    if set(tr_ids) & set(ts_ids):
        raise LeakageError(f"{task['task_id']}: isolate in both training and test")
    if not all(e["genome_verified"] for e in train + test):
        raise LeakageError(f"{task['task_id']}: isolate without a verified genome")
    assemblies = Counter(e["assembly_accession"] for e in train + test)
    if any(c > 1 for c in assemblies.values()):
        raise LeakageError(f"{task['task_id']}: one assembly used for two isolates")
    if strategy == "random_cluster":
        if {e["split_unit"] for e in train} & {e["split_unit"] for e in test}:
            raise LeakageError(f"{task['task_id']}: split unit on both sides")
    if strategy == "lineage":
        tr = [e for e in train if e["lineage_test_eligible"]]
        if not all(e["lineage_test_eligible"] for e in test):
            raise LeakageError(f"{task['task_id']}: isolate without a lineage component in test")
        if {e["lineage_component"] for e in tr} & {e["lineage_component"] for e in test}:
            raise LeakageError(f"{task['task_id']}: lineage component on both sides")
        if known(e["lineage_id"] for e in train) & known(e["lineage_id"] for e in test):
            raise LeakageError(f"{task['task_id']}: held-out ST in training")
        if known(e["snp_cluster"] for e in train) & known(e["snp_cluster"] for e in test):
            raise LeakageError(f"{task['task_id']}: SNP cluster on both sides")
    if strategy == "country":
        if any(e["country"] == NA for e in train) or {e["country"] for e in train} & {e["country"] for e in test}:
            raise LeakageError(f"{task['task_id']}: held-out country in training")
    if strategy == "project":
        if {e["project"] for e in train} & {e["project"] for e in test}:
            raise LeakageError(f"{task['task_id']}: held-out project in training")
    if strategy == "geography":
        ref = config["splits"]["strategies"]["geography"]["reference_country"]
        if any(e["country"] != ref for e in train) or any(e["country"] in (ref, NA) for e in test):
            raise LeakageError(f"{task['task_id']}: geography partition mixes countries")
    if strategy == "temporal":
        if any(e["year"] is None for e in train + test) or \
                max(e["year"] for e in train) >= min(e["year"] for e in test):
            raise LeakageError(f"{task['task_id']}: training years do not precede test years")


def validate_cv(tasks: list[dict], endpoints: list[dict]) -> None:
    """Folds of a CV strategy must cover each test-eligible endpoint exactly once."""
    groups = defaultdict(list)
    for t in tasks:
        if t["strategy"] in CV_STRATEGIES:
            groups[(t["cohort"], t["antimicrobial"], t["strategy"])].append(t)
    for (cohort, drug, strategy), ts in groups.items():
        seen = Counter(e["isolate_id"] for t in ts for e in t["test"])
        if any(c != 1 for c in seen.values()):
            raise LeakageError(f"{cohort} {drug} {strategy}: isolate tested in more than one fold")
        pool, _ = cohort_members([e for e in endpoints if e["antimicrobial"] == drug], cohort, strategy)
        expected = {e["isolate_id"] for e in pool if strategy != "lineage" or e["lineage_test_eligible"]}
        if set(seen) != expected:
            raise LeakageError(f"{cohort} {drug} {strategy}: folds do not reconcile with the cohort")
        for t in ts:
            if len(t["train"]) + len(t["test"]) != len(pool):
                raise LeakageError(f"{t['task_id']}: training and test do not add up to the cohort")


def share(test: list[dict], train: list[dict], key: str) -> tuple[int, int]:
    """(test isolates with a known value, of those how many share it with training)."""
    seen = known(e[key] for e in train)
    vals = [e[key] for e in test if e[key] not in (NA, "UNKNOWN", None)]
    return len(vals), sum(v in seen for v in vals)


def leakage_audit(task: dict) -> dict:
    train, test = task["train"], task["test"]
    out = {k: task[k] for k in ("task_id", "cohort", "antimicrobial", "strategy", "fold")}
    out["n_train"], out["n_test"] = len(train), len(test)
    out["isolate_overlap"] = len({e["isolate_id"] for e in train} & {e["isolate_id"] for e in test})
    for key, name in (("lineage_id", "st"), ("snp_cluster", "snp_cluster"), ("lineage_component", "component"),
                      ("project", "project"), ("country", "country"), ("period", "period")):
        n, s = share(test, train, key)
        out[f"test_with_{name}"], out[f"test_{name}_in_train"] = n, s
        out[f"fraction_{name}_in_train"] = s / n if n else None
    train_years = [e["year"] for e in train if e["year"] is not None]
    test_years = [e["year"] for e in test if e["year"] is not None]
    out["train_year_range"] = f"{min(train_years)}-{max(train_years)}" if train_years else NA
    out["test_year_range"] = f"{min(test_years)}-{max(test_years)}" if test_years else NA
    out["test_after_train_years"] = sum(y > max(train_years) for y in test_years) if train_years else NA
    return out


# ---------------------------------------------------------------- features per task and model

def feature_filter(X: sparse.csr_matrix, rows: np.ndarray, min_present: int, min_absent: int) -> np.ndarray:
    """Columns present in at least min_present and absent from at least min_absent of `rows`."""
    sub = X[rows]
    counts = np.bincount(sub.indices, minlength=X.shape[1])
    return np.flatnonzero((counts >= min_present) & (rows.shape[0] - counts >= min_absent))


def fit_model(X, y, C: float, model: dict, seed: int):
    from sklearn.linear_model import LogisticRegression
    if len(set(y.tolist())) < 2:
        raise SingleClassError("training data hold one phenotype")
    clf = LogisticRegression(penalty=model["penalty"], solver=model["solver"], C=C,
                             class_weight=model["class_weight"], max_iter=model["max_iter"], tol=model["tol"],
                             intercept_scaling=model["intercept_scaling"], random_state=seed)
    clf.fit(X, y)
    return clf


def select_C(X, y, units: list, model: dict, features: dict, seed: int) -> tuple[float, list]:
    """Choose C by grouped inner cross-validation on training rows only."""
    grid = model["C_grid"]
    if len(grid) == 1:
        return grid[0], []
    k = model["inner_folds"]
    per_unit = defaultdict(lambda: [0, 0])
    for u, label in zip(units, y.tolist()):
        per_unit[u][0 if label == 1 else 1] += 1
    folds = assign_folds({u: tuple(v) for u, v in per_unit.items()}, k, seed, "inner")
    fold_of = np.array([folds[u] for u in units])
    scores = defaultdict(list)
    for f in range(k):
        tr, va = np.flatnonzero(fold_of != f), np.flatnonzero(fold_of == f)
        if len(set(y[va].tolist())) < 2 or len(set(y[tr].tolist())) < 2:
            continue
        cols = feature_filter(X, tr, features["min_train_present"], features["min_train_absent"])
        Xtr, Xva = X[tr][:, cols], X[va][:, cols]
        for C in grid:
            clf = fit_model(Xtr, y[tr], C, model, seed)
            scores[C].append(auroc(y[va], clf.predict_proba(Xva)[:, 1]))
    if not scores:
        return min(grid), []
    mean = {C: float(np.mean(v)) for C, v in scores.items()}
    best = max(grid, key=lambda C: (mean.get(C, -1.0), -C))
    return best, [[C, mean.get(C)] for C in grid]


def run_task(task: dict, X, row_of: dict, config: dict) -> tuple[dict, list[dict]]:
    """Fit on the task's training rows, predict its test rows."""
    model, features = config["model"], config["features"]
    seed = config["seeds"]["model"]
    started = time.time()
    tr_rows = np.array([row_of[e["isolate_id"]] for e in task["train"]], dtype=np.int64)
    ts_rows = np.array([row_of[e["isolate_id"]] for e in task["test"]], dtype=np.int64)
    y = np.array([e["phenotype"] == "R" for e in task["train"]], dtype=np.int8)
    if len(set(y.tolist())) < 2:
        raise SingleClassError("training data hold one phenotype")
    Xtr_all = X[tr_rows]
    C, inner = select_C(Xtr_all, y, [e["split_unit"] for e in task["train"]], model, features,
                        config["seeds"]["inner_cv"])
    cols = feature_filter(X, tr_rows, features["min_train_present"], features["min_train_absent"])
    clf = fit_model(Xtr_all[:, cols], y, C, model, seed)
    prob = clf.predict_proba(X[ts_rows][:, cols])[:, 1]
    threshold = model["threshold"]
    preds = [{"isolate_id": e["isolate_id"], "phenotype": e["phenotype"], "predicted_probability": float(p),
              "predicted_class": "R" if p >= threshold else "S"} for e, p in zip(task["test"], prob)]
    status = {"task_id": task["task_id"], "status": "success", "C": C, "inner_cv_auroc": inner,
              "features_used": int(cols.shape[0]), "vocabulary_sha256": sha256_bytes(cols.astype(np.int64).tobytes()),
              "nonzero_coefficients": int(np.count_nonzero(clf.coef_)), "n_train": int(tr_rows.shape[0]),
              "n_test": int(ts_rows.shape[0]), "train_R": int(y.sum()), "threshold": threshold,
              "split_sha256": task_digest(task), "seconds": round(time.time() - started, 2)}
    return status, preds


# ---------------------------------------------------------------- metrics

def auroc(y, p):
    y, p = np.asarray(y).astype(bool), np.asarray(p, dtype=float)
    n1 = int(y.sum())
    n0 = y.shape[0] - n1
    if n1 == 0 or n0 == 0:
        return None
    ranks = rankdata(p)
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auprc(y, p):
    """Average precision: sum over thresholds of recall gain times precision."""
    y, p = np.asarray(y).astype(bool), np.asarray(p, dtype=float)
    n1 = int(y.sum())
    if n1 == 0 or n1 == y.shape[0]:
        return None
    order = np.argsort(-p, kind="mergesort")
    ys, ps = y[order], p[order]
    last = np.r_[np.flatnonzero(ps[1:] != ps[:-1]), ps.shape[0] - 1]
    tp = np.cumsum(ys)[last]
    fp = (last + 1) - tp
    precision = tp / (tp + fp)
    recall = tp / n1
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def div(a, b):
    return a / b if b else None


def confusion_metrics(y, pred) -> dict:
    y, pred = np.asarray(y).astype(bool), np.asarray(pred).astype(bool)
    tp, tn = int(np.sum(y & pred)), int(np.sum(~y & ~pred))
    fp, fn = int(np.sum(~y & pred)), int(np.sum(y & ~pred))
    sens, spec = div(tp, tp + fn), div(tn, tn + fp)
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn, "sensitivity": sens, "specificity": spec,
            "PPV": div(tp, tp + fp), "NPV": div(tn, tn + fn),
            "balanced_accuracy": (sens + spec) / 2 if sens is not None and spec is not None else None,
            "F1": div(2 * tp, 2 * tp + fp + fn), "MCC": (tp * tn - fp * fn) / denom if denom else None}


def calibration(y, p) -> tuple:
    """Calibration intercept and slope from a logistic fit of y on logit(p)."""
    from sklearn.linear_model import LogisticRegression
    y = np.asarray(y).astype(int)
    if y.min() == y.max():
        return None, None
    q = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    z = np.log(q / (1 - q)).reshape(-1, 1)
    fit = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000).fit(z, y)
    return float(fit.intercept_[0]), float(fit.coef_[0][0])


def metrics(y, p, threshold: float, with_calibration: bool = True) -> dict:
    y, p = np.asarray(y).astype(bool), np.asarray(p, dtype=float)
    out = {"N": int(y.shape[0]), "N_R": int(y.sum()), "N_S": int((~y).sum())}
    out.update(confusion_metrics(y, p >= threshold))
    out["AUROC"], out["AUPRC"] = auroc(y, p), auprc(y, p)
    out["brier"] = float(np.mean((p - y) ** 2)) if y.shape[0] else None
    if with_calibration:
        out["calibration_intercept"], out["calibration_slope"] = calibration(y, p) if y.shape[0] else (None, None)
    return out


CI_METRICS = ("AUROC", "AUPRC", "balanced_accuracy", "MCC", "sensitivity", "specificity")


def quick(y, p, threshold) -> dict:
    out = confusion_metrics(y, p >= threshold)
    out["AUROC"], out["AUPRC"] = auroc(y, p), auprc(y, p)
    return out


def bootstrap_indices(units: list, replicates: int, seed: int):
    """Yield row indices for each replicate, resampling whole units."""
    codes = {}
    unit_idx = np.array([codes.setdefault(u, len(codes)) for u in units], dtype=np.int64)
    n_units = len(codes)
    order = np.argsort(unit_idx, kind="stable")
    starts = np.searchsorted(unit_idx[order], np.arange(n_units))
    sizes = np.bincount(unit_idx, minlength=n_units)
    rng = np.random.default_rng(seed)
    for _ in range(replicates):
        pick = rng.integers(0, n_units, n_units)
        s = sizes[pick]
        offsets = np.arange(s.sum()) - np.repeat(np.cumsum(s) - s, s)
        yield order[np.repeat(starts[pick], s) + offsets]


def bootstrap_ci(y, p, units, threshold, config, q=None) -> dict:
    """Percentile intervals; with q, also for the paired difference p - q."""
    boot = config["bootstrap"]
    alpha = (1 - boot["level"]) / 2
    y, p = np.asarray(y).astype(bool), np.asarray(p, dtype=float)
    q = None if q is None else np.asarray(q, dtype=float)
    store = defaultdict(list)
    for idx in bootstrap_indices(units, boot["replicates"], config["seeds"]["bootstrap"]):
        a = quick(y[idx], p[idx], threshold)
        b = quick(y[idx], q[idx], threshold) if q is not None else None
        for m in CI_METRICS:
            if a[m] is not None:
                store[m].append(a[m])
            if b is not None and a[m] is not None and b[m] is not None:
                store["delta_" + m].append(b[m] - a[m])
    out = {}
    for m, values in store.items():
        out[m + "_ci"] = (float(np.quantile(values, alpha)), float(np.quantile(values, 1 - alpha)))
    return out


# ---------------------------------------------------------------- run stage

def load_inputs(config_path: Path, population: Path, baseline: Path, verification: Path) -> dict:
    config = load_config(config_path)
    up = upstream_hashes(config, population, baseline)
    root = Path(__file__).resolve().parents[1]
    pop_config = pop.load_config(root / config["endpoints"]["population_config"])
    if pop.sha256(root / config["endpoints"]["population_config"]) != config["endpoints"]["population_config_sha256"]:
        raise ValueError("population config differs from the benchmark config")
    status = {r["isolate_id"]: r for r in pop.read_tsv(verification)}
    verified = {i for i, r in status.items() if r["status"] == "verified"}
    endpoints, plan = load_endpoints(baseline, population, verified, pop_config)
    return {"config": config, "upstream": up, "plan": plan, "endpoints": endpoints, "verification": status}


def write_json(path: Path, value) -> None:
    pop.write_json_atomic(Path(path), value)


def _worker(args):
    task, matrix_dir, config = args
    X, row_of, _ = load_matrix(Path(matrix_dir))
    try:
        return run_task(task, X, row_of, config)
    except SingleClassError as error:
        return {"task_id": task["task_id"], "status": "failed", "reason": f"single_class_training: {error}"}, []
    except Exception as error:  # noqa: BLE001 - a failed task is recorded, not fatal
        return {"task_id": task["task_id"], "status": "failed", "reason": f"{type(error).__name__}: {error}"}, []


PRED_FIELDS = "isolate_id phenotype predicted_probability predicted_class".split()


def run_tasks(tasks: list[dict], matrix_dir: Path, outdir: Path, config: dict, workers: int = 1,
              only: str | None = None) -> Counter:
    """Run tasks without a valid checkpoint; returns status counts."""
    taskdir = Path(outdir) / "tasks"
    taskdir.mkdir(parents=True, exist_ok=True)
    todo = []
    counts = Counter()
    for t in tasks:
        if only and not re.search(only, t["task_id"]):
            continue
        st = taskdir / f"{t['task_id']}.json"
        if st.exists():
            prev = json.loads(st.read_text(encoding="utf-8"))
            if prev.get("status") == "success" and prev.get("split_sha256") == task_digest(t):
                counts["reused"] += 1
                continue
        todo.append(t)
    slim = lambda t: {**t, "train": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["train"]],  # noqa: E731
                      "test": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["test"]]}

    def save(status, preds):
        name = status["task_id"]
        if status["status"] == "success":
            write_tsv(taskdir / f"{name}.tsv", PRED_FIELDS, preds)
        write_json(taskdir / f"{name}.json", {**status, "finished_at": now()})
        counts[status["status"]] += 1
        print(f"{now()} {name} {status['status']} {status.get('seconds', '')}s", flush=True)

    if workers <= 1:
        X, row_of, _ = load_matrix(Path(matrix_dir))
        for t in todo:
            try:
                save(*run_task(slim(t), X, row_of, config))
            except SingleClassError as error:
                save({"task_id": t["task_id"], "status": "failed", "reason": f"single_class_training: {error}"}, [])
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_worker, (slim(t), str(matrix_dir), config)) for t in todo]
            for future in as_completed(futures):
                save(*future.result())
    return counts


# ---------------------------------------------------------------- report stage

def novelty(task: dict) -> dict:
    """Per test isolate: is its ST, SNP cluster, country seen in training; is it later than training?"""
    train = task["train"]
    sts, clusters, countries = known(e["lineage_id"] for e in train), known(e["snp_cluster"] for e in train), \
        known(e["country"] for e in train)
    years = [e["year"] for e in train if e["year"] is not None]
    last = max(years) if years else None

    def flag(value, seen):
        return NA if value in (NA, "UNKNOWN", None) else ("yes" if value in seen else "no")

    return {e["isolate_id"]: {"lineage_seen_in_training": flag(e["lineage_id"], sts),
                              "snp_cluster_seen_in_training": flag(e["snp_cluster"], clusters),
                              "country_seen_in_training": flag(e["country"], countries),
                              "year_after_training_range": NA if e["year"] is None or last is None
                              else ("yes" if e["year"] > last else "no")} for e in task["test"]}


def det_pred(e) -> str:
    return {"positive": "R", "negative": "S"}.get(e["determinant_status"], NA)


def metric_row(base: dict, y, p, units, threshold, config, ci: bool, support_ok=None) -> dict:
    m = metrics(y, p, threshold)
    row = {**base, **m}
    if ci and m["N_R"] and m["N_S"]:
        for name, (lo, hi) in bootstrap_ci(y, p, units, threshold, config).items():
            row[name.replace("_ci", "_ci_lower")], row[name.replace("_ci", "_ci_upper")] = lo, hi
    sup = config["support"]
    row["supported"] = "yes" if m["N_R"] >= sup["metric_min_R"] and m["N_S"] >= sup["metric_min_S"] else "no"
    return row


METRIC_FIELDS = ("N N_R N_S TP FP TN FN AUROC AUPRC sensitivity specificity PPV NPV balanced_accuracy F1 MCC brier "
                 "calibration_intercept calibration_slope supported").split() + \
    [f"{m}_ci_{s}" for m in CI_METRICS for s in ("lower", "upper")]


def summarise(values) -> dict:
    v = [x for x in values if x is not None]
    if not v:
        return {"folds_defined": 0, "median": None, "q1": None, "q3": None, "min": None, "max": None}
    return {"folds_defined": len(v), "median": float(np.median(v)), "q1": float(np.quantile(v, 0.25)),
            "q3": float(np.quantile(v, 0.75)), "min": float(min(v)), "max": float(max(v))}


def distribution(values) -> dict:
    c = Counter(values)
    total = sum(c.values()) or 1
    return {k: round(v / total, 4) for k, v in sorted(c.items())}


def tvd(a: dict, b: dict) -> float:
    return 0.5 * sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


def report(inputs: dict, tasks: list[dict], skipped: list[dict], outdir: Path, manifest: dict) -> dict:
    config, outdir = inputs["config"], Path(outdir)
    threshold = config["model"]["threshold"]
    taskdir = outdir / "tasks"
    status, preds = {}, {}
    for t in tasks:
        path = taskdir / f"{t['task_id']}.json"
        st = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"status": "missing"}
        if st["status"] == "success" and st["split_sha256"] != task_digest(t):
            st = {"status": "stale", "reason": "split changed since the task was run"}
        status[t["task_id"]] = st
        if st["status"] == "success":
            preds[t["task_id"]] = {r["isolate_id"]: float(r["predicted_probability"])
                                   for r in pop.read_tsv(taskdir / f"{t['task_id']}.tsv")}
    by_key = defaultdict(list)
    for t in tasks:
        by_key[(t["cohort"], t["antimicrobial"], t["strategy"])].append(t)

    pred_rows, fold_rows, pooled_rows, drop_rows = [], [], [], []
    seen_rows, det_rows, unexplained_rows = [], [], []
    pooled_pred = {}  # (cohort, drug, strategy) -> {isolate: prob}
    for (cohort, drug, strategy), ts in sorted(by_key.items()):
        merged = {}
        for t in ts:
            if t["task_id"] not in preds:
                continue
            nov = novelty(t)
            for e in t["test"]:
                p = preds[t["task_id"]][e["isolate_id"]]
                merged[e["isolate_id"]] = p
                pred_rows.append({"benchmark_version": config["benchmark_version"], "cohort": cohort,
                                  "antimicrobial": drug, "split_strategy": strategy, "fold_or_holdout_id": t["fold"],
                                  "isolate_id": e["isolate_id"], "phenotype": e["phenotype"],
                                  "predicted_probability": p, "predicted_class": "R" if p >= threshold else "S",
                                  "correct": "yes" if (p >= threshold) == (e["phenotype"] == "R") else "no",
                                  "lineage": e["lineage_id"], "lineage_component": e["lineage_component"],
                                  "snp_cluster": e["snp_cluster"], "country": e["country"],
                                  "year": e["year"], "source": e["source"], "project": e["project"],
                                  "established_determinant_status": e["determinant_status"],
                                  "concordance_category": e["concordance_category"], **nov[e["isolate_id"]]})
            y = [e["phenotype"] == "R" for e in t["test"]]
            p = [preds[t["task_id"]][e["isolate_id"]] for e in t["test"]]
            base = {"cohort": cohort, "antimicrobial": drug, "split_strategy": strategy, "fold_or_holdout_id": t["fold"],
                    "n_train": len(t["train"]), "train_R": sum(e["phenotype"] == "R" for e in t["train"]),
                    "C": status[t["task_id"]].get("C"), "features_used": status[t["task_id"]].get("features_used"),
                    "nonzero_coefficients": status[t["task_id"]].get("nonzero_coefficients")}
            fold_rows.append(metric_row(base, y, p, [e["split_unit"] for e in t["test"]], threshold, config,
                                        ci=strategy not in CV_STRATEGIES))
        pooled_pred[(cohort, drug, strategy)] = merged
        complete = all(t["task_id"] in preds for t in ts)
        test_all = [e for t in ts for e in t["test"] if e["isolate_id"] in merged]
        y = [e["phenotype"] == "R" for e in test_all]
        p = [merged[e["isolate_id"]] for e in test_all]
        fold_auc = summarise([r["AUROC"] for r in fold_rows if (r["cohort"], r["antimicrobial"], r["split_strategy"])
                              == (cohort, drug, strategy)])
        base = {"cohort": cohort, "antimicrobial": drug, "split_strategy": strategy, "folds": len(ts),
                "folds_completed": sum(t["task_id"] in preds for t in ts), "complete": "yes" if complete else "no",
                **{f"fold_auroc_{k}": v for k, v in fold_auc.items()}}
        if test_all:
            pooled_rows.append(metric_row(base, y, p, [e["split_unit"] for e in test_all], threshold, config, ci=True))
        # seen versus unseen lineage within this strategy
        nov_all = {}
        for t in ts:
            nov_all.update(novelty(t))
        for indicator in ("lineage_seen_in_training", "snp_cluster_seen_in_training"):
            for value in ("yes", "no"):
                sub = [e for e in test_all if nov_all[e["isolate_id"]][indicator] == value]
                if sub:
                    seen_rows.append(metric_row({"cohort": cohort, "antimicrobial": drug, "split_strategy": strategy,
                                                 "indicator": indicator, "value": value},
                                                [e["phenotype"] == "R" for e in sub],
                                                [merged[e["isolate_id"]] for e in sub],
                                                [e["split_unit"] for e in sub], threshold, config, ci=False))
        # determinant baseline on the same endpoints
        sub = [e for e in test_all if det_pred(e) != NA]
        if sub:
            y_s = np.array([e["phenotype"] == "R" for e in sub])
            ml = np.array([merged[e["isolate_id"]] for e in sub])
            dp = np.array([1.0 if det_pred(e) == "R" else 0.0 for e in sub])
            mm, dm = metrics(y_s, ml, threshold, False), metrics(y_s, dp, 0.5, False)
            row = {"cohort": cohort, "antimicrobial": drug, "split_strategy": strategy, "N": len(sub),
                   "N_R": int(y_s.sum()), "N_S": int((~y_s).sum())}
            for m in ("sensitivity", "specificity", "balanced_accuracy", "MCC", "AUROC", "PPV", "NPV"):
                row[f"determinant_{m}"], row[f"ml_{m}"] = dm[m], mm[m]
                row[f"delta_{m}"] = None if dm[m] is None or mm[m] is None else mm[m] - dm[m]
            cats = Counter(e["concordance_category"] for e in sub)
            ml_r = {e["isolate_id"]: merged[e["isolate_id"]] >= threshold for e in sub}
            for cat in ("R_DET_POS", "R_DET_NEG", "S_DET_POS", "S_DET_NEG"):
                row[f"n_{cat}"] = cats[cat]
                row[f"ml_R_{cat}"] = sum(ml_r[e["isolate_id"]] for e in sub if e["concordance_category"] == cat)
            if row["N_R"] and row["N_S"]:
                ci = bootstrap_ci(y_s, dp, [e["split_unit"] for e in sub], threshold, config, q=ml)
                for m in ("balanced_accuracy", "MCC", "sensitivity", "specificity"):
                    if "delta_" + m + "_ci" in ci:
                        row[f"delta_{m}_ci_lower"], row[f"delta_{m}_ci_upper"] = ci["delta_" + m + "_ci"]
            det_rows.append(row)
            for lineage_value in ("all", "yes", "no"):
                rdn = [e for e in sub if e["concordance_category"] == "R_DET_NEG" and
                       (lineage_value == "all" or nov_all[e["isolate_id"]]["lineage_seen_in_training"] == lineage_value)]
                sdp = [e for e in sub if e["concordance_category"] == "S_DET_POS" and
                       (lineage_value == "all" or nov_all[e["isolate_id"]]["lineage_seen_in_training"] == lineage_value)]
                unexplained_rows.append({
                    "cohort": cohort, "antimicrobial": drug, "split_strategy": strategy,
                    "lineage_seen_in_training": lineage_value, "n_R_DET_NEG": len(rdn),
                    "ml_predicted_R": sum(ml_r[e["isolate_id"]] for e in rdn),
                    "ml_sensitivity_R_DET_NEG": div(sum(ml_r[e["isolate_id"]] for e in rdn), len(rdn)),
                    "n_S_DET_POS": len(sdp), "ml_predicted_S": sum(not ml_r[e["isolate_id"]] for e in sdp),
                    "ml_specificity_S_DET_POS": div(sum(not ml_r[e["isolate_id"]] for e in sdp), len(sdp))})

    # generalisation drops: holdout against random predictions on the same isolates
    drop_metrics = ("AUROC", "AUPRC", "balanced_accuracy", "MCC", "sensitivity", "specificity")
    for (cohort, drug, strategy), ts in sorted(by_key.items()):
        if strategy in ("random", "random_cluster"):
            continue
        hold = pooled_pred[(cohort, drug, strategy)]
        for comparator in ("random", "random_cluster"):
            ref = pooled_pred.get((cohort, drug, comparator), {})
            levels = [("pooled", "all", [e for t in ts for e in t["test"]])] + \
                     [("fold", t["fold"], t["test"]) for t in ts] if len(ts) > 1 else \
                [("pooled", "all", [e for t in ts for e in t["test"]])]
            for level, fold, test in levels:
                common = [e for e in test if e["isolate_id"] in hold and e["isolate_id"] in ref]
                if not common:
                    continue
                y = np.array([e["phenotype"] == "R" for e in common])
                pr = np.array([ref[e["isolate_id"]] for e in common])
                ph = np.array([hold[e["isolate_id"]] for e in common])
                a, b = metrics(y, pr, threshold, False), metrics(y, ph, threshold, False)
                row = {"cohort": cohort, "antimicrobial": drug, "holdout_strategy": strategy, "comparator": comparator,
                       "comparison": "same_test_isolates", "level": level, "fold_or_holdout_id": fold,
                       "N": len(common), "N_R": int(y.sum()), "N_S": int((~y).sum())}
                for m in drop_metrics:
                    row[f"{m}_comparator"], row[f"{m}_holdout"] = a[m], b[m]
                    row[f"delta_{m}"] = None if a[m] is None or b[m] is None else a[m] - b[m]
                if level == "pooled" and row["N_R"] and row["N_S"]:
                    # Bootstrap of comparator minus holdout, resampling split units
                    ci = bootstrap_ci(y, ph, [e["split_unit"] for e in common], threshold, config, q=pr)
                    for m in drop_metrics:
                        if "delta_" + m + "_ci" in ci:
                            row[f"delta_{m}_ci_lower"], row[f"delta_{m}_ci_upper"] = ci["delta_" + m + "_ci"]
                drop_rows.append(row)

    # random-split optimism (wide), common cohort first
    pooled_index = {(r["cohort"], r["antimicrobial"], r["split_strategy"]): r for r in pooled_rows}
    drop_index = {(r["cohort"], r["antimicrobial"], r["holdout_strategy"], r["comparator"]): r
                  for r in drop_rows if r["level"] == "pooled"}
    optimism = []
    for cohort in COHORTS:
        for drug in sorted({k[1] for k in pooled_index if k[0] == cohort}):
            for m in ("AUROC", "AUPRC", "balanced_accuracy", "MCC"):
                row = {"cohort": cohort, "antimicrobial": drug, "metric": m}
                rnd = pooled_index.get((cohort, drug, "random"), {}).get(m)
                row["random_all_test"] = rnd
                row["random_cluster_all_test"] = pooled_index.get((cohort, drug, "random_cluster"), {}).get(m)
                for s in STRATEGIES[2:]:
                    h = pooled_index.get((cohort, drug, s), {}).get(m)
                    row[s] = h
                    row[f"random_all_minus_{s}"] = None if rnd is None or h is None else rnd - h
                    d = drop_index.get((cohort, drug, s, "random"), {})
                    row[f"random_same_isolates_minus_{s}"] = d.get(f"delta_{m}")
                    row[f"random_same_isolates_minus_{s}_ci_lower"] = d.get(f"delta_{m}_ci_lower")
                    row[f"random_same_isolates_minus_{s}_ci_upper"] = d.get(f"delta_{m}_ci_upper")
                optimism.append(row)

    # shift audits
    country_rows, temporal_rows, audit_rows = [], [], []
    for t in tasks:
        audit_rows.append(leakage_audit(t))
        train, test = t["train"], t["test"]
        if t["strategy"] in ("country", "geography", "project"):
            n, s = share(test, train, "lineage_id")
            sts_test = known(e["lineage_id"] for e in test)
            country_rows.append({**{k: t[k] for k in ("cohort", "antimicrobial", "strategy", "fold")},
                                 "test_N": len(test), "test_R": rs(test)[0], "test_S": rs(test)[1],
                                 "test_lineages": len(sts_test),
                                 "test_lineages_unseen": len(sts_test - known(e["lineage_id"] for e in train)),
                                 "fraction_test_isolates_lineage_unseen": div(n - s, n),
                                 "fraction_test_clustered_isolates_cluster_seen": div(*reversed(share(test, train, "snp_cluster"))),
                                 "source_composition": pop.canonical(distribution(e["source"] for e in test)),
                                 "year_distribution": pop.canonical(distribution(str(e["year"]) for e in test)),
                                 "project_composition": pop.canonical(distribution(e["project"] for e in test)),
                                 "country_composition": pop.canonical(distribution(e["country"] for e in test)),
                                 "train_R_fraction": div(rs(train)[0], len(train)),
                                 "test_R_fraction": div(rs(test)[0], len(test))})
        if t["strategy"] == "temporal":
            n, s = share(test, train, "lineage_id")
            nc, sc = share(test, train, "snp_cluster")
            ctr, cts = distribution(e["country"] for e in train), distribution(e["country"] for e in test)
            str_, sts_ = distribution(e["source"] for e in train), distribution(e["source"] for e in test)
            temporal_rows.append({**{k: t[k] for k in ("cohort", "antimicrobial", "strategy", "fold")},
                                  "test_N": len(test), "test_R": rs(test)[0], "test_S": rs(test)[1],
                                  "fraction_test_isolates_lineage_unseen": div(n - s, n),
                                  "fraction_test_clustered_isolates_cluster_unseen": div(nc - sc, nc),
                                  "train_R_fraction": div(rs(train)[0], len(train)),
                                  "test_R_fraction": div(rs(test)[0], len(test)),
                                  "country_tvd": tvd(ctr, cts), "source_tvd": tvd(str_, sts_),
                                  "test_country_composition": pop.canonical(cts),
                                  "test_source_composition": pop.canonical(sts_),
                                  "test_project_composition": pop.canonical(distribution(e["project"] for e in test))})

    # manuscript table: common cohort, pooled, drop against random on the same isolates
    manuscript = []
    for r in pooled_rows:
        if r["cohort"] != "common":
            continue
        d = drop_index.get(("common", r["antimicrobial"], r["split_strategy"], "random"), {})
        manuscript.append({"antimicrobial": r["antimicrobial"], "split_strategy": r["split_strategy"],
                           "test_N": r["N"], "R": r["N_R"], "S": r["N_S"], "AUROC": r["AUROC"], "AUPRC": r["AUPRC"],
                           "balanced_accuracy": r["balanced_accuracy"], "MCC": r["MCC"],
                           "sensitivity": r["sensitivity"], "specificity": r["specificity"],
                           "AUROC_CI_lower": r.get("AUROC_ci_lower"), "AUROC_CI_upper": r.get("AUROC_ci_upper"),
                           "delta_AUROC_vs_random_same_isolates": d.get("delta_AUROC"),
                           "delta_CI_lower": d.get("delta_AUROC_ci_lower"), "delta_CI_upper": d.get("delta_AUROC_ci_upper"),
                           "delta_balanced_accuracy_vs_random_same_isolates": d.get("delta_balanced_accuracy"),
                           "folds": r["folds"], "complete": r["complete"]})

    outdir.mkdir(parents=True, exist_ok=True)
    fold_fields = ["cohort", "antimicrobial", "split_strategy", "fold_or_holdout_id", "n_train", "train_R", "C",
                   "features_used", "nonzero_coefficients"] + METRIC_FIELDS
    pooled_fields = ["cohort", "antimicrobial", "split_strategy", "folds", "folds_completed", "complete"] + \
        METRIC_FIELDS + [f"fold_auroc_{k}" for k in ("folds_defined", "median", "q1", "q3", "min", "max")]
    write_tsv(outdir / "fold_metrics.tsv", fold_fields, fold_rows)
    write_tsv(outdir / "benchmark_metrics.tsv", pooled_fields, pooled_rows)
    drop_fields = ["cohort", "antimicrobial", "holdout_strategy", "comparator", "comparison", "level",
                   "fold_or_holdout_id", "N", "N_R", "N_S"]
    for m in drop_metrics:
        drop_fields += [f"{m}_comparator", f"{m}_holdout", f"delta_{m}", f"delta_{m}_ci_lower", f"delta_{m}_ci_upper"]
    write_tsv(outdir / "generalisation_drop.tsv", drop_fields, drop_rows)
    opt_fields = ["cohort", "antimicrobial", "metric", "random_all_test", "random_cluster_all_test"]
    for s in STRATEGIES[2:]:
        opt_fields += [s, f"random_all_minus_{s}", f"random_same_isolates_minus_{s}",
                       f"random_same_isolates_minus_{s}_ci_lower", f"random_same_isolates_minus_{s}_ci_upper"]
    write_tsv(outdir / "random_split_optimism.tsv", opt_fields, optimism)
    pred_fields = ("benchmark_version cohort antimicrobial split_strategy fold_or_holdout_id isolate_id phenotype "
                   "predicted_probability predicted_class correct lineage lineage_component snp_cluster country year "
                   "source project established_determinant_status concordance_category lineage_seen_in_training "
                   "snp_cluster_seen_in_training country_seen_in_training year_after_training_range").split()
    write_tsv(outdir / "predictions.tsv.gz", pred_fields, pred_rows)
    write_tsv(outdir / "unseen_lineage_summary.tsv",
              ["cohort", "antimicrobial", "split_strategy", "indicator", "value"] + METRIC_FIELDS, seen_rows)
    if country_rows:
        write_tsv(outdir / "country_shift_summary.tsv", list(country_rows[0]), country_rows)
    if temporal_rows:
        write_tsv(outdir / "temporal_shift_summary.tsv", list(temporal_rows[0]), temporal_rows)
    write_tsv(outdir / "split_leakage_audit.tsv", list(audit_rows[0]), audit_rows)
    det_fields = ["cohort", "antimicrobial", "split_strategy", "N", "N_R", "N_S"]
    for m in ("sensitivity", "specificity", "balanced_accuracy", "MCC", "AUROC", "PPV", "NPV"):
        det_fields += [f"determinant_{m}", f"ml_{m}", f"delta_{m}"]
    for m in ("balanced_accuracy", "MCC", "sensitivity", "specificity"):
        det_fields += [f"delta_{m}_ci_lower", f"delta_{m}_ci_upper"]
    for cat in ("R_DET_POS", "R_DET_NEG", "S_DET_POS", "S_DET_NEG"):
        det_fields += [f"n_{cat}", f"ml_R_{cat}"]
    write_tsv(outdir / "determinant_vs_ml.tsv", det_fields, det_rows)
    write_tsv(outdir / "unexplained_resistance_ml.tsv",
              ["cohort", "antimicrobial", "split_strategy", "lineage_seen_in_training", "n_R_DET_NEG",
               "ml_predicted_R", "ml_sensitivity_R_DET_NEG", "n_S_DET_POS", "ml_predicted_S",
               "ml_specificity_S_DET_POS"], unexplained_rows)
    write_tsv(outdir / "manuscript_table.tsv", list(manuscript[0]) if manuscript else ["antimicrobial"], manuscript)
    write_tsv(outdir / "not_feasible.tsv", ["cohort", "antimicrobial", "strategy", "fold", "status", "reason"],
              skipped)
    figure_data(outdir, pooled_rows, drop_rows, fold_rows, seen_rows)

    failed = {k: v for k, v in status.items() if v["status"] != "success"}
    summary = {"benchmark_version": config["benchmark_version"], "generated_at": now(),
               "tasks": len(tasks), "tasks_succeeded": len(tasks) - len(failed),
               "tasks_not_succeeded": {k: v.get("reason", v["status"]) for k, v in sorted(failed.items())},
               "not_feasible": len(skipped), "prediction_rows": len(pred_rows),
               "manifest_sha256": text_sha(manifest)}
    write_json(outdir / "benchmark_summary.json", summary)
    return summary


def figure_data(outdir: Path, pooled, drops, folds, seen) -> None:
    fig = outdir / "figure_data"
    fig.mkdir(exist_ok=True)
    keep = ("cohort", "antimicrobial", "split_strategy", "N", "N_R", "N_S", "AUROC", "AUROC_ci_lower",
            "AUROC_ci_upper", "balanced_accuracy", "balanced_accuracy_ci_lower", "balanced_accuracy_ci_upper", "MCC")
    write_tsv(fig / "fig_a_performance_by_strategy.tsv", list(keep), pooled)
    dkeep = ("cohort", "antimicrobial", "holdout_strategy", "comparator", "N", "delta_AUROC", "delta_AUROC_ci_lower",
             "delta_AUROC_ci_upper", "delta_balanced_accuracy", "delta_balanced_accuracy_ci_lower",
             "delta_balanced_accuracy_ci_upper")
    write_tsv(fig / "fig_b_drop_vs_random.tsv", list(dkeep), [d for d in drops if d["level"] == "pooled"])
    fkeep = ("cohort", "antimicrobial", "split_strategy", "fold_or_holdout_id", "N", "N_R", "N_S", "AUROC",
             "balanced_accuracy", "MCC", "supported")
    for name, strategy in (("fig_c_lineage_folds", "lineage"), ("fig_d_country_holdouts", "country"),
                           ("fig_e_temporal", "temporal")):
        write_tsv(fig / f"{name}.tsv", list(fkeep), [f for f in folds if f["split_strategy"] == strategy])
    skeep = ("cohort", "antimicrobial", "split_strategy", "indicator", "value", "N", "N_R", "N_S", "AUROC",
             "balanced_accuracy", "sensitivity", "specificity", "supported")
    write_tsv(fig / "fig_f_seen_vs_unseen_lineage.tsv", list(skeep), seen)


# ---------------------------------------------------------------- freeze stage

def eligibility_table(endpoints: list[dict], plan: dict) -> list[dict]:
    rows = []
    by = defaultdict(list)
    for e in endpoints:
        by[e["antimicrobial"]].append(e)
    for drug, eps in sorted(by.items()):
        verified = [e for e in eps if e["genome_verified"]]
        common, reasons = cohort_members(eps, "common", "random")
        r, s = rs(verified)
        row = {"antimicrobial": drug, "tier": plan[drug]["tier"], "endpoints": len(eps), "genome_verified": len(verified),
               "R": r, "S": s, "R_fraction": div(r, r + s), "common_cohort_N": len(common),
               "common_cohort_R": rs(common)[0], "common_cohort_S": rs(common)[1],
               "common_exclusions": pop.canonical(reasons)}
        for strategy, name in (("random", "random"), ("lineage", "lineage"), ("country", "country"),
                               ("temporal", "temporal")):
            pool, _ = cohort_members(eps, "maximal", strategy)
            if strategy == "lineage":
                pool = [e for e in pool if e["lineage_test_eligible"]]
            row[f"{name}_eligible_N"] = len(pool)
        row["lineages"] = len(known(e["lineage_id"] for e in verified))
        row["lineage_components"] = len({e["lineage_component"] for e in verified if e["lineage_test_eligible"]})
        row["countries"] = len(known(e["country"] for e in verified))
        years = [e["year"] for e in verified if e["year"] is not None]
        row["year_range"] = f"{min(years)}-{max(years)}" if years else NA
        for key in ("project_holdout", "geographic_transfer", "temporal"):
            row[f"plan_{key}"] = plan[drug][key]
        rows.append(row)
    return rows


def split_summary(task: dict) -> dict:
    tr, ts = rs(task["train"]), rs(task["test"])
    return {**{k: task[k] for k in ("task_id", "cohort", "antimicrobial", "strategy", "fold")},
            "train_N": len(task["train"]), "train_R": tr[0], "train_S": tr[1], "test_N": len(task["test"]),
            "test_R": ts[0], "test_S": ts[1], "split_sha256": task_digest(task), "status": "planned"}


SPLIT_FIELDS = ("benchmark_version cohort antimicrobial isolate_id split_strategy fold_or_holdout_id partition "
                "phenotype lineage lineage_component snp_cluster split_unit country year source project").split()


def freeze(inputs: dict, outdir: Path, config_path: Path, feature_manifest: dict) -> dict:
    config, outdir = inputs["config"], Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tasks, skipped = build_tasks(inputs["endpoints"], inputs["plan"], config)
    for t in tasks:
        validate_task(t, config)
    validate_cv(tasks, inputs["endpoints"])
    splitdir = outdir / "splits"
    splitdir.mkdir(exist_ok=True)
    for strategy in STRATEGIES:
        rows = []
        for t in tasks:
            if t["strategy"] != strategy:
                continue
            for part in ("train", "test"):
                for e in t[part]:
                    rows.append({"benchmark_version": config["benchmark_version"], "cohort": t["cohort"],
                                 "antimicrobial": t["antimicrobial"], "isolate_id": e["isolate_id"],
                                 "split_strategy": strategy, "fold_or_holdout_id": t["fold"], "partition": part,
                                 "phenotype": e["phenotype"], "lineage": e["lineage_id"],
                                 "lineage_component": e["lineage_component"], "snp_cluster": e["snp_cluster"],
                                 "split_unit": e["split_unit"], "country": e["country"], "year": e["year"],
                                 "source": e["source"], "project": e["project"]})
        write_tsv(splitdir / f"{strategy}.tsv.gz", SPLIT_FIELDS, rows)
    summaries = [split_summary(t) for t in tasks]
    write_tsv(outdir / "split_summary.tsv", list(summaries[0]), summaries)
    write_tsv(outdir / "not_feasible.tsv", ["cohort", "antimicrobial", "strategy", "fold", "status", "reason"], skipped)
    elig = eligibility_table(inputs["endpoints"], inputs["plan"])
    write_tsv(outdir / "cohort_eligibility.tsv", list(elig[0]), elig)
    plan_rows = sorted(inputs["plan"].values(), key=lambda p: p["antimicrobial"])
    write_tsv(outdir / "drug_plan_corrected.tsv", list(plan_rows[0]), plan_rows)
    audit = [leakage_audit(t) for t in tasks]
    write_tsv(outdir / "split_leakage_audit.tsv", list(audit[0]), audit)
    endpoints = inputs["endpoints"]
    manifest = {
        "benchmark_version": config["benchmark_version"], "builder_version": VERSION,
        "builder_sha256": pop.sha256(Path(__file__)), "config_sha256": config_sha(config_path),
        "config_path": str(Path(config_path).name), "decision_record": config["decision_record"],
        "upstream": inputs["upstream"], "genome_manifest_sha256": config["upstream"]["genome_manifest_sha256"],
        "genome_verification": dict(Counter(r["status"] for r in inputs["verification"].values())),
        "feature_manifest": {k: feature_manifest[k] for k in sorted(feature_manifest) if k != "created_at"},
        "phenotype_endpoint_definition": config["endpoints"]["definition"],
        "antimicrobials": {d: p["tier"] for d, p in sorted(inputs["plan"].items()) if p["tier"] != "excluded"},
        "endpoints": len(endpoints), "endpoints_genome_verified": sum(e["genome_verified"] for e in endpoints),
        "isolates": len({e["isolate_id"] for e in endpoints}),
        "endpoint_list_sha256": text_sha(sorted([e["antimicrobial"], e["isolate_id"], e["phenotype"]]
                                               for e in endpoints)),
        "seeds": config["seeds"], "splits": config["splits"], "cohorts": config["cohorts"],
        "model": config["model"], "features": config["features"],
        "tasks": len(tasks), "not_feasible": len(skipped),
        "task_split_sha256": {t["task_id"]: task_digest(t) for t in tasks},
        "frozen_at": now(),
    }
    write_json(outdir / "benchmark_manifest.json", manifest)
    return manifest


def check_frozen(tasks: list[dict], manifest: dict, config_path: Path) -> None:
    if manifest["config_sha256"] != config_sha(config_path):
        raise ValueError("benchmark config changed after freezing; start a new benchmark version")
    expected = manifest["task_split_sha256"]
    got = {t["task_id"]: task_digest(t) for t in tasks}
    if got != expected:
        changed = sorted(set(got.items()) ^ set(expected.items()))[:5]
        raise LeakageError(f"rebuilt splits differ from the frozen manifest: {changed}")


# ---------------------------------------------------------------- command line

def cohort_assemblies(baseline: Path) -> dict:
    out = {}
    with (Path(baseline) / "genotype_phenotype_endpoints.tsv").open(encoding="utf-8") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            out[r["isolate_id"]] = r["assembly_accession"]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("sketch", "matrix", "freeze", "run", "report"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--population", type=Path, required=True, help="PR #6 population build directory")
    parser.add_argument("--baseline", type=Path, required=True, help="PR #7 baseline build directory")
    parser.add_argument("--cache", type=Path, help="genome cache (sketch)")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--only", help="regular expression on task IDs (run)")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    outdir = args.outdir
    if args.command == "sketch":
        if args.cache is None:
            parser.error("sketch needs --cache")
        upstream_hashes(config, args.population, args.baseline)
        summary = sketch_all(config, args.cache, outdir, cohort_assemblies(args.baseline), args.workers, args.limit)
        print(json.dumps(summary["status_counts"]))
        return 0 if set(summary["status_counts"]) <= {"verified"} else 1
    verification = outdir / "genome_verification.tsv"
    if args.command == "matrix":
        rows = [(r["isolate_id"], r["assembly_accession"]) for r in pop.read_tsv(verification)
                if r["status"] == "verified"]
        f = config["features"]
        built = build_matrix(outdir / "sketches", rows, f["min_train_present"], f["min_train_absent"],
                             f["deduplicate_identical_columns"], scaled=f.get("matrix_scaled"))
        meta = {"params": sketch_params(config), "genome_manifest_sha256": config["upstream"]["genome_manifest_sha256"],
                "genome_verification_sha256": pop.sha256(verification), "global_prune": {
                    "min_present": f["min_train_present"], "min_absent": f["min_train_absent"],
                    "note": "drops only features every per-task training filter would drop"},
                "deduplicated": f["deduplicate_identical_columns"], "builder_version": VERSION}
        manifest = save_matrix(outdir / "matrix", built, rows, meta)
        print(json.dumps({k: manifest[k] for k in ("rows", "observed_features", "kept_features", "patterns", "nnz")}))
        return 0
    inputs = load_inputs(args.config, args.population, args.baseline, verification)
    feature_manifest = json.loads((outdir / "matrix" / "feature_manifest.json").read_text(encoding="utf-8"))
    if feature_manifest["genome_verification_sha256"] != pop.sha256(verification):
        raise ValueError("feature matrix was built from a different genome verification")
    if args.command == "freeze":
        manifest = freeze(inputs, outdir, args.config, feature_manifest)
        print(json.dumps({"tasks": manifest["tasks"], "not_feasible": manifest["not_feasible"]}))
        return 0
    manifest = json.loads((outdir / "benchmark_manifest.json").read_text(encoding="utf-8"))
    tasks, skipped = build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
    check_frozen(tasks, manifest, args.config)
    if args.command == "run":
        counts = run_tasks(tasks, outdir / "matrix", outdir, inputs["config"], args.workers, args.only)
        print(json.dumps(dict(counts)))
        return 0 if not counts.get("failed") else 1
    summary = report(inputs, tasks, skipped, outdir, manifest)
    print(json.dumps({k: summary[k] for k in ("tasks", "tasks_succeeded", "prediction_rows")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
