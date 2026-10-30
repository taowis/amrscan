#!/usr/bin/env python3
"""Representation sensitivity of the E. coli generalisation benchmark (representation_sensitivity_v1).

Refits the frozen generalisation_v1 tasks with a denser FracMinHash
representation (R1), and the error_decomposition_v1 mutation-augmented task
set with R1 plus the same known-mutation indicators (R3, R3_permuted), and
compares them with the original features (R0) and the reused augmentation
(R2, R2_permuted) on identical endpoints.

Commands, run in this order:
  check    verify every frozen upstream file and write the analysis manifest
  sketch   verify every genome against the PR #8 manifest, sketch it once at
           the densest candidate and check that downsampling reproduces the
           PR #14 sketch hash for hash (resumable)
  pilot    resources and k-mer capture for every candidate density; reads no
           prediction and no test phenotype
  fitpilot time and memory of model fits on the largest tasks at one density;
           the predictions are discarded unscored
  select   apply the frozen selection rule to the pilot and a budget, write
           selection.json
  matrix   build the R1 feature matrix
  run      refit the frozen tasks for R1, R3 or R3_permuted (resumable)
  report   benchmark tables for R1, then the paired comparison tables
  pairwise direct paired comparisons between the new representations
           (R2 and R3 against R1, R3 against R3_permuted and R2)

Nothing here changes a generalisation_v1 or error_decomposition_v1 output.
"""

import argparse
import csv
import gzip
import hashlib
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import json
import os
from pathlib import Path
import re
import resource
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import error_decomposition as ed  # noqa: E402
import generalisation_benchmark as gb  # noqa: E402
import population_structure as pop  # noqa: E402

VERSION = "1.0.0"
NA = "NA"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1.json"
REPRESENTATIONS = ("R0", "R1", "R2", "R2_permuted", "R3", "R3_permuted")
NEW_RUNS = ("R1", "R3", "R3_permuted")
GENOME_STATUSES = ("genome_verified", "genome_missing", "genome_hash_mismatch", "genome_unreadable",
                   "genome_version_mismatch")

write_tsv, now, sha256_bytes = gb.write_tsv, gb.now, gb.sha256_bytes


class UpstreamError(ValueError):
    pass


class IdentityError(ValueError):
    pass


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def load_config(path: Path) -> dict:
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("analysis_version", "upstream", "representations", "density", "benchmark", "mutation_panel",
                "capture", "analysis", "bootstrap", "support", "decision_record"):
        if key not in config:
            raise ValueError(f"analysis config lacks {key}")
    d = config["density"]
    if d["sketch_scaled"] != min(d["candidates"]):
        raise ValueError("the sketch must be taken at the densest candidate")
    if any(d["reference_scaled"] % s for s in d["candidates"]):
        raise ValueError("every candidate must divide the reference scaled value")
    if sorted(config["representations"]) != sorted(REPRESENTATIONS):
        raise ValueError("unexpected representation set")
    return config


def max_rss() -> int:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(r if sys.platform == "darwin" else r * 1024)


def paths(outdir: Path) -> dict:
    o = Path(outdir)
    return {"root": o, "manifest": o / "representation_sensitivity_manifest.json", "sketches": o / "sketches",
            "verification": o / "genome_verification.tsv", "pilot": o / "pilot", "selection": o / "selection.json",
            "R1": o / "R1", "R3": o / "R3", "R3_permuted": o / "R3_permuted"}


# ---------------------------------------------------------------- upstream

def check_upstream(config: dict, benchmark_dir: Path, decomposition_dir: Path) -> dict:
    """Stop unless every frozen input matches the hash in the config."""
    up = config["upstream"]
    outputs = resolve(up["error_decomposition_outputs"])
    checks = {
        "generalisation_config": resolve(up["generalisation_config"]),
        "benchmark_manifest": Path(benchmark_dir) / "benchmark_manifest.json",
        "predictions": Path(benchmark_dir) / "predictions.tsv.gz",
        "feature_manifest": Path(benchmark_dir) / "matrix" / "feature_manifest.json",
        "genome_manifest": resolve(up["genome_manifest"]),
        "error_decomposition_config": resolve(up["error_decomposition_config"]),
        "error_decomposition_sha256sums": outputs / "SHA256SUMS",
        "augmentation_plan": outputs / "augmentation_plan.json",
        "mechanism_catalogue": resolve(up["mechanism_catalogue"]),
    }
    got = {}
    for name, path in checks.items():
        if not path.is_file():
            raise UpstreamError(f"{name} missing: {path}")
        digest = pop.sha256(path)
        if digest != up[f"{name}_sha256"]:
            raise UpstreamError(f"{name} differs from the frozen input: {digest}")
        got[f"{name}_sha256"] = digest
    # every committed error_decomposition_v1 output still matches its SHA256SUMS
    for line in (outputs / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, name = line.split(maxsplit=1)
        if pop.sha256(outputs / name.strip()) != digest:
            raise UpstreamError(f"error_decomposition_v1 output {name} changed")
    plan = json.loads((outputs / "augmentation_plan.json").read_text(encoding="utf-8"))
    if pop.sha256(Path(decomposition_dir) / "augmented" / "augmentation_plan.json") != up["augmentation_plan_sha256"]:
        raise UpstreamError("the augmentation plan beside the reused R2 runs differs from the committed one")
    if len(plan["indicators"]) != config["mutation_panel"]["indicators"]:
        raise UpstreamError("mutation panel size differs")
    manifest = json.loads(checks["benchmark_manifest"].read_text(encoding="utf-8"))
    if manifest["tasks"] != config["benchmark"]["tasks"]:
        raise UpstreamError(f"benchmark has {manifest['tasks']} tasks, config expects {config['benchmark']['tasks']}")
    got.update({"benchmark_version": manifest["benchmark_version"], "tasks": manifest["tasks"],
                "mutation_panel_sha256": mutation_panel_sha(plan), "augmentation_drugs": plan["drugs"]})
    return got


def mutation_panel_sha(plan: dict) -> str:
    return gb.text_sha(sorted(plan["indicators"]))


# ---------------------------------------------------------------- sketches

def downsample(hashes: np.ndarray, scaled: int) -> np.ndarray:
    """FracMinHash downsample of a sorted sketch: the hashes at or below the coarser threshold."""
    return hashes[:int(np.searchsorted(hashes, np.uint64(gb.max_hash(scaled)), side="right"))]


def sketch_params(config: dict, scaled: int) -> dict:
    d = config["density"]
    return {"method": "FracMinHash k-mer presence", "tool": d["tool"], "tool_version": d["tool_version"],
            "ksize": d["ksize"], "scaled": scaled, "seed": d["hash_seed"]}


def sketch_one(cache: str, outdir: str, r0_dir: str, row: dict, expected_assembly: str | None, params: dict,
               reference_scaled: int) -> dict:
    """Verify and sketch one genome at the densest candidate; compare its downsample with the R0 sketch."""
    gca = row["assembly_accession"]
    meta_path, npy = Path(outdir) / f"{gca}.json", Path(outdir) / f"{gca}.npy"
    base = {"isolate_id": row["isolate_id"], "assembly_accession": gca, "genome_sha256": row["genome_sha256"]}
    hashes = None
    if meta_path.exists() and npy.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("params") == params and meta.get("genome_sha256") == row["genome_sha256"] and \
                (expected_assembly is None or expected_assembly == gca):
            hashes = np.load(npy)
            if hashes.shape[0] != meta.get("hash_count"):
                hashes = None
    action = "reused"
    if hashes is None:
        status, detail, data = gb.verify_genome(Path(cache), row, expected_assembly)
        if status != "verified":
            return {**base, "status": status, "detail": detail or NA, "hash_count": 0, "r0_reproduced": NA}
        try:
            hashes = gb.sketch_genome(data, params)
        except (ValueError, UnicodeDecodeError, OSError) as error:
            return {**base, "status": "genome_unreadable", "detail": str(error), "hash_count": 0,
                    "r0_reproduced": NA}
        tmp = npy.with_name(f".{npy.name}.{os.getpid()}.tmp.npy")
        np.save(tmp, hashes)
        tmp.replace(npy)
        pop.write_json_atomic(meta_path, {**base, "params": params, "hash_count": int(hashes.shape[0]),
                                          "created_at": now()})
        action = "sketched"
    r0_path = Path(r0_dir) / f"{gca}.npy"
    r0 = np.load(r0_path) if r0_path.exists() else None
    same = r0 is not None and np.array_equal(downsample(hashes, reference_scaled), r0)
    return {**base, "status": "genome_verified", "detail": action, "hash_count": int(hashes.shape[0]),
            "r0_reproduced": "yes" if same else "no"}


VERIFY_FIELDS = "isolate_id assembly_accession genome_sha256 status detail hash_count r0_reproduced".split()


def stage_sketch(config: dict, outdir: Path, benchmark_dir: Path, cache: Path, baseline: Path, workers: int,
                 limit: int | None = None) -> dict:
    p = paths(outdir)
    p["sketches"].mkdir(parents=True, exist_ok=True)
    up = config["upstream"]
    manifest = gb.read_genome_manifest(resolve(up["genome_manifest"]), up["genome_manifest_sha256"],
                                       up["genome_manifest_rows"])
    d = config["density"]
    params = sketch_params(config, d["sketch_scaled"])
    assemblies = gb.cohort_assemblies(baseline)
    rows = sorted(manifest.values(), key=lambda r: r["isolate_id"])[:limit]
    args = [(str(cache), str(p["sketches"]), str(Path(benchmark_dir) / "sketches"), r,
             assemblies.get(r["isolate_id"]), params, d["reference_scaled"]) for r in rows]
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, res in enumerate(pool.map(_sketch_star, args, chunksize=4), 1):
            results.append(res)
            if i % 500 == 0:
                print(f"{now()} sketched {i}/{len(rows)}", flush=True)
    results.sort(key=lambda r: r["isolate_id"])
    write_tsv(p["verification"], VERIFY_FIELDS, results)
    status = Counter(r["status"] for r in results)
    return {"genomes": len(results), "status": dict(status),
            "r0_reproduced": dict(Counter(r["r0_reproduced"] for r in results)),
            "actions": dict(Counter(r["detail"] for r in results if r["status"] == "genome_verified"))}


def _sketch_star(a):
    return sketch_one(*a)


def verified_rows(outdir: Path) -> list[tuple[str, str]]:
    rows = pop.read_tsv(paths(outdir)["verification"])
    bad = [r for r in rows if r["status"] != "genome_verified" or r["r0_reproduced"] != "yes"]
    if bad:
        raise UpstreamError(f"{len(bad)} genomes not verified or not reproducing R0, e.g. {bad[0]}")
    return [(r["isolate_id"], r["assembly_accession"]) for r in rows]



# ---------------------------------------------------------------- pilot: feature space per candidate

def npy_slices(path: Path, edges: np.ndarray, top: int) -> tuple[int, np.ndarray]:
    """Data offset of a sorted uint64 .npy file, the positions of the bucket `edges` in it and, last, the end of
    the hashes at or below `top`."""
    with open(path, "rb") as handle:
        major, _ = np.lib.format.read_magic(handle)
        read = np.lib.format.read_array_header_1_0 if major == 1 else np.lib.format.read_array_header_2_0
        shape, _, dtype = read(handle)
        offset = handle.tell()
    if dtype != np.uint64 or len(shape) != 1:
        raise ValueError(f"{path}: not a 1-d uint64 sketch")
    h = np.load(path)
    return offset, np.r_[np.searchsorted(h, edges, side="left"),
                         np.searchsorted(h, np.uint64(top), side="right")].astype(np.int64)


def candidate_feature_spaces(sketch_dir: Path, isolates: list[tuple[str, str]], scaled_values: list[int],
                             min_present: int, min_absent: int, buckets: int = 256) -> dict:
    """One pass over the densest sketches giving, for every candidate, what bin/generalisation_benchmark.py
    build_matrix would build from its downsampled sketches: observed and kept hashes, merged patterns and
    non-zeros, plus the kept hashes and their pattern ids.

    Hashes are visited in ascending order and patterns numbered at first sight, so the kept hashes of a
    candidate are a prefix of the densest candidate's and its pattern ids are those a standalone build
    assigns."""
    finest = min(scaled_values)
    top = gb.max_hash(finest)
    limits = {s: gb.max_hash(s) for s in scaled_values}
    edges = np.array([i * (top + 1) // buckets for i in range(buckets)], dtype=np.uint64)
    n = len(isolates)
    files = [Path(sketch_dir) / f"{gca}.npy" for _, gca in isolates]
    offsets, pos = [], np.empty((n, buckets + 1), dtype=np.int64)
    for r, f in enumerate(files):
        off, p = npy_slices(f, edges, top)
        offsets.append(off)
        pos[r] = p
    signatures, sizes, kept_hash, kept_pattern = {}, [], [], []
    observed = Counter()
    for b in range(buckets):
        parts, rows = [], []
        for r, f in enumerate(files):
            start, stop = pos[r, b], pos[r, b + 1]
            if stop > start:
                parts.append(np.fromfile(f, dtype=np.uint64, count=int(stop - start),
                                         offset=offsets[r] + 8 * int(start)))
                rows.append(np.full(int(stop - start), r, dtype=np.int32))
        if not parts:
            continue
        hv, rv = np.concatenate(parts), np.concatenate(rows)
        del parts, rows
        order = np.argsort(hv, kind="stable")
        hv, rv = hv[order], rv[order]
        del order
        starts = np.flatnonzero(np.r_[True, hv[1:] != hv[:-1]])
        counts = np.diff(np.r_[starts, hv.shape[0]])
        uniq = hv[starts]
        for s_, lim in limits.items():
            observed[s_] += int(np.count_nonzero(uniq <= np.uint64(lim)))
        keep = (counts >= min_present) & (n - counts >= min_absent)
        for s, c, h in zip(starts[keep], counts[keep], uniq[keep]):
            members = rv[s:s + c]
            key = hashlib.blake2b(members.tobytes(), digest_size=16).digest()
            pid = signatures.get(key)
            if pid is None:
                pid = signatures[key] = len(sizes)
                sizes.append(int(c))
            kept_hash.append(h)
            kept_pattern.append(pid)
    hashes = np.array(kept_hash, dtype=np.uint64)
    pattern = np.array(kept_pattern, dtype=np.int64)
    cum = np.r_[0, np.cumsum(np.array(sizes, dtype=np.int64))]
    out = {}
    for s_, lim in sorted(limits.items(), reverse=True):
        k = int(np.searchsorted(hashes, np.uint64(lim), side="right"))
        pat = pattern[:k]
        npat = int(pat.max()) + 1 if k else 0
        if k and npat != len(np.unique(pat)):
            raise AssertionError("pattern ids of a candidate are not a prefix")
        out[s_] = {"scaled": s_, "observed_features": observed[s_], "kept_features": k, "patterns": npat,
                   "nnz": int(cum[npat]), "hashes": hashes[:k], "hash_pattern": pat}
    return out



# ---------------------------------------------------------------- pilot: capture per candidate

def audit_isolate(row: dict, calls: list[dict], cache: Path, refs: dict, spaces: dict, cfg: dict) -> dict:
    """error_decomposition_v1 capture of every called gene and point mutation, under each feature space."""
    status, detail, data = gb.verify_genome(Path(cache), row)
    out = {"isolate_id": row["isolate_id"], "genome_status": status, "detail": detail, "genes": [],
           "mutations": []}
    if status != "verified":
        return out
    contigs = ed.read_fasta_records(data)
    for c in calls:
        if c["element_subtype"] == "AMR":
            res = {s: ed.gene_capture(contigs, c, sp)["status"] for s, sp in spaces.items()}
            out["genes"].append({"determinant": c["element_symbol"], "status": res})
        elif c["element_subtype"] == "POINT":
            res = {s: ed.mutation_capture(contigs, c, refs, sp, cfg)["status"] for s, sp in spaces.items()}
            out["mutations"].append({"mutation": c["element_symbol"], "method": c["method"], "status": res})
    return out


_AUDIT = {}


def _audit_init(cache, db, taxgroup, loci, space_files, params, cfg):
    refs = ed.reference_sequences(Path(db), taxgroup, loci)
    spaces = {int(s): ed.FeatureSpace(np.load(h), np.load(p), int(s), params["ksize"], params["seed"])
              for s, (h, p) in space_files.items()}
    _AUDIT.update(cache=Path(cache), refs=refs, spaces=spaces, cfg=cfg)


def _audit_one(args):
    row, calls = args
    return audit_isolate(row, calls, _AUDIT["cache"], _AUDIT["refs"], _AUDIT["spaces"], _AUDIT["cfg"])


def run_capture_audit(config: dict, ed_config: dict, decomposition_dir: Path, cache: Path, space_files: dict,
                      ckpt: Path, workers: int) -> dict:
    up = config["upstream"]
    manifest = gb.read_genome_manifest(resolve(up["genome_manifest"]), up["genome_manifest_sha256"])
    calls = defaultdict(list)
    for c in ed.read_tsv(ed.paths(decomposition_dir)["calls"]):
        calls[c["isolate_id"]].append(c)
    done = {}
    if ckpt.exists():
        for line in ckpt.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            done[rec["isolate_id"]] = rec
    eup = ed_config["upstream"]
    db = resolve(eup["amrfinder_database"])
    for name, expected in eup["amrfinder_database_files_sha256"].items():
        if pop.sha256(db / name) != expected:
            raise UpstreamError(f"AMRFinderPlus database file {name} differs")
    loci = ed.mutation_loci(db, eup["amrfinder_taxgroup"])
    params = sketch_params(config, config["density"]["reference_scaled"])
    todo = [(manifest[i], calls.get(i, [])) for i in sorted(manifest) if i not in done]
    init = (str(cache), str(db), eup["amrfinder_taxgroup"], loci,
            {str(s): (str(h), str(p)) for s, (h, p) in space_files.items()}, params, ed_config["capture"])
    with ckpt.open("a", encoding="utf-8") as handle:
        if workers <= 1:
            _audit_init(*init)
            results = map(_audit_one, todo)
            pool = None
        else:
            pool = ProcessPoolExecutor(max_workers=workers, initializer=_audit_init, initargs=init)
            results = pool.map(_audit_one, todo, chunksize=8)
        for n, rec in enumerate(results, 1):
            handle.write(json.dumps(rec, sort_keys=True) + "\n")
            handle.flush()
            done[rec["isolate_id"]] = rec
            if n % 500 == 0:
                print(f"{now()} audited {len(done)}/{len(manifest)}", flush=True)
        if pool is not None:
            pool.shutdown()
    return done


def capture_tables(audit: dict, endpoints: list[dict], catalogue_by_drug: dict, groups: dict, label) -> dict:
    """Per-call and per-drug capture under each representation. `endpoints` are the common-cohort
    endpoints (antimicrobial, isolate_id and their catalogue genes and mutations); `label` maps a scaled
    value to its representation name."""
    order = ed.CAPTURE_ORDER
    by_mut, by_gene = defaultdict(Counter), defaultdict(Counter)
    carriers = defaultdict(set)
    gcap, mcap = {}, {}
    for rec in audit.values():
        for g in rec["genes"]:
            for sc, st in g["status"].items():
                by_gene[(label(sc), g["determinant"])][st] += 1
                key = (rec["isolate_id"], g["determinant"], label(sc))
                gcap[key] = ed.best_capture([st, gcap.get(key)])  # several copies: the best, as in PR #15
        for m in rec["mutations"]:
            for sc, st in m["status"].items():
                by_mut[(label(sc), m["mutation"])][st] += 1
                key = (rec["isolate_id"], m["mutation"], label(sc))
                mcap[key] = ed.best_capture([st, mcap.get(key)])
            carriers[m["mutation"]].add(rec["isolate_id"])
    reps = sorted({k[0] for k in by_mut} | {k[0] for k in by_gene}, key=rep_order)
    mut_rows = [{"representation": r, "mutation": m, "calls": sum(c.values()), "carriers": len(carriers[m]),
                 **{st: c[st] for st in order}, "captured_fraction": ed.frac(c["captured"], sum(c.values()))}
                for (r, m), c in sorted(by_mut.items(), key=lambda kv: (kv[0][1], rep_order(kv[0][0])))]
    gene_rows = [{"representation": r, "determinant": d, "calls": sum(c.values()), **{st: c[st] for st in order},
                  "captured_fraction": ed.frac(c["captured"], sum(c.values()))}
                 for (r, d), c in sorted(by_gene.items(), key=lambda kv: (kv[0][1], rep_order(kv[0][0])))]
    agg = defaultdict(Counter)
    for e in endpoints:
        for kind, field in (("gene", "gene_determinants"), ("point_mutation", "catalogue_mutations")):
            if e[field] == NA:
                continue
            for d in e[field].split(";"):
                for r in reps:
                    st = (gcap if kind == "gene" else mcap).get((e["isolate_id"], d, r)) or "representation_unknown"
                    agg[(e["antimicrobial"], kind, d, r)][st] += 1
                    agg[(e["antimicrobial"], kind, "ALL", r)][st] += 1
    drug_rows = []
    for (drug, kind, det, r), c in sorted(agg.items(), key=lambda kv: (*kv[0][:3], rep_order(kv[0][3]))):
        n = sum(c.values())
        cls = catalogue_by_drug.get(drug, {}).get(det, {}).get("mechanism_class", NA) if det != "ALL" else "ALL"
        drug_rows.append({"antimicrobial": drug, "drug_group": groups.get(drug, NA), "determinant_kind": kind,
                          "determinant": det, "mechanism_class": cls, "representation": r, "endpoints": n,
                          **{st: c[st] for st in order}, "captured_fraction": ed.frac(c["captured"], n)})
    return {"mutations": mut_rows, "genes": gene_rows, "drugs": drug_rows, "gcap": gcap, "mcap": mcap}


def rep_order(name: str) -> tuple:
    if name.startswith("scaled_"):
        return (1, -int(name.split("_")[1]))
    return (0, REPRESENTATIONS.index(name)) if name in REPRESENTATIONS else (2, name)


def common_endpoints(decomposition_dir: Path) -> list[dict]:
    """Each common-cohort endpoint once (random strategy), with only its determinant columns: no phenotype,
    prediction or error field is kept."""
    keep = ("antimicrobial", "isolate_id", "gene_determinants", "catalogue_mutations")
    out = []
    with gzip.open(ed.paths(decomposition_dir)["context"], "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        fields = next(reader)
        idx = {f: fields.index(f) for f in ("cohort", "split_strategy", *keep)}
        for v in reader:
            if v[idx["cohort"]] == "common" and v[idx["split_strategy"]] == "random":
                out.append({k: sys.intern(v[idx[k]]) for k in keep})
    return out


RESOURCE_FIELDS = ("representation scaled ksize seed rows observed_features kept_features patterns nnz "
                   "matrix_density matrix_bytes sketch_hashes sketch_bytes build_seconds build_max_rss_bytes "
                   "projected_task_seconds_total projected_task_peak_rss_bytes projection_basis").split()


def built_feature_spaces(matrices: dict, isolates: list[tuple[str, str]]) -> dict:
    """Feature spaces read back from matrices built by stage_matrix, after checking they cover the verified
    genomes in order and match their own manifest."""
    want = [iid for iid, _ in isolates]
    out = {}
    for sc, m in matrices.items():
        meta = json.loads((m / "feature_manifest.json").read_text(encoding="utf-8"))
        rows = [line.split("\t")[0] for line in (m / "rows.tsv").read_text(encoding="utf-8").splitlines()[1:]]
        if rows != want or meta["params"]["scaled"] != sc:
            raise UpstreamError(f"{m}: built for other genomes or another density")
        hashes, pattern = np.load(m / "feature_hashes.npy"), np.load(m / "feature_pattern.npy")
        if (sha256_bytes(hashes.tobytes()) != meta["feature_hashes_sha256"]
                or sha256_bytes(pattern.tobytes()) != meta["feature_pattern_sha256"]):
            raise UpstreamError(f"{m}: feature files differ from their manifest")
        out[sc] = {"hashes": hashes, "hash_pattern": pattern,
                   **{k: meta[k] for k in ("observed_features", "kept_features", "patterns", "nnz")},
                   "build_seconds": meta.get("build_seconds"), "build_max_rss_bytes": meta.get("build_max_rss_bytes")}
    return out


def stage_pilot(config: dict, outdir: Path, benchmark_dir: Path, decomposition_dir: Path, cache: Path,
                workers: int) -> dict:
    p = paths(outdir)
    p["pilot"].mkdir(parents=True, exist_ok=True)
    isolates = verified_rows(outdir)
    d = config["density"]
    gconfig = json.loads(resolve(config["upstream"]["generalisation_config"]).read_text(encoding="utf-8"))
    f = gconfig["features"]
    scaled_values = sorted({d["reference_scaled"], *d["candidates"]}, reverse=True)
    label = lambda sc: "R0" if int(sc) == d["reference_scaled"] else f"scaled_{int(sc)}"  # noqa: E731
    built = {sc: p["root"] / f"pilot_scaled_{sc}" / "matrix" for sc in scaled_values}
    started = time.time()
    if all((m / "feature_manifest.json").exists() for m in built.values()):
        # the fit pilot already built every candidate matrix (on the cluster); its feature spaces are the
        # ones candidate_feature_spaces gives, and reading them avoids one random-access pass per bucket
        spaces, source = built_feature_spaces(built, isolates), "pilot_scaled_*/matrix"
        build_seconds = {sc: spaces[sc]["build_seconds"] for sc in scaled_values}
        build_rss = {sc: spaces[sc]["build_max_rss_bytes"] for sc in scaled_values}
    else:
        spaces = candidate_feature_spaces(p["sketches"], isolates, scaled_values, f["min_train_present"],
                                          f["min_train_absent"])
        source = "candidate_feature_spaces"
        build_seconds = dict.fromkeys(scaled_values, round(time.time() - started, 1))
        build_rss = dict.fromkeys(scaled_values, max_rss())
    # R0 must be the PR #14 feature space, hash for hash and pattern for pattern
    r0 = spaces[d["reference_scaled"]]
    fm = json.loads((Path(benchmark_dir) / "matrix" / "feature_manifest.json").read_text(encoding="utf-8"))
    r0_check = {
        "feature_hashes": sha256_bytes(r0["hashes"].tobytes()) == fm["feature_hashes_sha256"],
        "feature_pattern": sha256_bytes(r0["hash_pattern"].tobytes()) == fm["feature_pattern_sha256"],
        "observed_features": r0["observed_features"] == fm["observed_features"],
        "patterns": r0["patterns"] == fm["patterns"], "nnz": r0["nnz"] == fm["nnz"]}
    if not all(r0_check.values()):
        raise UpstreamError(f"downsampled sketches do not rebuild the PR #14 feature space: {r0_check}")
    space_files = {}
    for sc, sp in spaces.items():
        h, pt = p["pilot"] / f"scaled_{sc}_feature_hashes.npy", p["pilot"] / f"scaled_{sc}_feature_pattern.npy"
        np.save(h, sp["hashes"])
        np.save(pt, sp["hash_pattern"])
        space_files[sc] = (h, pt)
    sketch_hashes = {sc: 0 for sc in scaled_values}
    for _, gca in isolates:
        h = np.load(p["sketches"] / f"{gca}.npy", mmap_mode="r")
        for sc in scaled_values:
            sketch_hashes[sc] += int(np.searchsorted(h, np.uint64(gb.max_hash(sc)), side="right"))
    # resources: generalisation_v1 measured R0; denser candidates are projected by their non-zeros
    run = json.loads(resolve("benchmark/cohorts/ecoli_research_v1/generalisation_v1/run_record.json")
                     .read_text(encoding="utf-8"))
    r0_secs, r0_rss = run["stages"]["run"]["task_seconds_total"], run["stages"]["run"]["max_rss_bytes"]
    rows = []
    for sc in scaled_values:
        sp = spaces[sc]
        ratio = sp["nnz"] / r0["nnz"]
        rows.append({"representation": label(sc), "scaled": sc, "ksize": d["ksize"], "seed": d["hash_seed"],
                     "rows": len(isolates), "observed_features": sp["observed_features"],
                     "kept_features": sp["kept_features"], "patterns": sp["patterns"], "nnz": sp["nnz"],
                     "matrix_density": sp["nnz"] / (len(isolates) * sp["patterns"]),
                     "matrix_bytes": 5 * sp["nnz"] + 8 * (len(isolates) + 1),
                     "sketch_hashes": sketch_hashes[sc], "sketch_bytes": 8 * sketch_hashes[sc],
                     "build_seconds": build_seconds[sc], "build_max_rss_bytes": build_rss[sc],
                     "projected_task_seconds_total": round(r0_secs * ratio),
                     "projected_task_peak_rss_bytes": round(r0_rss * ratio),
                     "projection_basis": "generalisation_v1 measured R0 run (one worker, 16 GB laptop) scaled by "
                                         "non-zeros; replaced by measured fits where available"})
    write_tsv(p["pilot"] / "representation_resource_pilot.tsv", list(RESOURCE_FIELDS), rows)
    # capture of every called determinant under every candidate
    audit = run_capture_audit(config, ed.load_config(resolve(config["upstream"]["error_decomposition_config"])),
                              decomposition_dir, cache, space_files, p["pilot"] / "capture_checkpoint.jsonl",
                              workers)
    statuses = Counter(r["genome_status"] for r in audit.values())
    if set(statuses) != {"verified"}:
        raise UpstreamError(f"capture audit genome statuses {dict(statuses)}")
    ed_manifest = json.loads((resolve(config["upstream"]["error_decomposition_outputs"]) /
                              "error_decomposition_manifest.json").read_text(encoding="utf-8"))
    groups = {k: v["group"] for k, v in ed_manifest["drug_groups"].items()}
    catalogue = ed.catalogue_index(ed.read_tsv(resolve(config["upstream"]["mechanism_catalogue"])))
    tables = capture_tables(audit, common_endpoints(decomposition_dir), catalogue, groups, label)
    mismatch = r0_capture_mismatches(audit, d["reference_scaled"], decomposition_dir)
    order = list(ed.CAPTURE_ORDER)
    write_tsv(p["pilot"] / "mutation_capture_by_candidate.tsv",
              ["representation", "mutation", "calls", "carriers", *order, "captured_fraction"], tables["mutations"])
    write_tsv(p["pilot"] / "determinant_capture_by_candidate.tsv",
              ["representation", "determinant", "calls", *order, "captured_fraction"], tables["genes"])
    write_tsv(p["pilot"] / "drug_capture_by_candidate.tsv",
              ["antimicrobial", "drug_group", "determinant_kind", "determinant", "mechanism_class", "representation",
               "endpoints", *order, "captured_fraction"], tables["drugs"])
    summary = {"rows": len(isolates), "feature_space_source": source, "r0_feature_space_check": r0_check,
               "r0_capture_mismatches": mismatch,
               "resources": rows, "written_at": now()}
    gb.write_json(p["pilot"] / "pilot_summary.json", summary)
    return {"r0_check": r0_check, "r0_capture_mismatches": mismatch,
            "candidates": {r["representation"]: {k: r[k] for k in ("kept_features", "patterns", "nnz")} for r in rows}}


def r0_capture_mismatches(audit: dict, reference_scaled: int, decomposition_dir: Path) -> dict:
    """(isolate, determinant) pairs whose R0 call statuses differ from error_decomposition_v1. An isolate can
    carry several copies of one determinant, so the statuses of all its calls are compared as a multiset."""
    p = ed.paths(decomposition_dir)
    sc = str(reference_scaled)
    new, old = defaultdict(list), defaultdict(list)
    for rec in audit.values():
        for kind, field, key in (("gene", "genes", "determinant"), ("mutation", "mutations", "mutation")):
            for c in rec[field]:
                status = {str(k): v for k, v in c["status"].items()}  # int keys in memory, str from the checkpoint
                new[(kind, rec["isolate_id"], c[key])].append(status[sc])
    for kind, path, key in (("gene", p["gene_capture"], "determinant"), ("mutation", p["mut_capture"], "mutation")):
        for r in ed.read_tsv(path):
            old[(kind, r["isolate_id"], r[key])].append(r["status"])
    bad = Counter(k[0] for k in set(new) | set(old) if sorted(new.get(k, [])) != sorted(old.get(k, [])))
    return dict(bad)


# ---------------------------------------------------------------- feature matrix

def build_matrix_streaming(sketch_dir: Path, isolates: list[tuple[str, str]], scaled: int, min_present: int,
                           min_absent: int, spill: Path, buckets: int = 256) -> dict:
    """The matrix bin/generalisation_benchmark.py build_matrix builds from sketches downsampled to `scaled`,
    without holding the sketches in memory: pattern columns are spilled to disk in pattern order and turned
    into CSR rows in chunks."""
    top = gb.max_hash(scaled)
    edges = np.array([i * (top + 1) // buckets for i in range(buckets)], dtype=np.uint64)
    n = len(isolates)
    files = [Path(sketch_dir) / f"{gca}.npy" for _, gca in isolates]
    offsets, pos = [], np.empty((n, buckets + 1), dtype=np.int64)
    for r, f in enumerate(files):
        off, p = npy_slices(f, edges, top)
        offsets.append(off)
        pos[r] = p
    signatures, sizes, kept_hash, kept_pattern = {}, [], [], []
    observed = 0
    spill = Path(spill)
    with spill.open("wb") as out:
        for b in range(buckets):
            parts, rows = [], []
            for r, f in enumerate(files):
                start, stop = pos[r, b], pos[r, b + 1]
                if stop > start:
                    parts.append(np.fromfile(f, dtype=np.uint64, count=int(stop - start),
                                             offset=offsets[r] + 8 * int(start)))
                    rows.append(np.full(int(stop - start), r, dtype=np.int32))
            if not parts:
                continue
            hv, rv = np.concatenate(parts), np.concatenate(rows)
            del parts, rows
            order = np.argsort(hv, kind="stable")
            hv, rv = hv[order], rv[order]
            del order
            starts = np.flatnonzero(np.r_[True, hv[1:] != hv[:-1]])
            counts = np.diff(np.r_[starts, hv.shape[0]])
            observed += starts.shape[0]
            keep = (counts >= min_present) & (n - counts >= min_absent)
            for s, c in zip(starts[keep], counts[keep]):
                members = rv[s:s + c]
                key = hashlib.blake2b(members.tobytes(), digest_size=16).digest()
                pid = signatures.get(key)
                if pid is None:
                    pid = signatures[key] = len(sizes)
                    sizes.append(int(c))
                    out.write(members.tobytes())
                kept_hash.append(hv[s])
                kept_pattern.append(pid)
    del signatures
    p = len(sizes)
    if p == 0:
        raise ValueError("no informative features")
    sizes = np.array(sizes, dtype=np.int64)
    col_ptr = np.r_[0, np.cumsum(sizes)]
    nnz = int(col_ptr[-1])
    csc_indices = np.memmap(spill, dtype=np.int32, mode="r", shape=(nnz,))
    row_counts = np.zeros(n, dtype=np.int64)
    step = 1 << 26
    for a in range(0, nnz, step):
        row_counts += np.bincount(csc_indices[a:a + step], minlength=n)
    indptr = np.r_[0, np.cumsum(row_counts)].astype(gb.index_dtype(nnz))
    indices = np.empty(nnz, dtype=np.int32)
    fill = indptr[:-1].astype(np.int64).copy()
    # columns in ascending order, chunk by chunk, so every row's indices come out sorted
    col_step = max(1, int(np.searchsorted(col_ptr, step)))
    for c0 in range(0, p, col_step):
        c1 = min(p, c0 + col_step)
        a, z = int(col_ptr[c0]), int(col_ptr[c1])
        r = np.asarray(csc_indices[a:z])
        c = np.repeat(np.arange(c0, c1, dtype=np.int32), sizes[c0:c1])
        order = np.argsort(r, kind="stable")
        r, c = r[order], c[order]
        first = np.r_[0, np.flatnonzero(r[1:] != r[:-1]) + 1]
        group = np.repeat(np.arange(first.shape[0]), np.diff(np.r_[first, r.shape[0]]))
        rank = np.arange(r.shape[0]) - first[group]
        indices[fill[r] + rank] = c
        np.add.at(fill, r[first], np.diff(np.r_[first, r.shape[0]]))
    del csc_indices
    spill.unlink()
    from scipy import sparse
    csr = sparse.csr_matrix((np.ones(nnz, dtype=np.int8), indices, indptr), shape=(n, p))
    return {"matrix": csr, "hashes": np.array(kept_hash, dtype=np.uint64),
            "hash_pattern": np.array(kept_pattern, dtype=np.int64), "observed_features": int(observed),
            "kept_features": len(kept_hash), "patterns": p, "nnz": nnz}


def stage_matrix(config: dict, outdir: Path, scaled: int, name: str) -> dict:
    p = paths(outdir)
    isolates = verified_rows(outdir)
    gconfig = json.loads(resolve(config["upstream"]["generalisation_config"]).read_text(encoding="utf-8"))
    f = gconfig["features"]
    target = p["root"] / name / "matrix"
    target.mkdir(parents=True, exist_ok=True)
    started = time.time()
    built = build_matrix_streaming(p["sketches"], isolates, scaled, f["min_train_present"], f["min_train_absent"],
                                   target / ".csc_spill.bin")
    meta = {"params": sketch_params(config, scaled), "representation": name,
            "genome_manifest_sha256": config["upstream"]["genome_manifest_sha256"],
            "genome_verification_sha256": pop.sha256(p["verification"]),
            "global_prune": {"min_present": f["min_train_present"], "min_absent": f["min_train_absent"],
                             "note": "drops only features every per-task training filter would drop"},
            "deduplicated": True, "builder": f"bin/representation_sensitivity.py {VERSION}",
            "build_seconds": round(time.time() - started, 1), "build_max_rss_bytes": max_rss()}
    manifest = gb.save_matrix(target, built, isolates, meta)
    return {k: manifest[k] for k in ("rows", "observed_features", "kept_features", "patterns", "nnz",
                                     "build_seconds", "build_max_rss_bytes")}


# ---------------------------------------------------------------- frozen tasks and their identity

def load_benchmark(config: dict, benchmark_dir: Path, population: Path, baseline: Path):
    """Rebuild the generalisation_v1 tasks from the frozen inputs; stops unless every split hash matches."""
    gpath = resolve(config["upstream"]["generalisation_config"])
    inputs = gb.load_inputs(gpath, population, baseline, Path(benchmark_dir) / "genome_verification.tsv")
    manifest = json.loads((Path(benchmark_dir) / "benchmark_manifest.json").read_text(encoding="utf-8"))
    tasks, skipped = gb.build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
    gb.check_frozen(tasks, manifest, gpath)
    return inputs, tasks, skipped, manifest


def r2_task_set(tasks: list[dict], plan: dict) -> list[dict]:
    """The tasks error_decomposition_v1 refitted with indicators: common cohort, eligible drugs."""
    return [t for t in tasks if t["cohort"] == "common" and t["antimicrobial"] in plan["drugs"]]


IDENTITY_FIELDS = ("task_id cohort antimicrobial split_strategy fold_or_holdout_id n_train n_test train_R test_R "
                   "split_sha256 r0_split_sha256 r0_match r2_split_sha256 r2_match permuted_split_sha256 "
                   "permuted_match status").split()


def task_identity(tasks: list[dict], manifest: dict, benchmark_dir: Path, decomposition_dir: Path,
                  plan: dict) -> list[dict]:
    """Every rebuilt task against the manifest, the R0 task record and the R2 task records."""
    expected = manifest["task_split_sha256"]
    got = {t["task_id"]: t for t in tasks}
    if set(got) != set(expected):
        missing, extra = sorted(set(expected) - set(got)), sorted(set(got) - set(expected))
        raise IdentityError(f"task set differs from the manifest: missing {missing[:3]}, unexpected {extra[:3]}")
    r2_ids = {t["task_id"] for t in r2_task_set(tasks, plan)}
    aug = Path(decomposition_dir) / "augmented"
    rows = []
    for tid in sorted(got):
        t = got[tid]
        digest = gb.task_digest(t)
        r0_path = Path(benchmark_dir) / "tasks" / f"{tid}.json"
        r0 = json.loads(r0_path.read_text(encoding="utf-8"))
        r0_test = {r["isolate_id"]: r["phenotype"] for r in ed.read_tsv(r0_path.with_suffix(".tsv"))}
        phenotypes_match = r0_test == {e["isolate_id"]: e["phenotype"] for e in t["test"]} and \
            r0.get("train_R") == sum(e["phenotype"] == "R" for e in t["train"])
        rec = {"task_id": tid, "cohort": t["cohort"], "antimicrobial": t["antimicrobial"],
               "split_strategy": t["strategy"], "fold_or_holdout_id": t["fold"], "n_train": len(t["train"]),
               "n_test": len(t["test"]), "train_R": sum(e["phenotype"] == "R" for e in t["train"]),
               "test_R": sum(e["phenotype"] == "R" for e in t["test"]), "split_sha256": digest,
               "r0_split_sha256": r0.get("split_sha256"), "r0_match": "yes" if r0.get("split_sha256") == digest
               and digest == expected[tid] and r0.get("n_train") == len(t["train"])
               and r0.get("n_test") == len(t["test"]) and phenotypes_match else "no"}
        for variant, col in (("mutation_augmented", "r2"), ("permuted_control", "permuted")):
            path = aug / variant / "tasks" / f"{tid}.json"
            if tid in r2_ids:
                st = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
                rec[f"{col}_split_sha256"] = st.get("split_sha256", NA)
                rec[f"{col}_match"] = "yes" if st.get("split_sha256") == digest and st.get("status") == "success" \
                    else "no"
            else:
                rec[f"{col}_split_sha256"], rec[f"{col}_match"] = NA, "not_applicable"
        rec["status"] = "identical" if rec["r0_match"] == "yes" and "no" not in (rec["r2_match"],
                                                                                   rec["permuted_match"]) \
            else "mismatch"
        rows.append(rec)
    bad = [r["task_id"] for r in rows if r["status"] != "identical"]
    if bad:
        raise IdentityError(f"{len(bad)} tasks differ from their R0 or R2 records, e.g. {bad[:3]}")
    return rows


# ---------------------------------------------------------------- runs

def indicator_matrix(rep_matrix: Path, decomposition_dir: Path, plan: dict, control_seed: int, permuted: bool,
                     target: Path) -> dict:
    """R1 plus the error_decomposition_v1 indicator block, saved as a matrix directory."""
    X, row_of, meta = gb.load_matrix(rep_matrix)
    isolates = sorted(row_of, key=row_of.get)
    calls = defaultdict(list)
    for c in ed.read_tsv(ed.paths(decomposition_dir)["calls"]):
        calls[c["isolate_id"]].append(c)
    block = ed.indicator_block(isolates, calls, plan["indicators"])
    prevalence = {s: int(block[:, j].sum()) for j, s in enumerate(plan["indicators"])}
    if prevalence != plan["indicator_prevalence"]:
        raise UpstreamError("indicator block differs from the error_decomposition_v1 panel")
    if permuted:
        block = ed.permute_rows(block, control_seed)
    Xa = ed.append_columns(X.tocsr(), block)
    del X
    target.mkdir(parents=True, exist_ok=True)
    np.save(target / "indptr.npy", Xa.indptr)
    np.save(target / "indices.npy", Xa.indices.astype(np.int32))
    (target / "rows.tsv").write_text((rep_matrix / "rows.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    m = {**{k: meta[k] for k in ("params", "genome_manifest_sha256", "genome_verification_sha256")},
         "base_matrix_indices_sha256": meta["indices_sha256"], "indicators": len(plan["indicators"]),
         "indicator_rows_permuted": permuted, "control_seed": control_seed if permuted else None,
         "rows": Xa.shape[0], "patterns": Xa.shape[1], "nnz": int(Xa.nnz),
         "indices_sha256": sha256_bytes(Xa.indices.astype(np.int32).tobytes()),
         "indptr_sha256": sha256_bytes(Xa.indptr.tobytes()), "created_at": now()}
    gb.write_json(target / "feature_manifest.json", m)
    return m


def _packed_worker(args):
    status, preds = gb._worker(args)
    return {**status, "max_rss_bytes": max_rss()}, preds


def run_tasks_packed(tasks: list[dict], matrix_dir: Path, outdir: Path, config: dict, workers: int,
                     memory_cap_bytes: float, bytes_per_train: float, only: str | None = None) -> Counter:
    """gb.run_tasks with the same checkpoints, but a task starts only while the memory estimated for the
    running tasks (bytes_per_train times training size, from the fit pilot) stays under the cap. Largest
    tasks first; one task per child process, so each checkpoint records that task's peak memory."""
    taskdir = Path(outdir) / "tasks"
    taskdir.mkdir(parents=True, exist_ok=True)
    counts, todo = Counter(), []
    for t in tasks:
        if only and not re.search(only, t["task_id"]):
            continue
        st = taskdir / f"{t['task_id']}.json"
        if st.exists():
            prev = json.loads(st.read_text(encoding="utf-8"))
            if prev.get("status") == "success" and prev.get("split_sha256") == gb.task_digest(t):
                counts["reused"] += 1
                continue
        todo.append(t)
    todo.sort(key=lambda t: (-len(t["train"]), t["task_id"]))
    need = lambda t: bytes_per_train * len(t["train"])  # noqa: E731
    slim = lambda t: {**t, "train": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")}  # noqa: E731
                                     for e in t["train"]],
                      "test": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["test"]]}
    running = {}
    with ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=1) as pool:
        while todo or running:
            used = sum(running.values())
            for t in list(todo):
                if len(running) >= workers:
                    break
                if used + need(t) <= memory_cap_bytes or not running:
                    todo.remove(t)
                    running[pool.submit(_packed_worker, (slim(t), str(matrix_dir), config))] = need(t)
                    used += need(t)
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for future in done:
                del running[future]
                status, preds = future.result()
                if status["status"] == "success":
                    write_tsv(taskdir / f"{status['task_id']}.tsv", gb.PRED_FIELDS, preds)
                gb.write_json(taskdir / f"{status['task_id']}.json", {**status, "finished_at": now()})
                counts[status["status"]] += 1
                print(f"{now()} {status['task_id']} {status['status']} {status.get('seconds', '')}s "
                      f"{status['max_rss_bytes'] / 2 ** 30:.1f}GB", flush=True)
    return counts


def packing_rate(config: dict, outdir: Path, scaled: int) -> float:
    """Peak bytes per training isolate at `scaled`: the largest ratio in the fit pilot, plus 10%."""
    fit = json.loads((paths(outdir)["pilot"] / f"fit_pilot_scaled_{scaled}.json").read_text(encoding="utf-8"))
    return 1.1 * max(f["peak_rss_bytes"] / f["n_train"] for f in fit["fits"])


def stage_run(config: dict, outdir: Path, benchmark_dir: Path, decomposition_dir: Path, population: Path,
              baseline: Path, rep: str, workers: int, only: str | None, memory_gb: float | None = None) -> dict:
    p = paths(outdir)
    selection = json.loads(p["selection"].read_text(encoding="utf-8"))
    inputs, tasks, _, manifest = load_benchmark(config, benchmark_dir, population, baseline)
    plan = json.loads((resolve(config["upstream"]["error_decomposition_outputs"]) /
                       "augmentation_plan.json").read_text(encoding="utf-8"))
    identity = task_identity(tasks, manifest, benchmark_dir, decomposition_dir, plan)
    write_tsv(p["root"] / "task_identity_validation.tsv", list(IDENTITY_FIELDS), identity)
    r1_matrix = p["R1"] / "matrix"
    fm = json.loads((r1_matrix / "feature_manifest.json").read_text(encoding="utf-8"))
    if fm["params"]["scaled"] != selection["R1_scaled"]:
        raise ValueError("R1 matrix was built at a scaled value other than the selected one")
    if rep == "R1":
        todo, matrix_dir = tasks, r1_matrix
    else:
        ctrl = json.loads(resolve(config["upstream"]["error_decomposition_config"]).read_text(encoding="utf-8"))
        matrix_dir = p[rep] / "matrix"
        if not (matrix_dir / "feature_manifest.json").exists():
            indicator_matrix(r1_matrix, decomposition_dir, plan, ctrl["augmentation"]["control_seed"],
                             rep == "R3_permuted", matrix_dir)
        todo = r2_task_set(tasks, plan)
    if memory_gb:
        # a quarter of the job's memory stays free, as the selection rule requires of the largest task
        counts = run_tasks_packed(todo, matrix_dir, p[rep], inputs["config"], workers, 0.75 * memory_gb * 2 ** 30,
                                  packing_rate(config, outdir, selection["R1_scaled"]), only)
    else:
        counts = gb.run_tasks(todo, matrix_dir, p[rep], inputs["config"], workers, only)
    return {"representation": rep, "tasks": len(todo), "counts": dict(counts)}




# ---------------------------------------------------------------- pilot: model-fit resources

def fit_pilot_tasks(tasks: list[dict], n_largest: int = 2) -> list[dict]:
    """The largest tasks by training size and the median one; chosen from sizes only."""
    by_size = sorted(tasks, key=lambda t: (-len(t["train"]), t["task_id"]))
    median = by_size[len(by_size) // 2]
    return by_size[:n_largest] + ([median] if median not in by_size[:n_largest] else [])


def _fit_only(args):
    """Fit one task in a fresh process and return its time and peak memory; the predictions are discarded
    without being scored."""
    task, matrix_dir, gconfig = args
    started = time.time()
    X, row_of, _ = gb.load_matrix(Path(matrix_dir))
    status, _ = gb.run_task(task, X, row_of, gconfig)
    return {"task_id": task["task_id"], "n_train": len(task["train"]), "seconds": round(time.time() - started, 2),
            "fit_seconds": status["seconds"], "features_used": status["features_used"],
            "peak_rss_bytes": max_rss()}


def stage_fit_pilot(config: dict, outdir: Path, benchmark_dir: Path, population: Path, baseline: Path,
                    scaled: int) -> dict:
    p = paths(outdir)
    d = config["density"]
    if scaled not in (d["reference_scaled"], *d["candidates"]):
        raise ValueError(f"scaled {scaled} is not a candidate")
    inputs, tasks, _, _ = load_benchmark(config, benchmark_dir, population, baseline)
    name = f"pilot_scaled_{scaled}"
    mdir = p["root"] / name / "matrix"
    if not (mdir / "feature_manifest.json").exists():
        stage_matrix(config, outdir, scaled, name)
    chosen = fit_pilot_tasks(tasks)
    slim = lambda t: {**t, "train": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")}  # noqa: E731
                                     for e in t["train"]],
                      "test": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["test"]]}
    rows = []
    for t in chosen:
        with ProcessPoolExecutor(max_workers=1) as pool:
            rows.append(pool.submit(_fit_only, (slim(t), str(mdir), inputs["config"])).result())
    fm = json.loads((mdir / "feature_manifest.json").read_text(encoding="utf-8"))
    out = {"scaled": scaled, "nnz": fm["nnz"], "patterns": fm["patterns"], "build_seconds": fm.get("build_seconds"),
           "build_max_rss_bytes": fm.get("build_max_rss_bytes"), "fits": rows, "machine": os.uname().nodename,
           "cpus": os.cpu_count(), "measured_at": now()}
    print(json.dumps(out), flush=True)  # in the job log too, in case the write below fails
    p["pilot"].mkdir(parents=True, exist_ok=True)
    gb.write_json(p["pilot"] / f"fit_pilot_scaled_{scaled}.json", out)
    return out


# ---------------------------------------------------------------- selection

def select_density(config: dict, resources: list[dict], budget: dict) -> dict:
    """The frozen rule: the densest candidate whose largest task fits one machine with 25% free and whose
    projected R1, R3 and R3_permuted compute fits the owner's budget. Only pilot resources are read."""
    d = config["density"]
    extra = budget["r3_tasks"] * 2 / config["benchmark"]["tasks"]
    rows = []
    for r in resources:
        sc = int(r["scaled"])
        if sc not in d["candidates"]:
            continue
        measured = budget.get("measured", {}).get(str(sc), {})
        peak = measured.get("peak_rss_bytes") or float(r["projected_task_peak_rss_bytes"])
        hours = (measured.get("task_hours_total") or float(r["projected_task_seconds_total"]) / 3600 *
                 budget["cpu_speed_factor"]) * (1 + extra)
        cores = max(1.0, peak / 2 ** 30 / budget["memory_gb_per_core"])
        su = budget["su_per_core_hour"] * cores * hours
        fits = peak <= 0.75 * budget["node_memory_gb"] * 2 ** 30
        rows.append({"scaled": sc, "peak_rss_bytes": peak, "task_hours": hours, "cores_charged": cores,
                     "projected_su": su, "memory_fits": fits, "within_budget": su <= budget["budget_su"],
                     "basis": "measured" if measured else "projected"})
    ok = [r for r in rows if r["memory_fits"] and r["within_budget"]]
    if not ok:
        raise ValueError("no candidate is feasible under the budget")
    chosen = min(ok, key=lambda r: r["scaled"])
    return {"R1_scaled": chosen["scaled"], "rule": d["selection_rule"], "budget": budget, "candidates": rows,
            "selected_at": now()}


# ---------------------------------------------------------------- report

HOLDOUTS = ed.HOLDOUTS
DELTA_METRICS = ("AUROC", "AUPRC", "balanced_accuracy", "MCC", "sensitivity", "specificity")
FOUR_WAY = ("both_correct", "determinant_only_correct", "ml_only_correct", "both_wrong")


def read_endpoint_context(decomposition_dir: Path, cohort: str) -> dict:
    """(drug, strategy) -> isolate -> phenotype, split unit, determinant call and the R0 probability, from the
    error_decomposition_v1 context (which carries the frozen generalisation_v1 predictions)."""
    keep = ("antimicrobial", "split_strategy", "task_id", "isolate_id", "phenotype", "predicted_probability",
            "determinant_prediction", "split_unit")
    out = defaultdict(dict)
    with gzip.open(ed.paths(decomposition_dir)["context"], "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        fields = next(reader)
        idx = {f: fields.index(f) for f in ("cohort", *keep)}
        for v in reader:
            if v[idx["cohort"]] != cohort:
                continue
            r = {k: sys.intern(v[idx[k]]) for k in keep}
            out[(r["antimicrobial"], r["split_strategy"])][r["isolate_id"]] = r
    return out


def task_predictions(taskdir: Path, tasks: list[dict]) -> tuple[dict, dict]:
    """Pooled predictions per (drug, strategy) from task checkpoints whose split hash matches, and the status
    of every task (completed, failed, resource_failed or not_run)."""
    preds, status = defaultdict(dict), {}
    for t in tasks:
        js = Path(taskdir) / f"{t['task_id']}.json"
        if not js.exists():
            status[t["task_id"]] = {"status": "not_run", "reason": "no checkpoint"}
            continue
        st = json.loads(js.read_text(encoding="utf-8"))
        if st.get("status") == "success" and st.get("split_sha256") == gb.task_digest(t):
            status[t["task_id"]] = {"status": "completed", "reason": NA, "seconds": st.get("seconds"),
                                    "features_used": st.get("features_used"),
                                    "nonzero_coefficients": st.get("nonzero_coefficients"), "C": st.get("C")}
            for r in ed.read_tsv(js.with_suffix(".tsv")):
                preds[(t["antimicrobial"], t["strategy"])][r["isolate_id"]] = float(r["predicted_probability"])
        else:
            reason = st.get("reason", "split hash differs")
            kind = "resource_failed" if re.search(r"MemoryError|Killed|resource", str(reason)) else "failed"
            status[t["task_id"]] = {"status": kind, "reason": reason}
    return preds, status


def representation_predictions(config: dict, outdir: Path, benchmark_dir: Path, decomposition_dir: Path,
                               tasks: list[dict], plan: dict, ctx: dict) -> tuple[dict, dict, dict]:
    """Predictions, task statuses and task scope of every representation, common cohort."""
    p = paths(outdir)
    common = [t for t in tasks if t["cohort"] == "common"]
    r2 = r2_task_set(tasks, plan)
    aug = Path(decomposition_dir) / "augmented"
    preds, status, scope = {}, {}, {}
    preds["R0"] = {k: {i: float(r["predicted_probability"]) for i, r in v.items()} for k, v in ctx.items()}
    status["R0"] = {t["task_id"]: {"status": "completed", "reason": "reused generalisation_v1"} for t in common}
    scope["R0"] = common
    for rep, taskdir, ts in (("R1", p["R1"] / "tasks", common), ("R2", aug / "mutation_augmented" / "tasks", r2),
                             ("R2_permuted", aug / "permuted_control" / "tasks", r2),
                             ("R3", p["R3"] / "tasks", r2), ("R3_permuted", p["R3_permuted"] / "tasks", r2)):
        preds[rep], status[rep] = task_predictions(taskdir, ts)
        scope[rep] = ts
    return preds, status, scope


def keys_complete(status: dict, scope: list[dict]) -> set:
    by = defaultdict(list)
    for t in scope:
        by[(t["antimicrobial"], t["strategy"])].append(status[t["task_id"]]["status"] == "completed")
    return {k for k, v in by.items() if all(v)}


def perf_row(y, p, threshold, sup) -> dict:
    m = gb.metrics(np.asarray(y), np.asarray(p), threshold, True)
    m["supported"] = "yes" if m["N_R"] >= sup["metric_min_R"] and m["N_S"] >= sup["metric_min_S"] else "no"
    return m


def fnum(v):
    return None if v in (None, NA, "") else float(v)


def comparison_tables(config: dict, gconfig: dict, ctx: dict, preds: dict, status: dict, scope: dict,
                      explain: list[dict], groups: dict) -> tuple[dict, dict, dict]:
    thr = gconfig["model"]["threshold"]
    sup = config["support"]
    boot = config["bootstrap"]
    done = {rep: keys_complete(status[rep], scope[rep]) for rep in preds}
    reps = [r for r in REPRESENTATIONS if r in preds]
    out = defaultdict(list)
    gap = {}

    def arrays(key, iso, rep):
        rows = ctx[key]
        y = np.array([rows[i]["phenotype"] == "R" for i in iso])
        return y, np.array([preds[rep][key][i] for i in iso]), [rows[i]["split_unit"] for i in iso]

    keys = sorted(ctx)
    for key in keys:
        drug, strat = key
        iso_all = sorted(ctx[key])
        for rep in reps:
            if key not in done[rep]:
                out["metrics"].append({"representation": rep, "antimicrobial": drug, "split_strategy": strat,
                                       "status": "not_applicable" if not any(
                                           (t["antimicrobial"], t["strategy"]) == key for t in scope[rep])
                                       else "incomplete"})
                continue
            if set(preds[rep][key]) != set(iso_all):
                raise IdentityError(f"{rep} {key}: test endpoints differ from R0")
            y, pr, _ = arrays(key, iso_all, rep)
            out["metrics"].append({"representation": rep, "antimicrobial": drug, "split_strategy": strat,
                                   "status": "completed", **perf_row(y, pr, thr, sup)})
            if rep != "R0" and key in done["R0"]:
                y, p0, units = arrays(key, iso_all, "R0")
                _, px, _ = arrays(key, iso_all, rep)
                a, b = gb.metrics(y, p0, thr, False), gb.metrics(y, px, thr, False)
                row = {"representation": rep, "antimicrobial": drug, "split_strategy": strat, "N": len(y),
                       "N_R": int(y.sum())}
                ci = gb.bootstrap_ci(y, p0, units, thr, gconfig, q=px) if y.any() and not y.all() else {}
                for m in DELTA_METRICS:
                    row[f"R0_{m}"], row[f"Rx_{m}"] = a[m], b[m]
                    row[f"delta_{m}"] = None if a[m] is None or b[m] is None else b[m] - a[m]
                    lo, hi = ci.get(f"delta_{m}_ci", (None, None))
                    row[f"delta_{m}_ci_lower"], row[f"delta_{m}_ci_upper"] = lo, hi
                out["deltas"].append(row)
    # generalisation gap per representation: random minus holdout on the holdout's test isolates
    for (drug, strat) in keys:
        if strat not in HOLDOUTS:
            continue
        for rep in reps:
            k_h, k_r = (drug, strat), (drug, "random")
            if k_h not in done[rep] or k_r not in done[rep]:
                continue
            iso = sorted(i for i in ctx[k_h] if i in ctx[k_r])
            y = np.array([ctx[k_h][i]["phenotype"] == "R" for i in iso])
            units = [ctx[k_h][i]["split_unit"] for i in iso]
            ph = np.array([preds[rep][k_h][i] for i in iso])
            prd = np.array([preds[rep][k_r][i] for i in iso])
            a, b = gb.metrics(y, prd, thr, False), gb.metrics(y, ph, thr, False)
            ci = gb.bootstrap_ci(y, ph, units, thr, gconfig, q=prd) if y.any() and not y.all() else {}
            row = {"representation": rep, "antimicrobial": drug, "holdout_strategy": strat,
                   "drug_group": groups.get(drug, NA), "N": len(iso), "N_R": int(y.sum())}
            for m in DELTA_METRICS:
                row[f"{m}_random"], row[f"{m}_holdout"] = a[m], b[m]
                row[f"gap_{m}"] = None if a[m] is None or b[m] is None else a[m] - b[m]
                lo, hi = ci.get(f"delta_{m}_ci", (None, None))
                row[f"gap_{m}_ci_lower"], row[f"gap_{m}_ci_upper"] = lo, hi
            out["gaps"].append(row)
            gap[(rep, drug, strat)] = row
            if rep == "R0":
                continue
            if (drug, strat) not in done["R0"] or (drug, "random") not in done["R0"]:
                continue
            p0h = np.array([preds["R0"][k_h][i] for i in iso])
            p0r = np.array([preds["R0"][k_r][i] for i in iso])
            g0 = gap[("R0", drug, strat)]["gap_AUROC"]
            rescue_ci = ed.gap_rescue_ci(y, units, p0r, p0h, prd, ph, boot["replicates"], boot["seed"],
                                         boot["level"]) if y.any() and not y.all() else {"lower": None,
                                                                                         "upper": None}
            out["rescue"].append({"representation": rep, "antimicrobial": drug, "holdout_strategy": strat,
                                  "drug_group": groups.get(drug, NA), "N": len(iso), "N_R": int(y.sum()),
                                  "R0_gap_AUROC": g0, "Rx_gap_AUROC": row["gap_AUROC"],
                                  "gap_rescue": None if g0 is None or row["gap_AUROC"] is None
                                  else g0 - row["gap_AUROC"],
                                  "gap_rescue_ci_lower": rescue_ci["lower"], "gap_rescue_ci_upper": rescue_ci["upper"],
                                  "R0_AUROC_holdout": gap[("R0", drug, strat)]["AUROC_holdout"],
                                  "Rx_AUROC_holdout": row["AUROC_holdout"]})
    # determinant rule against each representation, and the four-way table
    for (drug, strat) in keys:
        rows = ctx[(drug, strat)]
        iso = sorted(i for i, r in rows.items() if r["determinant_prediction"] != NA)
        if not iso:
            continue
        y = np.array([rows[i]["phenotype"] == "R" for i in iso])
        dp = np.array([1.0 if rows[i]["determinant_prediction"] == "R" else 0.0 for i in iso])
        dm = gb.metrics(y, dp, 0.5, False)
        for rep in reps:
            if (drug, strat) not in done[rep]:
                continue
            px = np.array([preds[rep][(drug, strat)][i] for i in iso])
            mm = gb.metrics(y, px, thr, False)
            row = {"representation": rep, "antimicrobial": drug, "split_strategy": strat,
                   "drug_group": groups.get(drug, NA), "N": len(iso), "N_R": int(y.sum())}
            for m in ("AUROC", "sensitivity", "specificity", "balanced_accuracy", "MCC"):
                row[f"determinant_{m}"], row[f"ml_{m}"] = dm[m], mm[m]
                row[f"ml_minus_determinant_{m}"] = None if dm[m] is None or mm[m] is None else mm[m] - dm[m]
            c = Counter(ed.four_way(rows[i]["phenotype"], "R" if px[j] >= thr else "S",
                                    rows[i]["determinant_prediction"]) for j, i in enumerate(iso))
            out["determinant"].append(row)
            out["four_way"].append({"representation": rep, "antimicrobial": drug, "split_strategy": strat,
                                    "N": len(iso), **{k: c[k] for k in FOUR_WAY}})
    # pooled four-way per strategy, over all drugs and over the R2 drug set
    r2_drugs = {t["antimicrobial"] for t in scope["R2"]}
    pooled = defaultdict(Counter)
    for r in out["four_way"]:
        for drug_set in ("all", "r2_drugs"):
            if drug_set == "r2_drugs" and r["antimicrobial"] not in r2_drugs:
                continue
            k = (r["representation"], drug_set, r["split_strategy"])
            pooled[k].update({f: r[f] for f in FOUR_WAY})
            pooled[k]["N"] += r["N"]
            pooled[k]["drugs"] += 1
    for (rep, drug_set, strat), c in sorted(pooled.items(), key=lambda kv: (rep_order(kv[0][0]), kv[0][1:])):
        expected = len({d for d, s in ctx if s == strat} if drug_set == "all" else
                       {d for d, s in ctx if s == strat and d in r2_drugs})
        out["four_way_pooled"].append({"representation": rep, "drug_set": drug_set, "split_strategy": strat,
                                       "drugs": c["drugs"], "drugs_expected": expected, "N": c["N"],
                                       **{f: c[f] for f in FOUR_WAY}})
    # all-drug robustness and the correlations with the error_decomposition_v1 shift measures
    ex = {(r["antimicrobial"], r["holdout_strategy"]): r for r in explain}
    for rep in reps:
        for drug_set in ("all", "r2_drugs"):
            pairs = [k for k in sorted(ex) if (rep, *k) in gap and (drug_set == "all" or k[0] in r2_drugs)]
            if not pairs:
                continue
            g = [gap[(rep, *k)]["gap_AUROC"] for k in pairs]
            lo = [gap[(rep, *k)]["gap_AUROC_ci_lower"] for k in pairs]
            out["robustness"].append({
                "representation": rep, "drug_set": drug_set, "pairs": len(pairs), "pairs_expected": len(
                    [k for k in ex if drug_set == "all" or k[0] in r2_drugs]),
                "median_gap_AUROC": float(np.median([v for v in g if v is not None])),
                "pairs_gap_positive": sum(v is not None and v > 0 for v in g),
                "pairs_gap_ci_above_zero": sum(v is not None and v > 0 for v in lo),
                "largest_gap_AUROC": max(v for v in g if v is not None)})
            for field in ("unseen_mechanism_fraction", "lineage_shift", "mechanism_shift",
                          "phenotype_prevalence_shift", "mutation_associated_fraction"):
                s = ed._spearman([ex[k][field] for k in pairs], [gap[(rep, *k)]["gap_AUROC"] for k in pairs])
                out["correlations"].append({"representation": rep, "drug_set": drug_set, "measure": field,
                                            "pairs": s["n"], "spearman_rho": s["rho"]})
    return dict(out), gap, done


# Direct comparisons between the new representations (reference, other). Added after the report: decision
# 0011 says R3 shows whether the indicators still help once k-mers are denser, but its tables compare each
# representation with R0 only.
PAIRS = (("R1", "R2"), ("R1", "R3"), ("R3_permuted", "R3"), ("R2", "R3"))


def paired_auroc_ci(y, units, p, q, replicates: int, seed: int, level: float) -> dict:
    """Paired bootstrap over split units of AUROC(q) - AUROC(p)."""
    vals = []
    for idx in gb.bootstrap_indices(units, replicates, seed):
        a, b = gb.auroc(y[idx], p[idx]), gb.auroc(y[idx], q[idx])
        if a is not None and b is not None:
            vals.append(b - a)
    if not vals:
        return {"lower": None, "upper": None}
    alpha = (1 - level) / 2
    return {"lower": float(np.quantile(vals, alpha)), "upper": float(np.quantile(vals, 1 - alpha))}


def pairwise_tables(config: dict, ctx: dict, preds: dict, status: dict, scope: dict, groups: dict) -> list[dict]:
    boot = config["bootstrap"]
    done = {rep: keys_complete(status[rep], scope[rep]) for rep in preds}
    rows = []
    for ref, other in PAIRS:
        if ref not in preds or other not in preds:
            continue
        for drug, strat in sorted(done[ref] & done[other]):
            key = (drug, strat)
            iso = sorted(ctx[key])
            y = np.array([ctx[key][i]["phenotype"] == "R" for i in iso])
            if not y.any() or y.all():
                continue
            units = [ctx[key][i]["split_unit"] for i in iso]
            p = np.array([preds[ref][key][i] for i in iso])
            q = np.array([preds[other][key][i] for i in iso])
            a, b = gb.auroc(y, p), gb.auroc(y, q)
            ci = paired_auroc_ci(y, units, p, q, boot["replicates"], boot["seed"], boot["level"])
            row = {"reference": ref, "representation": other, "antimicrobial": drug, "split_strategy": strat,
                   "drug_group": groups.get(drug, NA), "N": len(iso), "N_R": int(y.sum()),
                   "reference_AUROC": a, "representation_AUROC": b, "delta_AUROC": b - a,
                   "delta_AUROC_ci_lower": ci["lower"], "delta_AUROC_ci_upper": ci["upper"],
                   "reference_gap": None, "representation_gap": None, "gap_reduction": None,
                   "gap_reduction_ci_lower": None, "gap_reduction_ci_upper": None}
            k_r = (drug, "random")
            if strat in HOLDOUTS and k_r in done[ref] and k_r in done[other]:
                hiso = [i for i in iso if i in ctx[k_r]]
                yh = np.array([ctx[key][i]["phenotype"] == "R" for i in hiso])
                hu = [ctx[key][i]["split_unit"] for i in hiso]
                arr = {(rep, s): np.array([preds[rep][(drug, s)][i] for i in hiso])
                       for rep in (ref, other) for s in (strat, "random")}
                g = {rep: (gb.auroc(yh, arr[(rep, "random")]), gb.auroc(yh, arr[(rep, strat)]))
                     for rep in (ref, other)}
                if None not in g[ref] + g[other]:
                    rc = ed.gap_rescue_ci(yh, hu, arr[(ref, "random")], arr[(ref, strat)], arr[(other, "random")],
                                          arr[(other, strat)], boot["replicates"], boot["seed"], boot["level"])
                    row.update({"reference_gap": g[ref][0] - g[ref][1], "representation_gap": g[other][0] - g[other][1],
                                "gap_reduction": (g[ref][0] - g[ref][1]) - (g[other][0] - g[other][1]),
                                "gap_reduction_ci_lower": rc["lower"], "gap_reduction_ci_upper": rc["upper"]})
            rows.append(row)
    return rows


def stage_pairwise(config: dict, outdir: Path, benchmark_dir: Path, decomposition_dir: Path, population: Path,
                   baseline: Path) -> dict:
    p = paths(outdir)
    inputs, tasks, _, _ = load_benchmark(config, benchmark_dir, population, baseline)
    plan = json.loads((resolve(config["upstream"]["error_decomposition_outputs"]) /
                       "augmentation_plan.json").read_text(encoding="utf-8"))
    ctx = read_endpoint_context(decomposition_dir, config["analysis"]["cohort"])
    preds, status, scope = representation_predictions(config, outdir, benchmark_dir, decomposition_dir, tasks,
                                                       plan, ctx)
    eout = resolve(config["upstream"]["error_decomposition_outputs"])
    groups = {k: v["group"] for k, v in json.loads((eout / "error_decomposition_manifest.json")
                                                   .read_text(encoding="utf-8"))["drug_groups"].items()}
    rows = pairwise_tables(config, ctx, preds, status, scope, groups)
    o = p["root"] / "tables"
    o.mkdir(parents=True, exist_ok=True)
    write_tsv(o / "representation_pairwise.tsv", list(rows[0].keys()), rows)
    return {"rows": len(rows), "pairs": sorted({f"{r['representation']}-{r['reference']}" for r in rows})}


def drug_capture(capture_rows: list[dict], r1_label: str) -> dict:
    """(drug, kind) -> {'R0': fraction, 'R1': fraction, 'endpoints': n} for the ALL determinant rows."""
    out = defaultdict(dict)
    for r in capture_rows:
        if r["determinant"] != "ALL" or r["representation"] not in ("R0", r1_label):
            continue
        name = "R0" if r["representation"] == "R0" else "R1"
        out[(r["antimicrobial"], r["determinant_kind"])][name] = fnum(r["captured_fraction"])
        out[(r["antimicrobial"], r["determinant_kind"])]["endpoints"] = int(r["endpoints"])
    return out


def drug_view(drugs: list[str], ctx: dict, metrics: dict, gap: dict, rescue: dict, cap: dict, det: dict,
              groups: dict, role: str) -> list[dict]:
    """One row per drug and strategy with every representation side by side."""
    rows = []
    for drug in drugs:
        for strat in ("random", *HOLDOUTS):
            if (drug, strat) not in ctx:
                continue
            r = {"antimicrobial": drug, "role": role, "drug_group": groups.get(drug, NA), "split_strategy": strat,
                 "N": len(ctx[(drug, strat)]),
                 "mutation_carrier_endpoints": cap.get((drug, "point_mutation"), {}).get("endpoints", 0),
                 "R0_mutation_capture": cap.get((drug, "point_mutation"), {}).get("R0"),
                 "R1_mutation_capture": cap.get((drug, "point_mutation"), {}).get("R1"),
                 "R0_gene_capture": cap.get((drug, "gene"), {}).get("R0"),
                 "R1_gene_capture": cap.get((drug, "gene"), {}).get("R1")}
            for rep in REPRESENTATIONS:
                m = metrics.get((rep, drug, strat), {})
                r[f"{rep}_AUROC"], r[f"{rep}_sensitivity"] = m.get("AUROC"), m.get("sensitivity")
                g = gap.get((rep, drug, strat), {})
                r[f"{rep}_gap"] = g.get("gap_AUROC")
                if rep != "R0":
                    r[f"{rep}_gap_rescue"] = rescue.get((rep, drug, strat), {}).get("gap_rescue")
            r["determinant_AUROC"] = det.get((drug, strat))
            rows.append(r)
    return rows


def stage_report(config: dict, outdir: Path, benchmark_dir: Path, decomposition_dir: Path, population: Path,
                 baseline: Path) -> dict:
    p = paths(outdir)
    o = p["root"] / "tables"
    o.mkdir(parents=True, exist_ok=True)
    selection = json.loads(p["selection"].read_text(encoding="utf-8"))
    r1_label = f"scaled_{selection['R1_scaled']}"
    inputs, tasks, skipped, manifest = load_benchmark(config, benchmark_dir, population, baseline)
    gconfig = inputs["config"]
    # R1 benchmark tables, written by the generalisation_v1 report code into R1/
    if (p["R1"] / "tasks").exists():
        gb.report(inputs, tasks, skipped, p["R1"], manifest)
    plan = json.loads((resolve(config["upstream"]["error_decomposition_outputs"]) /
                       "augmentation_plan.json").read_text(encoding="utf-8"))
    ctx = read_endpoint_context(decomposition_dir, config["analysis"]["cohort"])
    preds, status, scope = representation_predictions(config, outdir, benchmark_dir, decomposition_dir, tasks,
                                                       plan, ctx)
    eout = resolve(config["upstream"]["error_decomposition_outputs"])
    explain = ed.read_tsv(eout / "explain_generalisation_drop.tsv")
    groups = {k: v["group"] for k, v in json.loads((eout / "error_decomposition_manifest.json")
                                                   .read_text(encoding="utf-8"))["drug_groups"].items()}
    tables, gap, done = comparison_tables(config, gconfig, ctx, preds, status, scope, explain, groups)
    capture_rows = ed.read_tsv(p["pilot"] / "drug_capture_by_candidate.tsv")
    cap = drug_capture(capture_rows, r1_label)
    metrics = {(r["representation"], r["antimicrobial"], r["split_strategy"]): r for r in tables["metrics"]
               if r["status"] == "completed"}
    rescue = {(r["representation"], r["antimicrobial"], r["holdout_strategy"]): r for r in tables["rescue"]}
    det = {(r["antimicrobial"], r["split_strategy"]): r["determinant_AUROC"] for r in tables["determinant"]
           if r["representation"] == "R0"}
    mfields = ["representation", "antimicrobial", "split_strategy", "status", *gb.METRIC_FIELDS[:19], "supported"]
    write_tsv(o / "representation_metrics.tsv", mfields, tables["metrics"])
    for name, key in (("representation_metric_deltas.tsv", "deltas"),
                      ("generalisation_gap_by_representation.tsv", "gaps"),
                      ("generalisation_gap_rescue.tsv", "rescue"), ("determinant_vs_representation.tsv", "determinant"),
                      ("rule_ml_discordance.tsv", "four_way"), ("rule_ml_discordance_pooled.tsv", "four_way_pooled"),
                      ("all_drug_robustness.tsv", "robustness"), ("correlation_robustness.tsv", "correlations")):
        if tables.get(key):
            write_tsv(o / name, list(tables[key][0].keys()), tables[key])
    # task status of every representation
    st_rows = []
    for rep in REPRESENTATIONS:
        for t in scope[rep]:
            s = status[rep][t["task_id"]]
            st_rows.append({"representation": rep, "task_id": t["task_id"], "antimicrobial": t["antimicrobial"],
                            "split_strategy": t["strategy"], "fold_or_holdout_id": t["fold"], "status": s["status"],
                            "reason": s.get("reason", NA), "seconds": s.get("seconds", NA),
                            "features_used": s.get("features_used", NA), "C": s.get("C", NA)})
    write_tsv(o / "task_status.tsv", list(st_rows[0].keys()), st_rows)
    # drug groups the hypotheses name
    c = config["capture"]
    views = {"quinolone_representation_analysis.tsv": (c["quinolones"], "quinolone"),
             "stable_drug_controls.tsv": (c["stable_controls"], "stable_control"),
             "cefazolin_control.tsv": ([c["mechanism_shift_control"]], "mechanism_shift_control"),
             "unresolved_drug_controls.tsv": (c["unresolved_controls"], "unresolved"),
             "levofloxacin_counterexample.tsv": ([c["counterexample"]], "counterexample")}
    for name, (drugs, role) in views.items():
        rows = drug_view(drugs, ctx, metrics, gap, rescue, cap, det, groups, role)
        write_tsv(o / name, list(rows[0].keys()), rows)
    all_rows = drug_view(sorted({d for d, _ in ctx}), ctx, metrics, gap, rescue, cap, det, groups, "all")
    primary = [{"antimicrobial": r["antimicrobial"], "split_strategy": r["split_strategy"],
                "R0_AUROC": r["R0_AUROC"], "R1_AUROC": r["R1_AUROC"], "R2_AUROC": r["R2_AUROC"],
                "R3_AUROC": r["R3_AUROC"], "R0_generalisation_gap": r["R0_gap"], "R1_generalisation_gap": r["R1_gap"],
                "R2_generalisation_gap": r["R2_gap"], "R1_gap_rescue": r["R1_gap_rescue"],
                "R2_gap_rescue": r["R2_gap_rescue"], "R0_known_mutation_capture": r["R0_mutation_capture"],
                "R1_known_mutation_capture": r["R1_mutation_capture"],
                "determinant_baseline_AUROC": r["determinant_AUROC"]} for r in all_rows]
    write_tsv(o / "manuscript_primary_table.tsv", list(primary[0].keys()), primary)
    quin = [{"antimicrobial": r["antimicrobial"], "holdout": r["split_strategy"],
             "known_mutation_carriers": r["mutation_carrier_endpoints"],
             "R0_mutation_capture": r["R0_mutation_capture"], "R1_mutation_capture": r["R1_mutation_capture"],
             "R0_AUROC": r["R0_AUROC"], "R1_AUROC": r["R1_AUROC"], "R2_AUROC": r["R2_AUROC"],
             "permutation_AUROC": r["R2_permuted_AUROC"], "R3_AUROC": r["R3_AUROC"],
             "R3_permutation_AUROC": r["R3_permuted_AUROC"], "R0_sensitivity": r["R0_sensitivity"],
             "R1_sensitivity": r["R1_sensitivity"], "R2_sensitivity": r["R2_sensitivity"], "R0_gap": r["R0_gap"],
             "R1_gap": r["R1_gap"], "R2_gap": r["R2_gap"], "R3_gap": r["R3_gap"]}
            for r in all_rows if r["antimicrobial"] in c["quinolones"]]
    write_tsv(o / "manuscript_quinolone_table.tsv", list(quin[0].keys()), quin)
    figure_data(o / "figure_data", p, capture_rows, tables, all_rows)
    summary = {"analysis_version": config["analysis_version"], "R1_scaled": selection["R1_scaled"],
               "tasks": {rep: dict(Counter(s["status"] for s in status[rep].values())) for rep in REPRESENTATIONS},
               "complete_drug_strategy_pairs": {rep: len(done[rep]) for rep in REPRESENTATIONS},
               "correlations": [r for r in tables["correlations"]
                                if r["measure"] in ("unseen_mechanism_fraction", "lineage_shift")],
               "four_way_pooled": tables["four_way_pooled"], "robustness": tables["robustness"],
               "generated_at": now()}
    gb.write_json(o / "representation_sensitivity_summary.json", summary)
    return {"tables": sorted(x.name for x in o.glob("*.tsv")), "tasks": summary["tasks"]}


def figure_data(fig: Path, p: dict, capture_rows: list[dict], tables: dict, all_rows: list[dict]) -> None:
    fig.mkdir(parents=True, exist_ok=True)
    resources = ed.read_tsv(p["pilot"] / "representation_resource_pilot.tsv")
    mut = ed.read_tsv(p["pilot"] / "mutation_capture_by_candidate.tsv")
    a = [{**{k: r[k] for k in ("representation", "mutation", "calls", "captured", "partially_represented",
                               "not_represented", "captured_fraction")},
          "scaled": 200 if r["representation"] == "R0" else int(r["representation"].split("_")[1])} for r in mut]
    write_tsv(fig / "fig_a_mutation_capture_by_density.tsv", list(a[0].keys()), a)
    b = [r for r in all_rows if r["split_strategy"] in HOLDOUTS]
    write_tsv(fig / "fig_b_R0_vs_R1_holdout_auroc.tsv", ["antimicrobial", "split_strategy", "drug_group", "R0_AUROC",
                                                          "R1_AUROC"], b)
    write_tsv(fig / "fig_c_quinolone_rescue.tsv",
              ["antimicrobial", "split_strategy", "R1_gap_rescue", "R2_gap_rescue", "R2_permuted_gap_rescue",
               "R3_gap_rescue", "R3_permuted_gap_rescue"], [r for r in b if r["role"] == "all" and
                                                            r["antimicrobial"] in ("nalidixic acid", "ciprofloxacin",
                                                                                   "levofloxacin")])
    write_tsv(fig / "fig_d_gap_by_representation.tsv", list(tables["gaps"][0].keys()), tables["gaps"])
    write_tsv(fig / "fig_e_determinant_vs_representation.tsv", list(tables["determinant"][0].keys()),
              tables["determinant"])
    write_tsv(fig / "fig_f_gain_vs_cost.tsv", list(resources[0].keys()), resources)


# ---------------------------------------------------------------- command line

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("check", "sketch", "pilot", "fitpilot", "select", "matrix", "run",
                                            "report", "pairwise"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--benchmark", type=Path, required=True, help="generalisation_v1 output directory")
    parser.add_argument("--decomposition", type=Path, required=True, help="error_decomposition_v1 output directory")
    parser.add_argument("--population", type=Path, help="PR #6 population build directory")
    parser.add_argument("--baseline", type=Path, help="PR #7 baseline build directory")
    parser.add_argument("--cache", type=Path, help="genome cache (sketch, pilot)")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--representation", choices=NEW_RUNS, default="R1")
    parser.add_argument("--only", help="regular expression on task IDs (run)")
    parser.add_argument("--scaled", type=int, help="candidate scaled value (fitpilot)")
    parser.add_argument("--memory-gb", type=float,
                        help="job memory (run): start tasks only while their estimated memory fits 75%% of it")
    parser.add_argument("--budget", type=Path, help="JSON with the owner's compute budget and cost model (select)")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    up = check_upstream(config, args.benchmark, args.decomposition)
    p = paths(args.outdir)
    p["root"].mkdir(parents=True, exist_ok=True)
    if args.command == "check":
        manifest = {"analysis_version": config["analysis_version"], "builder_version": VERSION,
                    "config_sha256": pop.sha256(args.config), "upstream": up, "checked_at": now()}
        gb.write_json(p["manifest"], manifest)
        print(json.dumps({k: up[k] for k in ("tasks", "mutation_panel_sha256")}))
        return 0
    if pop.sha256(args.config) != json.loads(p["manifest"].read_text(encoding="utf-8"))["config_sha256"]:
        raise ValueError("analysis config changed after check; start a new analysis version")
    if args.command == "sketch":
        if args.cache is None or args.baseline is None:
            parser.error("sketch needs --cache and --baseline")
        s = stage_sketch(config, args.outdir, args.benchmark, args.cache, args.baseline, args.workers, args.limit)
        print(json.dumps(s))
        return 0 if set(s["status"]) == {"genome_verified"} and set(s["r0_reproduced"]) == {"yes"} else 1
    if args.command == "pilot":
        if args.cache is None:
            parser.error("pilot needs --cache")
        print(json.dumps(stage_pilot(config, args.outdir, args.benchmark, args.decomposition, args.cache,
                                     args.workers)))
        return 0
    if args.command == "matrix":
        sel = json.loads(p["selection"].read_text(encoding="utf-8"))
        print(json.dumps(stage_matrix(config, args.outdir, sel["R1_scaled"], "R1")))
        return 0
    if args.command == "run":
        if args.population is None or args.baseline is None:
            parser.error("run needs --population and --baseline")
        s = stage_run(config, args.outdir, args.benchmark, args.decomposition, args.population, args.baseline,
                      args.representation, args.workers, args.only, args.memory_gb)
        print(json.dumps(s))
        return 0 if not s["counts"].get("failed") else 1
    if args.command == "fitpilot":
        if args.population is None or args.baseline is None or args.scaled is None:
            parser.error("fitpilot needs --population, --baseline and --scaled")
        r = stage_fit_pilot(config, args.outdir, args.benchmark, args.population, args.baseline, args.scaled)
        print(json.dumps({"scaled": r["scaled"], "nnz": r["nnz"], "fits": r["fits"]}))
        return 0
    if args.command == "select":
        if args.budget is None:
            parser.error("select needs --budget")
        if p["selection"].exists():
            raise ValueError("R1 is already selected; the selection is not repeated")
        resources = ed.read_tsv(p["pilot"] / "representation_resource_pilot.tsv")
        sel = select_density(config, resources, json.loads(args.budget.read_text(encoding="utf-8")))
        gb.write_json(p["selection"], sel)
        print(json.dumps({"R1_scaled": sel["R1_scaled"]}))
        return 0
    if args.population is None or args.baseline is None:
        parser.error(f"{args.command} needs --population and --baseline")
    if args.command == "pairwise":
        print(json.dumps(stage_pairwise(config, args.outdir, args.benchmark, args.decomposition, args.population,
                                        args.baseline)))
        return 0
    print(json.dumps(stage_report(config, args.outdir, args.benchmark, args.decomposition, args.population,
                                  args.baseline)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
