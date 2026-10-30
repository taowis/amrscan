#!/usr/bin/env python3
"""Error decomposition of the E. coli generalisation benchmark (error_decomposition_v1).

Reads the frozen generalisation_v1 predictions and asks which errors come with
established determinants that the k-mer features did or did not encode, with
determinants not seen in training, or with no established determinant.

Commands, run in this order:
  catalogue  build the mechanism catalogue from determinant mapping v1 and the
             AMRFinderPlus database (committed file)
  freeze     check upstream hashes, fix the drug groups from the PR #14 drops
             and write the analysis manifest; nothing downstream is read yet
  calls      collect AMRFinderPlus AMR calls (genes and point mutations) with
             their coordinates, checked against each call's provenance
  audit      verify every genome against the PR #8 manifest, locate each called
             determinant and test whether the PR #14 features encode it
             (resumable)
  context    join every PR #14 test prediction to determinants, capture status,
             mechanism and population novelty, phenotype and genome QC
  augment    refit the frozen PR #14 tasks with explicit point-mutation
             indicators, and with the same indicators permuted (sensitivity
             analysis; resumable)
  report     decomposition, shift, novelty, QC, augmentation and summary tables

Nothing here changes a PR #14 split, endpoint, label or prediction.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
import math
from pathlib import Path
import re
import sys
import time

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generalisation_benchmark as gb  # noqa: E402
import population_structure as pop  # noqa: E402

VERSION = "1.0.0"
NA = "NA"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "benchmark/cohorts/ecoli_research_v1/error_decomposition_v1.json"
HOLDOUTS = ("lineage", "country", "project", "geography", "temporal")
CAPTURE_ORDER = ("captured", "partially_represented", "not_represented", "representation_unknown")
MECHANISM_CLASSES = ("acquired_gene", "intrinsic_gene", "target_site_mutation", "promoter_mutation",
                     "regulatory_mutation", "porin_mutation", "other_mutation")
MUTATION_RE = re.compile(r"^(?P<gene>[A-Za-z0-9]+?)_(?P<ref>[A-Z]+)(?P<pos>-?\d+)(?P<alt>[A-Z]+|Ter)$")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
PARTIAL_METHODS = ("PARTIALP", "PARTIALX", "PARTIAL_CONTIG_ENDP", "PARTIAL_CONTIG_ENDX", "INTERNAL_STOP")

write_tsv, fmt, now, sha256_bytes = gb.write_tsv, gb.fmt, gb.now, gb.sha256_bytes


class UpstreamError(ValueError):
    pass


class CatalogueError(ValueError):
    pass


def read_tsv(path: Path) -> list[dict]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path: Path, value) -> None:
    pop.write_json_atomic(Path(path), value)


class ContextWriter:
    """Streams context rows to a gzip TSV, renamed into place on close."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.tmp = self.path.with_name(f".{self.path.name}.tmp")
        self.handle = gzip.open(self.tmp, "wt", encoding="utf-8", newline="")
        self.writer = csv.DictWriter(self.handle, fieldnames=CONTEXT_FIELDS, delimiter="\t", lineterminator="\n",
                                     extrasaction="ignore")
        self.writer.writeheader()

    def write(self, row: dict) -> None:
        self.writer.writerow({k: fmt(row.get(k, NA)) for k in CONTEXT_FIELDS})

    def close(self) -> None:
        self.handle.close()
        self.tmp.replace(self.path)


class Row:
    """One context row; slots keep a million of them small."""
    __slots__ = ()

    def __getitem__(self, key):
        return getattr(self, key)

    def get(self, key, default=None):
        return getattr(self, key, default)


def read_context(path: Path, keep=None) -> list:
    """Context rows as slotted objects with interned values; `keep` filters on the raw dict."""
    fields = None
    out = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        fields = next(reader)
        cls = type("ContextRow", (Row,), {"__slots__": tuple(fields)})
        for values in reader:
            if keep is not None and not keep(dict(zip(fields, values))):
                continue
            obj = cls()
            for f, v in zip(fields, values):
                setattr(obj, f, sys.intern(v))
            out.append(obj)
    return out


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


# ---------------------------------------------------------------- config and upstream

def load_config(path: Path) -> dict:
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("analysis_version", "upstream", "catalogue", "mechanism_classes", "drug_groups", "capture",
                "novelty", "augmentation", "bootstrap", "support", "decision_record"):
        if key not in config:
            raise ValueError(f"analysis config lacks {key}")
    g = config["drug_groups"]
    if not 0 <= g["stable_max_delta"] < g["high_loss_min_delta"]:
        raise ValueError("stable_max_delta must lie below high_loss_min_delta")
    return config


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def check_upstream(config: dict, benchmark_dir: Path) -> dict:
    """Stop unless every frozen input matches the hash in the config."""
    up = config["upstream"]
    checks = {
        "generalisation_config": resolve(up["generalisation_config"]),
        "benchmark_manifest": Path(benchmark_dir) / "benchmark_manifest.json",
        "predictions": Path(benchmark_dir) / "predictions.tsv.gz",
        "feature_manifest": Path(benchmark_dir) / "matrix" / "feature_manifest.json",
        "mapping_rules": resolve(up["mapping_rules"]),
        "genome_manifest": resolve(up["genome_manifest"]),
    }
    got = {}
    for name, path in checks.items():
        if not path.is_file():
            raise UpstreamError(f"{name} missing: {path}")
        digest = pop.sha256(path)
        if digest != up[f"{name}_sha256"]:
            raise UpstreamError(f"{name} differs from the frozen input: {digest}")
        got[f"{name}_sha256"] = digest
    db = resolve(up["amrfinder_database"])
    for name, expected in up["amrfinder_database_files_sha256"].items():
        if pop.sha256(db / name) != expected:
            raise UpstreamError(f"AMRFinderPlus database file {name} differs")
    manifest = json.loads(checks["benchmark_manifest"].read_text(encoding="utf-8"))
    if manifest["benchmark_version"] != up["benchmark_version"]:
        raise UpstreamError("benchmark version differs")
    if manifest["genome_manifest_sha256"] != up["genome_manifest_sha256"]:
        raise UpstreamError("the benchmark used a different genome manifest")
    got.update({"benchmark_version": manifest["benchmark_version"], "benchmark_upstream": manifest["upstream"],
                "amrfinder_database_files_sha256": up["amrfinder_database_files_sha256"]})
    return got


# ---------------------------------------------------------------- mechanism catalogue

CATALOGUE_FIELDS = ("catalogue_version organism antimicrobial determinant determinant_type mechanism_class "
                    "gene_family gene_function interpretation requires combination mapping_rule_id mapping_version "
                    "evidence_source evidence_reference review_status locus_type locus_reference locus_position "
                    "ref_allele alt_allele amrfinder_class amrfinder_subclass").split()


def parse_mutation(symbol: str) -> dict | None:
    m = MUTATION_RE.match(symbol)
    if not m:
        return None
    d = m.groupdict()
    d["pos"] = int(d["pos"])
    if d["alt"] == "Ter":
        d["alt"] = "*"
    return d


def mutation_loci(db: Path, taxgroup: str) -> dict:
    """symbol -> locus from the AMRFinderPlus mutation tables (protein and DNA)."""
    loci = {}
    for r in read_tsv(Path(db) / "AMRProt-mutation.tsv"):
        if r["#taxgroup"] == taxgroup:
            loci[r["standard_mutation_symbol"]] = {
                "locus_type": "protein", "locus_reference": r["accession_version"],
                "locus_position": int(r["mutation_position"]), "amrfinder_class": r["class"],
                "amrfinder_subclass": r["subclass"]}
    for r in read_tsv(Path(db) / f"AMR_DNA-{taxgroup}.tsv"):
        loci[r["standard_mutation_symbol"]] = {
            "locus_type": "dna", "locus_reference": r["#accesion_version@gene_name@gene_symbol:start-stop"],
            "locus_position": int(r["mutation_position"]), "amrfinder_class": r["class"],
            "amrfinder_subclass": r["subclass"]}
    return loci


def gene_functions(db: Path) -> dict:
    out = {}
    for r in read_tsv(Path(db) / "fam.tsv"):
        out[r["#node_id"]] = "efflux" if "efflux" in r["family_name"].lower() else "other"
    return out


def mechanism_class(row: dict, classes: dict) -> str:
    if row["determinant_type"] == "gene":
        return "intrinsic_gene" if row["gene_family"] in classes["intrinsic_gene_families"] else "acquired_gene"
    gene = row["determinant"].split("_")[0]
    return classes["point_mutation_genes"].get(gene, "unclassified")


def build_catalogue(mapping_rows: list[dict], loci: dict, functions: dict, config: dict, study: set) -> list[dict]:
    cat = config["catalogue"]
    out = []
    for r in mapping_rows:
        if r["antimicrobial"] not in study:
            continue
        row = {"catalogue_version": cat["version"], "organism": r["organism"], "antimicrobial": r["antimicrobial"],
               "determinant": r["determinant"], "determinant_type": r["determinant_type"],
               "gene_family": r["gene_family"], "interpretation": r["interpretation"], "requires": r["requires"],
               "combination": "yes" if r["requires"] != NA else "no", "mapping_rule_id": r["rule_id"],
               "mapping_version": r["mapping_version"], "evidence_source": r["evidence_source"],
               "evidence_reference": r["evidence_reference"], "review_status": r["review_status"]}
        row["mechanism_class"] = mechanism_class(row, config["mechanism_classes"])
        row["gene_function"] = functions.get(r["determinant"], functions.get(r["gene_family"], NA)) \
            if r["determinant_type"] == "gene" else NA
        locus = loci.get(r["determinant"]) if r["determinant_type"] == "point_mutation" else None
        parsed = parse_mutation(r["determinant"]) if locus else None
        row.update(locus or {"locus_type": NA, "locus_reference": NA, "locus_position": NA,
                             "amrfinder_class": NA, "amrfinder_subclass": NA})
        row["ref_allele"] = parsed["ref"] if parsed else NA
        row["alt_allele"] = parsed["alt"] if parsed else NA
        out.append(row)
    out.sort(key=lambda x: (x["antimicrobial"], x["determinant_type"], x["determinant"]))
    validate_catalogue(out, study, config)
    return out


def validate_catalogue(rows: list[dict], study: set, config: dict) -> None:
    organism = config["catalogue"]["organism"]
    seen = {}
    for r in rows:
        key = (r["antimicrobial"], r["determinant"])
        if r["organism"] != organism:
            raise CatalogueError(f"{key}: organism {r['organism']}")
        if r["antimicrobial"] not in study:
            raise CatalogueError(f"{key}: not a study antimicrobial")
        if r["determinant_type"] not in ("gene", "point_mutation"):
            raise CatalogueError(f"{key}: determinant type {r['determinant_type']}")
        if r["mechanism_class"] not in MECHANISM_CLASSES:
            raise CatalogueError(f"{key}: no mechanism class for {r['determinant']}")
        if r["determinant_type"] == "point_mutation" and (r["locus_type"] == NA or r["ref_allele"] == NA):
            raise CatalogueError(f"{key}: point mutation without an AMRFinderPlus locus")
        if r["catalogue_version"] != config["catalogue"]["version"]:
            raise CatalogueError(f"{key}: catalogue version {r['catalogue_version']}")
        if key in seen:
            if seen[key]["interpretation"] != r["interpretation"] or seen[key]["requires"] != r["requires"]:
                raise CatalogueError(f"{key}: conflicting definitions")
            raise CatalogueError(f"{key}: duplicate rule")
        seen[key] = r
    order = [(r["antimicrobial"], r["determinant_type"], r["determinant"]) for r in rows]
    if order != sorted(order):
        raise CatalogueError("catalogue rows are not in canonical order")


def load_catalogue(config: dict) -> list[dict]:
    path = resolve(config["catalogue"]["path"])
    if pop.sha256(path) != config["catalogue"]["sha256"]:
        raise CatalogueError("mechanism catalogue differs from its recorded checksum")
    return read_tsv(path)


def catalogue_index(catalogue: list[dict]) -> dict:
    """antimicrobial -> determinant -> row; plus the point-mutation symbols used as indicators."""
    by_drug = defaultdict(dict)
    for r in catalogue:
        by_drug[r["antimicrobial"]][r["determinant"]] = r
    return by_drug


def indicator_symbols(catalogue: list[dict], interpretations: list[str]) -> list[str]:
    return sorted({r["determinant"] for r in catalogue
                   if r["determinant_type"] == "point_mutation" and r["interpretation"] in interpretations})


# ---------------------------------------------------------------- drug groups

def drug_groups(drops: list[dict], rule: dict) -> dict:
    """Group drugs from the PR #14 pooled drops; fixed before any mechanism table is read."""
    best = {}
    for r in drops:
        if (r["cohort"], r["comparator"], r["comparison"], r["level"]) != \
                (rule["cohort"], rule["comparator"], "same_test_isolates", "pooled"):
            continue
        if r["holdout_strategy"] not in HOLDOUTS:
            continue
        delta = r[f"delta_{rule['metric']}"]
        if delta == NA:
            continue
        cand = (float(delta), r["holdout_strategy"], float(r[f"delta_{rule['metric']}_ci_lower"]))
        d = r["antimicrobial"]
        if d not in best or cand[0] > best[d][0]:
            best[d] = cand
    groups = {}
    for d, (delta, strategy, lower) in sorted(best.items()):
        if delta >= rule["high_loss_min_delta"] and (lower > 0 or not rule["high_loss_requires_ci_above_zero"]):
            g = "high_loss"
        elif delta < rule["stable_max_delta"]:
            g = "stable"
        else:
            g = "intermediate"
        groups[d] = {"group": g, "max_delta": delta, "max_delta_strategy": strategy, "max_delta_ci_lower": lower}
    return groups


# ---------------------------------------------------------------- AMRFinderPlus calls

CALL_FIELDS = ("isolate_id assembly_accession element_symbol element_type element_subtype method contig start "
               "stop strand amrfinder_class amrfinder_subclass target_length reference_length coverage identity "
               "closest_reference").split()


def read_calls(call_dir: Path, iid: str, manifest_row: dict, config: dict) -> tuple[str, list[dict]]:
    d = Path(call_dir) / iid
    prov_path, report = d / "provenance.json", d / "amrfinder.tsv"
    if not prov_path.is_file() or not report.is_file():
        return "calls_missing", []
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    up = config["upstream"]
    if prov.get("caller_version") != up["amrfinder_caller_version"]:
        return "caller_version_mismatch", []
    files = prov.get("database_files_sha256", {})
    if any(files.get(k) != v for k, v in up["amrfinder_database_files_sha256"].items()):
        return "database_mismatch", []
    if prov.get("input_sha256") != manifest_row["genome_sha256"]:
        return "genome_hash_mismatch", []
    if prov.get("status") != "success":
        return "call_not_successful", []
    if prov.get("raw_sha256") != pop.sha256(report):
        return "report_hash_mismatch", []
    rows = []
    for r in read_tsv(report):
        if r["Type"] != "AMR":
            continue
        rows.append({"isolate_id": iid, "assembly_accession": manifest_row["assembly_accession"],
                     "element_symbol": r["Element symbol"], "element_type": r["Type"],
                     "element_subtype": r["Subtype"], "method": r["Method"], "contig": r["Contig id"],
                     "start": r["Start"], "stop": r["Stop"], "strand": r["Strand"], "amrfinder_class": r["Class"],
                     "amrfinder_subclass": r["Subclass"], "target_length": r["Target length"],
                     "reference_length": r["Reference sequence length"],
                     "coverage": r["% Coverage of reference"], "identity": r["% Identity to reference"],
                     "closest_reference": r["Closest reference accession"]})
    return "success", rows


def known_mutation_rows(calls: list[dict], cat_by_drug: dict, config: dict) -> list[dict]:
    up = config["upstream"]
    by_symbol = defaultdict(list)
    for drug, dets in cat_by_drug.items():
        for sym, r in dets.items():
            if r["determinant_type"] == "point_mutation":
                by_symbol[sym].append(r)
    out = []
    for c in calls:
        if c["element_subtype"] != "POINT":
            continue
        base = {"isolate_id": c["isolate_id"], "gene": c["element_symbol"].split("_")[0],
                "mutation": c["element_symbol"], "method": c["method"], "caller": "AMRFinderPlus",
                "caller_version": up["amrfinder_caller_version"], "database_version": up["amrfinder_database_version"]}
        rules = by_symbol.get(c["element_symbol"])
        if not rules:
            out.append({**base, "antimicrobial": NA, "mutation_class": NA, "evidence_source": "AMRFinderPlus",
                        "status": "not_in_catalogue"})
        for r in rules or []:
            out.append({**base, "antimicrobial": r["antimicrobial"], "mutation_class": r["mechanism_class"],
                        "evidence_source": f"{r['evidence_source']}; {r['evidence_reference']}",
                        "status": r["interpretation"]})
    out.sort(key=lambda x: (x["isolate_id"], x["mutation"], x["antimicrobial"]))
    return out


# ---------------------------------------------------------------- representation audit

class FeatureSpace:
    """Membership tests against the PR #14 k-mer universe and kept features."""

    def __init__(self, hashes: np.ndarray, pattern: np.ndarray, scaled: int, ksize: int, seed: int):
        self.hashes, self.pattern = np.asarray(hashes, dtype=np.uint64), np.asarray(pattern)
        self.max_hash, self.ksize, self.seed = gb.max_hash(scaled), ksize, seed
        import sourmash
        self._mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, seed=seed)

    def kmer_hash(self, kmer: str) -> int | None:
        if len(kmer) != self.ksize or set(kmer) - set("ACGT"):
            return None
        return int(self._mh.seq_to_hashes(kmer, force=True)[0])

    def column(self, h: int) -> int | None:
        i = int(np.searchsorted(self.hashes, np.uint64(h)))
        return int(self.pattern[i]) if i < self.hashes.shape[0] and int(self.hashes[i]) == h else None

    def classify(self, kmers: list[str]) -> dict:
        valid, sampled, cols = 0, 0, set()
        for k in kmers:
            h = self.kmer_hash(k)
            if h is None:
                continue
            valid += 1
            if h < self.max_hash:
                sampled += 1
                c = self.column(h)
                if c is not None:
                    cols.add(c)
        return {"kmers": valid, "sampled": sampled, "kept": len(cols), "columns": sorted(cols)}


def overlapping_kmers(seq: str, start: int, end: int, k: int) -> list[str]:
    """k-mers of seq overlapping [start, end); an empty span means the junction before `start`."""
    if end > start:
        first, last = start - k + 1, end - 1
    else:
        first, last = start - k + 1, start - 1
    return [seq[j:j + k] for j in range(max(0, first), min(last, len(seq) - k) + 1)]


def gene_capture(contigs: dict, call: dict, space: FeatureSpace) -> dict:
    seq = contigs.get(call["contig"])
    try:
        start, stop = int(call["start"]), int(call["stop"])
    except ValueError:
        seq = None
    if seq is None or not 1 <= start <= stop <= len(seq):
        return {"status": "representation_unknown", "detail": "locus_not_found", "kmers": 0, "sampled": 0,
                "kept": 0, "columns": []}
    region = seq[start - 1:stop]
    k = space.ksize
    res = space.classify([region[j:j + k] for j in range(len(region) - k + 1)])
    if res["kmers"] == 0:
        status = "representation_unknown"
    elif res["kept"]:
        status = "captured"
    elif res["sampled"]:
        status = "partially_represented"
    else:
        status = "not_represented"
    return {**res, "status": status, "detail": ""}


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def find_window(text: str, pattern: str, fixed: range, max_mm: int) -> int | None:
    """Unique offset where pattern matches text, `fixed` positions exactly, others with <= max_mm mismatches."""
    best, hits = None, []
    for o in range(len(text) - len(pattern) + 1):
        if any(text[o + i] != pattern[i] for i in fixed):
            continue
        mm = sum(1 for i in range(len(pattern)) if i not in fixed and text[o + i] != pattern[i])
        if mm > max_mm:
            continue
        if best is None or mm < best:
            best, hits = mm, [o]
        elif mm == best:
            hits.append(o)
    return hits[0] if len(hits) == 1 else None


def oriented_locus(seq: str, start: int, stop: int, strand: str, flank: int) -> tuple[str, int, int]:
    """Locus plus flanks in coding orientation; returns (sequence, offset of locus start, locus length)."""
    a, b = max(0, start - 1 - flank), min(len(seq), stop + flank)
    chunk = seq[a:b]
    if strand == "-":
        return revcomp(chunk), b - stop, stop - start + 1
    return chunk, start - 1 - a, stop - start + 1


def allele_status(mut: dict, wt: dict) -> str:
    if mut["kmers"] == 0 and wt["kmers"] == 0:
        return "representation_unknown"
    if mut["kept"]:
        return "captured"
    if mut["sampled"] or wt["kept"]:
        return "partially_represented"
    return "not_represented"


def protein_site(locus: str, offset: int, length: int, ref_protein: str, pos: int, ref_aa: str, alt_aa: str,
                 cfg: dict) -> dict | None:
    coding = locus[offset:offset + length]
    prot = translate(coding)
    w = cfg["protein_window"]
    lo, hi = max(0, pos - 1 - w), min(len(ref_protein), pos + w)
    if not 0 < pos <= len(ref_protein) or ref_protein[pos - 1] != ref_aa:
        return None
    pattern = ref_protein[lo:pos - 1] + alt_aa + ref_protein[pos:hi]
    centre = pos - 1 - lo
    o = find_window(prot, pattern, range(centre, centre + 1), cfg["protein_max_mismatch"])
    if o is None:
        return None
    codon_at = offset + 3 * (o + centre)
    return {"codon_start": codon_at, "codon": locus[codon_at:codon_at + 3]}


def protein_mutation_capture(locus: str, site: dict, ref_aa: str, space: FeatureSpace) -> dict:
    codon, at = site["codon"], site["codon_start"]
    wts = [c for c in CODON if CODON[c] == ref_aa]
    if not wts or CODON.get(codon) is None:
        return {"status": "representation_unknown", "detail": "codon_unresolved"}
    d = min(hamming(codon, c) for c in wts)
    nearest = [c for c in wts if hamming(codon, c) == d]
    changed = sorted({i for c in nearest for i in range(3) if c[i] != codon[i]})
    if not changed:
        return {"status": "representation_unknown", "detail": "no_changed_base"}
    s, e = at + changed[0], at + changed[-1] + 1
    k = space.ksize
    mut = space.classify(overlapping_kmers(locus, s, e, k))
    wt_kmers = []
    for c in nearest:
        w = locus[:at] + c + locus[at + 3:]
        wt_kmers += overlapping_kmers(w, s, e, k)
    wt = space.classify(wt_kmers)
    return {"status": allele_status(mut, wt), "detail": "wt_codon_ambiguous" if len(nearest) > 1 else "",
            "mutant": mut, "wild_type": wt}


def dna_mutation_capture(seq: str, call: dict, ref_dna: str, pos: int, ref: str, alt: str, space: FeatureSpace,
                         cfg: dict) -> dict:
    start, stop = int(call["start"]), int(call["stop"])
    margin = cfg["dna_search_margin"]
    a, b = max(0, start - 1 - margin), min(len(seq), stop + margin)
    w = cfg["dna_window"]
    if not 0 < pos <= len(ref_dna) or ref_dna[pos - 1:pos - 1 + len(ref)] != ref:
        return {"status": "representation_unknown", "detail": "reference_mismatch"}
    left, right = ref_dna[max(0, pos - 1 - w):pos - 1], ref_dna[pos - 1 + len(ref):pos - 1 + len(ref) + w]
    pattern = left + alt + right
    fixed = range(len(left), len(left) + len(alt))
    for strand, text in (("+", seq[a:b]), ("-", revcomp(seq[a:b]))):
        o = find_window(text, pattern, fixed, cfg["dna_max_mismatch"])
        if o is not None:
            break
    else:
        return {"status": "representation_unknown", "detail": "site_not_found"}
    # shared prefix of ref and alt is not a changed base
    common = 0
    while common < min(len(ref), len(alt)) and ref[common] == alt[common]:
        common += 1
    site = o + len(left)
    mut_span = (site + common, site + len(alt))
    wt_text = text[:site] + ref + text[site + len(alt):]
    wt_span = (site + common, site + len(ref))
    k = space.ksize
    mut = space.classify(overlapping_kmers(text, *mut_span, k))
    wt = space.classify(overlapping_kmers(wt_text, *wt_span, k))
    return {"status": allele_status(mut, wt), "detail": "", "mutant": mut, "wild_type": wt}


def mutation_capture(contigs: dict, call: dict, refs: dict, space: FeatureSpace, cfg: dict) -> dict:
    sym = call["element_symbol"]
    parsed, locus = parse_mutation(sym), refs["loci"].get(sym)
    seq = contigs.get(call["contig"])
    if parsed is None or locus is None:
        return {"status": "representation_unknown", "detail": "not_a_listed_substitution"}
    if seq is None:
        return {"status": "representation_unknown", "detail": "contig_not_found"}
    try:
        start, stop = int(call["start"]), int(call["stop"])
    except ValueError:
        return {"status": "representation_unknown", "detail": "no_coordinates"}
    if locus["locus_type"] == "protein":
        if len(parsed["ref"]) != 1 or len(parsed["alt"]) != 1:
            return {"status": "representation_unknown", "detail": "not_a_substitution"}
        ref_protein = refs["proteins"].get(locus["locus_reference"])
        if ref_protein is None:
            return {"status": "representation_unknown", "detail": "reference_protein_missing"}
        text, offset, length = oriented_locus(seq, start, stop, call["strand"], cfg["flank"])
        site = protein_site(text, offset, length, ref_protein, locus["locus_position"], parsed["ref"],
                            parsed["alt"], cfg)
        if site is None:
            return {"status": "representation_unknown", "detail": "site_not_found"}
        return protein_mutation_capture(text, site, parsed["ref"], space)
    ref_dna = refs["dna"].get(locus["locus_reference"])
    if ref_dna is None:
        return {"status": "representation_unknown", "detail": "reference_dna_missing"}
    ref_allele, alt_allele = parsed["ref"], parsed["alt"]
    return dna_mutation_capture(seq, call, ref_dna, locus["locus_position"], ref_allele, alt_allele, space, cfg)


def read_fasta_records(data: bytes) -> dict:
    out, name, chunks = {}, None, []
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks).upper()
    return out


def reference_sequences(db: Path, taxgroup: str, loci: dict) -> dict:
    proteins, dna = {}, {}
    wanted = {l["locus_reference"] for l in loci.values() if l["locus_type"] == "protein"}
    data = (Path(db) / "AMRProt.fa").read_bytes().decode("ascii")
    for block in data.split(">")[1:]:
        head, _, body = block.partition("\n")
        acc = head.split("|")[0]
        if acc in wanted:
            proteins[acc] = body.replace("\n", "").strip().upper()
    for name, seq in read_fasta_records((Path(db) / f"AMR_DNA-{taxgroup}.fa").read_bytes()).items():
        dna[name] = seq
    return {"loci": loci, "proteins": proteins, "dna": dna}


GENE_CAPTURE_FIELDS = ("isolate_id assembly_accession determinant method contig start stop kmers sampled kept "
                       "columns status detail").split()
MUT_CAPTURE_FIELDS = ("isolate_id assembly_accession mutation locus_type method contig start stop mutant_kmers "
                      "mutant_sampled mutant_kept mutant_columns wt_kmers wt_sampled wt_kept wt_columns status "
                      "detail").split()


def audit_isolate(row: dict, calls: list[dict], cache: Path, refs: dict, space: FeatureSpace, cfg: dict) -> dict:
    status, detail, data = gb.verify_genome(Path(cache), row)
    out = {"isolate_id": row["isolate_id"], "assembly_accession": row["assembly_accession"],
           "genome_status": status, "detail": detail, "genes": [], "mutations": []}
    if status != "verified":
        return out
    contigs = read_fasta_records(data)
    for c in calls:
        base = {"isolate_id": row["isolate_id"], "assembly_accession": row["assembly_accession"],
                "method": c["method"], "contig": c["contig"], "start": c["start"], "stop": c["stop"]}
        if c["element_subtype"] == "AMR":
            r = gene_capture(contigs, c, space)
            out["genes"].append({**base, "determinant": c["element_symbol"], "kmers": r["kmers"],
                                 "sampled": r["sampled"], "kept": r["kept"],
                                 "columns": ";".join(map(str, r["columns"])) or NA,
                                 "status": r["status"], "detail": r["detail"]})
        elif c["element_subtype"] == "POINT":
            r = mutation_capture(contigs, c, refs, space, cfg)
            m, w = r.get("mutant", {}), r.get("wild_type", {})
            locus = refs["loci"].get(c["element_symbol"], {})
            out["mutations"].append({
                **base, "mutation": c["element_symbol"], "locus_type": locus.get("locus_type", NA),
                "mutant_kmers": m.get("kmers", NA), "mutant_sampled": m.get("sampled", NA),
                "mutant_kept": m.get("kept", NA), "mutant_columns": ";".join(map(str, m.get("columns", []))) or NA,
                "wt_kmers": w.get("kmers", NA), "wt_sampled": w.get("sampled", NA), "wt_kept": w.get("kept", NA),
                "wt_columns": ";".join(map(str, w.get("columns", []))) or NA,
                "status": r["status"], "detail": r["detail"]})
    return out


_AUDIT = {}


def _audit_init(cache, db, taxgroup, loci, matrix_dir, params, cfg):
    refs = reference_sequences(Path(db), taxgroup, loci)
    space = FeatureSpace(np.load(Path(matrix_dir) / "feature_hashes.npy"),
                         np.load(Path(matrix_dir) / "feature_pattern.npy"),
                         params["scaled"], params["ksize"], params["seed"])
    _AUDIT.update(cache=Path(cache), refs=refs, space=space, cfg=cfg)


def _audit_one(args):
    row, calls = args
    return audit_isolate(row, calls, _AUDIT["cache"], _AUDIT["refs"], _AUDIT["space"], _AUDIT["cfg"])


# ---------------------------------------------------------------- error context

def best_capture(statuses) -> str:
    statuses = [s for s in statuses if s in CAPTURE_ORDER]
    return min(statuses, key=CAPTURE_ORDER.index) if statuses else NA


def split_determinants(ids: str, types: str) -> tuple[list[str], list[str]]:
    if ids in (NA, ""):
        return [], []
    genes, muts = [], []
    for d, t in zip(ids.split(";"), types.split(";")):
        (muts if t == "point_mutation" else genes).append(d)
    return sorted(genes), sorted(muts)


def endpoint_status(r: dict) -> str:
    """PR #7 determinant status of an endpoint. PR #14 kept it only for mapping-eligible endpoints and wrote NA
    for the mapping-uncertain ones, which do have calls; those are uncertain here, not genotype-unavailable."""
    if r["genotype_status"] in ("", NA, "genotype_unavailable", "calls_missing"):
        return NA
    return r["determinant_status"]


def mechanism_group(status: str, genes: list, muts: list) -> str:
    if status == "positive":
        return "both" if genes and muts else ("mutation" if muts else "acquired_gene")
    if status == "negative":
        return "no_established_determinant"
    return "unresolved"


def training_exposure(train: list[dict], dets: dict) -> tuple[Counter, Counter]:
    """Carrier counts of each determinant and of each exact combination, over all training endpoints."""
    single, combo = Counter(), Counter()
    for e in train:
        g, m = dets.get((e["isolate_id"], e["antimicrobial"]), ([], []))
        for d in g + m:
            single[d] += 1
        if g or m:
            combo[tuple(sorted(g + m))] += 1
    return single, combo


def seen_flag(items: list[str], single: Counter, minimum: int) -> str:
    if not items:
        return NA
    return "yes" if all(single[d] >= minimum for d in items) else "no"


def mechanism_novelty(status: str, genes: list, muts: list, single: Counter, minimum: int) -> str:
    if status == "negative":
        return "no_established_determinant"
    if status != "positive":
        return "unresolved"
    return "seen" if seen_flag(genes + muts, single, minimum) == "yes" else "unseen"


def error_category(phenotype: str, predicted: str, status: str, genes: list, muts: list, gene_cap: str,
                   mut_cap: str) -> str:
    if phenotype == predicted:
        return "correct"
    if phenotype == "R":
        if status == NA:
            return "R_GENOTYPE_UNAVAILABLE"
        if status == "negative":
            return "R_NO_ESTABLISHED_DET"
        if status != "positive":
            return "R_UNRESOLVED"
        caps = ([gene_cap] if genes else []) + ([mut_cap] if muts else [])
        if all(c in (NA, "representation_unknown") for c in caps):
            return "R_UNRESOLVED"
        rep = "REPRESENTED" if "captured" in caps else "NOT_REPRESENTED"
        kind = "DET_AND_MUT" if genes and muts else ("MUT" if muts else "DET")
        return f"R_KNOWN_{kind}_{rep}"
    if status == "positive":
        return "S_KNOWN_DET_PRESENT"
    if status == "negative":
        return "S_NO_KNOWN_DET"
    return "S_UNRESOLVED"


def four_way(phenotype: str, predicted: str, det_pred: str) -> str:
    if det_pred == NA:
        return NA
    ml, det = predicted == phenotype, det_pred == phenotype
    return {(True, True): "both_correct", (False, True): "determinant_only_correct",
            (True, False): "ml_only_correct", (False, False): "both_wrong"}[(ml, det)]


def usable_columns(sub: sparse.csr_matrix, rows: np.ndarray, min_present: int, min_absent: int) -> np.ndarray:
    counts = np.asarray(sub[rows].sum(axis=0)).ravel()
    return (counts >= min_present) & (rows.shape[0] - counts >= min_absent)


# ---------------------------------------------------------------- shift measures

def distribution(counter: Counter) -> dict:
    total = sum(counter.values())
    return {k: v / total for k, v in counter.items()} if total else {}


def tvd(a: Counter, b: Counter) -> float | None:
    p, q = distribution(a), distribution(b)
    if not p or not q:
        return None
    return 0.5 * sum(abs(p.get(k, 0) - q.get(k, 0)) for k in set(p) | set(q))


def jsd(a: Counter, b: Counter) -> float | None:
    """Jensen-Shannon divergence in bits (0 to 1)."""
    p, q = distribution(a), distribution(b)
    if not p or not q:
        return None
    out = 0.0
    for k in set(p) | set(q):
        pk, qk = p.get(k, 0), q.get(k, 0)
        mk = (pk + qk) / 2
        if pk:
            out += 0.5 * pk * math.log2(pk / mk)
        if qk:
            out += 0.5 * qk * math.log2(qk / mk)
    return out


# ---------------------------------------------------------------- augmentation

def indicator_block(isolates: list[str], calls_by_isolate: dict, symbols: list[str]) -> sparse.csr_matrix:
    col = {s: j for j, s in enumerate(symbols)}
    rows, cols = [], []
    for i, iid in enumerate(isolates):
        present = {c["element_symbol"] for c in calls_by_isolate.get(iid, [])
                   if c["element_subtype"] == "POINT" and c["method"] not in PARTIAL_METHODS}
        for s in sorted(present & col.keys()):
            rows.append(i)
            cols.append(col[s])
    return sparse.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(len(isolates), len(symbols)))


def permute_rows(block: sparse.csr_matrix, seed: int) -> sparse.csr_matrix:
    order = np.random.default_rng(seed).permutation(block.shape[0])
    return block[order].tocsr()


def append_columns(X: sparse.csr_matrix, block: sparse.csr_matrix) -> sparse.csr_matrix:
    """[X | block] as CSR without going through COO for the large part."""
    if X.shape[0] != block.shape[0]:
        raise ValueError("row counts differ")
    block = block.tocsr()
    x_len, b_len = np.diff(X.indptr), np.diff(block.indptr)
    indptr = np.zeros(X.shape[0] + 1, dtype=np.int64)
    np.cumsum(x_len + b_len, out=indptr[1:])
    indices = np.empty(indptr[-1], dtype=np.int32)
    x_pos = np.arange(X.indices.shape[0]) + np.repeat(indptr[:-1] - X.indptr[:-1], x_len)
    indices[x_pos] = X.indices
    b_pos = np.arange(block.indices.shape[0]) + np.repeat(indptr[:-1] + x_len - block.indptr[:-1], b_len)
    indices[b_pos] = block.indices + X.shape[1]
    data = np.ones(indices.shape[0], dtype=np.int8)
    if indptr[-1] < 2 ** 31:
        indptr = indptr.astype(np.int32)
    return sparse.csr_matrix((data, indices, indptr), shape=(X.shape[0], X.shape[1] + block.shape[1]))


def augmentation_eligibility(context_rows: list[dict], rule: dict, catalogue_by_drug: dict) -> dict:
    """Per drug, from the common-cohort random predictions (each endpoint once): carriers of a catalogue
    point mutation for that drug and the share whose mutations are captured. No phenotype is read."""
    out = {}
    for drug, dets in sorted(catalogue_by_drug.items()):
        has_rule = any(r["determinant_type"] == "point_mutation" and r["interpretation"] in rule["interpretations"]
                       for r in dets.values())
        rows = [r for r in context_rows if r["antimicrobial"] == drug]
        carriers = [r for r in rows if r["catalogue_mutations"] != NA]
        captured = sum(r["catalogue_mutation_capture"] == "captured" for r in carriers)
        frac = captured / len(carriers) if carriers else None
        if not has_rule:
            reason, ok = "no point-mutation rule in the catalogue", False
        elif len(carriers) < rule["min_carriers"]:
            reason, ok = f"{len(carriers)} carriers, below {rule['min_carriers']}", False
        elif frac >= rule["max_captured_fraction"]:
            reason, ok = f"captured fraction {frac:.3f} at or above {rule['max_captured_fraction']}", False
        else:
            reason, ok = "", True
        out[drug] = {"eligible": ok, "carriers": len(carriers), "captured_fraction": frac, "reason": reason}
    return out


def gap_rescue_ci(y, units, p_rand, p_hold, q_rand, q_hold, replicates: int, seed: int, level: float) -> dict:
    """Paired bootstrap over split units of (original gap - augmented gap) in AUROC."""
    y = np.asarray(y).astype(bool)
    vals = []
    for idx in gb.bootstrap_indices(units, replicates, seed):
        a = [gb.auroc(y[idx], np.asarray(p)[idx]) for p in (p_rand, p_hold, q_rand, q_hold)]
        if None in a:
            continue
        vals.append((a[0] - a[1]) - (a[2] - a[3]))
    if not vals:
        return {"lower": None, "upper": None}
    alpha = (1 - level) / 2
    return {"lower": float(np.quantile(vals, alpha)), "upper": float(np.quantile(vals, 1 - alpha))}


# ---------------------------------------------------------------- stages

def paths(outdir: Path) -> dict:
    o = Path(outdir)
    return {"manifest": o / "error_decomposition_manifest.json", "calls": o / "amr_calls.tsv.gz",
            "call_status": o / "call_status.tsv", "mutations": o / "known_mutations.tsv.gz",
            "audit_ckpt": o / "audit_checkpoint.jsonl", "gene_capture": o / "determinant_capture.tsv.gz",
            "mut_capture": o / "mutation_capture.tsv.gz", "genome_verification": o / "genome_verification.tsv",
            "context": o / "prediction_error_context.tsv.gz", "augment": o / "augmented"}


def stage_freeze(config: dict, config_path: Path, benchmark_dir: Path, outdir: Path) -> dict:
    up = check_upstream(config, benchmark_dir)
    load_catalogue(config)
    drops = read_tsv(Path(benchmark_dir) / "generalisation_drop.tsv")
    groups = drug_groups(drops, config["drug_groups"])
    case = sorted(d for d, g in groups.items() if g["group"] in ("high_loss", "stable")
                  and d in config["drug_groups"]["case_study_tiers_from"])
    manifest = {"analysis_version": config["analysis_version"], "builder_version": VERSION,
                "builder_sha256": pop.sha256(Path(__file__)), "config_sha256": pop.sha256(config_path),
                "catalogue_sha256": config["catalogue"]["sha256"], "upstream": up, "drug_groups": groups,
                "case_study_drugs": case, "frozen_at": now()}
    Path(outdir).mkdir(parents=True, exist_ok=True)
    write_json(paths(outdir)["manifest"], manifest)
    return manifest


def check_frozen(config: dict, config_path: Path, outdir: Path) -> dict:
    manifest = json.loads(paths(outdir)["manifest"].read_text(encoding="utf-8"))
    if manifest["config_sha256"] != pop.sha256(config_path):
        raise ValueError("analysis config changed after freezing; start a new analysis version")
    if manifest["catalogue_sha256"] != config["catalogue"]["sha256"]:
        raise ValueError("mechanism catalogue changed after freezing")
    return manifest


def stage_calls(config: dict, outdir: Path, call_dir: Path) -> dict:
    manifest = gb.read_genome_manifest(resolve(config["upstream"]["genome_manifest"]),
                                       config["upstream"]["genome_manifest_sha256"])
    cat = catalogue_index(load_catalogue(config))
    status_rows, calls = [], []
    for iid in sorted(manifest):
        st, rows = read_calls(call_dir, iid, manifest[iid], config)
        status_rows.append({"isolate_id": iid, "status": st, "calls": len(rows)})
        calls += rows
    p = paths(outdir)
    write_tsv(p["calls"], CALL_FIELDS, calls)
    write_tsv(p["call_status"], ["isolate_id", "status", "calls"], status_rows)
    muts = known_mutation_rows(calls, cat, config)
    write_tsv(p["mutations"], "isolate_id antimicrobial gene mutation mutation_class method caller caller_version "
              "database_version evidence_source status".split(), muts)
    return {"isolates": len(status_rows), "status": dict(Counter(r["status"] for r in status_rows)),
            "calls": len(calls), "known_mutation_rows": len(muts)}


def stage_audit(config: dict, outdir: Path, benchmark_dir: Path, cache: Path, workers: int) -> dict:
    p = paths(outdir)
    manifest = gb.read_genome_manifest(resolve(config["upstream"]["genome_manifest"]),
                                       config["upstream"]["genome_manifest_sha256"])
    calls = defaultdict(list)
    for c in read_tsv(p["calls"]):
        calls[c["isolate_id"]].append(c)
    done = {}
    if p["audit_ckpt"].exists():
        for line in p["audit_ckpt"].read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            done[rec["isolate_id"]] = rec
    db = resolve(config["upstream"]["amrfinder_database"])
    taxgroup = config["upstream"]["amrfinder_taxgroup"]
    loci = mutation_loci(db, taxgroup)
    gcfg = json.loads((Path(benchmark_dir) / "matrix" / "feature_manifest.json").read_text(encoding="utf-8"))
    params = gcfg["params"]
    todo = [(manifest[i], calls.get(i, [])) for i in sorted(manifest) if i not in done]
    init = (str(cache), str(db), taxgroup, loci, str(Path(benchmark_dir) / "matrix"), params, config["capture"])
    with p["audit_ckpt"].open("a", encoding="utf-8") as ckpt:
        if workers <= 1:
            _audit_init(*init)
            results = map(_audit_one, todo)
        else:
            pool = ProcessPoolExecutor(max_workers=workers, initializer=_audit_init, initargs=init)
            results = pool.map(_audit_one, todo, chunksize=8)
        for n, rec in enumerate(results, 1):
            ckpt.write(json.dumps(rec, sort_keys=True) + "\n")
            ckpt.flush()
            done[rec["isolate_id"]] = rec
            if n % 250 == 0:
                print(f"{now()} audited {len(done)}/{len(manifest)}", flush=True)
        if workers > 1:
            pool.shutdown()
    recs = [done[i] for i in sorted(done)]

    def na_detail(rows):
        return [{**r, "detail": r["detail"] or NA} for r in rows]

    write_tsv(p["genome_verification"], "isolate_id assembly_accession genome_status detail".split(), na_detail(recs))
    write_tsv(p["gene_capture"], GENE_CAPTURE_FIELDS, na_detail(g for r in recs for g in r["genes"]))
    write_tsv(p["mut_capture"], MUT_CAPTURE_FIELDS, na_detail(m for r in recs for m in r["mutations"]))
    return {"isolates": len(recs), "genome_status": dict(Counter(r["genome_status"] for r in recs)),
            "gene_calls": sum(len(r["genes"]) for r in recs), "mutation_calls": sum(len(r["mutations"]) for r in recs),
            "gene_status": dict(Counter(g["status"] for r in recs for g in r["genes"])),
            "mutation_status": dict(Counter(m["status"] for r in recs for m in r["mutations"]))}


def load_benchmark(config: dict, benchmark_dir: Path, population: Path, baseline: Path):
    gconfig = resolve(config["upstream"]["generalisation_config"])
    inputs = gb.load_inputs(gconfig, population, baseline, Path(benchmark_dir) / "genome_verification.tsv")
    manifest = json.loads((Path(benchmark_dir) / "benchmark_manifest.json").read_text(encoding="utf-8"))
    tasks, skipped = gb.build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
    gb.check_frozen(tasks, manifest, gconfig)
    return inputs, tasks


def phenotype_provenance(ast_dir: Path, keys: set) -> dict:
    """Per (isolate, antimicrobial): observation count, conflict, methods and breakpoint standards."""
    out = {}
    for r in read_tsv(Path(ast_dir) / "isolate_drug_endpoints.tsv"):
        k = (r["isolate_id"], r["antimicrobial"])
        if k in keys:
            out[k] = {"n_ast_observations": r["n_ast_observations"], "duplicate_status": r["duplicate_status"],
                      "conflict_status": r["conflict_status"], "methods": set(), "standards": set(),
                      "mic_reported": False}
    with open(Path(ast_dir) / "ast_observations.tsv", encoding="utf-8", newline="") as handle:
        csv.field_size_limit(sys.maxsize)
        reader = csv.reader(handle, delimiter="\t")
        head = next(reader)
        ix = {n: head.index(n) for n in ("isolate_id", "antimicrobial", "testing_method", "breakpoint_standard",
                                         "breakpoint_version", "mic")}
        for row in reader:
            k = (row[ix["isolate_id"]], row[ix["antimicrobial"]])
            rec = out.get(k)
            if rec is None:
                continue
            rec["methods"].add(row[ix["testing_method"]] or NA)
            std = row[ix["breakpoint_standard"]] or NA
            ver = row[ix["breakpoint_version"]]
            rec["standards"].add(f"{std} {ver}".strip() if ver not in ("", NA) else std)
            rec["mic_reported"] |= row[ix["mic"]] not in ("", NA)
    return out


def phenotype_flags(rec: dict | None, corrected: bool) -> dict:
    if rec is None:
        return {"phenotype_flags": "provenance_missing", "testing_method": NA, "breakpoint_standard": NA,
                "n_ast_observations": NA}
    flags = []
    if corrected:
        flags.append("label_corrected")
    if rec["n_ast_observations"] not in ("1", NA):
        flags.append("multiple_observations")
    if rec["conflict_status"] not in ("none", NA, ""):
        flags.append("conflict")
    methods = sorted(m for m in rec["methods"])
    standards = sorted(s for s in rec["standards"])
    if not methods or methods == [NA]:
        flags.append("method_missing")
    if not standards or standards == [NA]:
        flags.append("breakpoint_standard_missing")
    if len(standards) > 1:
        flags.append("mixed_breakpoint_standards")
    return {"phenotype_flags": ";".join(flags) or "none", "testing_method": "|".join(methods) or NA,
            "breakpoint_standard": "|".join(standards) or NA, "n_ast_observations": rec["n_ast_observations"]}


CONTEXT_FIELDS = (
    "benchmark_version cohort antimicrobial split_strategy fold_or_holdout_id task_id isolate_id assembly_accession "
    "genome_sha256 feature_version phenotype predicted_probability predicted_class correct error_type "
    "error_category determinant_prediction four_way established_determinant_status determinant_ids "
    "gene_determinants mutation_determinants mechanism_group gene_capture mutation_capture gene_usable_in_task "
    "mutation_usable_in_task catalogue_mutations catalogue_mutation_capture determinant_seen_in_training "
    "mutation_seen_in_training determinant_combination_seen_in_training mechanism_novelty "
    "lineage_seen_in_training snp_cluster_seen_in_training country_seen_in_training project_seen_in_training "
    "year_after_training_range lineage snp_cluster lineage_component country project year source split_unit "
    "phenotype_flags testing_method breakpoint_standard n_ast_observations genome_qc_status genome_qc_flags "
    "assembly_length contig_count n50 ambiguous_fraction partial_amrfinder_calls").split()


def stage_context(config: dict, outdir: Path, benchmark_dir: Path, population: Path, baseline: Path,
                  ast_dir: Path) -> dict:
    p = paths(outdir)
    inputs, tasks = load_benchmark(config, benchmark_dir, population, baseline)
    endpoints = {(e["isolate_id"], e["antimicrobial"]): e for e in inputs["endpoints"]}
    cat_by_drug = catalogue_index(load_catalogue(config))
    gmanifest = gb.read_genome_manifest(resolve(config["upstream"]["genome_manifest"]),
                                        config["upstream"]["genome_manifest_sha256"])
    fman = json.loads((Path(benchmark_dir) / "matrix" / "feature_manifest.json").read_text(encoding="utf-8"))
    feature_version = f"{config['upstream']['benchmark_version']}:{fman['feature_hashes_sha256'][:12]}"
    dets, det_status = {}, {}
    for r in read_tsv(Path(baseline) / "genotype_phenotype_endpoints.tsv"):
        k = (r["isolate_id"], r["antimicrobial"])
        if k in endpoints:
            dets[k] = split_determinants(r["determinant_ids"], r["determinant_types"])
            det_status[k] = endpoint_status(r)
    meta = {r["isolate_id"]: r for r in pop.read_tsv(Path(population) / "population_metadata.tsv")}
    gene_cap, mut_cap, cols_of = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    for g in read_tsv(p["gene_capture"]):
        prev = gene_cap[g["isolate_id"]].get(g["determinant"])
        gene_cap[g["isolate_id"]][g["determinant"]] = best_capture([g["status"], prev])
        if g["columns"] != NA:
            cols_of[g["isolate_id"]].setdefault(g["determinant"], set()).update(map(int, g["columns"].split(";")))
    for m in read_tsv(p["mut_capture"]):
        prev = mut_cap[m["isolate_id"]].get(m["mutation"])
        mut_cap[m["isolate_id"]][m["mutation"]] = best_capture([m["status"], prev])
        if m["mutant_columns"] != NA:
            cols_of[m["isolate_id"]].setdefault(m["mutation"], set()).update(map(int, m["mutant_columns"].split(";")))
    partial = defaultdict(set)
    point_calls = defaultdict(set)
    for c in read_tsv(p["calls"]):
        if c["method"] in PARTIAL_METHODS:
            partial[c["isolate_id"]].add(c["element_symbol"])
        if c["element_subtype"] == "POINT" and c["method"] not in PARTIAL_METHODS:
            point_calls[c["isolate_id"]].add(c["element_symbol"])
    keys = set(endpoints)
    prov = phenotype_provenance(ast_dir, keys)
    # columns of interest for task-level usability
    interest = sorted({c for per in cols_of.values() for cs in per.values() for c in cs})
    col_index = {c: j for j, c in enumerate(interest)}
    X, row_of, _ = gb.load_matrix(Path(benchmark_dir) / "matrix")
    sub = X[:, interest].tocsr() if interest else None
    del X
    fcfg = inputs["config"]["features"]
    preds = load_predictions(Path(benchmark_dir) / "predictions.tsv.gz")
    reconcile(tasks, preds)
    mincarriers = config["novelty"]["min_training_carriers"]
    n_pred, categories = 0, Counter()
    writer = ContextWriter(p["context"])
    for t in tasks:
        key = (t["cohort"], t["antimicrobial"], t["strategy"], t["fold"])
        tp = preds[key]
        single, combo = training_exposure(t["train"], dets)
        train_projects = {e["project"] for e in t["train"]}
        tr_rows = np.array([row_of[e["isolate_id"]] for e in t["train"]], dtype=np.int64)
        usable = usable_columns(sub, tr_rows, fcfg["min_train_present"], fcfg["min_train_absent"]) \
            if sub is not None else None
        drug = t["antimicrobial"]
        cat = cat_by_drug.get(drug, {})
        for e in sorted(t["test"], key=lambda x: x["isolate_id"]):
            iid = e["isolate_id"]
            pr = tp[iid]
            n_pred += 1
            k = (iid, drug)
            status = det_status.get(k, NA)
            genes, muts = dets.get(k, ([], []))
            gcap = best_capture(gene_cap[iid].get(d, "representation_unknown") for d in genes) if genes else NA
            mcap = best_capture(mut_cap[iid].get(d, "representation_unknown") for d in muts) if muts else NA

            def use(items):
                if not items:
                    return NA
                cs = [col_index[c] for d in items for c in cols_of[iid].get(d, ())]
                return "yes" if cs and usable is not None and bool(usable[cs].any()) else "no"

            catmuts = sorted(s for s in point_calls[iid] if s in cat and cat[s]["determinant_type"] == "point_mutation"
                             and cat[s]["interpretation"] in config["augmentation"]["interpretations"])
            det_pred = {"positive": "R", "negative": "S"}.get(status, NA)
            m = meta[iid]
            flags = phenotype_flags(prov.get(k), e["phenotype_correction"] not in (NA, "", "none"))
            row = {
                "benchmark_version": config["upstream"]["benchmark_version"], "cohort": t["cohort"],
                "antimicrobial": drug, "split_strategy": t["strategy"], "fold_or_holdout_id": t["fold"],
                "task_id": t["task_id"], "isolate_id": iid, "assembly_accession": e["assembly_accession"],
                "genome_sha256": gmanifest[iid]["genome_sha256"], "feature_version": feature_version,
                "phenotype": e["phenotype"], "predicted_probability": pr["predicted_probability"],
                "predicted_class": pr["predicted_class"], "correct": pr["correct"],
                "error_type": {("R", "R"): "TP", ("S", "S"): "TN", ("R", "S"): "FN", ("S", "R"): "FP"}[
                    (e["phenotype"], pr["predicted_class"])],
                "error_category": error_category(e["phenotype"], pr["predicted_class"], status, genes, muts, gcap,
                                                 mcap),
                "determinant_prediction": det_pred, "four_way": four_way(e["phenotype"], pr["predicted_class"],
                                                                         det_pred),
                "established_determinant_status": status, "determinant_ids": ";".join(genes + muts) or NA,
                "gene_determinants": ";".join(genes) or NA, "mutation_determinants": ";".join(muts) or NA,
                "mechanism_group": mechanism_group(status, genes, muts), "gene_capture": gcap,
                "mutation_capture": mcap, "gene_usable_in_task": use(genes), "mutation_usable_in_task": use(muts),
                "catalogue_mutations": ";".join(catmuts) or NA,
                "catalogue_mutation_capture": best_capture(mut_cap[iid].get(s, "representation_unknown")
                                                           for s in catmuts) if catmuts else NA,
                "determinant_seen_in_training": seen_flag(genes + muts, single, mincarriers),
                "mutation_seen_in_training": seen_flag(muts, single, mincarriers),
                "determinant_combination_seen_in_training": NA if not (genes or muts) else
                ("yes" if combo[tuple(sorted(genes + muts))] >= config["novelty"]["combination_min_training_carriers"]
                 else "no"),
                "mechanism_novelty": mechanism_novelty(status, genes, muts, single, mincarriers),
                "lineage_seen_in_training": pr["lineage_seen_in_training"],
                "snp_cluster_seen_in_training": pr["snp_cluster_seen_in_training"],
                "country_seen_in_training": pr["country_seen_in_training"],
                "project_seen_in_training": "yes" if e["project"] in train_projects else "no",
                "year_after_training_range": pr["year_after_training_range"], "lineage": e["lineage_id"],
                "snp_cluster": e["snp_cluster"], "lineage_component": e["lineage_component"],
                "country": e["country"], "project": e["project"], "year": e["year"], "source": e["source"],
                "split_unit": e["split_unit"], **flags,
                "genome_qc_status": m["genome_qc_status"], "genome_qc_flags": m["genome_qc_flags"] or NA,
                "assembly_length": m["assembly_length"], "contig_count": m["contig_count"], "n50": m["n50"],
                "ambiguous_fraction": m["ambiguous_fraction"],
                "partial_amrfinder_calls": "yes" if partial[iid] & cat.keys() else "no"}
            writer.write(row)
            categories[row["error_category"]] += 1
    writer.close()
    return {"prediction_rows": n_pred, "tasks": len(tasks), "error_categories": dict(categories)}


PRED_KEEP = ("predicted_probability", "predicted_class", "correct", "lineage_seen_in_training",
             "snp_cluster_seen_in_training", "country_seen_in_training", "year_after_training_range")


def load_predictions(path: Path) -> dict:
    """PR #14 predictions by task key and isolate; a repeated (task, isolate) row stops the analysis."""
    preds = defaultdict(dict)
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            key = (r["cohort"], r["antimicrobial"], r["split_strategy"], r["fold_or_holdout_id"])
            if r["isolate_id"] in preds[key]:
                raise UpstreamError(f"duplicate prediction for {r['isolate_id']} in {key}")
            preds[key][r["isolate_id"]] = {k: sys.intern(r[k]) for k in PRED_KEEP}
    return preds


def reconcile(tasks: list[dict], preds: dict) -> None:
    """Every frozen test endpoint has exactly one prediction and no prediction lies outside the tasks."""
    drugs = {t["antimicrobial"] for t in tasks}
    extra = sorted({k[1] for k in preds} - drugs)
    if extra:
        raise UpstreamError(f"predictions for antimicrobials outside the benchmark: {extra}")
    keys = set()
    for t in tasks:
        key = (t["cohort"], t["antimicrobial"], t["strategy"], t["fold"])
        keys.add(key)
        got, want = set(preds.get(key, {})), {e["isolate_id"] for e in t["test"]}
        if got != want:
            missing, unexpected = sorted(want - got), sorted(got - want)
            raise UpstreamError(f"{t['task_id']}: {len(missing)} test isolates without a prediction "
                                f"{missing[:3]}, {len(unexpected)} predictions outside the test set {unexpected[:3]}")
    stray = sorted(set(preds) - keys)
    if stray:
        raise UpstreamError(f"predictions for tasks that were not frozen: {stray[:3]}")


def stage_augment(config: dict, outdir: Path, benchmark_dir: Path, population: Path, baseline: Path,
                  variants: list[str], only: str | None) -> dict:
    p = paths(outdir)
    manifest = check_frozen(config, Path(config["_path"]), outdir)
    inputs, tasks = load_benchmark(config, benchmark_dir, population, baseline)
    catalogue = load_catalogue(config)
    aug = config["augmentation"]
    ctx = read_context(p["context"], keep=lambda r: r["cohort"] == aug["cohort"] and r["split_strategy"] == "random")
    elig = augmentation_eligibility(ctx, aug, catalogue_index(catalogue))
    drugs = sorted(d for d, v in elig.items() if v["eligible"])
    symbols = indicator_symbols(catalogue, aug["interpretations"])
    calls = defaultdict(list)
    for c in read_tsv(p["calls"]):
        calls[c["isolate_id"]].append(c)
    X, row_of, _ = gb.load_matrix(Path(benchmark_dir) / "matrix")
    isolates = sorted(row_of, key=row_of.get)
    block = indicator_block(isolates, calls, symbols)
    gconfig = inputs["config"]
    counts = Counter()
    p["augment"].mkdir(parents=True, exist_ok=True)
    write_json(p["augment"] / "augmentation_plan.json", {
        "eligibility": elig, "drugs": drugs, "indicators": symbols, "variants": variants,
        "indicator_prevalence": {s: int(block[:, j].sum()) for j, s in enumerate(symbols)},
        "analysis_manifest_frozen_at": manifest["frozen_at"], "written_at": now()})
    for variant in variants:
        b = block if variant == "mutation_augmented" else permute_rows(block, aug["control_seed"])
        Xa = append_columns(X.tocsr(), b)
        todo = [t for t in tasks if t["cohort"] == aug["cohort"] and t["antimicrobial"] in drugs
                and (only is None or re.search(only, t["task_id"]))]
        vdir = p["augment"] / variant / "tasks"
        vdir.mkdir(parents=True, exist_ok=True)
        for t in todo:
            st = vdir / f"{t['task_id']}.json"
            if st.exists():
                prev = json.loads(st.read_text(encoding="utf-8"))
                if prev.get("status") == "success" and prev.get("split_sha256") == gb.task_digest(t):
                    counts[f"{variant}:reused"] += 1
                    continue
            slim = {**t, "train": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["train"]],
                    "test": [{k: e[k] for k in ("isolate_id", "phenotype", "split_unit")} for e in t["test"]]}
            try:
                status, preds = gb.run_task(slim, Xa, row_of, gconfig)
            except gb.SingleClassError as error:
                status, preds = {"task_id": t["task_id"], "status": "failed", "reason": str(error)}, []
            status.update(variant=variant, antimicrobial=t["antimicrobial"], split_strategy=t["strategy"],
                          fold_or_holdout_id=t["fold"])
            status["indicator_columns"] = len(symbols)
            if status["status"] == "success":
                write_tsv(vdir / f"{t['task_id']}.tsv", gb.PRED_FIELDS, preds)
            write_json(st, {**status, "finished_at": now()})
            counts[f"{variant}:{status['status']}"] += 1
            print(f"{now()} {variant} {t['task_id']} {status['status']} {status.get('seconds', '')}s", flush=True)
        del Xa
    return {"drugs": drugs, "counts": dict(counts)}


# ---------------------------------------------------------------- report

def frac(n, d):
    return n / d if d else None


def support_ok(y, config) -> bool:
    s = config["support"]
    y = np.asarray(y)
    return int(y.sum()) >= s["min_R"] and int((~y.astype(bool)).sum()) >= s["min_S"]


def perf(rows: list[dict], config: dict, prob_key="predicted_probability", class_key="predicted_class") -> dict:
    y = np.array([r["phenotype"] == "R" for r in rows], dtype=bool)
    p = np.array([float(r[prob_key]) for r in rows])
    pred = np.array([r[class_key] == "R" for r in rows], dtype=bool)
    out = {"N": len(y), "N_R": int(y.sum()), "N_S": int((~y).sum()), **gb.confusion_metrics(y, pred)}
    out["AUROC"], out["AUPRC"] = gb.auroc(y, p), gb.auprc(y, p)
    out["supported"] = "yes" if support_ok(y, config) else "no"
    return out


PERF_FIELDS = "N N_R N_S TP FP TN FN AUROC AUPRC sensitivity specificity balanced_accuracy MCC supported".split()


def stage_report(config: dict, outdir: Path, benchmark_dir: Path, population: Path, baseline: Path) -> dict:
    p = paths(outdir)
    manifest = check_frozen(config, Path(config["_path"]), outdir)
    ctx = read_context(p["context"])
    catalogue = load_catalogue(config)
    cat_by_drug = catalogue_index(catalogue)
    groups = manifest["drug_groups"]
    o = Path(outdir)
    out = {}

    write_tsv(o / "mechanism_catalogue_snapshot.tsv", CATALOGUE_FIELDS,
              [r for r in catalogue if r["interpretation"] in ("predicts_resistance", "component", "uncertain")])

    # representation audit by determinant and antimicrobial, over common-cohort endpoints (each once)
    base = [r for r in ctx if r["cohort"] == "common" and r["split_strategy"] == "random"]
    gcap = {(g["isolate_id"], g["determinant"]): g["status"] for g in read_tsv(p["gene_capture"])}
    mcap = {(m["isolate_id"], m["mutation"]): m for m in read_tsv(p["mut_capture"])}
    audit = []
    for r in base:
        for kind, field in (("gene", "gene_determinants"), ("point_mutation", "catalogue_mutations")):
            if r[field] == NA:
                continue
            for d in r[field].split(";"):
                st = gcap.get((r["isolate_id"], d)) if kind == "gene" else \
                    (mcap.get((r["isolate_id"], d)) or {}).get("status")
                audit.append({"antimicrobial": r["antimicrobial"], "determinant": d, "kind": kind,
                              "status": st or "representation_unknown", "phenotype": r["phenotype"]})
    agg = defaultdict(Counter)
    for a in audit:
        agg[(a["antimicrobial"], a["kind"], a["determinant"])][a["status"]] += 1
        agg[(a["antimicrobial"], a["kind"], "ALL")][a["status"]] += 1
    rep_rows = []
    for (drug, kind, det), c in sorted(agg.items()):
        n = sum(c.values())
        cls = cat_by_drug.get(drug, {}).get(det, {}).get("mechanism_class", NA) if det != "ALL" else "ALL"
        rep_rows.append({"antimicrobial": drug, "drug_group": groups.get(drug, {}).get("group", NA),
                         "determinant_kind": kind, "determinant": det, "mechanism_class": cls, "endpoints": n,
                         **{s: c[s] for s in CAPTURE_ORDER},
                         "captured_fraction": frac(c["captured"], n)})
    write_tsv(o / "representation_audit.tsv", ["antimicrobial", "drug_group", "determinant_kind", "determinant",
                                               "mechanism_class", "endpoints", *CAPTURE_ORDER, "captured_fraction"],
              rep_rows)

    # mutation and determinant capture, one row per call (the k-mer level evidence)
    capture_cfg = config["capture"]
    write_tsv(o / "mutation_capture.tsv", MUT_CAPTURE_FIELDS, read_tsv(p["mut_capture"]))
    write_tsv(o / "determinant_capture_summary.tsv", ["determinant", "calls", *CAPTURE_ORDER],
              [{"determinant": d, "calls": sum(c.values()), **{s: c[s] for s in CAPTURE_ORDER}}
               for d, c in sorted(_count_by(read_tsv(p["gene_capture"]), "determinant", "status").items())])

    # error decomposition
    errors = [r for r in ctx if r["error_type"] in ("FN", "FP")]
    write_tsv(o / "error_decomposition.tsv.gz", CONTEXT_FIELDS, errors)
    cats = sorted({r["error_category"] for r in errors})
    summ = []
    for key, rows in sorted(_group(ctx, ("cohort", "antimicrobial", "split_strategy")).items()) + \
            sorted(_group(ctx, ("cohort", "split_strategy"), all_drugs=True).items()):
        err = [r for r in rows if r["error_type"] in ("FN", "FP")]
        fn = [r for r in err if r["error_type"] == "FN"]
        c = Counter(r["error_category"] for r in err)
        summ.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2],
                     "drug_group": groups.get(key[1], {}).get("group", "ALL" if key[1] == "ALL" else NA),
                     "N": len(rows), "R": sum(r["phenotype"] == "R" for r in rows), "errors": len(err),
                     "FN": len(fn), "FP": len(err) - len(fn), **{k: c[k] for k in cats},
                     **{f"{k}_fraction_of_FN": frac(c[k], len(fn)) for k in cats if k.startswith("R_")}})
    fields = ["cohort", "antimicrobial", "split_strategy", "drug_group", "N", "R", "errors", "FN", "FP", *cats,
              *[f"{k}_fraction_of_FN" for k in cats if k.startswith("R_")]]
    write_tsv(o / "error_decomposition_summary.tsv", fields, summ)

    # four-way comparison and the informative subsets
    fw = []
    for key, rows in sorted(_group(ctx, ("cohort", "antimicrobial", "split_strategy")).items()) + \
            sorted(_group(ctx, ("cohort", "split_strategy"), all_drugs=True).items()):
        c = Counter(r["four_way"] for r in rows if r["four_way"] != NA)
        n = sum(c.values())
        fw.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "N": n,
                   **{k: c[k] for k in ("both_correct", "determinant_only_correct", "ml_only_correct", "both_wrong")},
                   **{f"{k}_fraction": frac(c[k], n) for k in ("determinant_only_correct", "ml_only_correct",
                                                               "both_wrong")}})
    write_tsv(o / "four_way_comparison.tsv",
              ["cohort", "antimicrobial", "split_strategy", "N", "both_correct", "determinant_only_correct",
               "ml_only_correct", "both_wrong", "determinant_only_correct_fraction", "ml_only_correct_fraction",
               "both_wrong_fraction"], fw)
    write_tsv(o / "determinant_correct_ml_wrong.tsv.gz", CONTEXT_FIELDS,
              [r for r in ctx if r["four_way"] == "determinant_only_correct"])
    write_tsv(o / "ml_correct_determinant_wrong.tsv.gz", CONTEXT_FIELDS,
              [r for r in ctx if r["four_way"] == "ml_only_correct"])
    dc = []
    for key, rows in sorted(_group([r for r in ctx if r["four_way"] == "determinant_only_correct"],
                                   ("cohort", "antimicrobial", "split_strategy")).items()):
        c = Counter((r["phenotype"], r["mechanism_group"], r["gene_capture"], r["mutation_capture"],
                     r["determinant_seen_in_training"]) for r in rows)
        for (ph, mg, gc, mc, seen), n in sorted(c.items()):
            dc.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "phenotype": ph,
                       "mechanism_group": mg, "gene_capture": gc, "mutation_capture": mc,
                       "determinant_seen_in_training": seen, "endpoints": n})
    write_tsv(o / "determinant_correct_ml_wrong_summary.tsv",
              ["cohort", "antimicrobial", "split_strategy", "phenotype", "mechanism_group", "gene_capture",
               "mutation_capture", "determinant_seen_in_training", "endpoints"], dc)

    # population and mechanism novelty
    nov = []
    for key, rows in sorted(_group([r for r in ctx if r["established_determinant_status"] == "positive"],
                                   ("cohort", "antimicrobial", "split_strategy")).items()) + \
            sorted(_group([r for r in ctx if r["established_determinant_status"] == "positive"],
                          ("cohort", "split_strategy"), all_drugs=True).items()):
        for pop_axis in ("lineage_seen_in_training", "country_seen_in_training", "project_seen_in_training"):
            cells = defaultdict(list)
            for r in rows:
                if r[pop_axis] in ("yes", "no") and r["mechanism_novelty"] in ("seen", "unseen"):
                    cells[(r[pop_axis], r["mechanism_novelty"])].append(r)
            for (ps, ms), cr in sorted(cells.items()):
                nov.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2],
                            "population_axis": pop_axis.replace("_seen_in_training", ""), "population_seen": ps,
                            "mechanism": ms, **perf(cr, config)})
    write_tsv(o / "population_vs_mechanism_novelty.tsv",
              ["cohort", "antimicrobial", "split_strategy", "population_axis", "population_seen", "mechanism",
               *PERF_FIELDS], nov)

    # stratified by mechanism group (R endpoints) and by mechanism novelty
    strat = []
    for key, rows in sorted(_group(ctx, ("cohort", "antimicrobial", "split_strategy")).items()) + \
            sorted(_group(ctx, ("cohort", "split_strategy"), all_drugs=True).items()):
        for axis in ("mechanism_group", "mechanism_novelty", "mutation_capture", "gene_capture"):
            for value, cr in sorted(_group(rows, (axis,)).items()):
                rr = [r for r in cr if r["phenotype"] == "R"]
                strat.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "axis": axis,
                              "value": value[0], "N": len(cr), "R": len(rr),
                              "sensitivity": frac(sum(r["predicted_class"] == "R" for r in rr), len(rr)),
                              "determinant_sensitivity": frac(sum(r["determinant_prediction"] == "R" for r in rr),
                                                              len(rr)),
                              "S": len(cr) - len(rr),
                              "specificity": frac(sum(r["predicted_class"] == "S" for r in cr if r["phenotype"] == "S"),
                                                  len(cr) - len(rr))})
    write_tsv(o / "performance_by_mechanism.tsv",
              ["cohort", "antimicrobial", "split_strategy", "axis", "value", "N", "R", "sensitivity",
               "determinant_sensitivity", "S", "specificity"], strat)

    # mechanism train/test overlap and shift, per task
    inputs, tasks = load_benchmark(config, benchmark_dir, population, baseline)
    dets = {}
    for r in read_tsv(Path(baseline) / "genotype_phenotype_endpoints.tsv"):
        if r["eligible_binary_rs"] == "true" and r["population_eligible"] == "true":
            dets[(r["isolate_id"], r["antimicrobial"])] = split_determinants(r["determinant_ids"],
                                                                             r["determinant_types"])
    overlap, prevalence = [], []
    minc = config["novelty"]["min_training_carriers"]
    for t in tasks:
        tr, ts = t["train"], t["test"]

        def mech(es, resistant_only=False, mut_only=False):
            c = Counter()
            for e in es:
                if resistant_only and e["phenotype"] != "R":
                    continue
                g, m = dets.get((e["isolate_id"], e["antimicrobial"]), ([], []))
                items = m if mut_only else g + m
                c[";".join(items) if items else "none"] += 1
            return c

        def carriers(es, mut_only=False):
            c = Counter()
            for e in es:
                g, m = dets.get((e["isolate_id"], e["antimicrobial"]), ([], []))
                for d in (m if mut_only else g + m):
                    c[d] += 1
            return c

        trc, tsc = carriers(tr), carriers(ts)
        tr_mut, ts_mut = carriers(tr, True), carriers(ts, True)
        ts_r_pos = [e for e in ts if e["phenotype"] == "R" and e["determinant_status"] == "positive"]
        unseen = sum(1 for e in ts_r_pos if any(trc[d] < minc for d in sum(dets.get((e["isolate_id"],
                                                                                      e["antimicrobial"]),
                                                                                     ([], [])), [])))
        r_tr, r_ts = sum(e["phenotype"] == "R" for e in tr), sum(e["phenotype"] == "R" for e in ts)
        lin_tr, lin_ts = Counter(e["lineage_id"] for e in tr), Counter(e["lineage_id"] for e in ts)
        overlap.append({
            "task_id": t["task_id"], "cohort": t["cohort"], "antimicrobial": t["antimicrobial"],
            "split_strategy": t["strategy"], "fold_or_holdout_id": t["fold"], "n_train": len(tr), "n_test": len(ts),
            "train_R": r_tr, "test_R": r_ts, "train_R_fraction": frac(r_tr, len(tr)),
            "test_R_fraction": frac(r_ts, len(ts)),
            "phenotype_prevalence_shift": abs(frac(r_tr, len(tr)) - frac(r_ts, len(ts))),
            "determinants_in_train": len(trc), "determinants_in_test": len(tsc),
            "determinants_test_only": sum(1 for d in tsc if trc[d] == 0),
            "determinants_test_below_min_in_train": sum(1 for d in tsc if trc[d] < minc),
            "mutations_in_train": len(tr_mut), "mutations_in_test": len(ts_mut),
            "mutations_test_only": sum(1 for d in ts_mut if tr_mut[d] == 0),
            "combinations_in_train": len(mech(tr)) - ("none" in mech(tr)),
            "combinations_in_test": len(mech(ts)) - ("none" in mech(ts)),
            "combinations_test_only": sum(1 for k in mech(ts) if k != "none" and k not in mech(tr)),
            "test_R_det_positive": len(ts_r_pos), "test_R_mechanism_unseen": unseen,
            "test_R_mechanism_unseen_fraction": frac(unseen, len(ts_r_pos)),
            "R_mechanism_tvd": tvd(mech(tr, True), mech(ts, True)),
            "R_mechanism_jsd": jsd(mech(tr, True), mech(ts, True)),
            "R_mutation_tvd": tvd(mech(tr, True, True), mech(ts, True, True)),
            "lineage_tvd": tvd(lin_tr, lin_ts),
            "test_lineage_unseen_fraction": frac(sum(1 for e in ts if e["lineage_id"] not in lin_tr), len(ts))})
        if t["strategy"] in ("country", "project", "geography", "temporal") or t["strategy"] == "lineage":
            for d in sorted(set(trc) | set(tsc)):
                a, b = frac(trc[d], len(tr)), frac(tsc[d], len(ts))
                if max(a, b) < config["novelty"]["prevalence_report_min"]:
                    continue
                cls = cat_by_drug.get(t["antimicrobial"], {}).get(d, {}).get("mechanism_class", NA)
                prevalence.append({"task_id": t["task_id"], "cohort": t["cohort"],
                                   "antimicrobial": t["antimicrobial"], "split_strategy": t["strategy"],
                                   "fold_or_holdout_id": t["fold"], "determinant": d, "mechanism_class": cls,
                                   "train_prevalence": a, "test_prevalence": b, "difference": b - a})
    ov_fields = list(overlap[0].keys())
    write_tsv(o / "mechanism_train_test_overlap.tsv", ov_fields, overlap)
    write_tsv(o / "mechanism_prevalence_by_task.tsv",
              ["task_id", "cohort", "antimicrobial", "split_strategy", "fold_or_holdout_id", "determinant",
               "mechanism_class", "train_prevalence", "test_prevalence", "difference"], prevalence)
    shift = []
    for key, rows in sorted(_group(overlap, ("cohort", "antimicrobial", "split_strategy")).items()):
        w = [r["n_test"] for r in rows]

        def wmean(f):
            vals = [(r[f], n) for r, n in zip(rows, w) if r[f] is not None]
            return sum(v * n for v, n in vals) / sum(n for _, n in vals) if vals else None

        shift.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "tasks": len(rows),
                      "test_N": sum(w), **{f: wmean(f) for f in (
                          "R_mechanism_tvd", "R_mechanism_jsd", "R_mutation_tvd", "lineage_tvd",
                          "phenotype_prevalence_shift", "test_lineage_unseen_fraction")},
                      "test_R_det_positive": sum(r["test_R_det_positive"] for r in rows),
                      "test_R_mechanism_unseen": sum(r["test_R_mechanism_unseen"] for r in rows),
                      "test_R_mechanism_unseen_fraction": frac(sum(r["test_R_mechanism_unseen"] for r in rows),
                                                               sum(r["test_R_det_positive"] for r in rows)),
                      "determinants_test_only": sum(r["determinants_test_only"] for r in rows),
                      "combinations_test_only": sum(r["combinations_test_only"] for r in rows)})
    write_tsv(o / "mechanism_shift_summary.tsv", list(shift[0].keys()), shift)

    # phenotype and genome QC: errors against correct predictions
    pq, gq = [], []
    for key, rows in sorted(_group(ctx, ("cohort", "antimicrobial", "split_strategy")).items()) + \
            sorted(_group(ctx, ("cohort", "split_strategy"), all_drugs=True).items()):
        for outcome in ("error", "correct"):
            sel = [r for r in rows if (r["error_type"] in ("FN", "FP")) == (outcome == "error")]
            n = len(sel)
            flags = Counter(f for r in sel for f in r["phenotype_flags"].split(";") if f != "none")
            pq.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "outcome": outcome,
                       "N": n, **{f"{f}_fraction": frac(flags[f], n) for f in PHENO_FLAGS},
                       "method_mix": _top(Counter(r["testing_method"] for r in sel)),
                       "standard_mix": _top(Counter(r["breakpoint_standard"] for r in sel))})
            gq.append({"cohort": key[0], "antimicrobial": key[1], "split_strategy": key[2], "outcome": outcome,
                       "N": n, **{f"{f}_median": _median(sel, f) for f in ("assembly_length", "contig_count", "n50",
                                                                           "ambiguous_fraction")},
                       "qc_flagged_fraction": frac(sum(r["genome_qc_status"] != "pass" for r in sel), n),
                       "partial_amrfinder_call_fraction": frac(sum(r["partial_amrfinder_calls"] == "yes"
                                                                   for r in sel), n)})
    write_tsv(o / "phenotype_quality_audit.tsv", ["cohort", "antimicrobial", "split_strategy", "outcome", "N",
                                                  *[f"{f}_fraction" for f in PHENO_FLAGS], "method_mix",
                                                  "standard_mix"], pq)
    write_tsv(o / "genome_quality_audit.tsv", list(gq[0].keys()), gq)

    # augmentation
    aug_rows, rescue = augmentation_tables(config, outdir, benchmark_dir, ctx)
    if aug_rows:
        write_tsv(o / "mutation_augmented_metrics.tsv", list(aug_rows[0].keys()), aug_rows)
    if rescue:
        write_tsv(o / "generalisation_gap_rescue.tsv", list(rescue[0].keys()), rescue)

    # drug-level explanation table
    drops = {(r["antimicrobial"], r["holdout_strategy"]): r for r in read_tsv(Path(benchmark_dir) /
                                                                               "generalisation_drop.tsv")
             if r["cohort"] == "common" and r["comparator"] == "random" and r["level"] == "pooled"}
    shift_by = {(r["antimicrobial"], r["split_strategy"]): r for r in shift if r["cohort"] == "common"}
    resc_by = {(r["antimicrobial"], r["holdout_strategy"]): r for r in rescue
               if r["variant"] == "mutation_augmented"}
    det_auc = _determinant_auroc(ctx)
    explain = []
    for (drug, strat), d in sorted(drops.items()):
        rows = [r for r in ctx if r["cohort"] == "common" and r["antimicrobial"] == drug
                and r["split_strategy"] == strat]
        rr = [r for r in rows if r["phenotype"] == "R"]
        pos = [r for r in rr if r["established_determinant_status"] == "positive"]
        mut = [r for r in rr if r["mechanism_group"] in ("mutation", "both")]
        s = shift_by.get((drug, strat), {})
        rs = resc_by.get((drug, strat), {})
        explain.append({
            "antimicrobial": drug, "holdout_strategy": strat, "drug_group": groups.get(drug, {}).get("group", NA),
            "test_N": len(rows), "test_R": len(rr), "original_AUROC": d["AUROC_holdout"],
            "delta_AUROC": d["delta_AUROC"], "delta_AUROC_ci_lower": d["delta_AUROC_ci_lower"],
            "delta_AUROC_ci_upper": d["delta_AUROC_ci_upper"],
            "lineage_shift": s.get("lineage_tvd"), "lineage_unseen_fraction": s.get("test_lineage_unseen_fraction"),
            "mechanism_shift": s.get("R_mechanism_tvd"), "mutation_shift": s.get("R_mutation_tvd"),
            "phenotype_prevalence_shift": s.get("phenotype_prevalence_shift"),
            "known_mechanism_fraction": frac(len(pos), len(rr)),
            "known_mechanism_seen_fraction": frac(sum(r["mechanism_novelty"] == "seen" for r in pos), len(pos)),
            "unseen_mechanism_fraction": frac(sum(r["mechanism_novelty"] == "unseen" for r in pos), len(pos)),
            "mutation_associated_fraction": frac(len(mut), len(rr)),
            "mutation_not_captured_fraction": frac(sum(r["mutation_capture"] != "captured" for r in mut), len(mut)),
            "representation_failure_fraction": frac(sum(_not_captured(r) for r in pos), len(pos)),
            "determinant_baseline_AUROC": det_auc.get(("common", drug, strat)),
            "augmented_AUROC": rs.get("augmented_AUROC_holdout"),
            "augmented_delta_vs_random": rs.get("augmented_delta_AUROC"), "gap_rescue": rs.get("gap_rescue"),
            "gap_rescue_ci_lower": rs.get("gap_rescue_ci_lower"), "gap_rescue_ci_upper": rs.get("gap_rescue_ci_upper"),
            "control_gap_rescue": (next((x for x in rescue if x["variant"] == "permuted_control"
                                         and x["antimicrobial"] == drug and x["holdout_strategy"] == strat), {})
                                   .get("gap_rescue"))})
    write_tsv(o / "explain_generalisation_drop.tsv", list(explain[0].keys()), explain)
    manuscript = [{k: r[k] for k in (
        "antimicrobial", "holdout_strategy", "drug_group", "original_AUROC", "delta_AUROC",
        "known_mechanism_seen_fraction", "unseen_mechanism_fraction", "mutation_associated_fraction",
        "mutation_not_captured_fraction", "determinant_baseline_AUROC", "augmented_AUROC",
        "augmented_delta_vs_random", "gap_rescue")} for r in explain]
    write_tsv(o / "manuscript_table.tsv", list(manuscript[0].keys()), manuscript)

    out = {"analysis_version": config["analysis_version"], "generated_at": now(),
           "prediction_rows": len(ctx), "errors": len(errors),
           "error_categories": dict(Counter(r["error_category"] for r in errors)),
           "drug_groups": {d: g["group"] for d, g in groups.items()}, "case_study_drugs": manifest["case_study_drugs"],
           "associations": _associations(explain)}
    write_json(o / "error_decomposition_summary.json", out)
    return out


PHENO_FLAGS = ("label_corrected", "multiple_observations", "conflict", "method_missing",
               "breakpoint_standard_missing", "mixed_breakpoint_standards", "provenance_missing")


def _not_captured(r: dict) -> bool:
    caps = [c for c in (r["gene_capture"], r["mutation_capture"]) if c != NA]
    return bool(caps) and "captured" not in caps


def _group(rows, keys, all_drugs=False) -> dict:
    out = defaultdict(list)
    for r in rows:
        k = tuple(r[k] for k in keys)
        if all_drugs:
            k = (k[0], "ALL", *k[1:])
        out[k].append(r)
    return out


def _count_by(rows, key, value) -> dict:
    out = defaultdict(Counter)
    for r in rows:
        out[r[key]][r[value]] += 1
    return out


def _top(counter: Counter, n: int = 3) -> str:
    total = sum(counter.values())
    return "; ".join(f"{k} {v / total:.2f}" for k, v in counter.most_common(n)) if total else NA


def _median(rows, field):
    vals = [float(r[field]) for r in rows if r[field] not in (NA, "")]
    return float(np.median(vals)) if vals else None


def _determinant_auroc(ctx) -> dict:
    """AUROC of the 0/1 determinant call on the endpoints with an interpretable status."""
    out = {}
    for key, rows in _group(ctx, ("cohort", "antimicrobial", "split_strategy")).items():
        rows = [r for r in rows if r["determinant_prediction"] != NA]
        y = np.array([r["phenotype"] == "R" for r in rows], dtype=bool)
        s = np.array([r["determinant_prediction"] == "R" for r in rows], dtype=float)
        out[key] = gb.auroc(y, s) if rows else None
    return out


def _spearman(x, y):
    from scipy.stats import spearmanr
    pairs = [(a, b) for a, b in zip(x, y) if a not in (None, NA) and b not in (None, NA)]
    if len(pairs) < 5 or len({a for a, _ in pairs}) < 2 or len({b for _, b in pairs}) < 2:
        return {"n": len(pairs), "rho": None}
    rho = spearmanr([float(a) for a, _ in pairs], [float(b) for _, b in pairs]).statistic
    return {"n": len(pairs), "rho": float(rho)}


def _associations(explain: list[dict]) -> dict:
    y = [r["delta_AUROC"] for r in explain]
    return {f: _spearman([r[f] for r in explain], y) for f in (
        "lineage_shift", "mechanism_shift", "mutation_shift", "phenotype_prevalence_shift",
        "unseen_mechanism_fraction", "mutation_associated_fraction", "mutation_not_captured_fraction",
        "representation_failure_fraction")}


def augmentation_tables(config: dict, outdir: Path, benchmark_dir: Path, ctx: list[dict]) -> tuple[list, list]:
    adir = paths(outdir)["augment"]
    plan_path = adir / "augmentation_plan.json"
    if not plan_path.exists():
        return [], []
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    orig = defaultdict(dict)
    units = {}
    for r in ctx:
        if r["cohort"] == config["augmentation"]["cohort"] and r["antimicrobial"] in plan["drugs"]:
            orig[(r["antimicrobial"], r["split_strategy"])][r["isolate_id"]] = r
            units[r["isolate_id"]] = r["split_unit"]
    variant_preds = {}
    for variant in plan["variants"]:
        vp = defaultdict(dict)
        for f in sorted((adir / variant / "tasks").glob("*.tsv")):
            status = json.loads(f.with_suffix(".json").read_text(encoding="utf-8"))
            if status["status"] != "success":
                continue
            for r in read_tsv(f):
                vp[(status["antimicrobial"], status["split_strategy"])][r["isolate_id"]] = r
        variant_preds[variant] = vp
    metrics_rows, rescue = [], []
    boot = config["bootstrap"]
    for variant, vp in variant_preds.items():
        for (drug, strat), rows in sorted(orig.items()):
            aug = vp.get((drug, strat), {})
            if set(aug) != set(rows):
                continue
            iso = sorted(rows)
            o_rows = [rows[i] for i in iso]
            a_rows = [{"phenotype": rows[i]["phenotype"], "predicted_probability": aug[i]["predicted_probability"],
                       "predicted_class": aug[i]["predicted_class"]} for i in iso]
            po, pa = perf(o_rows, config), perf(a_rows, config)
            metrics_rows.append({"variant": variant, "antimicrobial": drug, "split_strategy": strat,
                                 **{f"original_{k}": po[k] for k in PERF_FIELDS},
                                 **{f"augmented_{k}": pa[k] for k in ("AUROC", "AUPRC", "balanced_accuracy", "MCC",
                                                                      "sensitivity", "specificity")}})
            if strat not in HOLDOUTS:
                continue
            rand_o, rand_a = orig.get((drug, "random"), {}), vp.get((drug, "random"), {})
            if not set(iso) <= set(rand_o) or not set(iso) <= set(rand_a):
                continue
            y = [rows[i]["phenotype"] == "R" for i in iso]
            pr = lambda d: [float(d[i]["predicted_probability"]) for i in iso]  # noqa: E731
            yo = np.array(y)
            a = {"o_rand": pr(rand_o), "o_hold": pr(rows), "a_rand": pr(rand_a), "a_hold": pr(aug)}
            au = {k: gb.auroc(yo, np.array(v)) for k, v in a.items()}
            ci = gap_rescue_ci(y, [units[i] for i in iso], a["o_rand"], a["o_hold"], a["a_rand"], a["a_hold"],
                               boot["replicates"], boot["seed"], boot["level"])
            rescue.append({"variant": variant, "antimicrobial": drug, "holdout_strategy": strat, "N": len(iso),
                           "R": int(sum(y)), "original_AUROC_random": au["o_rand"],
                           "original_AUROC_holdout": au["o_hold"],
                           "original_delta_AUROC": au["o_rand"] - au["o_hold"],
                           "augmented_AUROC_random": au["a_rand"], "augmented_AUROC_holdout": au["a_hold"],
                           "augmented_delta_AUROC": au["a_rand"] - au["a_hold"],
                           "gap_rescue": (au["o_rand"] - au["o_hold"]) - (au["a_rand"] - au["a_hold"]),
                           "gap_rescue_ci_lower": ci["lower"], "gap_rescue_ci_upper": ci["upper"],
                           "holdout_AUROC_change": au["a_hold"] - au["o_hold"]})
    return metrics_rows, rescue


# ---------------------------------------------------------------- command line

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("catalogue", "freeze", "calls", "audit", "context", "augment", "report"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--benchmark", type=Path, help="PR #14 generalisation_v1 output directory")
    parser.add_argument("--population", type=Path, help="PR #6 population build directory")
    parser.add_argument("--baseline", type=Path, help="PR #7 baseline build directory")
    parser.add_argument("--calls", type=Path, help="AMRFinderPlus call directory (calls)")
    parser.add_argument("--cache", type=Path, help="genome cache (audit)")
    parser.add_argument("--ast", type=Path, help="PR #5 processed AST directory (context)")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--variants", default="mutation_augmented,permuted_control")
    parser.add_argument("--only", help="regular expression on task IDs (augment)")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    config["_path"] = str(args.config)
    if args.command == "catalogue":
        up = config["upstream"]
        mapping = read_tsv(resolve(up["mapping_rules"]))
        db = resolve(up["amrfinder_database"])
        study = set(config["catalogue"]["antimicrobials"])
        rows = build_catalogue(mapping, mutation_loci(db, up["amrfinder_taxgroup"]), gene_functions(db), config,
                               study)
        path = resolve(config["catalogue"]["path"])
        write_tsv(path, CATALOGUE_FIELDS, rows)
        print(json.dumps({"rows": len(rows), "sha256": pop.sha256(path),
                          "point_mutations": len({r["determinant"] for r in rows
                                                  if r["determinant_type"] == "point_mutation"}),
                          "genes": len({r["determinant"] for r in rows if r["determinant_type"] == "gene"})}))
        return 0
    if args.outdir is None:
        parser.error("--outdir is required")
    need = {"freeze": ["benchmark"], "calls": ["calls"], "audit": ["benchmark", "cache"],
            "context": ["benchmark", "population", "baseline", "ast"],
            "augment": ["benchmark", "population", "baseline"], "report": ["benchmark", "population", "baseline"]}
    missing = [a for a in need[args.command] if getattr(args, a) is None]
    if missing:
        parser.error(f"{args.command} needs --{' --'.join(missing)}")
    if args.command == "freeze":
        m = stage_freeze(config, args.config, args.benchmark, args.outdir)
        print(json.dumps({"drug_groups": {d: g["group"] for d, g in m["drug_groups"].items()},
                          "case_study_drugs": m["case_study_drugs"]}))
        return 0
    check_frozen(config, args.config, args.outdir)
    if args.command == "calls":
        print(json.dumps(stage_calls(config, args.outdir, args.calls)))
    elif args.command == "audit":
        print(json.dumps(stage_audit(config, args.outdir, args.benchmark, args.cache, args.workers)))
    elif args.command == "context":
        check_upstream(config, args.benchmark)
        print(json.dumps(stage_context(config, args.outdir, args.benchmark, args.population, args.baseline,
                                       args.ast)))
    elif args.command == "augment":
        check_upstream(config, args.benchmark)
        print(json.dumps(stage_augment(config, args.outdir, args.benchmark, args.population, args.baseline,
                                       args.variants.split(","), args.only)))
    else:
        check_upstream(config, args.benchmark)
        s = stage_report(config, args.outdir, args.benchmark, args.population, args.baseline)
        print(json.dumps({k: s[k] for k in ("prediction_rows", "errors")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
