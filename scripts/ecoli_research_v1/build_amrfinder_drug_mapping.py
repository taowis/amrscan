#!/usr/bin/env python3
"""Generate the E. coli AMRFinderPlus determinant-to-antimicrobial mapping.

Usage:
  python3 scripts/ecoli_research_v1/build_amrfinder_drug_mapping.py \
      --amrfinder-db datasets/amrfinderplus-db/2026-08-07.1 \
      --resfinder-db datasets/mapping_sources/resfinder_db-4.7.2 \
      --pointfinder-db datasets/mapping_sources/pointfinder_db-4.7.2 \
      --evaluation-plan datasets/ecoli_research_v1/population_20261015/evaluation_plan.tsv \
      --outdir benchmark/mappings

Inputs are external reference tables only: the AMRFinderPlus symbol catalogue
(which symbols the caller can report), the ResFinder gene phenotype table and
the PointFinder E. coli mutation table (which antimicrobials each determinant
is documented to affect), plus the study antimicrobials fixed in decision 0001.
No caller output and no AST phenotype is read.

The drug-term translation, inference rules and overrides below are part of
mapping version 1 and were fixed before any genotype was compared with AST.
"""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

MAPPING_VERSION = "ecoli_amrfinder_drug_mapping_v1"
# Revision 2 (decision 0003) differs from v1 in one rule only: an allele with no
# ResFinder entry of its own whose AMRFinderPlus parent node is itself a
# ResFinder gene (rmtE1 under rmtE) takes that gene's phenotypes. It was
# written after the v1 baseline had been run and is exploratory.
REVISION = 1
ORGANISM = "Escherichia coli"

# Source phenotype term -> PR #5 antimicrobial name. Exact, case-insensitive
# string equality on the stripped term; nothing else is matched.
DIRECT_TERMS = {
    "ampicillin": "ampicillin",
    "amoxicillin+clavulanic acid": "amoxicillin-clavulanic acid",
    "cefoxitin": "cefoxitin",
    "ceftriaxone": "ceftriaxone",
    "ceftazidime": "ceftazidime",
    "cefepime": "cefepime",
    "ciprofloxacin": "ciprofloxacin",
    "nalidixic acid": "nalidixic acid",
    "gentamicin": "gentamicin",
    "gentamicin c": "gentamicin",
    "tetracycline": "tetracycline",
    "chloramphenicol": "chloramphenicol",
    "streptomycin": "streptomycin",
    "meropenem": "meropenem",
    "ertapenem": "ertapenem",
    "imipenem": "imipenem",
    "amikacin": "amikacin",
    "tobramycin": "tobramycin",
    "piperacillin+tazobactam": "piperacillin-tazobactam",
    "aztreonam": "aztreonam",
    "kanamycin": "kanamycin",
    "kanamycin a": "kanamycin",
    "ceftiofur": "ceftiofur",
    # components of the trimethoprim-sulfamethoxazole combination rule
    "trimethoprim": "trimethoprim-sulfamethoxazole#trimethoprim",
    "sulfamethoxazole": "trimethoprim-sulfamethoxazole#sulfonamide",
}

# Study antimicrobials without a source term of their own. Each is inferred
# from named source terms and flagged for review.
INFERRED_TERMS = {
    "levofloxacin": (("ciprofloxacin", "fluoroquinolone"),
                     "same target: QRDR mutations and PMQR determinants act on fluoroquinolones as a group "
                     "(Hooper & Jacoby 2015, PMID 26190223); source lists ciprofloxacin or fluoroquinolone"),
    "sulfisoxazole": (("sulfamethoxazole",),
                      "same target: sul genes and folP changes give an insensitive dihydropteroate synthase, "
                      "the target of all sulfonamides; source lists sulfamethoxazole"),
    "ceftriaxone": (("cefotaxime",),
                    "cefotaxime and ceftriaxone are oxyimino-cephalosporins with the same hydrolysis profile; "
                    "source lists cefotaxime"),
    "ceftiofur": (("cefotaxime", "ceftriaxone"),
                  "ceftiofur is an oxyimino-cephalosporin; source lists cefotaxime or ceftriaxone"),
    "cefazolin": (("cephalothin", "cefotaxime", "ceftriaxone", "ceftazidime"),
                  "first-generation cephalosporin: source lists cephalothin, or an oxyimino-cephalosporin whose "
                  "hydrolysis implies hydrolysis of narrow-spectrum cephalosporins (Bush & Jacoby 2010, "
                  "PMID 19995920)"),
}

# "Unknown <class>" source phenotypes give an uncertain row for the study
# antimicrobials of that class (class lists from ResFinder antibiotic_classes.txt
# restricted to study drugs).
UNKNOWN_CLASS = {
    "unknown beta-lactam": ("ampicillin", "amoxicillin-clavulanic acid", "cefoxitin", "cefazolin", "ceftriaxone",
                            "ceftazidime", "cefepime", "meropenem", "ertapenem", "imipenem",
                            "piperacillin-tazobactam", "aztreonam", "ceftiofur"),
    "unknown aminoglycoside": ("gentamicin", "streptomycin", "amikacin", "tobramycin", "kanamycin"),
    "unknown quinolone": ("ciprofloxacin", "levofloxacin", "nalidixic acid"),
}

LOW_LEVEL_FQ_NOTE = ("source phenotypes describe non-wild-type (ECOFF) resistance; one QRDR mutation or one PMQR "
                     "determinant often leaves the MIC below the clinical breakpoint (Hooper & Jacoby 2015, "
                     "PMID 26190223; Strahilevitz et al. 2009, PMID 19822894)")
PMQR_FAMILIES = ("qnr", "aac(6')-Ib-cr", "qepA", "oqxA", "oqxB")

BETA_LACTAM_STUDY = UNKNOWN_CLASS["unknown beta-lactam"]

REFERENCES = {
    "resfinder": "Bortolaia V et al. 2020. ResFinder 4.0 for predictions of phenotypes from genotypes. "
                 "J Antimicrob Chemother 75:3491-3500. doi:10.1093/jac/dkaa345. PMID 32780112",
    "pointfinder": "Zankari E et al. 2017. PointFinder. J Antimicrob Chemother 72:2764-2768. "
                   "doi:10.1093/jac/dkx217. PMID 29091202",
    "amrfinderplus": "Feldgarden M et al. 2021. AMRFinderPlus and the Reference Gene Catalog. Sci Rep 11:12728. "
                     "doi:10.1038/s41598-021-91456-0. PMID 34135355",
    "ampc": "Jacoby GA. 2009. AmpC beta-lactamases. Clin Microbiol Rev 22:161-182. doi:10.1128/CMR.00036-08. "
            "PMID 19136439",
    "dfra34": "Tagg KA et al. 2019. Novel trimethoprim resistance gene dfrA34 identified in Salmonella Heidelberg "
              "in the USA. J Antimicrob Chemother 74:38-41. doi:10.1093/jac/dky373. PMID 30202900",
    "quinolone": "Hooper DC, Jacoby GA. 2015. Mechanisms of drug resistance: quinolone resistance. "
                 "Ann N Y Acad Sci 1354:12-31. doi:10.1111/nyas.12830. PMID 26190223",
    "pmqr": "Strahilevitz J et al. 2009. Plasmid-mediated quinolone resistance: a multifaceted threat. "
            "Clin Microbiol Rev 22:664-689. doi:10.1128/CMR.00016-09. PMID 19822894",
    "bush": "Bush K, Jacoby GA. 2010. Updated functional classification of beta-lactamases. Antimicrob Agents "
            "Chemother 54:969-976. doi:10.1128/AAC.01009-09. PMID 19995920",
    "stoesser": "Stoesser N et al. 2013. Predicting antimicrobial susceptibilities for Escherichia coli and "
                "Klebsiella pneumoniae isolates using whole genomic sequence data. J Antimicrob Chemother "
                "68:2234-2244. doi:10.1093/jac/dkt180. PMID 23722448",
}

FIELDS = ("mapping_version organism rule_id determinant determinant_type gene_family allele_or_mutation "
          "antimicrobial interpretation requires evidence_source evidence_terms evidence_reference "
          "evidence_level review_status notes").split()

AA3 = re.compile(r"^[A-Z]{3}-(-?\d+)$")


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git_commit(path: Path) -> str:
    return subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], check=True, text=True,
                          capture_output=True).stdout.strip()


def read_rows(path: Path) -> list[list[str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t")]


# ---------------------------------------------------------------- AMRFinderPlus catalogue

def amrfinder_catalogue(db: Path) -> dict:
    """Reportable AMR symbols: alleles and gene/family nodes, plus E. coli point mutations."""
    parent, node_type, reportable = {}, {}, set()
    for row in read_rows(db / "fam.tsv")[1:]:
        if row and not row[0].startswith("#"):
            parent[row[0]] = row[1]
            node_type[row[0]] = row[13]
            if row[12] != "0":
                reportable.add(row[0])
    alleles = {}  # allele symbol -> node
    for line in (db / "AMRProt.fa").read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            fields = line[1:].split("|")
            symbol, node = fields[3], fields[4]
            if node_type.get(node) == "AMR" or node_type.get(symbol) == "AMR":
                alleles.setdefault(symbol, node)
    genes = {n for n, t in node_type.items() if t == "AMR" and n in reportable}
    mutations = {}
    for row in read_rows(db / "AMRProt-mutation.tsv")[1:]:
        if row[0] == "Escherichia":
            mutations[row[4]] = row[0]
    for row in read_rows(db / "AMR_DNA-Escherichia.tsv")[1:]:
        mutations[row[3]] = "Escherichia"
    return {"parent": parent, "alleles": alleles, "genes": genes, "mutations": mutations,
            "version": (db / "version.txt").read_text().strip()}


def ancestors(node: str, parent: dict) -> list[str]:
    chain, seen = [], set()
    while node and node not in seen and node in parent:
        chain.append(node)
        seen.add(node)
        node = parent[node]
    return chain


# ---------------------------------------------------------------- source phenotypes

def resfinder_genes(db: Path) -> dict:
    """Gene name -> list of (phenotype terms, PMIDs, entry id, notes)."""
    out = defaultdict(list)
    for row in read_rows(db / "phenotypes.txt")[1:]:
        if not row or row[0].startswith("#"):
            continue
        entry = row[0]
        name = re.sub(r"_\d+_[^_]+$", "", entry)
        terms = tuple(sorted({t.strip().lower() for t in row[2].split(",") if t.strip()}))
        out[name].append((terms, row[3].strip() if len(row) > 3 else "", entry,
                          row[5].strip() if len(row) > 5 else ""))
    return out


def pointfinder_mutations(db: Path, catalogue: dict) -> tuple[dict, list]:
    """AMRFinderPlus mutation symbol -> {'alone': terms, 'with': terms, 'requires': [...], 'pmids': ...}.

    Returns the mapping and the PointFinder entries with no AMRFinderPlus symbol."""
    gene_name = {"gyrA": "gyrA", "gyrB": "gyrB", "parC": "parC", "parE": "parE", "folP": "folP",
                 "pmrA": "basR", "pmrB": "pmrB", "rpoB": "rpoB", "23S": "23S", "16S-rrsB": "16S",
                 "16S-rrsC": "16S", "16S-rrsH": "16S", "ampC-promoter-size-53bp": "ampC"}
    # Promoter insertions use different coordinates in the two tables; the
    # PointFinder note states the -13/-14/-15/-16 GT insertions are one event.
    insertions = {"-14:GT": "ampC_T-14TGT"}
    out, unmatched = {}, []

    def symbols(gene, ref, pos, alts):
        return [f"{gene_name[gene]}_{ref}{pos}{alt}" for alt in alts]

    def required_symbols(text):
        found = []
        for part in [p.strip() for p in text.split(",") if p.strip()]:
            m = re.match(r"^(\w+?)_\d+_[A-Z0-9.]+_([A-Z])(-?\d+)([A-Z.]+)$", part)
            if not m:
                raise ValueError(f"unparsed PointFinder required mutation: {part}")
            gene, ref, pos, alts = m.groups()
            found += symbols(gene, ref, pos, alts.split("."))
        return sorted(found)

    for row in read_rows(db / "escherichia_coli" / "phenotypes.txt"):
        if not row or row[0].startswith("#") or len(row) < 9:
            continue
        gene = re.sub(r"_\d+_[^_]+$", "", row[0])
        if gene not in gene_name:
            raise ValueError(f"unknown PointFinder gene {gene}")
        terms = tuple(sorted({t.strip().lower() for t in row[8].split(",") if t.strip()}))
        required = row[12].strip() if len(row) > 12 else ""
        pmids = row[9].strip()
        if row[1] == "AA" or (row[1] == "NUC" and AA3.match(row[2])):
            pos = row[3]
            ref = row[5]
            syms = symbols(gene, ref, pos, [a.strip() for a in row[6].split(",") if a.strip()])
        else:  # promoter insertion rows
            key = f"{row[3]}:{row[6]}"
            syms = [insertions[key]] if key in insertions else []
            if not syms:
                unmatched.append({"pointfinder_entry": row[0], "mutation": f"{row[3]} ins {row[6]}",
                                  "reason": "no matching AMRFinderPlus symbol"})
                continue
        for sym in syms:
            if sym not in catalogue["mutations"]:
                unmatched.append({"pointfinder_entry": row[0], "mutation": sym,
                                  "reason": "not reported by AMRFinderPlus for Escherichia"})
                continue
            rec = out.setdefault(sym, {"alone": set(), "with": set(), "requires": set(), "pmids": set(),
                                       "entry": row[0]})
            rec["pmids"].add(pmids)
            if required:
                rec["with"].update(terms)
                rec["requires"].update(required_symbols(required))
            else:
                rec["alone"].update(terms)
    return out, unmatched


# ---------------------------------------------------------------- rule assembly

def translate(terms: set) -> dict:
    """Study antimicrobial -> (evidence level, matched terms)."""
    out = {}
    for term in sorted(terms):
        if term in DIRECT_TERMS:
            out.setdefault(DIRECT_TERMS[term], ("source_exact_term", set()))[1].add(term)
    for drug, (sources, _) in INFERRED_TERMS.items():
        hit = {t for t in terms if t in sources}
        if hit and drug not in out:
            out[drug] = ("inferred_same_target", hit)
    return out


def interpretations(entries: list[set], study: set) -> dict:
    """Per drug: predicts if every source entry supports it, uncertain if only some do.

    Entries that name only an unknown class say nothing about specific drugs, so
    they are left out when other entries are specific."""
    specific = [t for t in entries if not all(x.startswith("unknown ") for x in t)]
    if specific:
        entries = specific
    per_entry = []
    for terms in entries:
        drugs = {d: v for d, v in translate(set(terms)).items() if d.split("#")[0] in study}
        unknown = set()
        for term in terms:
            unknown.update(d for d in UNKNOWN_CLASS.get(term, ()) if d in study)
        per_entry.append((drugs, unknown))
    result = {}
    all_drugs = set().union(*[set(d) | u for d, u in per_entry]) if per_entry else set()
    for drug in sorted(all_drugs):
        supporting = [d[drug] for d, _ in per_entry if drug in d]
        if len(supporting) == len(per_entry):
            levels = {lvl for lvl, _ in supporting}
            level = "inferred_same_target" if "inferred_same_target" in levels else "source_exact_term"
            result[drug] = ("predicts_resistance", level, sorted(set().union(*[t for _, t in supporting])))
        else:
            terms = sorted(set().union(*[t for _, t in supporting])) if supporting else ["unknown class"]
            result[drug] = ("uncertain", "source_partial_support", terms)
    return result


def component(drug: str, interp: str) -> tuple[str, str | None]:
    """Split 'drug#component' keys of the combination rule."""
    if "#" in drug:
        base, part = drug.split("#")
        return base, part
    return drug, None


def build(args) -> dict:
    db = Path(args.amrfinder_db)
    cat = amrfinder_catalogue(db)
    rf = resfinder_genes(Path(args.resfinder_db))
    pf, pf_unmatched = pointfinder_mutations(Path(args.pointfinder_db), cat)
    plan = {r["antimicrobial"]: r["tier"] for r in csv.DictReader(
        Path(args.evaluation_plan).open(encoding="utf-8"), delimiter="\t")}
    study = {d for d, t in plan.items() if t in {"primary", "secondary"}}

    # descendants: node -> alleles below it
    below = defaultdict(set)
    for allele, node in cat["alleles"].items():
        for anc in ancestors(node, cat["parent"]):
            below[anc].add(allele)
        below[allele].add(allele)
    symbols = sorted(set(cat["alleles"]) | cat["genes"])
    rows, coverage = [], {"exact": 0, "family": 0, "none": 0}

    def add(sym, dtype, family, drug, interp, requires, source, terms, ref, level, review, notes):
        rows.append({"mapping_version": MAPPING_VERSION, "organism": ORGANISM, "rule_id": f"{sym}|{drug}",
                     "determinant": sym, "determinant_type": dtype, "gene_family": family,
                     "allele_or_mutation": sym if dtype == "point_mutation" or cat["alleles"].get(sym, sym) != sym
                     else "NA",
                     "antimicrobial": drug, "interpretation": interp, "requires": requires or "NA",
                     "evidence_source": source, "evidence_terms": ";".join(terms) or "NA",
                     "evidence_reference": ref or "NA", "evidence_level": level, "review_status": review,
                     "notes": notes or "NA"})

    for sym in symbols:
        family = cat["alleles"].get(sym, cat["parent"].get(sym, "NA")) or "NA"
        if sym == "blaEC" or family == "blaEC" or "blaEC" in ancestors(family, cat["parent"]):
            for drug in sorted(d for d in BETA_LACTAM_STUDY if d in study):
                add(sym, "gene", family, drug, "does_not_predict", None, "override", [], REFERENCES["ampc"],
                    "expert_rule", "reviewed_source",
                    "intrinsic chromosomal AmpC of E. coli; resistance needs promoter/attenuator changes, mapped "
                    "separately as ampC point mutations")
            continue
        # ResFinder splits some alleles into nucleotide variants with a letter
        # suffix (blaTEM-1A to -1D encode the TEM-1 protein).
        variants = sorted(n for n in rf if re.fullmatch(re.escape(sym) + r"[A-Z]", n)) if sym not in rf else []
        if sym in rf or variants:
            names = [sym] if sym in rf else variants
            entries = [e for n in names for e in rf[n]]
            level_note = "ResFinder entries for " + ", ".join(names)
            coverage["exact"] += 1
        else:
            # An allele missing from ResFinder borrows from its AMRFinderPlus node;
            # a gene or family node uses the alleles below it.
            node = family if cat["alleles"].get(sym, sym) != sym else sym
            if REVISION >= 2 and node != sym and node in rf:
                names = [node]
                level_note = f"ResFinder entries for parent gene {node}"
                coverage["parent_gene"] = coverage.get("parent_gene", 0) + 1
            else:
                names = sorted(a for a in below.get(node, set()) if a in rf)
                level_note = f"family consensus over {len(names)} ResFinder genes under {node}"
                coverage["family" if names else "none"] += 1
            entries = [e for a in names for e in rf[a]]
        if not entries:
            continue
        pmids = sorted({e[1] for e in entries if e[1]})
        ref = "ResFinder 4.7.2 phenotypes.txt; PMID " + ", ".join(pmids[:6]) + (" ..." if len(pmids) > 6 else "")
        notes = []
        for drug, (interp, level, terms) in interpretations([set(e[0]) for e in entries], study).items():
            if sym not in rf and not variants and interp == "predicts_resistance":
                level = "family_consensus" if level == "source_exact_term" else level
            review = "reviewed_source"
            note = level_note
            if level == "inferred_same_target":
                review = "HUMAN_REVIEW_REQUIRED"
                note += "; " + INFERRED_TERMS[drug][1]
            if drug in {"ciprofloxacin", "levofloxacin"} and any(sym.startswith(p) for p in PMQR_FAMILIES):
                review = "HUMAN_REVIEW_REQUIRED"
                note += "; " + LOW_LEVEL_FQ_NOTE
            if notes:
                note += "; source note: " + " | ".join(notes)[:300]
            drug, requires = component(drug, interp)
            if requires:
                interp = "component" if interp == "predicts_resistance" else interp
            add(sym, "gene", family, drug, interp, requires, "ResFinder", terms, ref, level, review, note)

    for sym in sorted(pf):
        rec = pf[sym]
        ref = "PointFinder 4.7.2 escherichia_coli/phenotypes.txt; PMID " + ", ".join(sorted(rec["pmids"]))
        alone = translate(rec["alone"])
        with_req = translate(rec["with"])
        for drug in sorted(d for d in set(alone) | set(with_req) if d.split("#")[0] in study):
            if drug in alone:
                level, terms = alone[drug]
                requires = None
            else:
                level, terms = with_req[drug]
                requires = "|".join(sorted(rec["requires"]))
            review, note = "reviewed_source", "PointFinder E. coli mutation"
            if level == "inferred_same_target":
                review, note = "HUMAN_REVIEW_REQUIRED", note + "; " + INFERRED_TERMS[drug][1]
            if drug in {"ciprofloxacin", "levofloxacin"} and requires is None:
                review, note = "HUMAN_REVIEW_REQUIRED", note + "; " + LOW_LEVEL_FQ_NOTE
            if requires:
                note += "; predicts only together with a required mutation listed in requires"
            interp = "predicts_resistance"
            drug, part = component(drug, interp)
            if part:
                interp, requires = "component", part
            add(sym, "point_mutation", sym.split("_")[0], drug, interp, requires, "PointFinder",
                sorted(terms), ref, "source_exact_term" if level == "source_exact_term" else level, review, note)
    # The -15 G insertion has no exact PointFinder counterpart; see the -13 to -16 note there.
    if "ampC_G-15GG" in cat["mutations"]:
        for drug in ("ampicillin", "amoxicillin-clavulanic acid", "cefoxitin", "ceftazidime"):
            add("ampC_G-15GG", "point_mutation", "ampC", drug, "uncertain", None, "PointFinder",
                ["ampicillin" if drug == "ampicillin" else drug], "PointFinder 4.7.2; PMID 21653764",
                "expert_rule", "HUMAN_REVIEW_REQUIRED",
                "promoter insertion near the PointFinder -13 to -16 insertions but not identical to them")

    rows.sort(key=lambda r: (r["determinant"], r["antimicrobial"]))
    ids = [r["rule_id"] for r in rows]
    if len(ids) != len(set(ids)):
        import collections; raise ValueError(f"duplicate rule ids: {[k for k, n in collections.Counter(ids).items() if n > 1][:20]}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / f"{MAPPING_VERSION}.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    study_list = sorted(study)
    meta = {
        "mapping_version": MAPPING_VERSION,
        "organism": ORGANISM,
        "status": "prespecified" if REVISION == 1 else "exploratory",
        "fixed_on": args.fixed_on,
        "rules_file": tsv.name,
        "rules_sha256": sha256(tsv),
        "rule_count": len(rows),
        "study_antimicrobials": {d: plan[d] for d in study_list},
        "not_mapped_antimicrobials": {d: plan[d] for d in sorted(set(plan) - study)},
        "decision_record": ("docs/decisions/0002-ecoli-determinant-mapping-v1.md" if REVISION == 1
                            else "docs/decisions/0003-ecoli-determinant-mapping-v2.md"),
        "evidence_policy": {
            "element_type": "AMR",
            "ignored_element_types": ["STRESS", "VIRULENCE"],
            "matching": "exact AMRFinderPlus element symbol; no class or text matching",
            "partial_methods_as_uncertain": ["INTERNAL_STOP", "PARTIALP", "PARTIALX"],
            "status_order": "positive if any predicts_resistance rule applies with its requires met; otherwise "
                            "uncertain if any uncertain rule applies or a requires is unmet; otherwise negative",
        },
        "combination_rules": {
            "trimethoprim-sulfamethoxazole": {
                "components": {"trimethoprim": ["trimethoprim"], "sulfonamide": ["sulfamethoxazole"]},
                "positive": "a determinant for each component",
                "partial": "uncertain",
                "reference": REFERENCES["dfra34"],
                "review_status": "HUMAN_REVIEW_REQUIRED",
            }
        },
        "rule_variants": {
            "primary": "uncertain endpoints excluded from metrics",
            "uncertain_as_positive": "uncertain endpoints counted as determinant-positive",
            "uncertain_as_negative": "uncertain endpoints counted as determinant-negative",
            "fq_two_mechanisms": "ciprofloxacin and levofloxacin positive only with at least two distinct "
                                 "predicting determinants; exactly one is uncertain",
        },
        "sources": {
            "amrfinderplus_database": {"version": cat["version"], "path": str(db),
                                       "files_sha256": {n: sha256(db / n) for n in (
                                           "fam.tsv", "AMRProt.fa", "AMRProt-mutation.tsv",
                                           "AMR_DNA-Escherichia.tsv")}},
            "resfinder_db": {"tag": "resfinder-4.7.2", "commit": git_commit(Path(args.resfinder_db)),
                             "url": "https://bitbucket.org/genomicepidemiology/resfinder_db",
                             "phenotypes_sha256": sha256(Path(args.resfinder_db) / "phenotypes.txt")},
            "pointfinder_db": {"tag": "resfinder-4.7.2", "commit": git_commit(Path(args.pointfinder_db)),
                               "url": "https://bitbucket.org/genomicepidemiology/pointfinder_db",
                               "phenotypes_sha256": sha256(Path(args.pointfinder_db) / "escherichia_coli"
                                                           / "phenotypes.txt")},
            "evaluation_plan_sha256": sha256(Path(args.evaluation_plan)),
        },
        "source_term_translation": {"direct": DIRECT_TERMS,
                                    "inferred": {d: list(v[0]) for d, v in INFERRED_TERMS.items()},
                                    "unknown_class": {k: list(v) for k, v in UNKNOWN_CLASS.items()}},
        "catalogue_coverage": coverage,
        "pointfinder_entries_without_amrfinder_symbol": pf_unmatched,
        "references": REFERENCES,
    }
    (outdir / f"{MAPPING_VERSION}.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n",
                                                    encoding="utf-8")
    return meta


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("amrfinder-db", "resfinder-db", "pointfinder-db", "evaluation-plan", "outdir"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--fixed-on", required=True, help="date the mapping was frozen (YYYY-MM-DD)")
    parser.add_argument("--revision", type=int, choices=(1, 2), default=1,
                        help="1: prespecified mapping v1; 2: exploratory v2 (decision 0003)")
    args = parser.parse_args()
    global MAPPING_VERSION, REVISION
    REVISION = args.revision
    MAPPING_VERSION = f"ecoli_amrfinder_drug_mapping_v{REVISION}"
    meta = build(args)
    print(json.dumps({k: meta[k] for k in ("rule_count", "rules_sha256", "catalogue_coverage")}, indent=2))


if __name__ == "__main__":
    main()
