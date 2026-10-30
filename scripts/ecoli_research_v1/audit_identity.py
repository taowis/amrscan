#!/usr/bin/env python3
"""Identity audit of an assembled census snapshot.

Usage: python3 scripts/ecoli_research_v1/audit_identity.py <census_dir> <candidate_list.tsv>
Writes <census_dir>/identity_audit.json; prints only nonzero issue groups.
"""
import csv, json, re, sys, collections
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "bin"))
import ast_dataset as ad
census, cand = Path(sys.argv[1]), Path(sys.argv[2])
snap = census / "snapshot"
man = json.loads((snap / "manifest.json").read_text())
ast = ad.read_table(snap / "ast.tsv", "ast")
meta = collections.defaultdict(list)
for s in man["sources"]:
    if s["role"] == "metadata":
        for r in ad.read_table(snap / s["filename"], "metadata"): meta[r["target_acc"]].append(r)
cl = collections.defaultdict(list)
with cand.open(encoding="utf-8-sig", newline="") as h:
    rd = csv.reader(h, delimiter="\t"); hdr = [x.lstrip("#") for x in next(rd)]
    for r in rd: d = dict(zip(hdr, r)); cl[d["target_acc"]].append(d)
out = collections.Counter(); issues = collections.defaultdict(list)
FIELDS = ("biosample_acc", "asm_acc", "scientific_name", "species_taxid", "taxid", "strain")
for pdt, rows in meta.items():
    out["metadata_rows_per_pdt=%d" % len(rows)] += 1
    for f in FIELDS:
        vals = {r.get(f, "") for r in rows}
        if len(vals) > 1: issues[f"metadata rows disagree on {f}"].append((pdt, sorted(vals)))
    r0 = rows[0]
    if not re.fullmatch(r"PDT\d+\.\d+", pdt): issues["unversioned PDT"].append(pdt)
    if not re.fullmatch(r"SAM[END]\d+", r0["biosample_acc"]): issues["invalid BioSample"].append((pdt, r0["biosample_acc"]))
    if not re.fullmatch(r"GCA_\d{9}\.\d+", r0["asm_acc"]): issues["not a versioned GCA"].append((pdt, r0["asm_acc"]))
    if r0["scientific_name"] != "Escherichia coli": issues["organism not E. coli"].append((pdt, r0["scientific_name"]))
    if r0.get("species_taxid") != "562": issues["species_taxid not 562"].append((pdt, r0.get("species_taxid")))
    old = {(c["biosample_acc"], c["asm_acc"], c["scientific_name"], c.get("strain", "")) for c in cl.get(pdt, [])}
    new = {(r["biosample_acc"], r["asm_acc"], r["scientific_name"], r.get("strain", "")) for r in rows}
    if old != new: issues["changed since 2026-10-11 candidate list"].append((pdt, sorted(old), sorted(new)))
bs_to_pdt = collections.defaultdict(set); asm_to_pdt = collections.defaultdict(set)
for pdt, rows in meta.items():
    bs_to_pdt[rows[0]["biosample_acc"]].add(pdt); asm_to_pdt[rows[0]["asm_acc"]].add(pdt)
issues["BioSample shared by >1 PDT"] = [(b, sorted(p)) for b, p in bs_to_pdt.items() if len(p) > 1]
issues["Assembly shared by >1 PDT"] = [(a, sorted(p)) for a, p in asm_to_pdt.items() if len(p) > 1]
ids = collections.Counter(r["id"] for r in ast)
issues["duplicate AST source id"] = [i for i, n in ids.items() if n > 1]
for r in ast:
    m = meta.get(r["target_acc"])
    if not m: issues["AST PDT without metadata"].append(r["id"]); continue
    if r["biosample_acc"] != m[0]["biosample_acc"]: issues["AST BioSample differs from metadata"].append((r["id"], r["biosample_acc"], m[0]["biosample_acc"]))
    if r["scientific_name"] != m[0]["scientific_name"]: issues["AST organism differs from metadata"].append((r["id"], r["scientific_name"]))
issues["metadata PDT without AST"] = sorted(set(meta) - {r["target_acc"] for r in ast})
res = {"ast_rows": len(ast), "ast_pdts": len({r["target_acc"] for r in ast}), "metadata_pdts": len(meta),
       "metadata_rows": sum(len(v) for v in meta.values()), "counts": dict(out),
       "issues": {k: {"n": len(v), "examples": v[:5]} for k, v in issues.items()}}
(census / "identity_audit.json").write_text(json.dumps(res, indent=1) + "\n")
print(json.dumps({k: v for k, v in res.items() if k != "issues"}), "\n", json.dumps({k: v["n"] for k, v in res["issues"].items()}, indent=1))
for k, v in res["issues"].items():
    if v["n"]: print(k, v["examples"][:3])
