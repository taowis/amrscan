#!/usr/bin/env python3
"""Post-build checks on the E. coli genotype-phenotype baseline.

Usage:
  python3 scripts/ecoli_research_v1/baseline_qc.py \
      --baseline datasets/ecoli_research_v1/genotype_baseline_20261018 \
      --processed datasets/ecoli_research_v1/ast_census_20261012/processed \
      --outdir benchmark/cohorts/ecoli_research_v1/genotype_baseline_20261018/qc \
      [--corrections benchmark/cohorts/ecoli_research_v1/phenotype_corrections_v1.tsv]

Reads a finished baseline build and the PR #5 AST observations. Nothing here
changes the baseline, the PR #5 cohort or the mapping.

With --corrections, the label check is also run on the corrected labels
(R and S swapped for the listed projects), and every MIC observation of a
corrected project is compared with the medians of its project, drug and
breakpoint standard: a corrected R at or below the median MIC of corrected S,
or a corrected S at or above the median MIC of corrected R, is listed as out
of order. Without breakpoints this is a coarse check.

1. Phenotype label check: for each submitting project and antimicrobial with
   at least 3 R and 3 S MIC observations, compare the median MIC of R and S
   observations. R with a lower median MIC than S means the labels run against
   the measurements. A project is flagged when this holds for at least half
   of its tested antimicrobials. Disk-only data are not checked.
2. Discordance concentration: for each study drug and concordance category,
   the largest lineage, country, project, source, period and (for DET+)
   determinant set, with its share.
3. Post hoc metrics: primary-variant metrics recomputed without the flagged
   projects. This is a descriptive sensitivity table chosen after seeing the
   data, not a replacement for the primary baseline.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
import genotype_baseline as gb  # noqa: E402

MIN_PER_CLASS = 3
CATEGORIES = ("R_DET_POS", "R_DET_NEG", "S_DET_POS", "S_DET_NEG")
CONCENTRATION_KEYS = ("lineage_id", "country", "project", "source_category", "period")


def label_check(observations, project_of):
    mics = defaultdict(lambda: {"R": [], "S": []})
    for o in observations:
        if o["phenotype"] not in ("R", "S") or o["isolate_id"] not in project_of:
            continue
        try:
            mic = float(o["mic"])
        except ValueError:
            continue
        mics[(project_of[o["isolate_id"]], o["antimicrobial"])][o["phenotype"]].append(mic)
    rows = []
    for (project, drug), v in sorted(mics.items()):
        if len(v["R"]) < MIN_PER_CLASS or len(v["S"]) < MIN_PER_CLASS:
            continue
        r, s = statistics.median(v["R"]), statistics.median(v["S"])
        rows.append({"project": project, "antimicrobial": drug, "mic_R": len(v["R"]), "mic_S": len(v["S"]),
                     "median_mic_R": r, "median_mic_S": s, "inverted": "true" if r < s else "false"})
    tested, inverted = Counter(), Counter()
    for row in rows:
        tested[row["project"]] += 1
        inverted[row["project"]] += row["inverted"] == "true"
    flagged = sorted(p for p in tested if inverted[p] * 2 >= tested[p] and inverted[p])
    return rows, tested, inverted, flagged


def corrected_observations(observations, project_of, corrections):
    out = []
    for o in observations:
        if project_of.get(o["isolate_id"]) in corrections and o["phenotype"] in gb.SWAP_RS:
            o = {**o, "phenotype": gb.SWAP_RS[o["phenotype"]]}
        out.append(o)
    return out


def out_of_order(observations, project_of, corrections):
    groups = defaultdict(list)
    for o in observations:
        if project_of.get(o["isolate_id"]) in corrections and o["phenotype"] in ("R", "S"):
            try:
                groups[(project_of[o["isolate_id"]], o["antimicrobial"], o["breakpoint_standard"])].append(
                    (o, float(o["mic"])))
            except ValueError:
                continue
    rows, totals = [], Counter()
    for (project, drug, standard), obs in sorted(groups.items()):
        r = [m for o, m in obs if o["phenotype"] == "R"]
        s = [m for o, m in obs if o["phenotype"] == "S"]
        totals["observations"] += len(obs)
        if not r or not s:
            continue
        med_r, med_s = statistics.median(r), statistics.median(s)
        for o, m in obs:
            if (o["phenotype"] == "R" and m <= med_s) or (o["phenotype"] == "S" and m >= med_r):
                rows.append({"project": project, "antimicrobial": drug, "breakpoint_standard": standard,
                             "isolate_id": o["isolate_id"], "record_id": o["record_id"],
                             "corrected_phenotype": o["phenotype"], "mic": o["mic_operator"] + o["mic"],
                             "median_mic_R": med_r, "median_mic_S": med_s})
    totals["out_of_order"] = len(rows)
    return rows, totals


def concentration(eligible, drugs):
    rows = []
    for drug in sorted(drugs):
        for cat in CATEGORIES:
            group = [e for e in eligible if e["antimicrobial"] == drug and e["concordance_category"] == cat]
            row = {"antimicrobial": drug, "study_tier": drugs[drug], "category": cat, "n": len(group)}
            keys = CONCENTRATION_KEYS + (("determinant_ids",) if cat.endswith("POS") else ())
            for key in CONCENTRATION_KEYS + ("determinant_ids",):
                if key in keys and group:
                    value, count = sorted(Counter(e[key] for e in group).items(), key=lambda x: (-x[1], x[0]))[0]
                    row[f"top_{key}"], row[f"top_{key}_fraction"] = value, f"{count / len(group):.4f}"
                else:
                    row[f"top_{key}"], row[f"top_{key}_fraction"] = gb.NA, gb.NA
            rows.append(row)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--processed", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--corrections", type=Path)
    args = ap.parse_args(argv)

    endpoints = gb.read_table(args.baseline / "genotype_phenotype_endpoints.tsv")
    observations = gb.read_table(args.processed / "ast_observations.tsv")
    project_of = {e["isolate_id"]: e["project"] for e in endpoints}
    eligible = [e for e in endpoints if e["endpoint_eligible"] == "true"]
    drugs = {e["antimicrobial"]: e["study_tier"] for e in eligible}

    checks, tested, inverted, flagged = label_check(observations, project_of)
    args.outdir.mkdir(parents=True, exist_ok=True)
    check_fields = "project antimicrobial mic_R mic_S median_mic_R median_mic_S inverted".split()
    gb.write_table(args.outdir / "phenotype_label_check.tsv", check_fields, checks)
    corrections = {}
    if args.corrections:
        corrections = {r["project"]: r["correction_id"] for r in gb.read_table(args.corrections, gb.CORRECTION_FIELDS)}
        fixed = corrected_observations(observations, project_of, corrections)
        c_checks, _, c_inverted, c_flagged = label_check(fixed, project_of)
        gb.write_table(args.outdir / "phenotype_label_check_corrected.tsv", check_fields, c_checks)
        order_rows, order_totals = out_of_order(fixed, project_of, corrections)
        gb.write_table(args.outdir / "corrected_mic_out_of_order.tsv",
                       "project antimicrobial breakpoint_standard isolate_id record_id corrected_phenotype mic "
                       "median_mic_R median_mic_S".split(), order_rows)

    conc = concentration(eligible, drugs)
    conc_fields = ["antimicrobial", "study_tier", "category", "n"]
    for key in CONCENTRATION_KEYS + ("determinant_ids",):
        conc_fields += [f"top_{key}", f"top_{key}_fraction"]
    gb.write_table(args.outdir / "discordance_concentration.tsv", conc_fields, conc)

    metric_rows = []
    for drug in sorted(drugs):
        for subset, keep in (("all_projects", lambda e: True),
                             ("without_flagged_projects", lambda e: e["project"] not in flagged)):
            pairs = [(e["phenotype"], e["determinant_status"]) for e in eligible
                     if e["antimicrobial"] == drug and keep(e)]
            m = gb.metrics(*gb.confusion(pairs))
            metric_rows.append({"antimicrobial": drug, "study_tier": drugs[drug], "subset": subset,
                                **{k: gb.fmt(m[k]) for k in gb.METRIC_FIELDS}})
    gb.write_table(args.outdir / "metrics_without_flagged_projects.tsv",
                   ["antimicrobial", "study_tier", "subset"] + gb.METRIC_FIELDS, metric_rows)

    flagged_isolates = sorted({e["isolate_id"] for e in endpoints if e["project"] in flagged})
    summary = {
        "baseline_manifest_sha256": hashlib.sha256((args.baseline / "baseline_manifest.json").read_bytes()).hexdigest(),
        "ast_observations_sha256": hashlib.sha256((args.processed / "ast_observations.tsv").read_bytes()).hexdigest(),
        "label_check_rule": f"median MIC of R below median MIC of S, at least {MIN_PER_CLASS} R and "
                            f"{MIN_PER_CLASS} S MIC observations; project flagged when at least half of its "
                            "tested antimicrobials are inverted",
        "projects_tested": len(tested),
        "projects_flagged": {p: {"antimicrobials_tested": tested[p], "inverted": inverted[p]} for p in flagged},
        "flagged_isolates": len(flagged_isolates),
        "flagged_metric_eligible_endpoints": sum(e["project"] in flagged for e in eligible),
        "baseline_phenotype_corrections": sorted({e["phenotype_correction"] for e in endpoints
                                                  if e.get("phenotype_correction", gb.NA) != gb.NA}),
        "status": "post hoc descriptive check; the PR #5 endpoints are unchanged; the label check on the "
                  "PR #5 labels is the evidence for decision 0004",
    }
    if args.corrections:
        summary["corrections_sha256"] = hashlib.sha256(args.corrections.read_bytes()).hexdigest()
        summary["corrections"] = corrections
        summary["projects_flagged_after_corrections"] = c_flagged
        summary["inverted_after_corrections"] = sum(c_inverted.values())
        summary["corrected_mic_observations"] = order_totals["observations"]
        summary["corrected_mic_out_of_order"] = order_totals["out_of_order"]
    (args.outdir / "qc_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
