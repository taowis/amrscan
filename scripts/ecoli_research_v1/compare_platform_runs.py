#!/usr/bin/env python3
"""Compare two runs of the same generalisation_v1 tasks task by task (noise baseline).

Used for the representation_sensitivity_v1 platform check: the PR #14 R0 fits
from the workstation against the same tasks refitted on Gadi (R0_gadi).

    .venv/bin/python scripts/ecoli_research_v1/compare_platform_runs.py \
        datasets/ecoli_research_v1/generalisation_v1 \
        datasets/ecoli_research_v1/representation_sensitivity_v1/R0_gadi \
        --out benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1/r0_platform_check.tsv
"""
import argparse
import csv
import json
from pathlib import Path

from sklearn.metrics import roc_auc_score


def read_predictions(path: Path) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["isolate_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def auroc(rows: dict) -> float:
    y = [r["phenotype"] == "R" for r in rows.values()]
    if len(set(y)) < 2:
        return float("nan")
    return roc_auc_score(y, [float(r["predicted_probability"]) for r in rows.values()])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", type=Path, help="run directory with tasks/ (e.g. generalisation_v1)")
    ap.add_argument("other", type=Path, help="second run directory with tasks/")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    rows = []
    for meta in sorted((a.other / "tasks").glob("*.json")):
        ref_meta = a.reference / "tasks" / meta.name
        m1, m2 = json.loads(ref_meta.read_text()), json.loads(meta.read_text())
        p1 = read_predictions(ref_meta.with_suffix(".tsv"))
        p2 = read_predictions(meta.with_suffix(".tsv"))
        same_isolates = p1.keys() == p2.keys()
        diffs = [abs(float(p1[k]["predicted_probability"]) - float(p2[k]["predicted_probability"]))
                 for k in p1 if k in p2]
        rows.append({
            "task_id": m2["task_id"],
            "same_split": m1.get("split_sha256") == m2.get("split_sha256"),
            "same_isolates": same_isolates,
            "same_C": m1.get("C") == m2.get("C"),
            "same_features_used": m1.get("features_used") == m2.get("features_used"),
            "n_test": len(p2),
            "class_disagreements": sum(p1[k]["predicted_class"] != p2[k]["predicted_class"] for k in p1 if k in p2),
            "max_abs_probability_diff": f"{max(diffs):.3g}" if diffs else "NA",
            "auroc_reference": f"{auroc(p1):.6f}",
            "auroc_other": f"{auroc(p2):.6f}",
            "auroc_diff": f"{auroc(p2) - auroc(p1):.3g}",
        })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    flags = ("same_split", "same_isolates", "same_C", "same_features_used")
    summary = {"tasks": len(rows), **{f: sum(r[f] for r in rows) for f in flags},
               "tasks_with_class_disagreement": sum(r["class_disagreements"] > 0 for r in rows),
               "class_disagreements": sum(r["class_disagreements"] for r in rows),
               "max_abs_probability_diff": max(float(r["max_abs_probability_diff"]) for r in rows),
               "max_abs_auroc_diff": max(abs(float(r["auroc_diff"])) for r in rows
                                         if r["auroc_diff"] != "nan")}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
