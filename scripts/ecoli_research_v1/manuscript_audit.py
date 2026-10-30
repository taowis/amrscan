#!/usr/bin/env python3
"""Descriptive checks for the Microbial Genomics manuscript, from existing outputs only.

No model is fitted and no frozen table is changed. The script reads the committed result tables of
generalisation_v1, error_decomposition_v1 and representation_sensitivity_v1 and, for the tables that need
prediction-level data, the uncommitted predictions, splits and error context under datasets/. Everything here
was added after those results had been read, while revising the manuscript, and is reported as such.

    .venv/bin/python scripts/ecoli_research_v1/manuscript_audit.py \
        --out benchmark/cohorts/ecoli_research_v1/manuscript_audit_v1
"""
import argparse
import csv
import gzip
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
COMMITTED = ROOT / "benchmark/cohorts/ecoli_research_v1"
DATASETS = ROOT / "datasets/ecoli_research_v1"
HOLDOUTS = ("lineage", "country", "project", "geography", "temporal")
STRATEGIES = ("random", "random_cluster") + HOLDOUTS
MIN_CARRIERS = 5
FLAGGED_PROJECT = "PRJNA1322038"


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def frac(a, b):
    return a / b if b else None


def auroc(y, p):
    return roc_auc_score(y, p) if len(set(y)) == 2 else None


def split_sizes(out: Path) -> None:
    rows = [r for r in read_tsv(COMMITTED / "generalisation_v1/split_summary.tsv") if r["cohort"] == "common"]
    by = defaultdict(list)
    for r in rows:
        by[r["strategy"]].append(r)
    summary = []
    for s in STRATEGIES:
        rs = by[s]
        tn = [int(r["train_N"]) for r in rs]
        summary.append({
            "strategy": s, "drugs": len({r["antimicrobial"] for r in rs}), "tasks": len(rs),
            "train_N_min": min(tn), "train_N_median": statistics.median(tn), "train_N_max": max(tn),
            "train_R_min": min(int(r["train_R"]) for r in rs), "train_S_min": min(int(r["train_S"]) for r in rs),
            "test_N_min": min(int(r["test_N"]) for r in rs),
            "test_N_median": statistics.median(int(r["test_N"]) for r in rs),
            "test_N_max": max(int(r["test_N"]) for r in rs),
            "test_R_min": min(int(r["test_R"]) for r in rs), "test_S_min": min(int(r["test_S"]) for r in rs),
            "tasks_train_N_below_500": sum(t < 500 for t in tn),
            "tasks_train_N_below_1000": sum(t < 1000 for t in tn)})
    write_tsv(out / "split_sizes_by_strategy.tsv", summary)
    small = sorted((r for r in rows if r["strategy"] in HOLDOUTS and int(r["train_N"]) < 1000),
                   key=lambda r: int(r["train_N"]))
    write_tsv(out / "holdout_tasks_train_below_1000.tsv",
              [{k: r[k] for k in ("antimicrobial", "strategy", "fold", "train_N", "train_R", "train_S", "test_N",
                                  "test_R", "test_S")} for r in small])


def pooled_vs_group(out: Path) -> None:
    """Pooled gap (held-out predictions of several models ranked together) against per-group gaps."""
    drops = [r for r in read_tsv(COMMITTED / "generalisation_v1/generalisation_drop.tsv")
             if r["cohort"] == "common" and r["comparator"] == "random" and r["comparison"] == "same_test_isolates"]
    pooled = {(r["antimicrobial"], r["holdout_strategy"]): r for r in drops if r["level"] == "pooled"}
    groups = defaultdict(list)
    for r in drops:
        if r["level"] == "fold" and r["delta_AUROC"] not in ("", "NA"):
            groups[(r["antimicrobial"], r["holdout_strategy"])].append(r)
    rows = []
    for (drug, s), p in sorted(pooled.items()):
        g = groups.get((drug, s), [])
        sup = [r for r in g if int(r["N_R"]) >= 10 and int(r["N_S"]) >= 10]
        n = [int(r["N"]) for r in sup]
        d = [float(r["delta_AUROC"]) for r in sup]
        h = [float(r["AUROC_holdout"]) for r in sup]
        rows.append({
            "antimicrobial": drug, "holdout_strategy": s, "groups": len(g), "groups_supported": len(sup),
            "pooled_AUROC_holdout": float(p["AUROC_holdout"]), "pooled_gap": float(p["delta_AUROC"]),
            "group_AUROC_holdout_mean": statistics.mean(h) if h else None,
            "group_AUROC_holdout_min": min(h) if h else None, "group_AUROC_holdout_max": max(h) if h else None,
            "group_gap_mean": statistics.mean(d) if d else None,
            "group_gap_weighted_mean": sum(x * w for x, w in zip(d, n)) / sum(n) if n else None,
            "group_gap_min": min(d) if d else None, "group_gap_max": max(d) if d else None})
    write_tsv(out / "pooled_vs_group_gap.tsv", rows)
    summ = []
    for s in HOLDOUTS:
        rs = [r for r in rows if r["holdout_strategy"] == s and r["group_gap_mean"] is not None]
        if not rs:
            continue
        diff = [r["pooled_gap"] - r["group_gap_mean"] for r in rs]
        summ.append({"holdout_strategy": s, "drugs": len(rs),
                     "median_pooled_gap": statistics.median(r["pooled_gap"] for r in rs),
                     "median_group_mean_gap": statistics.median(r["group_gap_mean"] for r in rs),
                     "median_group_weighted_gap": statistics.median(r["group_gap_weighted_mean"] for r in rs),
                     "median_pooled_minus_group_mean": statistics.median(diff),
                     "max_abs_pooled_minus_group_mean": max(abs(x) for x in diff),
                     "drugs_abs_difference_above_0.02": sum(abs(x) > 0.02 for x in diff)})
    write_tsv(out / "pooled_vs_group_gap_summary.tsv", summ)


def table4(out: Path) -> None:
    cats = ["R_KNOWN_DET_REPRESENTED", "R_KNOWN_MUT_NOT_REPRESENTED", "R_NO_ESTABLISHED_DET", "R_UNRESOLVED",
            "R_KNOWN_DET_AND_MUT_REPRESENTED", "R_KNOWN_DET_NOT_REPRESENTED", "R_KNOWN_MUT_REPRESENTED"]
    rows = []
    for r in read_tsv(COMMITTED / "error_decomposition_v1/error_decomposition_summary.tsv"):
        if r["cohort"] == "common" and r["antimicrobial"] == "ALL":
            fn = int(r["FN"])
            rows.append({"strategy": r["split_strategy"], "FN": fn, "categorised": sum(int(r[c]) for c in cats),
                         **{c: int(r[c]) for c in cats}, **{f"{c}_fraction": int(r[c]) / fn for c in cats}})
    write_tsv(out / "false_negative_categories_all.tsv", rows)


def mechanism_novelty(out: Path) -> None:
    """Carrier counts per determinant in each task's training set, and alternative novelty definitions."""
    ctx_path = DATASETS / "error_decomposition_v1/prediction_error_context.tsv.gz"
    dets, ctx = {}, []
    with gzip.open(ctx_path, "rt", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cohort"] != "common":
                continue
            d = [x for f in ("gene_determinants", "mutation_determinants")
                 for x in (r[f].split(";") if r[f] not in ("", "NA") else [])]
            dets[(r["antimicrobial"], r["isolate_id"])] = d
            if r["phenotype"] == "R" and r["established_determinant_status"] == "positive":
                ctx.append({k: r[k] for k in ("antimicrobial", "split_strategy", "fold_or_holdout_id", "isolate_id",
                                              "predicted_class", "mechanism_novelty")})
    train = defaultdict(Counter)
    for s in STRATEGIES:
        for r in read_tsv(DATASETS / f"generalisation_v1/splits/{s}.tsv.gz"):
            if r["cohort"] == "common" and r["partition"] == "train":
                for d in dets.get((r["antimicrobial"], r["isolate_id"]), []):
                    train[(r["antimicrobial"], s, r["fold_or_holdout_id"])][d] += 1

    def classes(e):
        c = [train[(e["antimicrobial"], e["split_strategy"], e["fold_or_holdout_id"])][d]
             for d in dets[(e["antimicrobial"], e["isolate_id"])]]
        lo = min(c)
        return {"any_below_min": lo < MIN_CARRIERS, "any_absent": lo == 0,
                "all_below_min": max(c) < MIN_CARRIERS,
                "band": "all_seen" if lo >= MIN_CARRIERS else ("min_0" if lo == 0 else "min_1_4"),
                "unseen_with_seen_determinant": lo < MIN_CARRIERS <= max(c)}

    check = Counter()
    cells = defaultdict(lambda: [0, 0])
    pair = defaultdict(lambda: Counter())
    for e in ctx:
        c = classes(e)
        check[(c["any_below_min"], e["mechanism_novelty"] == "unseen")] += 1
        tp = e["predicted_class"] == "R"
        s = e["split_strategy"]
        for label in (c["band"], "all_below_min" if c["all_below_min"] else None,
                      "unseen_with_seen_determinant" if c["unseen_with_seen_determinant"] else None):
            if label:
                cells[(s, label)][0] += 1
                cells[(s, label)][1] += tp
        p = pair[(e["antimicrobial"], s)]
        p["n"] += 1
        p["any_below_min"] += c["any_below_min"]
        p["any_absent"] += c["any_absent"]
        p["all_below_min"] += c["all_below_min"]
    rows = [{"strategy": s, "class": label, "R_det_positive_endpoints": n, "sensitivity": tp / n}
            for (s, label), (n, tp) in sorted(cells.items())]
    write_tsv(out / "mechanism_novelty_definitions.tsv", rows)
    agree = sum(v for (a, b), v in check.items() if a == b)
    explain = {(r["antimicrobial"], r["holdout_strategy"]): r for r in
               read_tsv(COMMITTED / "error_decomposition_v1/explain_generalisation_drop.tsv")}
    pr, corr = [], []
    for (drug, s), r in sorted(explain.items()):
        p = pair[(drug, s)]
        pr.append({"antimicrobial": drug, "holdout_strategy": s, "gap": float(r["delta_AUROC"]),
                   "R_det_positive": p["n"], "frac_any_below_min": frac(p["any_below_min"], p["n"]),
                   "frac_any_absent": frac(p["any_absent"], p["n"]),
                   "frac_all_below_min": frac(p["all_below_min"], p["n"])})
    write_tsv(out / "mechanism_novelty_by_pair.tsv", pr)
    for f in ("frac_any_below_min", "frac_any_absent", "frac_all_below_min"):
        xs = [(r["gap"], r[f]) for r in pr if r[f] is not None]
        corr.append({"measure": f, "pairs": len(xs),
                     "spearman_rho": spearmanr([a for a, _ in xs], [b for _, b in xs]).statistic})
    write_tsv(out / "mechanism_novelty_correlations.tsv", corr)
    (out / "mechanism_novelty_check.json").write_text(json.dumps(
        {"R_det_positive_endpoints": len(ctx), "agree_with_frozen_flag": agree,
         "disagree": len(ctx) - agree}, indent=1) + "\n")


def correlation_robustness(out: Path) -> None:
    ex = read_tsv(COMMITTED / "error_decomposition_v1/explain_generalisation_drop.tsv")
    rows = []

    def rho(rs, f):
        xs = [(float(r["delta_AUROC"]), float(r[f])) for r in rs if r[f] not in ("", "NA")]
        return spearmanr([a for a, _ in xs], [b for _, b in xs]).statistic, len(xs)

    for f in ("unseen_mechanism_fraction", "lineage_shift"):
        r_all, n = rho(ex, f)
        lodo = [rho([r for r in ex if r["antimicrobial"] != d], f)[0] for d in sorted({r["antimicrobial"] for r in ex})]
        rows.append({"measure": f, "subset": "all pairs", "pairs": n, "spearman_rho": r_all,
                     "lodo_min": min(lodo), "lodo_max": max(lodo)})
        for s in HOLDOUTS:
            v, n = rho([r for r in ex if r["holdout_strategy"] == s], f)
            rows.append({"measure": f, "subset": f"{s} only", "pairs": n, "spearman_rho": v,
                         "lodo_min": None, "lodo_max": None})
        for drop in (("nalidixic acid",), ("nalidixic acid", "cefazolin")):
            v, n = rho([r for r in ex if r["antimicrobial"] not in drop], f)
            rows.append({"measure": f, "subset": "without " + " and ".join(drop), "pairs": n, "spearman_rho": v,
                         "lodo_min": None, "lodo_max": None})
    write_tsv(out / "correlation_robustness.tsv", rows)


def without_flagged_project(out: Path) -> None:
    """Gaps recomputed with the test endpoints of the relabelled project left out (models unchanged)."""
    pred = defaultdict(dict)
    for r in read_tsv(DATASETS / "generalisation_v1/predictions.tsv.gz"):
        if r["cohort"] == "common":
            pred[(r["antimicrobial"], r["split_strategy"])][r["isolate_id"]] = r
    rows = []
    for (drug, s), hold in sorted(pred.items()):
        if s not in HOLDOUTS:
            continue
        rnd = pred[(drug, "random")]
        out_rows = {}
        for label, keep in (("all", lambda r: True), ("without", lambda r: r["project"] != FLAGGED_PROJECT)):
            ids = [i for i, r in hold.items() if keep(r)]
            y = [hold[i]["phenotype"] == "R" for i in ids]
            a_h = auroc(y, [float(hold[i]["predicted_probability"]) for i in ids])
            a_r = auroc(y, [float(rnd[i]["predicted_probability"]) for i in ids])
            out_rows[label] = (len(ids), a_r - a_h if a_h is not None and a_r is not None else None)
        rows.append({"antimicrobial": drug, "holdout_strategy": s, "N": out_rows["all"][0],
                     "N_without": out_rows["without"][0], "gap": out_rows["all"][1],
                     "gap_without": out_rows["without"][1],
                     "flagged_test_endpoints": out_rows["all"][0] - out_rows["without"][0]})
    write_tsv(out / "gap_without_relabelled_project.tsv", rows)
    summ = []
    for s in HOLDOUTS:
        rs = [r for r in rows if r["holdout_strategy"] == s and r["gap_without"] is not None]
        summ.append({"holdout_strategy": s, "drugs": len(rs),
                     "median_gap": statistics.median(r["gap"] for r in rs),
                     "median_gap_without": statistics.median(r["gap_without"] for r in rs),
                     "max_abs_change": max(abs(r["gap"] - r["gap_without"]) for r in rs),
                     "flagged_test_endpoints": sum(r["flagged_test_endpoints"] for r in rs)})
    write_tsv(out / "gap_without_relabelled_project_summary.tsv", summ)


def cohort_flow(out: Path) -> None:
    elig = read_tsv(COMMITTED / "generalisation_v1/cohort_eligibility.tsv")
    excl = Counter()
    for r in elig:
        excl.update(json.loads(r["common_exclusions"]))
    iso_all, iso_common = set(), set()
    for r in read_tsv(DATASETS / "generalisation_v1/splits/random.tsv.gz"):
        (iso_common if r["cohort"] == "common" else iso_all).add(r["isolate_id"])
    flow = {"benchmark_drugs": len(elig),
            "benchmark_endpoints": sum(int(r["genome_verified"]) for r in elig),
            "benchmark_isolates_maximal_random": len(iso_all),
            "common_cohort_endpoints": sum(int(r["common_cohort_N"]) for r in elig),
            "common_cohort_isolates": len(iso_common),
            "common_exclusions_endpoint_reasons": dict(excl)}
    (out / "cohort_flow.json").write_text(json.dumps(flow, indent=1) + "\n")


def representation_effects(out: Path) -> None:
    """Holdout AUROC change, random AUROC change on the same isolates and gap change, side by side."""
    rs = COMMITTED / "representation_sensitivity_v1"
    gaps = {(r["representation"], r["antimicrobial"], r["holdout_strategy"]): r
            for r in read_tsv(rs / "generalisation_gap_by_representation.tsv")}
    deltas = {(r["representation"], r["antimicrobial"], r["split_strategy"]): r
              for r in read_tsv(rs / "representation_metric_deltas.tsv")}
    rescue = {(r["representation"], r["antimicrobial"], r["holdout_strategy"]): r
              for r in read_tsv(rs / "generalisation_gap_rescue.tsv")}
    rows = []
    for (rep, drug, s), g in sorted(gaps.items()):
        if rep == "R0":
            continue
        g0 = gaps[("R0", drug, s)]
        d = deltas.get((rep, drug, s), {})
        rr = rescue.get((rep, drug, s), {})
        rows.append({
            "representation": rep, "antimicrobial": drug, "holdout_strategy": s, "N": g["N"], "N_R": g["N_R"],
            "R0_holdout_AUROC": float(g0["AUROC_holdout"]), "Rx_holdout_AUROC": float(g["AUROC_holdout"]),
            "holdout_AUROC_change": float(g["AUROC_holdout"]) - float(g0["AUROC_holdout"]),
            "holdout_change_ci_lower": d.get("delta_AUROC_ci_lower"), "holdout_change_ci_upper": d.get("delta_AUROC_ci_upper"),
            "R0_random_AUROC_same_isolates": float(g0["AUROC_random"]),
            "Rx_random_AUROC_same_isolates": float(g["AUROC_random"]),
            "random_AUROC_change": float(g["AUROC_random"]) - float(g0["AUROC_random"]),
            "R0_gap": float(g0["gap_AUROC"]), "Rx_gap": float(g["gap_AUROC"]),
            "gap_change": float(g["gap_AUROC"]) - float(g0["gap_AUROC"]),
            "rescue": rr.get("gap_rescue"), "rescue_ci_lower": rr.get("gap_rescue_ci_lower"),
            "rescue_ci_upper": rr.get("gap_rescue_ci_upper"),
            "R0_holdout_sensitivity": float(g0["sensitivity_holdout"]),
            "Rx_holdout_sensitivity": float(g["sensitivity_holdout"])})
    write_tsv(out / "representation_three_way.tsv", rows)


def baseline_exceptions(out: Path) -> None:
    rows = []
    for r in read_tsv(COMMITTED / "generalisation_v1/determinant_vs_ml.tsv"):
        if r["cohort"] != "common":
            continue
        lo, hi = float(r["delta_balanced_accuracy_ci_lower"]), float(r["delta_balanced_accuracy_ci_upper"])
        side = "model_better" if lo > 0 else "baseline_better" if hi < 0 else "interval_includes_zero"
        rows.append({"strategy": r["split_strategy"], "antimicrobial": r["antimicrobial"], "N": r["N"],
                     "baseline_balanced_accuracy": r["determinant_balanced_accuracy"],
                     "model_balanced_accuracy": r["ml_balanced_accuracy"],
                     "difference_model_minus_baseline": r["delta_balanced_accuracy"],
                     "ci_lower": lo, "ci_upper": hi, "classification": side})
    write_tsv(out / "baseline_vs_model_by_drug.tsv", sorted(rows, key=lambda r: (r["strategy"], r["antimicrobial"])))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    for step in (split_sizes, pooled_vs_group, table4, correlation_robustness, representation_effects,
                 baseline_exceptions, cohort_flow, without_flagged_project, mechanism_novelty):
        step(a.out)
        print("done", step.__name__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
