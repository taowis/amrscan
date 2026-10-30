"""Offline tests for the generalisation benchmark (synthetic endpoints, genomes and matrices)."""

from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
spec = importlib.util.spec_from_file_location("generalisation_benchmark", ROOT / "bin/generalisation_benchmark.py")
gb = importlib.util.module_from_spec(spec)
sys.modules["generalisation_benchmark"] = gb  # process-pool workers pickle its functions
spec.loader.exec_module(gb)
pop = gb.pop
NA = gb.NA
CONFIG = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/generalisation_v1.json").read_text())
POP_CONFIG = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/population_config.json").read_text())


def endpoint(n, phen="R", st="ST1", cluster=NA, country="USA", year=2015, project="P1", source="human",
             component=None, test_ok=True, verified=True, drug="ampicillin", det="positive", cat=None, tier="primary"):
    iid = f"PDT{n:09d}.1"
    return {"isolate_id": iid, "assembly_accession": f"GCA_{n:09d}.1", "antimicrobial": drug, "phenotype": phen,
            "lineage_id": st, "lineage_component": component or (st if test_ok else NA), "lineage_test_eligible": test_ok,
            "snp_cluster": cluster, "split_unit": cluster if cluster != NA else iid, "country": country,
            "year": year, "period": NA if year is None else ("recent" if year >= 2021 else "earlier"),
            "source": source, "project": project, "geography_group": NA, "genome_verified": verified,
            "determinant_status": det, "concordance_category": cat or f"{phen}_DET_{'POS' if det == 'positive' else 'NEG'}",
            "tier": tier}


def world(n=200, seed=1, drug="ampicillin"):
    """Synthetic endpoints over 10 STs, 3 countries, 4 projects and 2012-2024."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        st = f"ST{i % 10}"
        out.append(endpoint(i, phen="R" if rng.random() < 0.4 else "S", st=st,
                            cluster=f"PDS{i % 10}{(i // 10) % 3}" if i % 3 else NA, country=("USA", "Canada", "Netherlands")[i % 3],
                            year=2012 + i % 13, project=f"P{i % 4}", drug=drug))
    return with_components(out)


def with_components(eps):
    """Decision 0001 lineage components (STs joined through shared SNP clusters)."""
    comp = pop.lineage_components([{"isolate_id": e["isolate_id"], "eligible_population_analysis": "true",
                                    "lineage_id": e["lineage_id"], "ncbi_snp_cluster": e["snp_cluster"]} for e in eps])
    for e in eps:
        e["lineage_component"] = comp[e["isolate_id"]]
        e["lineage_test_eligible"] = comp[e["isolate_id"]] != NA
    return eps


def plan(drug="ampicillin", tier="primary", **flags):
    row = {"antimicrobial": drug, "tier": tier, "project_holdout": "yes", "geographic_transfer": "yes", "temporal": "yes"}
    row.update(flags)
    return row


def cfg(**changes):
    config = json.loads(json.dumps(CONFIG))
    config["model"]["C_grid"] = [0.3]
    for path, value in changes.items():
        node = config
        keys = path.split(".")
        for k in keys[:-1]:
            node = node[k]
        node[keys[-1]] = value
    return config


class SplitTests(unittest.TestCase):
    def setUp(self):
        self.eps = world()
        self.config = cfg()

    def tasks(self, strategy, cohort="common", eps=None, config=None, plan_row=None):
        return gb.strategy_tasks(eps or self.eps, "ampicillin", cohort, strategy, config or self.config,
                                 plan_row or plan())

    def test_random_split_is_deterministic_and_seeded(self):
        a, _ = self.tasks("random")
        b, _ = self.tasks("random")
        self.assertEqual([gb.task_digest(t) for t in a], [gb.task_digest(t) for t in b])
        c, _ = self.tasks("random", config=cfg(**{"seeds.split": 7}))
        self.assertNotEqual([gb.task_digest(t) for t in a], [gb.task_digest(t) for t in c])

    def test_random_split_is_isolate_level_and_stratified(self):
        tasks, _ = self.tasks("random")
        self.assertEqual(len(tasks), 5)
        tested = Counter(e["isolate_id"] for t in tasks for e in t["test"])
        self.assertEqual(set(tested.values()), {1})
        self.assertEqual(len(tested), len(self.eps))
        r = [gb.rs(t["test"])[0] for t in tasks]
        s = [gb.rs(t["test"])[1] for t in tasks]
        self.assertLessEqual(max(r) - min(r), 1)
        self.assertLessEqual(max(s) - min(s), 1)
        for t in tasks:
            gb.validate_task(t, self.config)
            self.assertFalse({e["isolate_id"] for e in t["train"]} & {e["isolate_id"] for e in t["test"]})
        gb.validate_cv(tasks, self.eps)

    def test_random_split_can_divide_snp_clusters_but_cluster_split_cannot(self):
        rnd, _ = self.tasks("random")
        audit = [gb.leakage_audit(t) for t in rnd]
        self.assertTrue(any(a["test_snp_cluster_in_train"] > 0 for a in audit))
        grouped, _ = self.tasks("random_cluster")
        for t in grouped:
            gb.validate_task(t, self.config)
            self.assertEqual(gb.leakage_audit(t)["test_snp_cluster_in_train"], 0)

    def test_lineage_holdout_excludes_test_lineages(self):
        tasks, _ = self.tasks("lineage")
        for t in tasks:
            gb.validate_task(t, self.config)
            self.assertFalse(gb.known(e["lineage_id"] for e in t["train"]) & gb.known(e["lineage_id"] for e in t["test"]))
            self.assertEqual(gb.leakage_audit(t)["fraction_st_in_train"], 0)
        gb.validate_cv(tasks, self.eps)

    def test_missing_lineage_is_training_only(self):
        eps = self.eps + [endpoint(900, st="UNKNOWN", test_ok=False)]
        common, reasons = gb.cohort_members(eps, "common", "lineage")
        self.assertEqual(reasons, {"no_lineage_component": 1})
        tasks, _ = self.tasks("lineage", cohort="maximal", eps=eps)
        self.assertTrue(all(any(e["isolate_id"] == "PDT000000900.1" for e in t["train"]) for t in tasks))
        self.assertFalse(any(e["isolate_id"] == "PDT000000900.1" for t in tasks for e in t["test"]))
        gb.validate_cv(tasks, eps)

    def test_country_holdout_excludes_test_country_and_unknown_country(self):
        eps = self.eps + [endpoint(901, country=NA)]
        tasks, _ = self.tasks("country", cohort="maximal", eps=eps)
        self.assertEqual(sorted(t["fold"] for t in tasks), ["Canada", "Netherlands", "USA"])
        for t in tasks:
            gb.validate_task(t, self.config)
            self.assertNotIn(t["fold"], {e["country"] for e in t["train"]})
            self.assertNotIn(NA, {e["country"] for e in t["train"]})
        _, reasons = gb.cohort_members(eps, "maximal", "country")
        self.assertEqual(reasons, {"missing_country": 1})

    def test_project_and_geography_holdouts(self):
        tasks, _ = self.tasks("project")
        self.assertEqual(len(tasks), 4)
        for t in tasks:
            gb.validate_task(t, self.config)
        geo, _ = self.tasks("geography")
        self.assertEqual(geo[0]["fold"], "USA_to_other")
        gb.validate_task(geo[0], self.config)
        self.assertEqual({e["country"] for e in geo[0]["train"]}, {"USA"})

    def test_temporal_training_precedes_test_and_missing_year(self):
        eps = self.eps + [endpoint(902, year=None)]
        tasks, _ = self.tasks("temporal", cohort="maximal", eps=eps)
        t = tasks[0]
        gb.validate_task(t, self.config)
        self.assertLess(max(e["year"] for e in t["train"]), min(e["year"] for e in t["test"]))
        self.assertEqual(min(e["year"] for e in t["test"]), 2021)
        _, reasons = gb.cohort_members(eps, "common", "random")
        self.assertEqual(reasons, {"missing_year": 1})

    def test_insufficient_r_or_s_is_not_feasible(self):
        few_r = [e for e in self.eps if not (e["year"] >= 2021 and e["phenotype"] == "R")][:-1]
        few_r += [endpoint(950 + i, year=2022) for i in range(3)]
        tasks, skipped = self.tasks("temporal", eps=few_r)
        self.assertEqual(tasks, [])
        self.assertEqual(skipped[0]["status"], "NOT_FEASIBLE")
        self.assertIn("test 3 R", skipped[0]["reason"])
        few_s = [e for e in self.eps if not (e["year"] < 2021 and e["phenotype"] == "S")]
        tasks, skipped = self.tasks("temporal", eps=few_s)
        self.assertIn("training", skipped[0]["reason"])

    def test_undersized_holdout_is_not_a_task(self):
        eps = self.eps + [endpoint(960 + i, country="Zambia", phen="R") for i in range(12)] + \
            [endpoint(980 + i, country="Zambia", phen="S") for i in range(3)]
        tasks, skipped = self.tasks("country", eps=eps)
        self.assertNotIn("Zambia", [t["fold"] for t in tasks])
        tasks, skipped = self.tasks("country", eps=self.eps[:20])
        self.assertEqual(tasks, [])
        self.assertEqual(skipped[0]["fold"], "all")

    def test_tier_and_plan_flags_gate_strategies(self):
        tasks, skipped = self.tasks("temporal", plan_row=plan(tier="secondary"))
        self.assertEqual(tasks, [])
        self.assertIn("secondary drugs", skipped[0]["reason"])
        tasks, skipped = self.tasks("temporal", plan_row=plan(temporal="no"))
        self.assertIn("plan temporal = no", skipped[0]["reason"])
        tasks, _ = self.tasks("lineage", plan_row=plan(tier="secondary"))
        self.assertEqual(len(tasks), 5)

    def test_build_tasks_covers_drugs_cohorts_and_strategies(self):
        eps = self.eps + world(150, seed=3, drug="meropenem")
        tasks, skipped = gb.build_tasks(eps, {"ampicillin": plan(), "meropenem": plan("meropenem", "secondary")},
                                        self.config)
        strategies = Counter((t["antimicrobial"], t["strategy"]) for t in tasks)
        self.assertEqual(strategies[("meropenem", "lineage")], 10)
        self.assertNotIn(("meropenem", "country"), strategies)
        self.assertEqual(strategies[("ampicillin", "temporal")], 2)
        self.assertTrue(all(s["antimicrobial"] == "meropenem" for s in skipped))

    def test_reconciliation_and_validation_failures(self):
        tasks, _ = self.tasks("random")
        broken = [dict(t) for t in tasks]
        broken[0] = {**broken[0], "test": broken[0]["test"] + broken[1]["test"][:1]}
        with self.assertRaisesRegex(gb.LeakageError, "more than one fold"):
            gb.validate_cv(broken, self.eps)
        short = [{**t, "test": t["test"][1:]} if i == 0 else t for i, t in enumerate(tasks)]
        with self.assertRaisesRegex(gb.LeakageError, "do not reconcile"):
            gb.validate_cv(short, self.eps)
        lost = [{**t, "train": t["train"][1:]} if i == 0 else t for i, t in enumerate(tasks)]
        with self.assertRaisesRegex(gb.LeakageError, "add up"):
            gb.validate_cv(lost, self.eps)

    def test_leakage_detection(self):
        t = self.tasks("random")[0][0]
        a, b = t["test"][0], t["train"][0]
        cases = [
            ({**t, "train": t["train"] + [a]}, "both training and test"),
            ({**t, "test": t["test"] + [a]}, "duplicate isolate"),
            ({**t, "train": t["train"] + [{**a, "isolate_id": "X", "genome_verified": False}]}, "verified genome"),
            ({**t, "train": t["train"] + [{**b, "isolate_id": "Y"}]}, "one assembly"),
        ]
        lineage = self.tasks("lineage")[0][0]
        cases += [({**lineage, "train": lineage["train"] + [{**lineage["test"][0], "isolate_id": "Z",
                                                              "assembly_accession": "GCA_Z"}]}, "lineage component"),
                  ({**lineage, "train": lineage["train"] + [{**lineage["train"][0], "isolate_id": "Z2",
                                                              "assembly_accession": "GCA_Z2",
                                                              "lineage_id": lineage["test"][0]["lineage_id"]}]},
                   "held-out ST"),
                  ({**lineage, "train": lineage["train"] + [{**lineage["train"][0], "isolate_id": "Z3",
                                                              "assembly_accession": "GCA_Z3", "lineage_id": "ST99",
                                                              "snp_cluster": "PDSX"}],
                    "test": lineage["test"] + [{**lineage["test"][0], "isolate_id": "Z4", "assembly_accession": "GCA_Z4",
                                                "snp_cluster": "PDSX"}]}, "SNP cluster"),
                  ({**lineage, "test": lineage["test"] + [{**lineage["test"][0], "isolate_id": "Z5",
                                                           "assembly_accession": "GCA_Z5", "lineage_test_eligible": False}]},
                   "without a lineage component")]
        country = self.tasks("country")[0][0]
        cases.append(({**country, "train": country["train"] + [{**country["test"][0], "isolate_id": "C",
                                                                 "assembly_accession": "GCA_C"}]}, "held-out country"))
        project = self.tasks("project")[0][0]
        cases.append(({**project, "train": project["train"] + [{**project["test"][0], "isolate_id": "P",
                                                                 "assembly_accession": "GCA_P"}]}, "held-out project"))
        geo = self.tasks("geography")[0][0]
        cases.append(({**geo, "train": geo["train"] + [{**geo["test"][0], "isolate_id": "G",
                                                         "assembly_accession": "GCA_G"}]}, "mixes countries"))
        temporal = self.tasks("temporal")[0][0]
        cases.append(({**temporal, "train": temporal["train"] + [{**temporal["test"][0], "isolate_id": "T",
                                                                   "assembly_accession": "GCA_T"}]}, "do not precede"))
        grouped = self.tasks("random_cluster")[0][0]
        clustered = next(e for e in grouped["test"] if e["snp_cluster"] != NA)
        cases.append(({**grouped, "train": grouped["train"] + [{**clustered, "isolate_id": "R",
                                                                 "assembly_accession": "GCA_R"}]}, "split unit"))
        for task, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(gb.LeakageError, message):
                gb.validate_task(task, self.config)

    def test_snp_cluster_overlap_audit(self):
        train = [endpoint(1, cluster="PDS1"), endpoint(2, cluster="PDS2", st="ST2")]
        test = [endpoint(3, cluster="PDS1"), endpoint(4, cluster="PDS3"), endpoint(5, st="UNKNOWN")]
        task = gb.make_task("common", "ampicillin", "random", "fold1", train, test)
        audit = gb.leakage_audit(task)
        self.assertEqual((audit["test_with_snp_cluster"], audit["test_snp_cluster_in_train"]), (2, 1))
        self.assertEqual(audit["fraction_snp_cluster_in_train"], 0.5)
        self.assertEqual((audit["test_with_st"], audit["test_st_in_train"]), (2, 2))
        nov = gb.novelty(task)
        self.assertEqual(nov["PDT000000004.1"]["snp_cluster_seen_in_training"], "no")
        self.assertEqual(nov["PDT000000005.1"]["lineage_seen_in_training"], NA)
        self.assertEqual(nov["PDT000000003.1"]["year_after_training_range"], "no")

    def test_fold_assignment_balances_groups(self):
        folds = gb.assign_folds({"a": (10, 0), "b": (0, 10), "c": (5, 5), "d": (1, 1), "e": (1, 0)}, 2, 1, "x")
        self.assertEqual(set(folds), set("abcde"))
        load = [[0, 0], [0, 0]]
        for u, f in folds.items():
            load[f][0] += {"a": 10, "b": 0, "c": 5, "d": 1, "e": 1}[u]
        self.assertLessEqual(abs(load[0][0] - load[1][0]), 10)


class CommonCohortTests(unittest.TestCase):
    def test_common_cohort_intersection_and_reasons(self):
        eps = world(60) + [endpoint(700, st="UNKNOWN", test_ok=False), endpoint(701, country=NA),
                           endpoint(702, year=None), endpoint(703, verified=False)]
        common, reasons = gb.cohort_members(eps, "common", "random")
        self.assertEqual(len(common), 60)
        self.assertEqual(reasons, {"no_lineage_component": 1, "missing_country": 1, "missing_year": 1,
                                   "genome_not_verified": 1})

    def test_maximal_cohort_differs_by_strategy(self):
        eps = world(60) + [endpoint(701, country=NA), endpoint(702, year=None)]
        sizes = {s: len(gb.cohort_members(eps, "maximal", s)[0]) for s in ("random", "country", "temporal", "project")}
        self.assertEqual(sizes, {"random": 62, "country": 61, "temporal": 61, "project": 62})

    def test_common_cohort_is_the_same_isolate_set_for_cv_strategies(self):
        eps = world(120) + [endpoint(701, country=NA)]
        config = cfg()
        sets = []
        for strategy in gb.CV_STRATEGIES:
            tasks, _ = gb.strategy_tasks(eps, "ampicillin", "common", strategy, config, plan())
            sets.append({e["isolate_id"] for t in tasks for e in t["test"]})
        self.assertEqual(sets[0], sets[1])
        self.assertEqual(sets[1], sets[2])
        self.assertNotIn("PDT000000701.1", sets[0])

    def test_eligibility_table_reconciles(self):
        eps = world(60) + [endpoint(701, country=NA), endpoint(703, verified=False)]
        row = gb.eligibility_table(eps, {"ampicillin": plan()})[0]
        self.assertEqual(row["endpoints"], 62)
        self.assertEqual(row["genome_verified"], 61)
        self.assertEqual(row["R"] + row["S"], 61)
        self.assertEqual(row["common_cohort_N"], 60)
        self.assertEqual(row["country_eligible_N"], 60)
        self.assertEqual(row["random_eligible_N"], 61)
        self.assertEqual(json.loads(row["common_exclusions"]), {"genome_not_verified": 1, "missing_country": 1})

    def test_cohort_labels_in_task_ids(self):
        config = cfg()
        eps = world(120)
        tasks, _ = gb.build_tasks(eps, {"ampicillin": plan()}, config)
        self.assertEqual({t["cohort"] for t in tasks}, {"common", "maximal"})
        self.assertTrue(all(t["task_id"].startswith(t["cohort"] + "__") for t in tasks))


class MetricTests(unittest.TestCase):
    y = np.array([1, 1, 0, 0, 1, 0])
    p = np.array([0.9, 0.4, 0.6, 0.2, 0.8, 0.1])

    def test_auroc_by_hand(self):
        # positives 0.9 0.4 0.8, negatives 0.6 0.2 0.1: 8 of 9 pairs ordered (0.4 < 0.6)
        self.assertAlmostEqual(gb.auroc(self.y, self.p), 8 / 9)
        self.assertAlmostEqual(gb.auroc([1, 0], [0.5, 0.5]), 0.5)

    def test_auprc_by_hand(self):
        # ranked: 0.9 R, 0.8 R, 0.6 S, 0.4 R -> AP = (1 + 1 + 3/4) / 3
        self.assertAlmostEqual(gb.auprc(self.y, self.p), (1 + 1 + 0.75) / 3)
        # tied scores form one threshold
        self.assertAlmostEqual(gb.auprc([1, 0, 1], [0.5, 0.5, 0.1]), 0.5 * 0.5 + 0.5 * (2 / 3))

    def test_confusion_metrics_by_hand(self):
        m = gb.confusion_metrics(self.y, self.p >= 0.5)  # TP 2, FN 1, FP 1, TN 2
        self.assertEqual((m["TP"], m["FN"], m["FP"], m["TN"]), (2, 1, 1, 2))
        self.assertAlmostEqual(m["sensitivity"], 2 / 3)
        self.assertAlmostEqual(m["specificity"], 2 / 3)
        self.assertAlmostEqual(m["PPV"], 2 / 3)
        self.assertAlmostEqual(m["NPV"], 2 / 3)
        self.assertAlmostEqual(m["balanced_accuracy"], 2 / 3)
        self.assertAlmostEqual(m["F1"], 4 / 6)
        self.assertAlmostEqual(m["MCC"], (4 - 1) / 9)

    def test_brier_and_calibration(self):
        m = gb.metrics(self.y, self.p, 0.5)
        self.assertAlmostEqual(m["brier"], float(np.mean((self.p - self.y) ** 2)))
        self.assertEqual(m["N"], 6)
        intercept, slope = gb.calibration([0, 1] * 50, [0.3, 0.7] * 50)
        self.assertGreater(slope, 0)
        self.assertEqual(gb.calibration([1, 1], [0.2, 0.3]), (None, None))

    def test_undefined_metrics(self):
        self.assertIsNone(gb.auroc([1, 1], [0.2, 0.9]))
        self.assertIsNone(gb.auprc([0, 0], [0.2, 0.9]))
        m = gb.confusion_metrics([1, 1], [True, True])
        self.assertIsNone(m["specificity"])
        self.assertIsNone(m["balanced_accuracy"])
        self.assertIsNone(m["MCC"])
        self.assertIsNone(gb.confusion_metrics([], [])["F1"])
        self.assertEqual(gb.fmt(None), NA)
        self.assertEqual(gb.fmt(float("nan")), NA)
        self.assertEqual(gb.fmt(0.5), "0.5")
        self.assertEqual(gb.fmt(2.0), "2")

    def test_bootstrap_is_reproducible_and_resamples_units(self):
        config = cfg(**{"bootstrap.replicates": 50})
        y = np.array([1, 0] * 20)
        p = np.linspace(0, 1, 40)
        units = [f"u{i // 4}" for i in range(40)]
        a = gb.bootstrap_ci(y, p, units, 0.5, config)
        b = gb.bootstrap_ci(y, p, units, 0.5, config)
        self.assertEqual(a, b)
        self.assertLessEqual(a["AUROC_ci"][0], a["AUROC_ci"][1])
        for idx in gb.bootstrap_indices(units, 3, 1):
            picked = Counter(units[i] for i in idx)
            self.assertTrue(all(c % 4 == 0 for c in picked.values()))
            self.assertEqual(len(idx) % 4, 0)
        c = gb.bootstrap_ci(y, p, units, 0.5, cfg(**{"bootstrap.replicates": 50, "seeds.bootstrap": 3}))
        self.assertNotEqual(a, c)

    def test_paired_delta(self):
        config = cfg(**{"bootstrap.replicates": 30})
        y = np.array([1, 0] * 10)
        good = np.where(y == 1, 0.9, 0.1)
        bad = np.full(20, 0.5)
        ci = gb.bootstrap_ci(y, bad, list(range(20)), 0.5, config, q=good)
        self.assertGreater(ci["delta_AUROC_ci"][0], 0)
        self.assertEqual(gb.summarise([None])["folds_defined"], 0)
        s = gb.summarise([0.6, 0.8, 0.7])
        self.assertEqual((s["median"], s["min"], s["max"]), (0.7, 0.6, 0.8))
        self.assertAlmostEqual(gb.tvd({"a": 1.0}, {"b": 1.0}), 1.0)


class FeatureAndModelTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.n = 80
        self.y = np.array([1] * 40 + [0] * 40, dtype=np.int8)
        signal = self.y.reshape(-1, 1).repeat(3, axis=1)
        noise = rng.random((self.n, 20)) < 0.4
        const = np.ones((self.n, 2), dtype=bool)
        dense = np.hstack([signal, noise, const]).astype(np.int8)
        self.X = sparse.csr_matrix(dense)
        self.model = cfg()["model"]
        self.features = cfg()["features"]

    def test_feature_filter_counts_training_rows_only(self):
        dense = np.zeros((10, 3), dtype=np.int8)
        dense[:6, 0] = 1      # present in 6 training rows
        dense[6:, 1] = 1      # present only in test rows
        dense[:, 2] = 1       # constant
        cols = gb.feature_filter(sparse.csr_matrix(dense), np.arange(6), 2, 0)
        self.assertEqual(cols.tolist(), [0, 2])
        cols = gb.feature_filter(sparse.csr_matrix(dense), np.arange(6), 2, 1)
        self.assertEqual(cols.tolist(), [])
        cols = gb.feature_filter(self.X, np.arange(self.n), 5, 5)
        self.assertNotIn(self.X.shape[1] - 1, cols.tolist())

    def test_model_trains_predicts_and_is_deterministic(self):
        a = gb.fit_model(self.X, self.y, 1.0, self.model, 1)
        b = gb.fit_model(self.X, self.y, 1.0, self.model, 1)
        pa = a.predict_proba(self.X)[:, 1]
        self.assertTrue(np.allclose(pa, b.predict_proba(self.X)[:, 1]))
        self.assertTrue(((pa >= 0) & (pa <= 1)).all())
        self.assertGreater(gb.auroc(self.y, pa), 0.95)

    def test_single_class_training_fails(self):
        with self.assertRaises(gb.SingleClassError):
            gb.fit_model(self.X, np.ones(self.n, dtype=np.int8), 1.0, self.model, 1)

    def test_select_c_uses_training_rows(self):
        model = {**self.model, "C_grid": [0.01, 1.0], "inner_folds": 3}
        units = [f"u{i // 2}" for i in range(self.n)]
        C, scores = gb.select_C(self.X, self.y, units, model, self.features, 1)
        self.assertIn(C, (0.01, 1.0))
        self.assertEqual([s[0] for s in scores], [0.01, 1.0])
        self.assertEqual(gb.select_C(self.X, self.y, units, {**model, "C_grid": [0.5]}, self.features, 1), (0.5, []))
        one = np.array([1] + [0] * (self.n - 1), dtype=np.int8)
        self.assertEqual(gb.select_C(self.X, one, units, model, self.features, 1)[0], 0.01)

    def task(self, flip_test=False):
        eps = [endpoint(i, phen="R" if self.y[i] else "S", cluster=f"PDS{i // 2}") for i in range(self.n)]
        test = eps[::4]
        if flip_test:
            test = [{**e, "phenotype": "S" if e["phenotype"] == "R" else "R"} for e in test]
        train = [e for i, e in enumerate(eps) if i % 4]
        return gb.make_task("common", "ampicillin", "random", "fold1", train, test)

    def test_run_task_ignores_test_labels_and_records_provenance(self):
        config = cfg()
        config["model"]["C_grid"] = [0.1, 1.0]
        row_of = {f"PDT{i:09d}.1": i for i in range(self.n)}
        status, preds = gb.run_task(self.task(), self.X, row_of, config)
        status2, preds2 = gb.run_task(self.task(flip_test=True), self.X, row_of, config)
        self.assertEqual([p["predicted_probability"] for p in preds], [p["predicted_probability"] for p in preds2])
        self.assertEqual(status["C"], status2["C"])
        self.assertEqual(status["vocabulary_sha256"], status2["vocabulary_sha256"])
        self.assertEqual(status["n_test"], 20)
        self.assertEqual(status["threshold"], 0.5)
        for p in preds:
            self.assertEqual(p["predicted_class"], "R" if p["predicted_probability"] >= 0.5 else "S")
        with self.assertRaises(gb.SingleClassError):
            t = self.task()
            gb.run_task({**t, "train": [{**e, "phenotype": "S"} for e in t["train"]]}, self.X, row_of, config)

    def test_missing_features_still_predict(self):
        row_of = {f"PDT{i:09d}.1": i for i in range(self.n)}
        dense = self.X.toarray()
        dense[0, :] = 0  # a test genome with none of the training features
        status, preds = gb.run_task(self.task(), sparse.csr_matrix(dense), row_of, cfg())
        self.assertTrue(0 < preds[0]["predicted_probability"] < 1)

    def test_build_matrix_prunes_and_merges_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            sk = {"A": [1, 2, 3, 10], "B": [1, 2, 3], "C": [1, 4], "D": [1, 2, 3, 4, 2 ** 63 + 5]}
            for name, h in sk.items():
                np.save(d / f"GCA_{name}.npy", np.array(sorted(h), dtype=np.uint64))
            rows = [(f"I{n}", f"GCA_{n}") for n in sk]
            built = gb.build_matrix(d, rows, 1, 1, dedup=True, buckets=4)
            # 1 is in every row (pruned); 2 and 3 share a pattern; 10 and 2^63+5 are single-row features
            self.assertEqual(built["observed_features"], 6)
            self.assertEqual(built["kept_features"], 5)
            self.assertEqual(built["patterns"], 4)
            m = built["matrix"].toarray()
            self.assertEqual(m.shape, (4, 4))
            self.assertEqual(sorted(m.sum(axis=0).tolist()), [1, 1, 2, 3])
            nodedup = gb.build_matrix(d, rows, 1, 1, dedup=False, buckets=4)
            self.assertEqual(nodedup["patterns"], 5)
            strict = gb.build_matrix(d, rows, 2, 1, dedup=True, buckets=4)
            self.assertEqual(strict["patterns"], 2)
            with self.assertRaises(ValueError):
                gb.build_matrix(d, rows, 4, 4, dedup=True, buckets=4)
            down = gb.build_matrix(d, rows, 1, 1, dedup=True, buckets=4, scaled=2)
            self.assertLess(down["kept_features"], built["kept_features"])
            out = d / "matrix"
            manifest = gb.save_matrix(out, built, rows, {"params": {"k": 1}})
            X, row_of, meta = gb.load_matrix(out, check=True)
            self.assertEqual(X.toarray().tolist(), m.tolist())
            self.assertEqual(row_of["IA"], 0)
            self.assertEqual(manifest["feature_hashes_sha256"], hashlib.sha256(built["hashes"].tobytes()).hexdigest())
            np.save(out / "indices.npy", np.zeros_like(np.load(out / "indices.npy")))
            with self.assertRaisesRegex(ValueError, "differs from its manifest"):
                gb.load_matrix(out, check=True)

    def test_sketch_genome_is_reproducible(self):
        params = gb.sketch_params(cfg(**{"features.scaled": 1, "features.ksize": 21}))
        data = (ROOT / "tests/fixtures/tiny.fna").read_bytes()
        a, b = gb.sketch_genome(data, params), gb.sketch_genome(data, params)
        self.assertTrue(np.array_equal(a, b))
        self.assertTrue((np.diff(a.astype(np.float64)) > 0).all())
        with self.assertRaises(RuntimeError):
            gb.sketch_genome(data, {**params, "tool_version": "0.0"})
        with self.assertRaises(ValueError):
            gb.parse_fasta(b">only\n")


class UpstreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def manifest(self, rows):
        path = self.root / "manifest.tsv"
        fields = ["isolate_id", "assembly_accession", "biosample_accession", "genome_sha256", "sequence_records",
                  "sequence_bases"]
        with path.open("w", newline="") as h:
            w = csv.DictWriter(h, fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        return path

    def row(self, n, digest=None):
        return {"isolate_id": f"I{n}", "assembly_accession": f"GCA_{n:09d}.1", "biosample_accession": f"S{n}",
                "genome_sha256": digest or hashlib.sha256(str(n).encode()).hexdigest(), "sequence_records": "1",
                "sequence_bases": "10"}

    def test_genome_manifest_checks(self):
        path = self.manifest([self.row(1), self.row(2)])
        self.assertEqual(len(gb.read_genome_manifest(path, pop.sha256(path), 2)), 2)
        with self.assertRaisesRegex(ValueError, "SHA-256 differs"):
            gb.read_genome_manifest(path, "0" * 64)
        with self.assertRaisesRegex(ValueError, "expected 3"):
            gb.read_genome_manifest(path, expected_rows=3)
        with self.assertRaisesRegex(ValueError, "twice"):
            gb.read_genome_manifest(self.manifest([self.row(1), {**self.row(2), "isolate_id": "I1"}]))
        with self.assertRaisesRegex(ValueError, "listed for"):
            gb.read_genome_manifest(self.manifest([self.row(1), {**self.row(2), "assembly_accession": "GCA_000000001.1"}]))
        with self.assertRaisesRegex(ValueError, "unversioned"):
            gb.read_genome_manifest(self.manifest([{**self.row(1), "assembly_accession": "GCA_000000001"}]))
        with self.assertRaisesRegex(ValueError, "lacks"):
            gb.read_genome_manifest(self._missing_field())
        with self.assertRaisesRegex(ValueError, "malformed"):
            bad = self.root / "bad.tsv"
            bad.write_text("isolate_id\tassembly_accession\nI1\n")
            gb.read_genome_manifest(bad)

    def _missing_field(self):
        path = self.root / "missing.tsv"
        path.write_text("isolate_id\tassembly_accession\nI1\tGCA_000000001.1\n")
        return path

    def test_duplicate_genome_detection(self):
        same = hashlib.sha256(b"x").hexdigest()
        manifest = gb.read_genome_manifest(self.manifest([self.row(1, same), self.row(2, same), self.row(3)]))
        self.assertEqual(gb.duplicate_genomes(manifest), [["I1", "I2"]])

    def test_genome_verification_classes(self):
        cache = self.root / "cache"
        (cache / "GCA_000000001.1").mkdir(parents=True)
        (cache / "GCA_000000001.1/genome.fna").write_bytes(b">a\nACGT\n")
        good = self.row(1, hashlib.sha256(b">a\nACGT\n").hexdigest())
        self.assertEqual(gb.verify_genome(cache, good)[0], "verified")
        self.assertEqual(gb.verify_genome(cache, {**good, "genome_sha256": "0" * 64})[0], "genome_hash_mismatch")
        self.assertEqual(gb.verify_genome(cache, self.row(2))[0], "genome_missing")
        self.assertEqual(gb.verify_genome(cache, good, "GCA_000000001.2")[0], "genome_version_mismatch")
        with patch.object(Path, "read_bytes", side_effect=OSError("io")):
            self.assertEqual(gb.verify_genome(cache, good)[0], "genome_unreadable")
        params = gb.sketch_params(cfg(**{"features.scaled": 1, "features.ksize": 3}))
        out = self.root / "sk"
        out.mkdir()
        first = gb.sketch_one(str(cache), str(out), {**good, "isolate_id": "I1"}, None, params)
        self.assertEqual((first["status"], first["action"]), ("verified", "sketched"))
        again = gb.sketch_one(str(cache), str(out), {**good, "isolate_id": "I1"}, None, params)
        self.assertEqual(again["action"], "reused")
        changed = gb.sketch_one(str(cache), str(out), {**good, "isolate_id": "I1"}, None, {**params, "scaled": 2})
        self.assertEqual(changed["action"], "sketched")
        wrong = gb.sketch_one(str(cache), str(out), {**good, "genome_sha256": "0" * 64}, None, params)
        self.assertEqual(wrong["status"], "genome_hash_mismatch")
        (cache / "GCA_000000001.1/genome.fna").write_bytes(b">a\n\xff\n")
        bad = gb.sketch_one(str(cache), str(out / ".."), self.row(1, hashlib.sha256(b">a\n\xff\n").hexdigest()),
                            None, params)
        self.assertEqual(bad["status"], "genome_unreadable")

    def test_config_validation(self):
        path = self.root / "c.json"
        for change, message in ((lambda c: c.pop("model"), "lacks model"),
                                (lambda c: c["model"].update(penalty="l2"), "only exact"),
                                (lambda c: c["model"].update(C_grid=[]), "C_grid"),
                                (lambda c: c["model"].update(threshold=1.0), "threshold"),
                                (lambda c: c["seeds"].update(split=None), "seed split"),
                                (lambda c: c["splits"].update(folds=1), "two folds")):
            config = json.loads(json.dumps(CONFIG))
            change(config)
            path.write_text(json.dumps(config))
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                gb.load_config(path)
        path.write_text(json.dumps(CONFIG))
        self.assertEqual(gb.load_config(path)["benchmark_version"], "generalisation_v1")
        self.assertEqual(gb.config_sha(path), pop.sha256(path))

    def test_check_sums(self):
        d = self.root / "build"
        d.mkdir()
        (d / "a.tsv").write_text("x\n")
        (d / "SHA256SUMS").write_text(f"{pop.sha256(d / 'a.tsv')}  a.tsv\n")
        gb.check_sums(d, ["a.tsv"])
        with self.assertRaisesRegex(ValueError, "not listed"):
            gb.check_sums(d, ["b.tsv"])
        (d / "a.tsv").write_text("y\n")
        with self.assertRaisesRegex(ValueError, "differs"):
            gb.check_sums(d, ["a.tsv"])


# ------------------------------------------------------------ end-to-end synthetic study

POP_FIELDS = ("isolate_id assembly_accession lineage_id ncbi_snp_cluster country collection_year source_category "
              "project eligible_temporal_holdout").split()
GROUP_FIELDS = "isolate_id lineage_component lineage_cv_test_eligible split_unit period geography_group".split()
EP_FIELDS = ("isolate_id assembly_accession antimicrobial phenotype phenotype_pr5 phenotype_correction "
             "eligible_binary_rs population_eligible endpoint_eligible determinant_status concordance_category").split()


def write_rows(path, fields, rows):
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


class EndToEndTests(unittest.TestCase):
    """A 150-isolate synthetic study run through every command."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        root = cls.root = Path(cls.tmp.name)
        rng = random.Random(5)
        base = "".join(rng.choice("ACGT") for _ in range(1500))
        gene = "".join(rng.choice("ACGT") for _ in range(400))
        markers = ["".join(rng.choice("ACGT") for _ in range(150)) for _ in range(8)]
        cache, pop_dir, base_dir = root / "cache", root / "population", root / "baseline"
        for d in (cache, pop_dir, base_dir):
            d.mkdir()
        meta, groups, eps, manifest = [], [], [], []
        countries = ("USA", "Canada", "Netherlands")
        for i in range(150):
            iid, gca = f"PDT{i:09d}.1", f"GCA_{i:09d}.1"
            st = f"ST{i % 8}" if i != 149 else "UNKNOWN"
            resistant = (i * 7) % 10 < 4
            noisy = i % 23 == 0
            seq = base + markers[i % 8] + (gene if resistant != noisy else "")
            fasta = f">{gca}\n{seq}\n".encode()
            (cache / gca).mkdir()
            (cache / gca / "genome.fna").write_bytes(fasta)
            # clusters stay within an ST except PDSX, which joins ST5 and ST6 into one component
            cluster = "PDSX" if i in (5, 6) else (f"PDS{i % 8}{(i // 8) % 4}" if i % 2 else NA)
            year = 2014 + i % 10
            meta.append({"isolate_id": iid, "assembly_accession": gca, "lineage_id": st,
                         "ncbi_snp_cluster": cluster if i != 149 else NA, "country": countries[i % 3],
                         "collection_year": str(year), "source_category": ("human", "animal")[i % 2],
                         "project": f"PRJ{i % 4}", "eligible_temporal_holdout": "true"})
            groups.append({"isolate_id": iid, "lineage_component": st if i != 149 else NA,
                           "lineage_cv_test_eligible": "false" if i == 149 else "true",
                           "split_unit": cluster if cluster != NA and i != 149 else iid,
                           "period": "recent" if year >= 2021 else "earlier",
                           "geography_group": "reference" if i % 3 == 0 else "other"})
            phen = "R" if resistant else "S"
            det = "positive" if resistant and i % 5 else "negative"
            for drug in ("ampicillin", "tetracycline"):
                p = phen if drug == "ampicillin" else ("R" if i % 5 < 2 else "S")
                eps.append({"isolate_id": iid, "assembly_accession": gca, "antimicrobial": drug, "phenotype": p,
                            "phenotype_pr5": p, "phenotype_correction": NA, "eligible_binary_rs": "true",
                            "population_eligible": "true", "endpoint_eligible": "true" if i % 11 else "false",
                            "determinant_status": det,
                            "concordance_category": f"{p}_DET_{'POS' if det == 'positive' else 'NEG'}"})
            manifest.append({"isolate_id": iid, "assembly_accession": gca, "biosample_accession": f"SAMN{i}",
                             "genome_sha256": hashlib.sha256(fasta).hexdigest(), "sequence_records": "1",
                             "sequence_bases": str(len(seq))})
        eps.append({**eps[0], "antimicrobial": "cefotaxime"})
        eps.append({**eps[0], "antimicrobial": "colistin", "eligible_binary_rs": "false", "phenotype": "ND"})
        comp = pop.lineage_components([{**m, "eligible_population_analysis": "true"} for m in meta])
        for g in groups:
            g["lineage_component"] = comp[g["isolate_id"]]
            g["lineage_cv_test_eligible"] = "false" if comp[g["isolate_id"]] == NA else "true"
        write_rows(pop_dir / "population_metadata.tsv", POP_FIELDS, meta)
        write_rows(pop_dir / "evaluation_groups.tsv", GROUP_FIELDS, groups)
        (pop_dir / "population_manifest.json").write_text(json.dumps({"pr5_cohort_manifest_sha256": "a" * 64}))
        write_rows(base_dir / "genotype_phenotype_endpoints.tsv", EP_FIELDS, eps)
        (base_dir / "baseline_manifest.json").write_text(json.dumps({
            "pr5_cohort_manifest_sha256": "a" * 64, "phenotype_corrections_sha256": "c" * 64,
            "mapping_version": "ecoli_amrfinder_drug_mapping_v1"}))
        for d, names in ((pop_dir, ("population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json")),
                         (base_dir, ("genotype_phenotype_endpoints.tsv", "baseline_manifest.json"))):
            (d / "SHA256SUMS").write_text("".join(f"{pop.sha256(d / n)}  {n}\n" for n in names))
        manifest_path = root / "genome_manifest.tsv"
        write_rows(manifest_path, ["isolate_id", "assembly_accession", "biosample_accession", "genome_sha256",
                                   "sequence_records", "sequence_bases"], manifest)
        popc = json.loads(json.dumps(POP_CONFIG))
        popc["evaluation_design"]["drug_tiers"].update({
            "primary": {"min_R": 40, "min_S": 40, "max_top_project_R_fraction": 0.6, "min_projects_with_both": 2,
                        "min_sts_with_both": 2},
            "secondary": {"min_R": 10, "min_S": 10}, "project_both_min": 5, "st_both_min": 2})
        popc["evaluation_design"]["temporal_min_per_period"] = 5
        popc["evaluation_design"]["geography"]["min_pooled"] = 10
        popc_path = root / "population_config.json"
        popc_path.write_text(json.dumps(popc))
        config = json.loads(json.dumps(CONFIG))
        config["upstream"].update({"pr5_cohort_manifest_sha256": "a" * 64,
                                   "pr6_population_manifest_sha256": pop.sha256(pop_dir / "population_manifest.json"),
                                   "pr7_baseline_manifest_sha256": pop.sha256(base_dir / "baseline_manifest.json"),
                                   "pr7_phenotype_corrections_sha256": "c" * 64,
                                   "genome_manifest": str(manifest_path),
                                   "genome_manifest_sha256": pop.sha256(manifest_path), "genome_manifest_rows": 150})
        config["endpoints"].update({"population_config": str(popc_path),
                                    "population_config_sha256": pop.sha256(popc_path)})
        config["features"].update({"scaled": 1, "ksize": 21})
        config["model"]["C_grid"] = [0.3, 1.0]
        config["bootstrap"]["replicates"] = 20
        for s in ("country", "project", "geography", "temporal"):
            config["splits"]["strategies"][s].update({"min_test_R": 5, "min_test_S": 5, "min_train_R": 5,
                                                      "min_train_S": 5})
        cls.config_path = root / "benchmark.json"
        cls.config_path.write_text(json.dumps(config))
        cls.config = config
        cls.out = root / "out"
        cls.args = ["--config", str(cls.config_path), "--population", str(pop_dir), "--baseline", str(base_dir),
                    "--outdir", str(cls.out)]
        cls.cache = cache
        with patch("sys.stdout"):
            cls.codes = [gb.main(["sketch", *cls.args, "--cache", str(cache), "--workers", "2"]),
                         gb.main(["matrix", *cls.args]), gb.main(["freeze", *cls.args]),
                         gb.main(["run", *cls.args, "--workers", "2", "--only", "common__ampicillin__random__"]),
                         gb.main(["run", *cls.args]), gb.main(["report", *cls.args])]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def read(self, name):
        path = self.out / name
        if name.endswith(".gz"):
            with gzip.open(path, "rt") as h:
                return list(csv.DictReader(h, delimiter="\t"))
        return pop.read_tsv(path)

    def test_commands_succeed(self):
        self.assertEqual(self.codes, [0, 0, 0, 0, 0, 0])
        summary = json.loads((self.out / "benchmark_summary.json").read_text())
        self.assertEqual(summary["tasks"], summary["tasks_succeeded"])
        self.assertEqual(summary["tasks_not_succeeded"], {})

    def test_genome_verification_and_matrix(self):
        ver = self.read("genome_verification.tsv")
        self.assertEqual({r["status"] for r in ver}, {"verified"})
        fm = json.loads((self.out / "matrix/feature_manifest.json").read_text())
        self.assertEqual(fm["rows"], 150)
        self.assertGreater(fm["kept_features"], fm["patterns"])

    def test_manifest_records_upstream_and_design(self):
        m = json.loads((self.out / "benchmark_manifest.json").read_text())
        self.assertEqual(m["config_sha256"], pop.sha256(self.config_path))
        self.assertEqual(m["upstream"]["pr7_baseline_manifest_sha256"], self.config["upstream"]["pr7_baseline_manifest_sha256"])
        self.assertEqual(m["antimicrobials"], {"ampicillin": "primary", "tetracycline": "primary"})
        self.assertEqual(m["genome_verification"], {"verified": 150})
        self.assertEqual(m["seeds"], self.config["seeds"])
        plan_rows = {r["antimicrobial"]: r for r in self.read("drug_plan_corrected.tsv")}
        self.assertEqual(plan_rows["cefotaxime"]["tier"], "excluded")

    def test_split_outputs(self):
        summary = self.read("split_summary.tsv")
        self.assertEqual(len(summary), json.loads((self.out / "benchmark_manifest.json").read_text())["tasks"])
        rows = self.read("splits/lineage.tsv.gz")
        self.assertEqual({r["partition"] for r in rows}, {"train", "test"})
        unknown = [r for r in rows if r["isolate_id"] == "PDT000000149.1"]
        self.assertTrue(unknown and all(r["partition"] == "train" for r in unknown if r["cohort"] == "maximal"))
        self.assertFalse(any(r["cohort"] == "common" for r in unknown))
        elig = {r["antimicrobial"]: r for r in self.read("cohort_eligibility.tsv")}
        self.assertEqual(elig["ampicillin"]["common_cohort_N"], "149")

    def test_predictions_and_metrics(self):
        preds = self.read("predictions.tsv.gz")
        self.assertTrue(preds)
        common_random = [p for p in preds if p["cohort"] == "common" and p["split_strategy"] == "random"
                         and p["antimicrobial"] == "ampicillin"]
        self.assertEqual(len(common_random), 149)
        self.assertEqual({p["correct"] for p in preds} <= {"yes", "no"}, True)
        metrics = {(r["cohort"], r["antimicrobial"], r["split_strategy"]): r for r in self.read("benchmark_metrics.tsv")}
        rnd = metrics[("common", "ampicillin", "random")]
        self.assertGreater(float(rnd["AUROC"]), 0.8)
        self.assertNotEqual(rnd["AUROC_ci_lower"], NA)
        self.assertEqual(rnd["complete"], "yes")
        folds = self.read("fold_metrics.tsv")
        self.assertEqual(len([f for f in folds if f["cohort"] == "common" and f["antimicrobial"] == "ampicillin"
                              and f["split_strategy"] == "lineage"]), 5)

    def test_drop_tables(self):
        drops = self.read("generalisation_drop.tsv")
        pooled = [d for d in drops if d["level"] == "pooled" and d["comparator"] == "random"]
        self.assertEqual({d["holdout_strategy"] for d in pooled},
                         {"lineage", "country", "project", "geography", "temporal"})
        for d in pooled:
            if d["delta_AUROC"] != NA:
                self.assertAlmostEqual(float(d["delta_AUROC"]),
                                       float(d["AUROC_comparator"]) - float(d["AUROC_holdout"]), places=5)
        opt = self.read("random_split_optimism.tsv")
        self.assertEqual({o["metric"] for o in opt}, {"AUROC", "AUPRC", "balanced_accuracy", "MCC"})
        table = self.read("manuscript_table.tsv")
        self.assertTrue(all(r["delta_AUROC_vs_random_same_isolates"] == NA for r in table
                            if r["split_strategy"] == "random"))

    def test_audits_and_comparisons(self):
        audit = self.read("split_leakage_audit.tsv")
        lineage = [a for a in audit if a["strategy"] == "lineage"]
        self.assertTrue(all(a["test_st_in_train"] == "0" for a in lineage))
        self.assertTrue(self.read("country_shift_summary.tsv"))
        self.assertTrue(self.read("temporal_shift_summary.tsv"))
        det = self.read("determinant_vs_ml.tsv")
        self.assertTrue(any(int(r["n_R_DET_NEG"]) > 0 for r in det))
        un = self.read("unexplained_resistance_ml.tsv")
        self.assertEqual({u["lineage_seen_in_training"] for u in un}, {"all", "yes", "no"})
        seen = self.read("unseen_lineage_summary.tsv")
        self.assertTrue(any(s["value"] == "no" for s in seen if s["split_strategy"] == "lineage"))
        for name in ("fig_a_performance_by_strategy", "fig_b_drop_vs_random", "fig_c_lineage_folds",
                     "fig_d_country_holdouts", "fig_e_temporal", "fig_f_seen_vs_unseen_lineage"):
            self.assertTrue((self.out / "figure_data" / f"{name}.tsv").exists())

    def test_rerun_reuses_checkpoints_and_detects_changes(self):
        with patch("sys.stdout"):
            self.assertEqual(gb.main(["run", *self.args]), 0)
        inputs = gb.load_inputs(self.config_path, Path(self.args[3]), Path(self.args[5]),
                                self.out / "genome_verification.tsv")
        tasks, _ = gb.build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
        counts = gb.run_tasks(tasks, self.out / "matrix", self.out, inputs["config"])
        self.assertEqual(counts, Counter({"reused": len(tasks)}))
        manifest = json.loads((self.out / "benchmark_manifest.json").read_text())
        gb.check_frozen(tasks, manifest, self.config_path)
        tampered = {**manifest, "task_split_sha256": {**manifest["task_split_sha256"], tasks[0]["task_id"]: "x"}}
        with self.assertRaises(gb.LeakageError):
            gb.check_frozen(tasks, tampered, self.config_path)
        with self.assertRaisesRegex(ValueError, "config changed"):
            gb.check_frozen(tasks, {**manifest, "config_sha256": "0"}, self.config_path)

    def test_freeze_is_deterministic(self):
        first = hashlib.sha256((self.out / "splits/random.tsv.gz").read_bytes()).hexdigest()
        other = self.root / "again"
        other.mkdir()
        (other / "matrix").mkdir()
        for name in ("genome_verification.tsv",):
            (other / name).write_bytes((self.out / name).read_bytes())
        (other / "matrix/feature_manifest.json").write_bytes((self.out / "matrix/feature_manifest.json").read_bytes())
        args = [a if a != str(self.out) else str(other) for a in self.args]
        with patch("sys.stdout"):
            gb.main(["freeze", *args])
        with gzip.open(self.out / "splits/random.tsv.gz", "rt") as a, gzip.open(other / "splits/random.tsv.gz", "rt") as b:
            self.assertEqual(a.read(), b.read())
        m1 = json.loads((self.out / "benchmark_manifest.json").read_text())
        m2 = json.loads((other / "benchmark_manifest.json").read_text())
        self.assertEqual(m1["task_split_sha256"], m2["task_split_sha256"])
        self.assertTrue(first)

    def test_failed_and_missing_tasks_are_reported(self):
        inputs = gb.load_inputs(self.config_path, Path(self.args[3]), Path(self.args[5]),
                                self.out / "genome_verification.tsv")
        tasks, skipped = gb.build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
        out = self.root / "partial"
        (out / "tasks").mkdir(parents=True)
        for t in tasks[:3]:
            for ext in ("json", "tsv"):
                (out / "tasks" / f"{t['task_id']}.{ext}").write_bytes((self.out / "tasks" / f"{t['task_id']}.{ext}").read_bytes())
        stale = json.loads((out / "tasks" / f"{tasks[2]['task_id']}.json").read_text())
        (out / "tasks" / f"{tasks[2]['task_id']}.json").write_text(json.dumps({**stale, "split_sha256": "x"}))
        manifest = json.loads((self.out / "benchmark_manifest.json").read_text())
        summary = gb.report(inputs, tasks, skipped, out, manifest)
        self.assertEqual(summary["tasks_succeeded"], 2)
        self.assertEqual(summary["tasks_not_succeeded"][tasks[2]["task_id"]], "split changed since the task was run")
        metrics = self.read_from(out, "benchmark_metrics.tsv")
        self.assertEqual({m["complete"] for m in metrics}, {"no"})

    def read_from(self, out, name):
        return pop.read_tsv(out / name)

    def test_worker_records_failures(self):
        inputs = gb.load_inputs(self.config_path, Path(self.args[3]), Path(self.args[5]),
                                self.out / "genome_verification.tsv")
        tasks, _ = gb.build_tasks(inputs["endpoints"], inputs["plan"], inputs["config"])
        t = tasks[0]
        single = {**t, "train": [{**e, "phenotype": "S"} for e in t["train"]]}
        status, preds = gb._worker((single, str(self.out / "matrix"), inputs["config"]))
        self.assertEqual(status["status"], "failed")
        self.assertIn("single_class_training", status["reason"])
        broken = {**t, "test": [{**t["test"][0], "isolate_id": "missing"}]}
        status, _ = gb._worker((broken, str(self.out / "matrix"), inputs["config"]))
        self.assertIn("KeyError", status["reason"])
        out = self.root / "fail"
        with patch("sys.stdout"):
            counts = gb.run_tasks([single], self.out / "matrix", out, inputs["config"])
        self.assertEqual(counts, Counter({"failed": 1}))

    def test_upstream_mismatches_fail(self):
        config = json.loads(json.dumps(self.config))
        pop_dir, base_dir = Path(self.args[3]), Path(self.args[5])
        for key, message in (("pr6_population_manifest_sha256", "PR #6"), ("pr7_baseline_manifest_sha256", "PR #7 baseline"),
                             ("pr5_cohort_manifest_sha256", "PR #5"), ("pr7_phenotype_corrections_sha256", "correction"),
                             ("pr7_mapping_version", "mapping")):
            changed = json.loads(json.dumps(config))
            changed["upstream"][key] = "0" * 64
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, message):
                gb.upstream_hashes(changed, pop_dir, base_dir)
        path = self.root / "changed.json"
        changed = json.loads(json.dumps(config))
        changed["endpoints"]["population_config_sha256"] = "0" * 64
        path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError, "population config differs"):
            gb.load_inputs(path, pop_dir, base_dir, self.out / "genome_verification.tsv")

    def test_sketch_failure_exit_code_and_cli_errors(self):
        out = self.root / "sk2"
        cache = self.root / "cache2"
        cache.mkdir()
        with patch("sys.stdout"):
            code = gb.main(["sketch", *[a if a != str(self.out) else str(out) for a in self.args],
                            "--cache", str(cache), "--limit", "2"])
        self.assertEqual(code, 1)
        self.assertEqual(json.loads((out / "genome_verification.json").read_text())["status_counts"],
                         {"genome_missing": 2})
        with self.assertRaises(SystemExit), patch("sys.stderr"):
            gb.main(["sketch", *self.args])

    def test_matrix_from_other_verification_is_refused(self):
        other = self.root / "mismatch"
        (other / "matrix").mkdir(parents=True)
        (other / "genome_verification.tsv").write_bytes((self.out / "genome_verification.tsv").read_bytes() + b"\n")
        fm = json.loads((self.out / "matrix/feature_manifest.json").read_text())
        (other / "matrix/feature_manifest.json").write_text(json.dumps(fm))
        with self.assertRaisesRegex(ValueError, "different genome verification"):
            gb.main(["freeze", *[a if a != str(self.out) else str(other) for a in self.args]])

    def test_load_endpoints_rejects_inconsistent_inputs(self):
        pop_dir, base_dir = Path(self.args[3]), Path(self.args[5])
        bad = self.root / "bad_base"
        bad.mkdir()
        rows = pop.read_tsv(base_dir / "genotype_phenotype_endpoints.tsv")
        write_rows(bad / "genotype_phenotype_endpoints.tsv", EP_FIELDS, [{**rows[0], "phenotype": "I"}])
        with self.assertRaisesRegex(ValueError, "binary endpoint with phenotype"):
            gb.load_endpoints(bad, pop_dir, set(), POP_CONFIG)
        write_rows(bad / "genotype_phenotype_endpoints.tsv", EP_FIELDS, [{**rows[0], "assembly_accession": "GCA_X"}])
        with self.assertRaisesRegex(ValueError, "assemblies differ"):
            gb.load_endpoints(bad, pop_dir, set(), POP_CONFIG)


if __name__ == "__main__":
    unittest.main()
