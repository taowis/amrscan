"""Offline tests for the representation sensitivity analysis (synthetic sketches, genomes, tasks and predictions)."""

from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
spec = importlib.util.spec_from_file_location("representation_sensitivity", ROOT / "bin/representation_sensitivity.py")
rs = importlib.util.module_from_spec(spec)
sys.modules["representation_sensitivity"] = rs
spec.loader.exec_module(rs)
ed, gb = rs.ed, rs.gb
NA = rs.NA
CONFIG_PATH = ROOT / "benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1.json"
CONFIG = json.loads(CONFIG_PATH.read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rnd_dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def write_sketches(d, n=12, seed=3, universe=4000, size=900):
    """Sorted uint64 sketches drawn from a shared universe below the scaled-2 threshold."""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    top = gb.max_hash(2)
    pool = np.unique(rng.integers(0, top, universe, dtype=np.uint64))
    rows = []
    for i in range(n):
        h = np.sort(rng.choice(pool, size, replace=False))
        gca = f"GCA_{i:09d}.1"
        np.save(d / f"{gca}.npy", h)
        rows.append((f"I{i:03d}", gca))
    return rows


class ConfigTests(unittest.TestCase):
    def test_committed_config_loads(self):
        c = rs.load_config(CONFIG_PATH)
        self.assertEqual(c["density"]["sketch_scaled"], 20)
        self.assertEqual(c["density"]["reference_scaled"], 200)
        self.assertEqual(sorted(c["representations"]), sorted(rs.REPRESENTATIONS))

    def test_config_rules(self):
        cases = [(lambda c: c.pop("density"), "lacks density"),
                 (lambda c: c["density"].update(sketch_scaled=50), "densest candidate"),
                 (lambda c: c["density"].update(candidates=[30, 20]), "divide"),
                 (lambda c: c["representations"].pop("R3"), "representation set")]
        for change, msg in cases:
            c = json.loads(json.dumps(CONFIG))
            change(c)
            with tempfile.TemporaryDirectory() as t, self.subTest(msg=msg):
                path = Path(t) / "c.json"
                path.write_text(json.dumps(c))
                with self.assertRaisesRegex(ValueError, msg):
                    rs.load_config(path)

    def test_same_feature_family_across_representations(self):
        g = json.loads((ROOT / CONFIG["upstream"]["generalisation_config"]).read_text())["features"]
        d = CONFIG["density"]
        self.assertEqual((d["ksize"], d["hash_seed"], d["tool"], d["tool_version"], d["reference_scaled"]),
                         (g["ksize"], g["hash_seed"], g["tool"], g["tool_version"], g["scaled"]))
        p200, p20 = rs.sketch_params(CONFIG, 200), rs.sketch_params(CONFIG, 20)
        self.assertEqual({k: v for k, v in p200.items() if k != "scaled"}, {k: v for k, v in p20.items()
                                                                              if k != "scaled"})
        self.assertEqual(p200, gb.sketch_params(json.loads((ROOT / CONFIG["upstream"]["generalisation_config"])
                                                           .read_text())))

    def test_rep_order(self):
        names = ["scaled_20", "R3", "R0", "scaled_100", "other"]
        self.assertEqual(sorted(names, key=rs.rep_order), ["R0", "R3", "scaled_100", "scaled_20", "other"])


class SketchTests(unittest.TestCase):
    def test_downsample_is_a_prefix(self):
        h = np.sort(np.random.default_rng(1).integers(0, gb.max_hash(2), 500, dtype=np.uint64))
        d4 = rs.downsample(h, 4)
        self.assertTrue(np.array_equal(d4, h[h <= np.uint64(gb.max_hash(4))]))
        self.assertTrue(np.array_equal(rs.downsample(d4, 8), rs.downsample(h, 8)))

    def test_sketch_reproduces_reference(self):
        params = {**rs.sketch_params(CONFIG, 1), "ksize": 21}
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            cache, out, r0 = t / "cache", t / "out", t / "r0"
            for d in (out, r0):
                d.mkdir()
            data = (">c1\n" + rnd_dna(3000, 5) + "\n").encode()
            (cache / "GCA_000000001.1").mkdir(parents=True)
            (cache / "GCA_000000001.1" / "genome.fna").write_bytes(data)
            row = {"isolate_id": "I1", "assembly_accession": "GCA_000000001.1",
                   "genome_sha256": hashlib.sha256(data).hexdigest()}
            full = gb.sketch_genome(data, params)
            np.save(r0 / "GCA_000000001.1.npy", rs.downsample(full, 3))
            res = rs.sketch_one(str(cache), str(out), str(r0), row, None, params, 3)
            self.assertEqual((res["status"], res["r0_reproduced"], res["detail"]), ("genome_verified", "yes",
                                                                                   "sketched"))
            again = rs.sketch_one(str(cache), str(out), str(r0), row, None, params, 3)
            self.assertEqual(again["detail"], "reused")
            np.save(r0 / "GCA_000000001.1.npy", rs.downsample(full, 3)[1:])
            self.assertEqual(rs.sketch_one(str(cache), str(out), str(r0), row, None, params, 3)["r0_reproduced"],
                             "no")
            # deterministic: a second sketch from scratch gives the same hashes
            self.assertTrue(np.array_equal(gb.sketch_genome(data, params), full))

    def test_genome_checks(self):
        params = {**rs.sketch_params(CONFIG, 1), "ksize": 21}
        with tempfile.TemporaryDirectory() as t:
            t = Path(t)
            (t / "out").mkdir()
            (t / "cache" / "GCA_000000002.1").mkdir(parents=True)
            (t / "cache" / "GCA_000000002.1" / "genome.fna").write_bytes(b">c\nACGT\n")
            row = {"isolate_id": "I2", "assembly_accession": "GCA_000000002.1", "genome_sha256": "0" * 64}
            args = (str(t / "cache"), str(t / "out"), str(t))
            self.assertEqual(rs.sketch_one(*args, row, None, params, 1)["status"], "genome_hash_mismatch")
            self.assertEqual(rs.sketch_one(*args, row, "GCA_000000002.2", params, 1)["status"],
                             "genome_version_mismatch")
            missing = {**row, "assembly_accession": "GCA_000000003.1"}
            self.assertEqual(rs.sketch_one(*args, missing, None, params, 1)["status"], "genome_missing")
            bad = b">c\n\xff\xfe\n"
            (t / "cache" / "GCA_000000002.1" / "genome.fna").write_bytes(bad)
            unreadable = {**row, "genome_sha256": hashlib.sha256(bad).hexdigest()}
            self.assertEqual(rs.sketch_one(*args, unreadable, None, params, 1)["status"], "genome_unreadable")

    def test_verified_rows_stop_on_failure(self):
        with tempfile.TemporaryDirectory() as t:
            fields = rs.VERIFY_FIELDS
            rows = [{"isolate_id": "I1", "assembly_accession": "G1", "status": "genome_verified",
                     "r0_reproduced": "yes"}]
            gb.write_tsv(Path(t) / "genome_verification.tsv", fields, rows)
            self.assertEqual(rs.verified_rows(Path(t)), [("I1", "G1")])
            gb.write_tsv(Path(t) / "genome_verification.tsv", fields, rows + [{**rows[0], "isolate_id": "I2",
                                                                               "r0_reproduced": "no"}])
            with self.assertRaises(rs.UpstreamError):
                rs.verified_rows(Path(t))


class FeatureSpaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name)
        cls.rows = write_sketches(cls.dir)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_one_pass_matches_build_matrix_at_every_density(self):
        spaces = rs.candidate_feature_spaces(self.dir, self.rows, [8, 4, 2], 2, 2, buckets=5)
        for s in (8, 4, 2):
            with self.subTest(scaled=s):
                b = gb.build_matrix(self.dir, self.rows, 2, 2, True, scaled=s)
                m = spaces[s]
                self.assertTrue(np.array_equal(b["hashes"], m["hashes"]))
                self.assertTrue(np.array_equal(b["hash_pattern"], m["hash_pattern"]))
                self.assertEqual((b["observed_features"], b["patterns"], b["nnz"]),
                                 (m["observed_features"], m["patterns"], m["nnz"]))
        self.assertLess(spaces[8]["kept_features"], spaces[4]["kept_features"])
        self.assertLess(spaces[4]["patterns"], spaces[2]["patterns"])

    def test_bucket_count_does_not_change_the_result(self):
        a = rs.candidate_feature_spaces(self.dir, self.rows, [4, 2], 2, 2, buckets=3)
        b = rs.candidate_feature_spaces(self.dir, self.rows, [4, 2], 2, 2, buckets=11)
        for s in (4, 2):
            self.assertTrue(np.array_equal(a[s]["hash_pattern"], b[s]["hash_pattern"]))
            self.assertEqual(a[s]["nnz"], b[s]["nnz"])

    def test_streaming_matrix_matches_build_matrix(self):
        for s in (4, 2):
            with self.subTest(scaled=s):
                b = gb.build_matrix(self.dir, self.rows, 2, 2, True, scaled=s)
                m = rs.build_matrix_streaming(self.dir, self.rows, s, 2, 2, self.dir / "spill.bin", buckets=7)
                self.assertFalse((self.dir / "spill.bin").exists())
                for k in ("indptr", "indices"):
                    self.assertTrue(np.array_equal(getattr(b["matrix"], k), getattr(m["matrix"], k)))
                self.assertTrue(np.array_equal(b["hashes"], m["hashes"]))
                self.assertEqual(b["matrix"].shape, m["matrix"].shape)
                # sparse rows sorted, presence only
                self.assertTrue(m["matrix"].has_sorted_indices)
                self.assertEqual(set(m["matrix"].data.tolist()), {1})

    def test_no_informative_feature(self):
        with self.assertRaisesRegex(ValueError, "no informative"):
            rs.build_matrix_streaming(self.dir, self.rows[:3], 2, 2, 2, self.dir / "s.bin", buckets=2)

    def test_vocabulary_hash_and_dimensions_saved(self):
        with tempfile.TemporaryDirectory() as t:
            built = rs.build_matrix_streaming(self.dir, self.rows, 2, 2, 2, Path(t) / "s.bin", buckets=3)
            m = gb.save_matrix(Path(t), built, self.rows, {"params": {"scaled": 2}})
            X, row_of, meta = gb.load_matrix(Path(t), check=True)
            self.assertEqual(X.shape, (len(self.rows), m["patterns"]))
            self.assertEqual(meta["feature_hashes_sha256"], gb.sha256_bytes(built["hashes"].tobytes()))
            self.assertEqual(sorted(row_of), [r[0] for r in self.rows])


class CaptureTests(unittest.TestCase):
    """A point mutation whose k-mers a sparse space never samples but a dense one keeps."""

    def setUp(self):
        self.k = 31
        self.locus = rnd_dna(200, 11)
        self.mh = ed.FeatureSpace(np.array([], dtype=np.uint64), np.array([], dtype=np.int64), 1, 31, 42)
        self.site = 100
        self.mut = self.locus[:self.site] + ("A" if self.locus[self.site] != "A" else "C") + self.locus[self.site + 1:]

    def space(self, kept, sample_all):
        hashes = np.array(sorted(kept), dtype=np.uint64)
        s = ed.FeatureSpace(hashes, np.arange(hashes.shape[0]), 200, 31, 42)
        s.max_hash = 2 ** 64 if sample_all else 0
        return s

    def mutant_hashes(self, seq):
        return [self.mh.kmer_hash(k) for k in ed.overlapping_kmers(seq, self.site, self.site + 1, self.k)]

    def test_missed_by_sparse_recovered_by_dense(self):
        mutant = self.mutant_hashes(self.mut)
        sparse_space = self.space([], sample_all=False)
        dense_space = self.space(mutant[:1], sample_all=True)
        k = ed.overlapping_kmers(self.mut, self.site, self.site + 1, self.k)
        self.assertEqual(sparse_space.classify(k)["kept"], 0)
        self.assertEqual(sparse_space.classify(k)["sampled"], 0)
        self.assertEqual(dense_space.classify(k)["kept"], 1)
        self.assertEqual(len(k), 31)

    def test_gene_captured_by_both_and_absent(self):
        contigs = {"c": rnd_dna(500, 2)}
        call = {"contig": "c", "start": "50", "stop": "250"}
        region = contigs["c"][49:250]
        h = self.mh.kmer_hash(region[10:41])
        dense, coarse = self.space([h], True), self.space([h], True)
        self.assertEqual(ed.gene_capture(contigs, call, dense)["status"], "captured")
        self.assertEqual(ed.gene_capture(contigs, call, coarse)["status"], "captured")
        self.assertEqual(ed.gene_capture(contigs, call, self.space([], False))["status"], "not_represented")
        self.assertEqual(ed.gene_capture({}, call, dense)["status"], "representation_unknown")

    def test_capture_tables_counts_and_fractions(self):
        audit = {
            "I1": {"isolate_id": "I1", "genes": [{"determinant": "tetA", "status": {200: "captured", 20: "captured"}}],
                   "mutations": [{"mutation": "gyrA_S83L", "method": "POINTX",
                                  "status": {200: "not_represented", 20: "captured"}},
                                 {"mutation": "parC_S80I", "method": "POINTX",
                                  "status": {200: "not_represented", 20: "partially_represented"}}]},
            "I2": {"isolate_id": "I2", "genes": [],
                   "mutations": [{"mutation": "gyrA_S83L", "method": "POINTX",
                                  "status": {200: "partially_represented", 20: "captured"}}]}}
        endpoints = [{"antimicrobial": "ciprofloxacin", "isolate_id": "I1", "gene_determinants": NA,
                      "catalogue_mutations": "gyrA_S83L;parC_S80I"},
                     {"antimicrobial": "ciprofloxacin", "isolate_id": "I2", "gene_determinants": NA,
                      "catalogue_mutations": "gyrA_S83L"},
                     {"antimicrobial": "tetracycline", "isolate_id": "I1", "gene_determinants": "tetA",
                      "catalogue_mutations": NA},
                     {"antimicrobial": "tetracycline", "isolate_id": "I3", "gene_determinants": "tetB",
                      "catalogue_mutations": NA}]
        label = lambda s: "R0" if s == 200 else f"scaled_{s}"  # noqa: E731
        t = rs.capture_tables(audit, endpoints, {}, {"ciprofloxacin": "high_loss"}, label)
        mut = {(r["representation"], r["mutation"]): r for r in t["mutations"]}
        self.assertEqual((mut[("R0", "gyrA_S83L")]["calls"], mut[("R0", "gyrA_S83L")]["captured"]), (2, 0))
        self.assertEqual(mut[("scaled_20", "gyrA_S83L")]["captured_fraction"], 1.0)
        self.assertEqual(mut[("R0", "gyrA_S83L")]["carriers"], 2)
        drugs = {(r["antimicrobial"], r["determinant"], r["representation"]): r for r in t["drugs"]}
        cip = drugs[("ciprofloxacin", "ALL", "scaled_20")]
        self.assertEqual((cip["endpoints"], cip["captured"], cip["partially_represented"]), (3, 2, 1))
        self.assertEqual(drugs[("ciprofloxacin", "ALL", "R0")]["captured"], 0)
        self.assertEqual(drugs[("ciprofloxacin", "ALL", "R0")]["drug_group"], "high_loss")
        tet = drugs[("tetracycline", "ALL", "R0")]
        self.assertEqual((tet["captured"], tet["representation_unknown"]), (1, 1))
        # deterministic
        self.assertEqual(rs.capture_tables(audit, endpoints, {}, {}, label)["drugs"][0], t["drugs"][0] | {
            "drug_group": NA})

    def test_several_copies_take_the_best_status_whatever_their_order(self):
        label = lambda s: "R0" if s == 200 else f"scaled_{s}"  # noqa: E731
        copies = [{"determinant": "dfrA12", "status": {200: "not_represented"}},
                  {"determinant": "dfrA12", "status": {200: "captured"}},
                  {"determinant": "dfrA12", "status": {200: "not_represented"}}]
        endpoints = [{"antimicrobial": "trimethoprim", "isolate_id": "I1", "gene_determinants": "dfrA12",
                      "catalogue_mutations": NA}]
        for genes in (copies, copies[::-1]):
            t = rs.capture_tables({"I1": {"isolate_id": "I1", "genes": genes, "mutations": []}}, endpoints, {}, {},
                                  label)
            self.assertEqual(t["gcap"][("I1", "dfrA12", "R0")], "captured")
            row = next(r for r in t["genes"] if r["determinant"] == "dfrA12")
            self.assertEqual((row["calls"], row["captured"], row["not_represented"]), (3, 1, 2))


class SelectionTests(unittest.TestCase):
    def resources(self):
        rows = []
        for sc, nnz in ((200, 1.0), (100, 2.0), (50, 3.5), (20, 7.0)):
            rows.append({"scaled": str(sc), "projected_task_peak_rss_bytes": str(7e9 * nnz),
                         "projected_task_seconds_total": str(56000 * nnz)})
        return rows

    def budget(self, su, node=256, measured=None):
        return {"budget_su": su, "node_memory_gb": node, "su_per_core_hour": 1.25, "memory_gb_per_core": 9.14,
                "cpu_speed_factor": 1.5, "r3_tasks": 262, "measured": measured or {}}

    def test_densest_feasible_candidate(self):
        sel = rs.select_density(CONFIG, self.resources(), self.budget(10 ** 6))
        self.assertEqual(sel["R1_scaled"], 20)
        self.assertEqual({r["scaled"] for r in sel["candidates"]}, {100, 50, 20})

    def test_memory_and_budget_limits(self):
        self.assertEqual(rs.select_density(CONFIG, self.resources(), self.budget(10 ** 6, node=50))["R1_scaled"], 50)
        tight = rs.select_density(CONFIG, self.resources(), self.budget(200))
        self.assertEqual(tight["R1_scaled"], 100)
        with self.assertRaisesRegex(ValueError, "no candidate"):
            rs.select_density(CONFIG, self.resources(), self.budget(1))

    def test_measured_values_replace_projections(self):
        measured = {"20": {"peak_rss_bytes": 20e9, "task_hours_total": 10}}
        sel = rs.select_density(CONFIG, self.resources(), self.budget(150, measured=measured))
        row = next(r for r in sel["candidates"] if r["scaled"] == 20)
        self.assertEqual(row["basis"], "measured")
        self.assertEqual(sel["R1_scaled"], 20)


def task(tid, drug, strat, fold, train, test, cohort="common"):
    e = lambda i, ph: {"isolate_id": i, "phenotype": ph, "split_unit": f"u{i}"}  # noqa: E731
    return {"task_id": tid, "cohort": cohort, "antimicrobial": drug, "strategy": strat, "fold": fold,
            "train": [e(*x) for x in train], "test": [e(*x) for x in test]}


class IdentityAndStatusTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = self.dir = Path(self.tmp.name)
        self.tasks = [task("common__cip__random__fold1", "cip", "random", "fold1", [("a", "R"), ("b", "S")],
                           [("c", "R"), ("d", "S")]),
                      task("common__amp__random__fold1", "amp", "random", "fold1", [("a", "R"), ("b", "S")],
                           [("c", "S"), ("d", "R")])]
        self.manifest = {"task_split_sha256": {x["task_id"]: gb.task_digest(x) for x in self.tasks}}
        (t / "bench" / "tasks").mkdir(parents=True)
        for x in self.tasks:
            (t / "bench" / "tasks" / f"{x['task_id']}.json").write_text(json.dumps(
                {"split_sha256": gb.task_digest(x), "n_train": 2, "n_test": 2, "train_R": 1}))
            gb.write_tsv(t / "bench" / "tasks" / f"{x['task_id']}.tsv", gb.PRED_FIELDS,
                         [{"isolate_id": e["isolate_id"], "phenotype": e["phenotype"]} for e in x["test"]])
        for v in ("mutation_augmented", "permuted_control"):
            d = t / "dec" / "augmented" / v / "tasks"
            d.mkdir(parents=True)
            x = self.tasks[0]
            (d / f"{x['task_id']}.json").write_text(json.dumps({"split_sha256": gb.task_digest(x),
                                                                "status": "success"}))
        self.plan = {"drugs": ["cip"]}

    def tearDown(self):
        self.tmp.cleanup()

    def test_identical_tasks(self):
        rows = rs.task_identity(self.tasks, self.manifest, self.dir / "bench", self.dir / "dec", self.plan)
        self.assertEqual(Counter(r["status"] for r in rows), Counter({"identical": 2}))
        by = {r["task_id"]: r for r in rows}
        self.assertEqual(by["common__amp__random__fold1"]["r2_match"], "not_applicable")
        self.assertEqual(by["common__cip__random__fold1"]["r2_match"], "yes")
        self.assertEqual(rs.r2_task_set(self.tasks, self.plan), [self.tasks[0]])

    def test_changed_test_isolates_fail(self):
        changed = json.loads(json.dumps(self.tasks))
        changed[0]["test"][0]["isolate_id"] = "z"
        with self.assertRaises(rs.IdentityError):
            rs.task_identity(changed, self.manifest, self.dir / "bench", self.dir / "dec", self.plan)

    def test_changed_phenotype_or_train_fail(self):
        for part, field, value in (("train", "phenotype", "S"), ("train", "isolate_id", "q"),
                                   ("test", "phenotype", "S")):
            changed = json.loads(json.dumps(self.tasks))
            changed[0][part][0][field] = value
            with self.subTest(part=part, field=field), self.assertRaises(rs.IdentityError):
                rs.task_identity(changed, self.manifest, self.dir / "bench", self.dir / "dec", self.plan)

    def test_missing_and_unexpected_tasks(self):
        with self.assertRaisesRegex(rs.IdentityError, "missing"):
            rs.task_identity(self.tasks[:1], self.manifest, self.dir / "bench", self.dir / "dec", self.plan)
        extra = self.tasks + [task("common__x__random__fold1", "x", "random", "fold1", [], [])]
        with self.assertRaisesRegex(rs.IdentityError, "unexpected"):
            rs.task_identity(extra, self.manifest, self.dir / "bench", self.dir / "dec", self.plan)

    def test_r2_record_missing_fails(self):
        (self.dir / "dec" / "augmented" / "permuted_control" / "tasks" / "common__cip__random__fold1.json").unlink()
        with self.assertRaises(rs.IdentityError):
            rs.task_identity(self.tasks, self.manifest, self.dir / "bench", self.dir / "dec", self.plan)

    def test_task_predictions_statuses(self):
        d = self.dir / "run"
        d.mkdir()
        a, b = self.tasks
        (d / f"{a['task_id']}.json").write_text(json.dumps({"status": "success", "split_sha256": gb.task_digest(a),
                                                            "seconds": 1, "features_used": 3}))
        gb.write_tsv(d / f"{a['task_id']}.tsv", gb.PRED_FIELDS,
                     [{"isolate_id": "c", "phenotype": "R", "predicted_probability": 0.9, "predicted_class": "R"},
                      {"isolate_id": "d", "phenotype": "S", "predicted_probability": 0.2, "predicted_class": "S"}])
        (d / f"{b['task_id']}.json").write_text(json.dumps({"status": "failed", "reason": "MemoryError: x"}))
        preds, status = rs.task_predictions(d, self.tasks + [task("t3", "z", "random", "f", [], [])])
        self.assertEqual(status[a["task_id"]]["status"], "completed")
        self.assertEqual(status[b["task_id"]]["status"], "resource_failed")
        self.assertEqual(status["t3"]["status"], "not_run")
        self.assertEqual(preds[("cip", "random")], {"c": 0.9, "d": 0.2})
        (d / f"{b['task_id']}.json").write_text(json.dumps({"status": "failed", "reason": "single_class"}))
        self.assertEqual(rs.task_predictions(d, self.tasks)[1][b["task_id"]]["status"], "failed")
        stale = {**json.loads((d / f"{a['task_id']}.json").read_text()), "split_sha256": "0"}
        (d / f"{a['task_id']}.json").write_text(json.dumps(stale))
        self.assertEqual(rs.task_predictions(d, self.tasks)[1][a["task_id"]]["status"], "failed")
        self.assertEqual(rs.keys_complete({"x": {"status": "completed"}, "y": {"status": "failed"}},
                                          [task("x", "d", "s", "1", [], []), task("y", "d", "s", "2", [], [])]),
                         set())


class ComparisonTests(unittest.TestCase):
    """Paired metrics, gaps, rescue, four-way and correlations on a small synthetic context."""

    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(7)
        n = 80
        iso = [f"I{i:03d}" for i in range(n)]
        y = rng.random(n) < 0.4
        ctx = {}
        for drug in ("cip", "tet"):
            for strat in ("random", "country", "lineage"):
                ctx[(drug, strat)] = {
                    i: {"phenotype": "R" if y[j] else "S", "split_unit": f"u{j // 2}", "isolate_id": i,
                        "determinant_prediction": ("R" if y[j] != (j % 9 == 0) else "S") if j % 11 else NA,
                        "predicted_probability": NA} for j, i in enumerate(iso)}
        cls.ctx = ctx
        noise = {r: rng.normal(0, s, n) for r, s in (("R0", 0.45), ("R1", 0.3), ("R2", 0.15), ("R2_permuted", 0.45),
                                                     ("R3", 0.1), ("R3_permuted", 0.4))}
        preds = {}
        for rep, e in noise.items():
            preds[rep] = {}
            for key in ctx:
                shift = 0.25 if key[1] != "random" else 0.0
                p = np.clip(y * 0.6 + 0.2 + e + rng.normal(0, shift, n), 0.001, 0.999)
                preds[rep][key] = dict(zip(iso, p.tolist()))
        cls.preds = preds
        tasks = [task(f"common__{d}__{s}__f", d, s, "f", [], []) for d, s in ctx]
        cls.scope = {"R0": tasks, "R1": tasks, **{r: [t for t in tasks if t["antimicrobial"] == "cip"]
                                                  for r in ("R2", "R2_permuted", "R3", "R3_permuted")}}
        cls.status = {r: {t["task_id"]: {"status": "completed"} for t in ts} for r, ts in cls.scope.items()}
        for r in ("R2", "R2_permuted", "R3", "R3_permuted"):
            for key in [k for k in preds[r] if k[0] == "tet"]:
                del preds[r][key]
        cls.explain = [{"antimicrobial": d, "holdout_strategy": s, "unseen_mechanism_fraction": str(0.1 * j),
                        "lineage_shift": str(1 - 0.1 * j), "mechanism_shift": "0.2", "phenotype_prevalence_shift":
                        "0.1", "mutation_associated_fraction": "0.3"}
                       for j, (d, s) in enumerate([(d, s) for d in ("cip", "tet") for s in ("country", "lineage")])]
        cfg = json.loads(json.dumps(CONFIG))
        cfg["bootstrap"]["replicates"] = 30
        cls.config = cfg
        gconfig = json.loads((ROOT / CONFIG["upstream"]["generalisation_config"]).read_text())
        gconfig["bootstrap"]["replicates"] = 30
        cls.gconfig = gconfig
        cls.tables, cls.gap, cls.done = rs.comparison_tables(cfg, gconfig, ctx, preds, cls.status, cls.scope,
                                                             cls.explain, {"cip": "high_loss", "tet": "stable"})

    def test_metrics_and_not_applicable(self):
        m = {(r["representation"], r["antimicrobial"], r["split_strategy"]): r for r in self.tables["metrics"]}
        self.assertEqual(m[("R2", "tet", "random")]["status"], "not_applicable")
        r = m[("R1", "cip", "random")]
        key = ("cip", "random")
        y = np.array([self.ctx[key][i]["phenotype"] == "R" for i in sorted(self.ctx[key])])
        p = np.array([self.preds["R1"][key][i] for i in sorted(self.ctx[key])])
        self.assertAlmostEqual(r["AUROC"], gb.auroc(y, p))
        self.assertIn("calibration_slope", r)
        self.assertIn(r["supported"], ("yes", "no"))

    def test_paired_deltas(self):
        d = next(r for r in self.tables["deltas"] if r["representation"] == "R1" and r["antimicrobial"] == "cip"
                 and r["split_strategy"] == "country")
        for m in rs.DELTA_METRICS:
            self.assertAlmostEqual(d[f"delta_{m}"], d[f"Rx_{m}"] - d[f"R0_{m}"])
        self.assertLessEqual(d["delta_AUROC_ci_lower"], d["delta_AUROC_ci_upper"])
        again = rs.comparison_tables(self.config, self.gconfig, self.ctx, self.preds, self.status, self.scope,
                                     self.explain, {})[0]
        self.assertEqual([r["delta_AUROC_ci_lower"] for r in again["deltas"]],
                         [r["delta_AUROC_ci_lower"] for r in self.tables["deltas"]])

    def test_gap_and_rescue(self):
        g = self.gap[("R0", "cip", "country")]
        self.assertAlmostEqual(g["gap_AUROC"], g["AUROC_random"] - g["AUROC_holdout"])
        res = next(r for r in self.tables["rescue"] if r["representation"] == "R3" and r["antimicrobial"] == "cip"
                   and r["holdout_strategy"] == "country")
        self.assertAlmostEqual(res["gap_rescue"], res["R0_gap_AUROC"] - res["Rx_gap_AUROC"])
        self.assertIsNotNone(res["gap_rescue_ci_lower"])
        self.assertFalse(any(r["holdout_strategy"] == "random" for r in self.tables["rescue"]))
        self.assertNotIn(("R2", "tet", "country"), self.gap)

    def test_pairwise(self):
        rows = rs.pairwise_tables(self.config, self.ctx, self.preds, self.status, self.scope, {"cip": "high_loss"})
        self.assertEqual({(r["reference"], r["representation"]) for r in rows}, set(rs.PAIRS))
        self.assertFalse(any(r["antimicrobial"] == "tet" for r in rows))
        r = next(r for r in rows if (r["reference"], r["representation"], r["split_strategy"]) ==
                 ("R1", "R3", "country"))
        self.assertAlmostEqual(r["delta_AUROC"], r["representation_AUROC"] - r["reference_AUROC"])
        self.assertLessEqual(r["delta_AUROC_ci_lower"], r["delta_AUROC_ci_upper"])
        # the gap reduction against R1 equals the difference of the two rescues against R0
        res = {x["representation"]: x["gap_rescue"] for x in self.tables["rescue"]
               if x["antimicrobial"] == "cip" and x["holdout_strategy"] == "country"}
        r1 = self.gap[("R1", "cip", "country")]["gap_AUROC"]
        r0 = self.gap[("R0", "cip", "country")]["gap_AUROC"]
        self.assertAlmostEqual(r["gap_reduction"], res["R3"] - (r0 - r1))
        self.assertIsNotNone(r["gap_reduction_ci_lower"])
        rnd = next(x for x in rows if x["split_strategy"] == "random")
        self.assertIsNone(rnd["gap_reduction"])
        self.assertEqual(rows, rs.pairwise_tables(self.config, self.ctx, self.preds, self.status, self.scope,
                                                  {"cip": "high_loss"}))

    def test_four_way_and_determinant(self):
        fw = [r for r in self.tables["four_way"] if r["representation"] == "R0"]
        for r in fw:
            self.assertEqual(sum(r[k] for k in rs.FOUR_WAY), r["N"])
        pooled = {(r["representation"], r["drug_set"], r["split_strategy"]): r for r in self.tables["four_way_pooled"]}
        self.assertEqual(pooled[("R2", "r2_drugs", "random")]["drugs"], 1)
        self.assertEqual(pooled[("R0", "all", "random")]["drugs"], 2)
        self.assertEqual(pooled[("R0", "all", "random")]["drugs_expected"], 2)
        det = next(r for r in self.tables["determinant"] if r["representation"] == "R1")
        self.assertAlmostEqual(det["ml_minus_determinant_AUROC"], det["ml_AUROC"] - det["determinant_AUROC"])

    def test_correlations_and_robustness(self):
        cor = {(r["representation"], r["drug_set"], r["measure"]): r for r in self.tables["correlations"]}
        self.assertEqual(cor[("R0", "all", "unseen_mechanism_fraction")]["pairs"], 4)
        self.assertEqual(cor[("R2", "r2_drugs", "lineage_shift")]["pairs"], 2)
        self.assertIsNone(cor[("R2", "r2_drugs", "lineage_shift")]["spearman_rho"])
        rob = {(r["representation"], r["drug_set"]): r for r in self.tables["robustness"]}
        self.assertEqual(rob[("R1", "all")]["pairs"], 4)
        self.assertEqual(rob[("R2", "r2_drugs")]["pairs_expected"], 2)
        gaps = [self.gap[("R1", d, s)]["gap_AUROC"] for d in ("cip", "tet") for s in ("country", "lineage")]
        self.assertAlmostEqual(rob[("R1", "all")]["median_gap_AUROC"], float(np.median(gaps)))

    def test_changed_test_endpoints_stop(self):
        preds = json.loads(json.dumps({r: {"|".join(k): v for k, v in d.items()} for r, d in self.preds.items()}))
        preds = {r: {tuple(k.split("|")): v for k, v in d.items()} for r, d in preds.items()}
        preds["R1"][("cip", "random")].pop("I000")
        with self.assertRaises(rs.IdentityError):
            rs.comparison_tables(self.config, self.gconfig, self.ctx, preds, self.status, self.scope, self.explain, {})

    def test_drug_views(self):
        metrics = {(r["representation"], r["antimicrobial"], r["split_strategy"]): r for r in self.tables["metrics"]
                   if r["status"] == "completed"}
        rescue = {(r["representation"], r["antimicrobial"], r["holdout_strategy"]): r for r in self.tables["rescue"]}
        cap = rs.drug_capture([{"antimicrobial": "cip", "determinant_kind": "point_mutation", "determinant": "ALL",
                                "representation": "R0", "captured_fraction": "0.01", "endpoints": "10"},
                               {"antimicrobial": "cip", "determinant_kind": "point_mutation", "determinant": "ALL",
                                "representation": "scaled_20", "captured_fraction": "0.9", "endpoints": "10"},
                               {"antimicrobial": "cip", "determinant_kind": "point_mutation",
                                "determinant": "gyrA_S83L", "representation": "R0", "captured_fraction": "0",
                                "endpoints": "4"}], "scaled_20")
        rows = rs.drug_view(["cip", "absent"], self.ctx, metrics, self.gap, rescue, cap, {}, {}, "quinolone")
        self.assertEqual(len(rows), 3)
        self.assertEqual((rows[0]["R0_mutation_capture"], rows[0]["R1_mutation_capture"]), (0.01, 0.9))
        country = next(r for r in rows if r["split_strategy"] == "country")
        self.assertIsNotNone(country["R3_gap_rescue"])
        self.assertIsNone(rows[0]["R2_gap"])


def load_ed_study():
    """The error_decomposition_v1 synthetic study (itself built on the generalisation_v1 one), run through every
    command."""
    spec = importlib.util.spec_from_file_location("_ed_study", ROOT / "tests/test_error_decomposition.py")
    mod = importlib.util.module_from_spec(spec)
    previous = {k: sys.modules.get(k) for k in ("generalisation_benchmark", "error_decomposition")}
    try:
        spec.loader.exec_module(mod)

        class Study(mod.EndToEndTests):
            pass

        Study.setUpClass()
    finally:
        for k, v in previous.items():
            if v is not None:
                sys.modules[k] = v
    assert Study.codes == [0] * 9, Study.codes
    return Study


class EndToEndTests(unittest.TestCase):
    """Every command on the synthetic study. The synthetic benchmark used scaled 1, so R1 is built at the same
    density: its predictions must equal R0's task for task, and R3's must equal R2's, which checks the whole
    refit path against the frozen runs."""

    @classmethod
    def setUpClass(cls):
        cls.study = study = load_ed_study()
        gstudy = study.study
        root = cls.root = study.root
        bench, dec = gstudy.out, study.out
        shutil.copy(dec / "augmented" / "augmentation_plan.json", dec / "augmentation_plan.json")
        names = sorted(p.name for p in dec.iterdir() if p.is_file() and p.suffix in (".tsv", ".json", ".gz")
                       and p.name != "SHA256SUMS")
        (dec / "SHA256SUMS").write_text("".join(f"{sha(dec / x)}  {x}\n" for x in names))
        ed_config = json.loads(study.config_path.read_text())
        plan = json.loads((dec / "augmentation_plan.json").read_text())
        config = json.loads(json.dumps(CONFIG))
        config["upstream"].update({
            "generalisation_config": str(gstudy.config_path), "generalisation_config_sha256": sha(gstudy.config_path),
            "benchmark_manifest_sha256": sha(bench / "benchmark_manifest.json"),
            "predictions_sha256": sha(bench / "predictions.tsv.gz"),
            "feature_manifest_sha256": sha(bench / "matrix" / "feature_manifest.json"),
            "genome_manifest": str(root / "genome_manifest.tsv"), "genome_manifest_sha256":
                sha(root / "genome_manifest.tsv"), "genome_manifest_rows": 150,
            "error_decomposition_config": str(study.config_path),
            "error_decomposition_config_sha256": sha(study.config_path), "error_decomposition_outputs": str(dec),
            "error_decomposition_sha256sums_sha256": sha(dec / "SHA256SUMS"),
            "augmentation_plan_sha256": sha(dec / "augmentation_plan.json"),
            "mechanism_catalogue": ed_config["catalogue"]["path"],
            "mechanism_catalogue_sha256": sha(ed_config["catalogue"]["path"])})
        config["density"].update({"ksize": 21, "reference_scaled": 1, "candidates": [1], "sketch_scaled": 1})
        config["benchmark"]["tasks"] = json.loads((bench / "benchmark_manifest.json").read_text())["tasks"]
        config["mutation_panel"]["indicators"] = len(plan["indicators"])
        config["capture"].update({"quinolones": ["tetracycline"], "stable_controls": ["ampicillin"],
                                  "mechanism_shift_control": "ampicillin", "unresolved_controls": ["ampicillin"],
                                  "counterexample": "tetracycline"})
        config["bootstrap"]["replicates"] = 20
        config["support"] = {"metric_min_R": 5, "metric_min_S": 5}
        cls.config_path = root / "representation.json"
        cls.config_path.write_text(json.dumps(config))
        cls.config = config
        out = cls.out = root / "representation"
        base = ["--config", str(cls.config_path), "--benchmark", str(bench), "--decomposition", str(dec),
                "--outdir", str(out)]
        inputs = ["--population", str(root / "population"), "--baseline", str(root / "baseline")]
        budget = root / "budget.json"
        budget.write_text(json.dumps({"budget_su": 10, "node_memory_gb": 64, "su_per_core_hour": 1.25,
                                      "memory_gb_per_core": 9.14, "cpu_speed_factor": 1.5, "r3_tasks": 2,
                                      "measured": {"1": {"peak_rss_bytes": 1e8, "task_hours_total": 0.01}}}))
        cls.args = {"base": base, "inputs": inputs}
        with patch("sys.stdout"):
            cls.codes = [rs.main(["check", *base]),
                         rs.main(["sketch", *base, "--cache", str(gstudy.cache), "--baseline", str(root / "baseline"),
                                  "--workers", "2"]),
                         rs.main(["pilot", *base, "--cache", str(gstudy.cache), "--workers", "2"]),
                         rs.main(["fitpilot", *base, *inputs, "--scaled", "1"]),
                         rs.main(["select", *base, "--budget", str(budget)]),
                         rs.main(["matrix", *base]),
                         rs.main(["run", *base, *inputs, "--representation", "R1", "--only", "__random__"]),
                         rs.main(["run", *base, *inputs, "--representation", "R1"]),
                         rs.main(["run", *base, *inputs, "--representation", "R3"]),
                         rs.main(["run", *base, *inputs, "--representation", "R3_permuted"]),
                         rs.main(["report", *base, *inputs])]
        cls.bench, cls.dec, cls.cache = bench, dec, gstudy.cache

    @classmethod
    def tearDownClass(cls):
        cls.study.tearDownClass()

    def table(self, name):
        return ed.read_tsv(self.out / "tables" / name)

    def test_commands_succeed(self):
        self.assertEqual(self.codes, [0] * 11)
        sel = json.loads((self.out / "selection.json").read_text())
        self.assertEqual(sel["R1_scaled"], 1)
        with self.assertRaisesRegex(ValueError, "already selected"), patch("sys.stdout"):
            rs.main(["select", *self.args["base"], "--budget", str(self.root / "budget.json")])

    def test_genomes_verified_and_r0_reproduced(self):
        rows = ed.read_tsv(self.out / "genome_verification.tsv")
        self.assertEqual(Counter(r["status"] for r in rows), Counter({"genome_verified": 150}))
        self.assertEqual({r["r0_reproduced"] for r in rows}, {"yes"})
        pilot = json.loads((self.out / "pilot" / "pilot_summary.json").read_text())
        self.assertTrue(all(pilot["r0_feature_space_check"].values()))
        self.assertEqual(pilot["r0_capture_mismatches"], {})

    def test_fit_pilot_records_resources_only(self):
        fit = json.loads((self.out / "pilot" / "fit_pilot_scaled_1.json").read_text())
        self.assertEqual(len(fit["fits"]), 3)
        self.assertTrue(all(f["peak_rss_bytes"] > 0 and f["seconds"] >= 0 for f in fit["fits"]))
        self.assertFalse(any("AUROC" in k or "predicted" in k for f in fit["fits"] for k in f))
        self.assertFalse((self.out / "pilot_scaled_1" / "tasks").exists())
        with self.assertRaisesRegex(ValueError, "not a candidate"), patch("sys.stdout"):
            rs.main(["fitpilot", *self.args["base"], *self.args["inputs"], "--scaled", "7"])

    def test_fit_pilot_runs_without_a_local_pilot(self):
        # on the cluster only the sketches and verification are copied, so pilot/ does not exist yet
        fresh = self.out.parent / "fresh_fitpilot"
        fresh.mkdir(exist_ok=True)
        for name in ("representation_sensitivity_manifest.json", "genome_verification.tsv"):
            shutil.copy(self.out / name, fresh / name)
        (fresh / "sketches").symlink_to(self.out / "sketches")
        base = [*self.args["base"][:-1], str(fresh)]
        with patch("sys.stdout"):
            self.assertEqual(rs.main(["fitpilot", *base, *self.args["inputs"], "--scaled", "1"]), 0)
        self.assertTrue((fresh / "pilot" / "fit_pilot_scaled_1.json").exists())

    def test_pilot_reads_feature_spaces_from_built_matrices(self):
        fresh = self.out.parent / "fresh_pilot"
        fresh.mkdir(exist_ok=True)
        for name in ("representation_sensitivity_manifest.json", "genome_verification.tsv"):
            shutil.copy(self.out / name, fresh / name)
        (fresh / "sketches").symlink_to(self.out / "sketches")
        base = [*self.args["base"][:-1], str(fresh)]
        with patch("sys.stdout"):
            rs.main(["fitpilot", *base, *self.args["inputs"], "--scaled", "1"])
            self.assertEqual(rs.main(["pilot", *base, "--cache", str(self.cache)]), 0)
        summary = json.loads((fresh / "pilot" / "pilot_summary.json").read_text())
        self.assertEqual(summary["feature_space_source"], "pilot_scaled_*/matrix")
        self.assertTrue(all(summary["r0_feature_space_check"].values()))
        for name in ("mutation_capture_by_candidate.tsv", "determinant_capture_by_candidate.tsv",
                     "drug_capture_by_candidate.tsv"):
            self.assertEqual((fresh / "pilot" / name).read_text(), (self.out / "pilot" / name).read_text())
        drop = {"build_seconds", "build_max_rss_bytes"}
        strip = lambda path: [{k: v for k, v in r.items() if k not in drop} for r in rs.ed.read_tsv(path)]  # noqa: E731
        self.assertEqual(strip(fresh / "pilot" / "representation_resource_pilot.tsv"),
                         strip(self.out / "pilot" / "representation_resource_pilot.tsv"))
        rows = fresh / "pilot_scaled_1" / "matrix" / "rows.tsv"
        rows.write_text(rows.read_text().replace("\n", "\nextra\t-\n", 1))
        with self.assertRaisesRegex(rs.UpstreamError, "other genomes"), patch("sys.stdout"):
            rs.main(["pilot", *base, "--cache", str(self.cache)])

    def test_packed_run_gives_the_same_predictions(self):
        fresh = self.out.parent / "packed_run"
        if fresh.exists():
            shutil.rmtree(fresh)
        fresh.mkdir()
        for name in ("representation_sensitivity_manifest.json", "genome_verification.tsv", "selection.json"):
            shutil.copy(self.out / name, fresh / name)
        for name in ("pilot", "R1/matrix"):
            (fresh / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(self.out / name, fresh / name)
        base = [*self.args["base"][:-1], str(fresh)]
        with patch("sys.stdout"):
            code = rs.main(["run", *base, *self.args["inputs"], "--representation", "R1", "--workers", "3",
                            "--memory-gb", "64"])
        self.assertEqual(code, 0)
        a, b = sorted((self.out / "R1" / "tasks").glob("*.tsv")), sorted((fresh / "R1" / "tasks").glob("*.tsv"))
        self.assertEqual([x.name for x in a], [x.name for x in b])
        self.assertTrue(all(x.read_bytes() == y.read_bytes() for x, y in zip(a, b)))
        rec = json.loads(next((fresh / "R1" / "tasks").glob("*.json")).read_text())
        self.assertGreater(rec["max_rss_bytes"], 0)
        with patch("sys.stdout"):
            rs.main(["run", *base, *self.args["inputs"], "--representation", "R1", "--memory-gb", "64"])
        self.assertEqual(json.loads((fresh / "R1" / "tasks" / b[0].name.replace(".tsv", ".json")).read_text())
                         ["finished_at"], json.loads(b[0].with_suffix(".json").read_text())["finished_at"])

    def test_packing_never_exceeds_the_memory_cap(self):
        e = {"isolate_id": "a", "phenotype": "R", "split_unit": "u"}
        tasks = [{"task_id": f"t{i}", "train": [e] * n, "test": []} for i, n in enumerate([8, 6, 5, 3, 2, 1, 1])]
        live, peak = [0], [0]

        class Future:
            def __init__(self, n):
                self.n = n

        class Pool:
            def __init__(self, *a, **k):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def submit(self, fn, args):
                live[0] += len(args[0]["train"])
                peak[0] = max(peak[0], live[0])
                f = Future(len(args[0]["train"]))
                f.result = lambda: ({"task_id": args[0]["task_id"], "status": "success", "max_rss_bytes": 1}, [])
                return f

        def wait(running, return_when):
            f = max(running, key=lambda f: f.n)
            live[0] -= f.n
            return {f}, set()

        with tempfile.TemporaryDirectory() as d, patch.object(rs, "ProcessPoolExecutor", Pool), \
                patch.object(rs, "wait", wait), patch.object(rs.gb, "task_digest", lambda t: "x"), \
                patch("sys.stdout"), patch.object(rs.gb, "write_json"), patch.object(rs, "write_tsv"):
            counts = rs.run_tasks_packed(tasks, Path(d), Path(d), {}, workers=4, memory_cap_bytes=10,
                                         bytes_per_train=1.0)
        self.assertEqual(counts["success"], 7)
        self.assertLessEqual(peak[0], 10)

    def test_r1_matrix_equals_r0_matrix(self):
        a = json.loads((self.bench / "matrix" / "feature_manifest.json").read_text())
        b = json.loads((self.out / "R1" / "matrix" / "feature_manifest.json").read_text())
        for k in ("indices_sha256", "indptr_sha256", "feature_hashes_sha256", "feature_pattern_sha256", "nnz"):
            self.assertEqual(a[k], b[k], k)

    def test_task_identity(self):
        rows = ed.read_tsv(self.out / "task_identity_validation.tsv")
        self.assertEqual(len(rows), self.config["benchmark"]["tasks"])
        self.assertEqual({r["status"] for r in rows}, {"identical"})

    def test_r1_reproduces_r0_and_r3_reproduces_r2(self):
        for rep, ref in (("R1", self.bench / "tasks"), ("R3", self.dec / "augmented" / "mutation_augmented" / "tasks"),
                         ("R3_permuted", self.dec / "augmented" / "permuted_control" / "tasks")):
            runs = sorted((self.out / rep / "tasks").glob("*.tsv"))
            self.assertTrue(runs)
            for f in runs:
                with self.subTest(rep=rep, task=f.name):
                    self.assertEqual(f.read_text(), (ref / f.name).read_text())

    def test_original_outputs_not_overwritten(self):
        self.assertEqual(sha(self.bench / "predictions.tsv.gz"), self.config["upstream"]["predictions_sha256"])
        for line in (self.dec / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split()
            self.assertEqual(sha(self.dec / name), digest)

    def test_tables_written(self):
        for name in ("representation_metrics.tsv", "representation_metric_deltas.tsv",
                     "generalisation_gap_by_representation.tsv", "generalisation_gap_rescue.tsv",
                     "determinant_vs_representation.tsv", "rule_ml_discordance.tsv", "rule_ml_discordance_pooled.tsv",
                     "all_drug_robustness.tsv", "correlation_robustness.tsv", "quinolone_representation_analysis.tsv",
                     "stable_drug_controls.tsv", "cefazolin_control.tsv", "unresolved_drug_controls.tsv",
                     "levofloxacin_counterexample.tsv", "manuscript_primary_table.tsv",
                     "manuscript_quinolone_table.tsv", "task_status.tsv"):
            with self.subTest(name=name):
                self.assertTrue(self.table(name))
        for name in ("fig_a_mutation_capture_by_density.tsv", "fig_b_R0_vs_R1_holdout_auroc.tsv",
                     "fig_d_gap_by_representation.tsv", "fig_e_determinant_vs_representation.tsv",
                     "fig_f_gain_vs_cost.tsv"):
            self.assertTrue(ed.read_tsv(self.out / "tables" / "figure_data" / name))
        self.assertTrue((self.out / "R1" / "benchmark_metrics.tsv").exists())

    def test_r0_tables_reproduce_upstream(self):
        drops = {(r["antimicrobial"], r["holdout_strategy"]): float(r["delta_AUROC"])
                 for r in ed.read_tsv(self.bench / "generalisation_drop.tsv")
                 if r["cohort"] == "common" and r["comparator"] == "random" and r["level"] == "pooled"
                 and r["delta_AUROC"] != NA}
        gaps = {(r["antimicrobial"], r["holdout_strategy"]): float(r["gap_AUROC"])
                for r in self.table("generalisation_gap_by_representation.tsv")
                if r["representation"] == "R0" and r["gap_AUROC"] != NA}
        self.assertTrue(drops)
        for k, v in drops.items():
            self.assertAlmostEqual(gaps[k], v, places=5)
        # R1 at the same density as R0: zero deltas, zero rescue
        for r in self.table("representation_metric_deltas.tsv"):
            if r["representation"] == "R1" and r["delta_AUROC"] != NA:
                self.assertAlmostEqual(float(r["delta_AUROC"]), 0.0, places=9)
        rescue = {(r["antimicrobial"], r["holdout_strategy"]): float(r["gap_rescue"])
                  for r in ed.read_tsv(self.dec / "generalisation_gap_rescue.tsv")
                  if r["variant"] == "mutation_augmented"}
        mine = {(r["antimicrobial"], r["holdout_strategy"]): float(r["gap_rescue"])
                for r in self.table("generalisation_gap_rescue.tsv") if r["representation"] == "R2"
                and r["gap_rescue"] != NA}
        self.assertTrue(rescue)
        for k, v in rescue.items():
            self.assertAlmostEqual(mine[k], v, places=5)
        fw = {(r["antimicrobial"], r["split_strategy"]): r for r in ed.read_tsv(self.dec / "four_way_comparison.tsv")
              if r["cohort"] == "common" and r["antimicrobial"] != "ALL"}
        for r in self.table("rule_ml_discordance.tsv"):
            if r["representation"] == "R0":
                ref = fw[(r["antimicrobial"], r["split_strategy"])]
                self.assertEqual([r[k] for k in rs.FOUR_WAY], [ref[k] for k in rs.FOUR_WAY])
        summary = json.loads((self.dec / "error_decomposition_summary.json").read_text())["associations"]
        cor = {r["measure"]: r for r in self.table("correlation_robustness.tsv")
               if r["representation"] == "R0" and r["drug_set"] == "all"}
        for m in ("unseen_mechanism_fraction", "lineage_shift"):
            ref = summary[m]["rho"]
            got = cor[m]["spearman_rho"]
            if ref is None:
                self.assertEqual(got, NA)
            else:
                self.assertAlmostEqual(float(got), ref, places=5)

    def test_report_is_deterministic(self):
        before = (self.out / "tables" / "generalisation_gap_rescue.tsv").read_bytes()
        with patch("sys.stdout"):
            rs.main(["report", *self.args["base"], *self.args["inputs"]])
        self.assertEqual((self.out / "tables" / "generalisation_gap_rescue.tsv").read_bytes(), before)

    def test_changed_config_or_inputs_stop(self):
        changed = json.loads(self.config_path.read_text())
        changed["analysis"]["cohort"] = "maximal"
        path = self.root / "changed.json"
        path.write_text(json.dumps(changed))
        args = ["--config", str(path), *self.args["base"][2:]]
        with self.assertRaisesRegex(ValueError, "changed after check"):
            rs.main(["report", *args, *self.args["inputs"]])
        bad = json.loads(self.config_path.read_text())
        bad["upstream"]["predictions_sha256"] = "0" * 64
        path.write_text(json.dumps(bad))
        with self.assertRaises(rs.UpstreamError):
            rs.main(["check", "--config", str(path), *self.args["base"][2:]])

    def test_cli_requires_inputs(self):
        with patch("sys.stderr"), self.assertRaises(SystemExit):
            rs.main(["run", *self.args["base"]])
        with patch("sys.stderr"), self.assertRaises(SystemExit):
            rs.main(["sketch", *self.args["base"]])


if __name__ == "__main__":
    unittest.main()
