"""Offline tests for the established-determinant genotype-phenotype baseline."""

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
spec = importlib.util.spec_from_file_location("genotype_baseline", ROOT / "bin/genotype_baseline.py")
gb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gb)
amrscan = gb.amrscan
NA = gb.NA
VERSION = "test_mapping_v1"
ORG = "Escherichia coli"
RAW_HEADER = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name", "Scope",
              "Type", "Subtype", "Class", "Subclass", "Method", "Target length", "Reference sequence length",
              "% Coverage of reference", "% Identity to reference", "Alignment length",
              "Closest reference accession", "Closest reference name", "HMM accession", "HMM description",
              "Hierarchy node"]


def rule(det, drug, interp="predicts_resistance", requires=NA, dtype="gene", version=VERSION, org=ORG):
    return {"mapping_version": version, "organism": org, "rule_id": f"{det}|{drug}", "determinant": det,
            "determinant_type": dtype, "gene_family": det.split("-")[0], "allele_or_mutation": det,
            "antimicrobial": drug, "interpretation": interp, "requires": requires, "evidence_source": "test",
            "evidence_terms": drug, "evidence_reference": "PMID 1", "evidence_level": "source_exact_term",
            "review_status": "reviewed_source", "notes": NA}


STANDARD_RULES = [
    rule("blaTEM-1", "ampicillin"),
    rule("blaCTX-M-15", "ampicillin"), rule("blaCTX-M-15", "ceftriaxone"),
    rule("blaEC", "ampicillin", "does_not_predict"),
    rule("blaTEM", "ceftriaxone", "uncertain"),
    rule("gyrA_S83L", "ciprofloxacin", dtype="point_mutation"),
    rule("parC_S80I", "ciprofloxacin", requires="gyrA_S83L|gyrA_S83W", dtype="point_mutation"),
    rule("qnrS1", "ciprofloxacin"),
    rule("dfrA17", "trimethoprim-sulfamethoxazole", "component", "trimethoprim"),
    rule("sul1", "trimethoprim-sulfamethoxazole", "component", "sulfonamide"),
]
STUDY = {"ampicillin": "primary", "ceftriaxone": "primary", "ciprofloxacin": "primary",
         "trimethoprim-sulfamethoxazole": "primary"}
COMBOS = {"trimethoprim-sulfamethoxazole": {"components": {"trimethoprim": ["trimethoprim"],
                                                           "sulfonamide": ["sulfamethoxazole"]}}}


def det(symbol, method="EXACTX", etype="AMR", subtype="AMR"):
    return {"element_type": etype, "gene": symbol, "evidence_type": method, "element_subtype": subtype,
            "subclass": "X"}


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write_mapping(self, rows, study=STUDY, combos=COMBOS, name="map", sha=None, version=VERSION, raw=None):
        tsv = self.root / f"{name}.tsv"
        if raw is not None:
            tsv.write_text(raw)
        else:
            with tsv.open("w", newline="") as h:
                w = csv.DictWriter(h, gb.MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        meta = {"mapping_version": version, "rules_sha256": sha or amrscan.sha256(tsv), "study_antimicrobials": study,
                "combination_rules": combos}
        (self.root / f"{name}.json").write_text(json.dumps(meta))
        return tsv, self.root / f"{name}.json"

    def mapping(self, rows=STANDARD_RULES, **kw):
        return gb.load_mapping(*self.write_mapping(rows, **kw))

    def status(self, symbols, drug, rows=STANDARD_RULES, methods=None):
        m = self.mapping(rows)
        dets = gb.amr_determinants([det(s, (methods or {}).get(s, "EXACTX")) for s in symbols])
        return gb.interpret(dets, m["rules"].get(drug, {}), m["combinations"].get(drug))


# ---------------------------------------------------------------- mapping

class MappingTests(Base):
    def test_exact_mapping(self):
        self.assertEqual(self.status(["blaTEM-1"], "ampicillin")["status"], "positive")

    def test_unmapped_determinant_is_negative(self):
        self.assertEqual(self.status(["mph(A)"], "ampicillin")["status"], "negative")

    def test_uncertain_mapping(self):
        r = self.status(["blaTEM"], "ceftriaxone")
        self.assertEqual(r["status"], "uncertain")
        self.assertEqual(r["uncertain"], ["blaTEM"])

    def test_wrong_antimicrobial(self):
        self.assertEqual(self.status(["blaTEM-1"], "ceftriaxone")["status"], "negative")
        with self.assertRaisesRegex(ValueError, "outside the study set"):
            self.mapping([rule("blaTEM-1", "colistin")])

    def test_wrong_organism(self):
        with self.assertRaisesRegex(ValueError, "organism"):
            self.mapping([rule("blaTEM-1", "ampicillin", org="Klebsiella pneumoniae")])

    def test_mutation_specific_mapping(self):
        self.assertEqual(self.status(["gyrA_S83L"], "ciprofloxacin")["status"], "positive")
        self.assertEqual(self.status(["gyrA_S83A"], "ciprofloxacin")["status"], "negative")

    def test_mutation_with_unmet_requirement_is_uncertain(self):
        self.assertEqual(self.status(["parC_S80I"], "ciprofloxacin")["status"], "uncertain")
        r = self.status(["parC_S80I", "gyrA_S83L"], "ciprofloxacin")
        self.assertEqual(r["positive"], ["gyrA_S83L", "parC_S80I"])

    def test_requirement_needs_full_evidence(self):
        r = self.status(["parC_S80I", "gyrA_S83L"], "ciprofloxacin", methods={"gyrA_S83L": "PARTIALX"})
        self.assertEqual(r["status"], "uncertain")

    def test_allele_specific_mapping(self):
        self.assertEqual(self.status(["blaCTX-M-15"], "ceftriaxone")["status"], "positive")
        self.assertEqual(self.status(["blaCTX-M-14"], "ceftriaxone")["status"], "negative")

    def test_multiple_determinants_preserved_in_order(self):
        r = self.status(["qnrS1", "gyrA_S83L"], "ciprofloxacin")
        self.assertEqual(r["positive"], ["gyrA_S83L", "qnrS1"])

    def test_one_determinant_several_drugs(self):
        m = self.mapping()
        dets = gb.amr_determinants([det("blaCTX-M-15")])
        for drug in ("ampicillin", "ceftriaxone"):
            self.assertEqual(gb.interpret(dets, m["rules"][drug])["status"], "positive")

    def test_stress_and_virulence_ignored(self):
        dets = gb.amr_determinants([det("blaTEM-1", etype="STRESS"), det("blaTEM-1", etype="VIRULENCE")])
        self.assertEqual(dets, [])

    def test_partial_hit_is_uncertain(self):
        for method in ("PARTIALX", "PARTIALP", "INTERNAL_STOP"):
            self.assertEqual(self.status(["blaTEM-1"], "ampicillin", methods={"blaTEM-1": method})["status"],
                             "uncertain")
        dets = gb.amr_determinants([det("blaTEM-1", "PARTIALX"), det("blaTEM-1", "EXACTX")])
        self.assertFalse(dets[0]["partial"])

    def test_does_not_predict_rule(self):
        self.assertEqual(self.status(["blaEC"], "ampicillin")["status"], "negative")

    def test_combination_rule(self):
        drug = "trimethoprim-sulfamethoxazole"
        self.assertEqual(self.status(["dfrA17", "sul1"], drug)["status"], "positive")
        self.assertEqual(self.status(["sul1"], drug)["status"], "uncertain")
        self.assertEqual(self.status(["dfrA17"], drug)["status"], "uncertain")
        self.assertEqual(self.status([], drug)["status"], "negative")

    def test_human_review_interpretation_counts_as_uncertain(self):
        rows = [rule("blaTEM-1", "ampicillin", "HUMAN_REVIEW_REQUIRED")]
        self.assertEqual(self.status(["blaTEM-1"], "ampicillin", rows)["status"], "uncertain")

    def test_mapping_version_validation(self):
        tsv, meta = self.write_mapping(STANDARD_RULES)
        with self.assertRaisesRegex(ValueError, "differs from expected"):
            gb.load_mapping(tsv, meta, "other_version")
        with self.assertRaisesRegex(ValueError, "mapping version"):
            self.mapping([rule("blaTEM-1", "ampicillin", version="v0")])
        (self.root / "nov.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "lacks mapping_version"):
            gb.load_mapping(tsv, self.root / "nov.json")

    def test_mapping_checksum_validation(self):
        with self.assertRaisesRegex(ValueError, "checksum"):
            gb.load_mapping(*self.write_mapping(STANDARD_RULES, sha="0" * 64))

    def test_duplicate_and_conflicting_rules(self):
        with self.assertRaisesRegex(ValueError, "duplicate mapping rule"):
            self.mapping(STANDARD_RULES + [rule("blaTEM-1", "ampicillin")])
        with self.assertRaisesRegex(ValueError, "conflicting"):
            self.mapping(STANDARD_RULES + [rule("blaTEM-1", "ampicillin", "uncertain")])

    def test_malformed_mapping_rows(self):
        bad = rule("blaTEM-1", "ampicillin", "maybe")
        with self.assertRaisesRegex(ValueError, "unknown interpretation"):
            self.mapping([bad])
        bad = dict(rule("blaTEM-1", "ampicillin"), rule_id="x")
        with self.assertRaisesRegex(ValueError, "rule id"):
            self.mapping([bad])
        with self.assertRaisesRegex(ValueError, "undefined component"):
            self.mapping([rule("dfrA17", "trimethoprim-sulfamethoxazole", "component", "folate")])
        with self.assertRaisesRegex(ValueError, "malformed row"):
            self.mapping(raw="\t".join(gb.MAPPING_FIELDS) + "\nonly\ttwo\n")

    def test_generated_mapping_file_is_valid(self):
        m = gb.load_mapping(ROOT / "benchmark/mappings/ecoli_amrfinder_drug_mapping_v1.tsv",
                            ROOT / "benchmark/mappings/ecoli_amrfinder_drug_mapping_v1.json",
                            "ecoli_amrfinder_drug_mapping_v1")
        self.assertEqual(len(m["study"]), 25)
        self.assertEqual(m["rules"]["ampicillin"]["blaEC"]["interpretation"], "does_not_predict")
        self.assertEqual(m["rules"]["ampicillin"]["blaTEM-1"]["interpretation"], "predicts_resistance")
        self.assertNotIn("blaTEM-1", m["rules"]["amoxicillin-clavulanic acid"])
        self.assertEqual(m["rules"]["trimethoprim-sulfamethoxazole"]["sul1"]["requires"], "sulfonamide")
        self.assertNotIn("rmtE1", m["rules"]["gentamicin"])

    def test_exploratory_mapping_v2_differs_only_for_parent_gene_alleles(self):
        base = ROOT / "benchmark/mappings/ecoli_amrfinder_drug_mapping_v"
        v1 = gb.load_mapping(Path(f"{base}1.tsv"), Path(f"{base}1.json"), "ecoli_amrfinder_drug_mapping_v1")
        v2 = gb.load_mapping(Path(f"{base}2.tsv"), Path(f"{base}2.json"), "ecoli_amrfinder_drug_mapping_v2")
        self.assertEqual(json.loads(Path(f"{base}2.json").read_text())["status"], "exploratory")
        changed = set()
        for drug in v1["study"]:
            a, b = v1["rules"].get(drug, {}), v2["rules"].get(drug, {})
            for det in set(a) | set(b):
                ra, rb = a.get(det), b.get(det)
                if (ra and (ra["interpretation"], ra["requires"])) != (rb and (rb["interpretation"], rb["requires"])):
                    changed.add(det)
        self.assertEqual(changed, {"rmtE1", "rmtE2", "rmtF1", "rmtF2", "rmtB1", "rmtB3", "rmtB4", "rmtD1",
                                   "blaTRU-1", "blaNPS-1"})
        self.assertEqual(v2["rules"]["gentamicin"]["rmtE1"]["interpretation"], "predicts_resistance")


# ---------------------------------------------------------------- variants and metrics

class MetricTests(unittest.TestCase):
    def test_counts_and_metrics(self):
        pairs = [("R", "positive")] * 8 + [("R", "negative")] * 2 + [("S", "positive")] * 1 + \
            [("S", "negative")] * 9 + [("R", "uncertain")]
        tp, fn, fp, tn = gb.confusion(pairs)
        self.assertEqual((tp, fn, fp, tn), (8, 2, 1, 9))
        m = gb.metrics(tp, fn, fp, tn)
        self.assertEqual((m["N"], m["N_R"], m["N_S"]), (20, 10, 10))
        self.assertAlmostEqual(m["sensitivity"], 0.8)
        self.assertAlmostEqual(m["specificity"], 0.9)
        self.assertAlmostEqual(m["PPV"], 8 / 9)
        self.assertAlmostEqual(m["NPV"], 9 / 11)
        self.assertAlmostEqual(m["accuracy"], 0.85)
        self.assertAlmostEqual(m["balanced_accuracy"], 0.85)
        self.assertAlmostEqual(m["F1"], 16 / 19)
        self.assertAlmostEqual(m["MCC"], (72 - 2) / math.sqrt(9 * 10 * 11 * 10))

    def test_wilson_interval(self):
        low, high = gb.wilson(8, 10)
        self.assertAlmostEqual(low, 0.4902, places=4)
        self.assertAlmostEqual(high, 0.9433, places=4)
        self.assertEqual(gb.wilson(0, 0), (None, None))
        self.assertEqual(gb.wilson(5, 5)[1], 1.0)
        self.assertEqual(gb.metrics(1, 1, 1, 1)["sensitivity_ci_low"], gb.metrics(1, 1, 1, 1)["sensitivity_ci_low"])

    def test_zero_denominators(self):
        m = gb.metrics(0, 0, 0, 0)
        for k in ("sensitivity", "specificity", "PPV", "NPV", "accuracy", "balanced_accuracy", "F1", "MCC"):
            self.assertIsNone(m[k])
        self.assertEqual(gb.fmt(m["MCC"]), "NA")

    def test_all_r_and_all_s(self):
        m = gb.metrics(5, 5, 0, 0)
        self.assertEqual(m["sensitivity"], 0.5)
        self.assertIsNone(m["specificity"])
        self.assertIsNone(m["balanced_accuracy"])
        self.assertIsNone(m["MCC"])
        m = gb.metrics(0, 0, 3, 7)
        self.assertAlmostEqual(m["specificity"], 0.7)
        self.assertIsNone(m["sensitivity"])
        self.assertEqual(m["PPV"], 0.0)

    def test_perfect_and_discordant(self):
        m = gb.metrics(4, 0, 0, 6)
        self.assertEqual((m["sensitivity"], m["specificity"], m["MCC"], m["F1"]), (1.0, 1.0, 1.0, 1.0))
        m = gb.metrics(0, 4, 6, 0)
        self.assertEqual((m["sensitivity"], m["specificity"], m["MCC"], m["F1"]), (0.0, 0.0, -1.0, 0.0))

    def test_rule_variants(self):
        unc = {"status": "uncertain", "positive": [], "uncertain": ["x"]}
        self.assertEqual(gb.variant_status(unc, "uncertain_as_positive", "ampicillin"), "positive")
        self.assertEqual(gb.variant_status(unc, "uncertain_as_negative", "ampicillin"), "negative")
        self.assertEqual(gb.variant_status(unc, "primary", "ampicillin"), "uncertain")
        one = {"status": "positive", "positive": ["gyrA_S83L"], "uncertain": []}
        two = {"status": "positive", "positive": ["gyrA_S83L", "qnrS1"], "uncertain": []}
        self.assertEqual(gb.variant_status(one, "fq_two_mechanisms", "ciprofloxacin"), "uncertain")
        self.assertEqual(gb.variant_status(two, "fq_two_mechanisms", "levofloxacin"), "positive")
        self.assertEqual(gb.variant_status(unc, "fq_two_mechanisms", "ciprofloxacin"), "uncertain")
        self.assertIsNone(gb.variant_status(two, "fq_two_mechanisms", "ampicillin"))


# ---------------------------------------------------------------- tables

class TableTests(Base):
    def write(self, text):
        path = self.root / "t.tsv"
        path.write_text(text)
        return path

    def test_strict_reader(self):
        with self.assertRaisesRegex(ValueError, "malformed"):
            gb.read_table(self.write("a\tb\n1\n"))
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            gb.read_table(self.write("a\ta\n1\t2\n"))
        with self.assertRaisesRegex(ValueError, "missing columns"):
            gb.read_table(self.write("a\tb\n1\t2\n"), ("c",))
        with self.assertRaisesRegex(ValueError, "missing header"):
            gb.read_table(self.write(""))
        self.assertEqual(gb.read_table(self.write("a\tb\n")), [])


# ---------------------------------------------------------------- caller layer and end-to-end build

class CohortFixture(Base):
    """Synthetic PR #5 cohort, PR #6 population build, genome cache and calls."""

    ISOLATES = [
        # isolate, ST, country, year, period, source, project, cluster
        ("PDT000000001.1", "ST131", "USA", "2015", "earlier", "human", "PRJNA1", "PDS000000001"),
        ("PDT000000002.1", "ST131", "USA", "2016", "earlier", "human", "PRJNA1", "PDS000000001"),
        ("PDT000000003.1", "ST10", "Canada", "2022", "recent", "food", "PRJNA2", "NA"),
        ("PDT000000004.1", "ST10", "Canada", "2023", "recent", "unknown", "PRJNA2", "NA"),
        ("PDT000000005.1", "UNKNOWN", "NA", "NA", "NA", "unknown", "PRJNA3", "NA"),
        ("PDT000000006.1", "ST69", "USA", "2019", "earlier", "human", "PRJNA1", "NA"),
    ]
    HITS = {
        "PDT000000001.1": [("blaTEM-1", "AMR", "EXACTX"), ("blaEC", "AMR", "ALLELEX"), ("gyrA_S83L", "AMR", "POINTX")],
        "PDT000000002.1": [("blaCTX-M-15", "AMR", "ALLELEX"), ("dfrA17", "AMR", "EXACTX"), ("sul1", "AMR", "EXACTX")],
        "PDT000000003.1": [("blaTEM", "AMR", "BLASTX"), ("asr", "STRESS", "EXACTX")],
        "PDT000000004.1": [],
        "PDT000000005.1": [("blaTEM-1", "AMR", "EXACTX")],
    }
    PHENO = {  # isolate -> {drug: phenotype}
        "PDT000000001.1": {"ampicillin": "R", "ceftriaxone": "S", "ciprofloxacin": "R", "colistin": "S"},
        "PDT000000002.1": {"ampicillin": "R", "ceftriaxone": "R", "trimethoprim-sulfamethoxazole": "R"},
        "PDT000000003.1": {"ampicillin": "S", "ceftriaxone": "R", "ciprofloxacin": "CONFLICT"},
        "PDT000000004.1": {"ampicillin": "R", "ceftriaxone": "S", "ciprofloxacin": "S"},
        "PDT000000005.1": {"ampicillin": "R"},
        "PDT000000006.1": {"ampicillin": "S"},
    }

    def setUp(self):
        super().setUp()
        self.processed = self.root / "pr5"
        self.population = self.root / "pr6"
        self.cache = self.root / "cache"
        self.calls = self.root / "calls"
        self.database = self.root / "db"
        self.database.mkdir()
        (self.database / "AMR.LIB").write_text("fixture database\n")
        self.database_id = gb.database_fingerprint(self.database)
        self.rules, self.meta = self.write_mapping(STANDARD_RULES)

    def write_rows(self, path, rows):
        with path.open("w", newline="") as h:
            w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    def sums(self, folder, names):
        (folder / "SHA256SUMS").write_text("".join(f"{amrscan.sha256(folder / n)}  {n}\n" for n in names))

    def build_inputs(self, endpoint_extra=None, population_changes=None):
        self.processed.mkdir(exist_ok=True)
        self.population.mkdir(exist_ok=True)
        cands, eps, meta, groups = [], [], [], []
        for n, (iso, st, country, year, period, source, project, cluster) in enumerate(self.ISOLATES, 1):
            asm, bs = f"GCA_00000000{n}.1", f"SAMN0000000{n}"
            cands.append({"cohort_id": "ecoli_research_v1", "isolate_id": iso, "pdt_accession": iso,
                          "biosample_accession": bs, "assembly_accession": asm})
            for drug, ph in sorted(self.PHENO[iso].items()):
                eps.append({"cohort_id": "ecoli_research_v1", "isolate_id": iso, "pdt_accession": iso,
                            "biosample_accession": bs, "assembly_accession": asm, "antimicrobial": drug,
                            "phenotype_endpoint": ph, "phenotype_status": "binary" if ph in "RS" else "conflict",
                            "eligible_binary_rs": "true" if ph in ("R", "S") else "false",
                            "country": country, "collection_year": year, "source_category": source})
            genome = self.cache / asm / "genome.fna"
            genome.parent.mkdir(parents=True, exist_ok=True)
            genome.write_text(f">c\nACGT{n}\n")
            qc = "fail" if iso == "PDT000000006.1" else "pass"
            meta.append({"isolate_id": iso, "pdt_accession": iso, "biosample_accession": bs, "assembly_accession": asm,
                         "genome_sha256": amrscan.sha256(genome), "genome_qc_status": qc, "genome_qc_flags": NA,
                         "mlst_st": st.replace("ST", ""), "lineage_id": st, "phylogroup": "B2",
                         "ncbi_snp_cluster": cluster, "project": project,
                         "eligible_population_analysis": "false" if qc == "fail" else "true"})
            groups.append({"isolate_id": iso, "lineage_component": st, "period": period,
                           "geography_group": "reference" if country == "USA" else "non_reference"})
        eps += endpoint_extra or []
        for row in meta:
            row.update((population_changes or {}).get(row["isolate_id"], {}))
        self.write_rows(self.processed / "candidate_isolates.tsv", cands)
        self.write_rows(self.processed / "isolate_drug_endpoints.tsv", eps)
        (self.processed / "cohort_manifest.json").write_text(json.dumps(
            {"cohort_id": "ecoli_research_v1", "candidate_isolate_count": len(cands)}))
        self.sums(self.processed, ("candidate_isolates.tsv", "isolate_drug_endpoints.tsv", "cohort_manifest.json"))
        self.write_rows(self.population / "population_metadata.tsv", meta)
        self.write_rows(self.population / "evaluation_groups.tsv", groups)
        (self.population / "population_manifest.json").write_text(json.dumps({
            "population_id": "test_population", "code_commit": "abc",
            "pr5_cohort_manifest_sha256": amrscan.sha256(self.processed / "cohort_manifest.json")}))
        self.sums(self.population, ("population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json"))

    def fake_caller(self, hits):
        def caller(args):
            rows = []
            for symbol, etype, method in hits.get(args.sample, []):
                rows.append(["NA", "c", "1", "10", "+", symbol, symbol, "core", etype, etype, "CLASS", "SUB", method,
                             "10", "10", "100.00", "100.00", "10", "WP_1", symbol, "NA", "NA", symbol])
            args.raw.write_text("\t".join(RAW_HEADER) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
            prov = {"schema_version": "1.0", "sample_id": args.sample, "tool": "AMRFinderPlus",
                    "caller_version": "4.2.7", "command": ["amrfinder"], "input_sha256": amrscan.sha256(args.input),
                    "database_files_sha256": {}, "database_id": self.database_id,
                    "thresholds": {"ident_min": args.identity, "coverage_min": args.coverage,
                                   "translation_table": 11, "report_all_equal": True,
                                   "organism": args.organism, "plus": args.plus},
                    "raw_sha256": amrscan.sha256(args.raw), "status": "success"}
            args.provenance.write_text(json.dumps(prov))
        return caller

    def run_calls(self, hits=None, **kw):
        return gb.call_cohort(self.population, self.cache, self.calls, self.database, workers=2,
                              caller=self.fake_caller(self.HITS if hits is None else hits),
                              version_probe=lambda: "4.2.7", **kw)

    def config(self, **changes):
        cfg = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/baseline_config.json").read_text())
        cfg["mapping"] = {"version": VERSION, "rules": str(self.rules), "metadata": str(self.meta)}
        cfg["expected_upstream"] = {"pr5_cohort_manifest_sha256": amrscan.sha256(self.processed / "cohort_manifest.json"),
                                    "pr6_population_manifest_sha256":
                                        amrscan.sha256(self.population / "population_manifest.json")}
        cfg["caller"]["database_id"] = self.database_id
        cfg["stratum_support"].update({"min_R": 1, "min_S": 1})
        cfg["mapping_audit"].update({"min_carriers": 1})
        for key, value in changes.items():
            cfg[key] = value
        path = self.root / "config.json"
        path.write_text(json.dumps(cfg))
        return path

    def build(self, name="out", config=None):
        return gb.build(self.processed, self.population, self.calls, self.database, config or self.config(),
                        self.root / name, workers=2, code_commit="test")

    def read(self, name, out="out"):
        return gb.read_table(self.root / out / name)


class CallerTests(CohortFixture):
    def test_calls_are_resumable_and_reused(self):
        self.build_inputs()
        first = self.run_calls()
        self.assertEqual(first["results"], {"success:called": 6})
        second = self.run_calls()
        self.assertEqual(second["results"], {"success:reused": 6})
        log = (self.calls / "call_log.jsonl").read_text().splitlines()
        self.assertEqual(len(log), 12)
        self.assertEqual(gb.call_times(self.calls), max(json.loads(l)["logged_at"] for l in log[:6]))

    def test_invalid_call_is_quarantined_and_rerun(self):
        self.build_inputs()
        self.run_calls(limit=1)
        raw = self.calls / "PDT000000001.1" / "amrfinder.tsv"
        raw.write_text(raw.read_text() + "tampered\n")
        status, detail = gb.call_status(raw.parent, "PDT000000001.1", "x", self.database_id)
        self.assertEqual(status, "provenance_mismatch")
        result = self.run_calls(limit=1)
        self.assertEqual(result["results"], {"success:called": 1})
        self.assertEqual(len(list((self.calls / ".quarantine").iterdir())), 1)

    def test_call_status_classes(self):
        self.build_inputs()
        self.run_calls(limit=1)
        d = self.calls / "PDT000000001.1"
        sha = json.loads((d / "provenance.json").read_text())["input_sha256"]
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id, "4.2.7"), ("success", ""))
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, "sha256:other")[0], "provenance_mismatch")
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id, "4.3.0")[0],
                         "provenance_mismatch")
        self.assertEqual(gb.call_status(d, "PDT000000009.1", sha, self.database_id)[0], "provenance_mismatch")
        self.assertEqual(gb.call_status(self.calls / "none", "x", sha, self.database_id)[0], "missing")
        (d / "harmonized.tsv").write_text((d / "harmonized.tsv").read_text() + "extra\n")
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id)[0], "harmonization_failed")
        (d / "harmonized.tsv").unlink()
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id)[0], "harmonization_failed")
        prov = json.loads((d / "provenance.json").read_text())
        prov["thresholds"]["plus"] = True
        (d / "provenance.json").write_text(json.dumps(prov))
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id)[0], "provenance_mismatch")
        prov["status"] = "failed"
        (d / "provenance.json").write_text(json.dumps(prov))
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id)[0], "caller_failed")
        (d / "provenance.json").write_text("{")
        self.assertEqual(gb.call_status(d, "PDT000000001.1", sha, self.database_id)[0], "provenance_mismatch")

    def test_caller_and_harmonizer_failures_are_recorded(self):
        self.build_inputs()

        def failing(args):
            raise subprocess.CalledProcessError(1, "amrfinder")
        result = gb.call_cohort(self.population, self.cache, self.calls, self.database, workers=1, caller=failing,
                                version_probe=lambda: "4.2.7", limit=1)
        self.assertEqual(result["results"], {"caller_failed:called": 1})

        def bad_harmonize(*a):
            raise ValueError("bad report")
        result = gb.call_cohort(self.population, self.cache, self.calls, self.database, workers=1,
                                caller=self.fake_caller(self.HITS), harmonizer=bad_harmonize,
                                version_probe=lambda: "4.2.7", limit=1)
        self.assertEqual(result["results"], {"harmonization_failed:called": 1})

    def test_missing_genome_is_not_called(self):
        self.build_inputs()
        (self.cache / "GCA_000000001.1" / "genome.fna").write_text(">c\nTTTT\n")
        result = self.run_calls(limit=1)
        self.assertEqual(result["results"], {"genome_unavailable:skipped": 1})

    def test_database_fingerprint_and_executable(self):
        empty = self.root / "empty"
        empty.mkdir()
        with self.assertRaisesRegex(ValueError, "empty"):
            gb.database_fingerprint(empty)
        self.build_inputs()
        with self.assertRaisesRegex(ValueError, "not found"):
            gb.call_cohort(self.population, self.cache, self.calls, self.database, executable="no-such-amrfinder")


class BuildTests(CohortFixture):
    def ready(self, **kw):
        self.build_inputs(**kw)
        self.run_calls()

    def endpoints(self, out="out"):
        return {(r["isolate_id"], r["antimicrobial"]): r for r in self.read("genotype_phenotype_endpoints.tsv", out)}

    def test_genotype_states_and_categories(self):
        self.ready()
        self.build()
        e = self.endpoints()
        a = e[("PDT000000001.1", "ampicillin")]
        self.assertEqual((a["genotype_status"], a["concordance_category"], a["determinant_ids"]),
                         ("determinant_positive", "R_DET_POS", "blaTEM-1"))
        self.assertEqual(e[("PDT000000001.1", "ceftriaxone")]["concordance_category"], "S_DET_NEG")
        self.assertEqual(e[("PDT000000001.1", "ciprofloxacin")]["determinant_types"], "point_mutation")
        self.assertEqual(e[("PDT000000001.1", "ciprofloxacin")]["fq_two_mechanisms_status"], "uncertain")
        u = e[("PDT000000003.1", "ceftriaxone")]
        self.assertEqual((u["genotype_status"], u["exclusion_reason"], u["uncertain_determinant_ids"]),
                         ("determinant_uncertain", "mapping_uncertain", "blaTEM"))
        self.assertEqual(e[("PDT000000003.1", "ampicillin")]["concordance_category"], "S_DET_NEG")
        # no-hit successful call: negative, not unavailable
        n = e[("PDT000000004.1", "ampicillin")]
        self.assertEqual((n["genotype_status"], n["concordance_category"]), ("determinant_negative", "R_DET_NEG"))
        sxt = e[("PDT000000002.1", "trimethoprim-sulfamethoxazole")]
        self.assertEqual(sxt["determinant_ids"], "dfrA17;sul1")
        self.assertEqual(e[("PDT000000001.1", "colistin")]["genotype_status"], "mapping_unavailable")
        self.assertEqual(e[("PDT000000003.1", "ciprofloxacin")]["exclusion_reason"], "non_binary_or_conflict")
        q = e[("PDT000000006.1", "ampicillin")]
        self.assertEqual((q["genotype_status"], q["exclusion_reason"]), ("genome_qc_failure", "population_ineligible"))

    def test_caller_and_provenance_failures_are_not_negative(self):
        self.build_inputs()
        self.run_calls()
        prov = self.calls / "PDT000000004.1" / "provenance.json"
        data = json.loads(prov.read_text())
        data["status"] = "failed"
        prov.write_text(json.dumps(data))
        prov = self.calls / "PDT000000002.1" / "provenance.json"
        data = json.loads(prov.read_text())
        data["input_sha256"] = "0" * 64
        prov.write_text(json.dumps(data))
        import shutil
        shutil.rmtree(self.calls / "PDT000000005.1")
        self.build()
        e = self.endpoints()
        self.assertEqual(e[("PDT000000004.1", "ampicillin")]["genotype_status"], "caller_failed")
        self.assertEqual(e[("PDT000000004.1", "ampicillin")]["concordance_category"], NA)
        self.assertEqual(e[("PDT000000004.1", "ampicillin")]["exclusion_reason"], "genotype_unavailable:caller_failed")
        self.assertEqual(e[("PDT000000002.1", "ampicillin")]["genotype_status"], "provenance_mismatch")
        self.assertEqual(e[("PDT000000005.1", "ampicillin")]["genotype_status"], "caller_missing")
        d = {r["antimicrobial"]: r for r in self.read("baseline_denominators.tsv")}["ampicillin"]
        self.assertEqual(d["excluded_genotype_unavailable"], "3")

    def test_irrelevant_and_non_amr_elements(self):
        self.ready()
        self.build()
        e = self.endpoints()
        # blaEC is does_not_predict; blaTEM-1 makes isolate 1 positive only for ampicillin
        self.assertEqual(e[("PDT000000001.1", "ceftriaxone")]["determinant_status"], "negative")
        self.assertNotIn("asr", e[("PDT000000003.1", "ampicillin")]["determinant_ids"])
        unexplained = self.read("unexplained_resistance.tsv")
        self.assertEqual([(r["isolate_id"], r["antimicrobial"]) for r in unexplained],
                         [("PDT000000004.1", "ampicillin")])
        s_pos = self.read("susceptible_with_determinant.tsv")
        self.assertEqual(s_pos, [])
        r_neg_other = [r for r in unexplained if r["other_amr_determinants"] != NA]
        self.assertEqual(r_neg_other, [])

    def test_unexplained_lists_other_determinants(self):
        hits = dict(self.HITS)
        hits["PDT000000004.1"] = [("tet(A)", "AMR", "EXACTX")]
        self.build_inputs()
        self.run_calls(hits)
        self.build()
        row = self.read("unexplained_resistance.tsv")[0]
        self.assertEqual(row["other_amr_determinants"], "tet(A)")
        self.assertEqual(self.read("unmapped_determinants.tsv")[0]["determinant"], "tet(A)")

    def test_metrics_and_denominators_reconcile(self):
        self.ready()
        self.build()
        metrics = {(r["antimicrobial"], r["rule_variant"]): r for r in self.read("baseline_metrics.tsv")}
        amp = metrics[("ampicillin", "primary")]
        # eligible: iso1 R/pos, iso2 R/pos (CTX-M-15), iso3 S/neg, iso4 R/neg, iso5 R/pos
        self.assertEqual((amp["TP"], amp["FN"], amp["FP"], amp["TN"]), ("3", "1", "0", "1"))
        self.assertEqual(amp["sensitivity"], "0.7500")
        self.assertEqual(amp["PPV"], "1.0000")
        cro = metrics[("ceftriaxone", "uncertain_as_positive")]
        self.assertEqual((cro["TP"], cro["N"]), ("2", "4"))
        self.assertEqual(metrics[("ceftriaxone", "primary")]["excluded_uncertain"], "1")
        self.assertIn(("ciprofloxacin", "fq_two_mechanisms"), metrics)
        self.assertNotIn(("ampicillin", "fq_two_mechanisms"), metrics)
        for d in self.read("baseline_denominators.tsv"):
            parts = sum(int(d[k]) for k in ("excluded_non_binary_or_conflict", "excluded_population_ineligible",
                                            "excluded_not_study_antimicrobial", "excluded_genotype_unavailable",
                                            "excluded_mapping_uncertain", "metric_eligible"))
            self.assertEqual(parts, int(d["total_endpoints"]))
            self.assertEqual(d["reconciled"], "true")

    def test_stratification(self):
        self.ready()
        self.build()
        lineage = {(r["dimension"], r["stratum"], r["antimicrobial"]): r for r in self.read("lineage_baseline.tsv")}
        st131 = lineage[("lineage", "ST131", "ampicillin")]
        self.assertEqual((st131["N"], st131["R"], st131["R_DET_POS"], st131["support_status"]),
                         ("2", "2", "2", "below_support: R 2, S 0; need R>=1 and S>=1"))
        self.assertEqual(st131["sensitivity"], NA)
        self.assertIn("only_R", st131["flags"])
        self.assertIn(("lineage", "UNKNOWN", "ampicillin"), lineage)
        country = {(r["dimension"], r["stratum"], r["antimicrobial"]): r for r in self.read("country_baseline.tsv")}
        self.assertIn(("country", "unknown", "ampicillin"), country)
        canada = country[("country", "Canada", "ampicillin")]
        self.assertEqual((canada["support_status"], canada["sensitivity"], canada["specificity"]),
                         ("reported", "0.0000", "1.0000"))
        temporal = {(r["dimension"], r["stratum"], r["antimicrobial"]) for r in self.read("temporal_baseline.tsv")}
        self.assertTrue({("period", "recent", "ampicillin"), ("period", "unknown", "ampicillin"),
                         ("collection_year", "unknown", "ampicillin")} <= temporal)
        source = {(r["stratum"], r["antimicrobial"]): r for r in self.read("source_baseline.tsv")}
        self.assertEqual(source[("unknown", "ampicillin")]["N"], "2")
        errors = {r["dimension"] for r in self.read("error_distribution.tsv")}
        self.assertEqual(errors, {"project", "snp_cluster"})
        # every eligible endpoint lands in exactly one stratum per dimension
        eligible = sum(r["endpoint_eligible"] == "true" for r in self.read("genotype_phenotype_endpoints.tsv"))
        self.assertEqual(sum(int(r["N"]) for r in self.read("source_baseline.tsv")), eligible)

    def test_mapping_audit(self):
        self.ready()
        self.build()
        audit = {(r["antimicrobial"], r["determinant"]): r for r in self.read("mapping_audit.tsv")}
        tem = audit[("ampicillin", "blaTEM-1")]
        self.assertEqual((tem["carriers"], tem["carriers_R"], tem["sole_support_R"]), ("2", "2", "2"))
        self.assertEqual(audit[("ampicillin", "blaEC")]["interpretation"], "does_not_predict")
        self.assertEqual(audit[("ampicillin", "blaEC")]["flag"], "HUMAN_REVIEW_REQUIRED: high R fraction for a "
                                                                 "non-predicting rule")
        self.assertEqual(audit[("ceftriaxone", "blaCTX-M-15")]["flag"], NA)

    def test_low_r_fraction_flag(self):
        self.PHENO = {k: dict(v) for k, v in self.PHENO.items()}
        self.PHENO["PDT000000001.1"]["ampicillin"] = "S"
        self.PHENO["PDT000000005.1"]["ampicillin"] = "S"
        self.ready()
        self.build()
        audit = {(r["antimicrobial"], r["determinant"]): r for r in self.read("mapping_audit.tsv")}
        self.assertTrue(audit[("ampicillin", "blaTEM-1")]["flag"].startswith("HUMAN_REVIEW_REQUIRED: low R"))
        self.assertEqual(len(self.read("susceptible_with_determinant.tsv")), 2)

    def test_manifest_summary_and_determinism(self):
        self.ready()
        first = self.build("a")
        self.build("b")
        self.assertEqual((self.root / "a" / "SHA256SUMS").read_text(), (self.root / "b" / "SHA256SUMS").read_text())
        m = first["manifest"]
        self.assertEqual(m["mapping_version"], VERSION)
        self.assertEqual(m["amrfinderplus_database_fingerprint"], self.database_id)
        self.assertEqual(m["cohort_isolate_count"], 6)
        self.assertEqual(m["caller_successful_isolates"], 5)
        self.assertEqual(m["generated_at"], gb.call_times(self.calls))
        s = first["summary"]
        self.assertEqual(s["caller_status"], {"genome_qc_failure": 1, "success": 5})
        self.assertEqual(s["unexplained_resistance_endpoints"], 1)
        rows = self.read("genotype_phenotype_endpoints.tsv", "a")
        self.assertEqual([(r["isolate_id"], r["antimicrobial"]) for r in rows],
                         sorted((r["isolate_id"], r["antimicrobial"]) for r in rows))
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.build("a")

    def test_join_validation(self):
        self.ready()
        bad = self.root / "bad"
        bad.mkdir()
        import shutil
        for name in ("candidate_isolates.tsv", "cohort_manifest.json"):
            shutil.copy(self.processed / name, bad / name)
        rows = gb.read_table(self.processed / "isolate_drug_endpoints.tsv")
        for field, label in (("assembly_accession", "assembly_accession mismatch"),
                             ("biosample_accession", "biosample_accession mismatch")):
            changed = [dict(r) for r in rows]
            changed[0][field] = "X"
            (bad / "isolate_drug_endpoints.tsv").unlink(missing_ok=True)
            self.write_rows(bad / "isolate_drug_endpoints.tsv", changed)
            self.sums(bad, ("candidate_isolates.tsv", "isolate_drug_endpoints.tsv", "cohort_manifest.json"))
            with self.assertRaisesRegex(ValueError, label):
                gb.build(bad, self.population, self.calls, self.database, self.config(), self.root / "o", 1)
        dup = rows + [rows[0]]
        (bad / "isolate_drug_endpoints.tsv").unlink()
        self.write_rows(bad / "isolate_drug_endpoints.tsv", dup)
        self.sums(bad, ("candidate_isolates.tsv", "isolate_drug_endpoints.tsv", "cohort_manifest.json"))
        with self.assertRaisesRegex(ValueError, "duplicate PR #5 endpoint"):
            gb.build(bad, self.population, self.calls, self.database, self.config(), self.root / "o", 1)

    def test_missing_phenotype_and_genotype(self):
        extra = [{"cohort_id": "ecoli_research_v1", "isolate_id": "PDT000000004.1",
                  "pdt_accession": "PDT000000004.1", "biosample_accession": "SAMN00000004",
                  "assembly_accession": "GCA_000000004.1", "antimicrobial": "trimethoprim-sulfamethoxazole",
                  "phenotype_endpoint": "NA", "phenotype_status": "missing", "eligible_binary_rs": "false",
                  "country": "Canada", "collection_year": "2023", "source_category": "unknown"}]
        self.build_inputs(endpoint_extra=extra, population_changes={"PDT000000005.1": {"genome_sha256": NA}})
        self.run_calls()
        self.build()
        e = self.endpoints()
        self.assertEqual(e[("PDT000000004.1", "trimethoprim-sulfamethoxazole")]["exclusion_reason"],
                         "non_binary_or_conflict")
        self.assertEqual(e[("PDT000000005.1", "ampicillin")]["genotype_status"], "genome_unavailable")

    def test_upstream_hash_checks(self):
        self.ready()
        cfg = self.config()
        data = json.loads(cfg.read_text())
        data["expected_upstream"]["pr5_cohort_manifest_sha256"] = "0" * 64
        cfg.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "PR #5 cohort manifest hash"):
            self.build(config=cfg)
        data = json.loads(self.config().read_text())
        data["expected_upstream"]["pr6_population_manifest_sha256"] = "0" * 64
        cfg.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "PR #6 population manifest hash"):
            self.build(config=cfg)
        data = json.loads(self.config().read_text())
        data["caller"]["database_id"] = "sha256:other"
        cfg.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "database fingerprint"):
            self.build(config=cfg)
        manifest = json.loads((self.population / "population_manifest.json").read_text())
        manifest["pr5_cohort_manifest_sha256"] = "0" * 64
        (self.population / "population_manifest.json").write_text(json.dumps(manifest))
        self.sums(self.population, ("population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json"))
        with self.assertRaisesRegex(ValueError, "different PR #5 cohort"):
            gb.load_population(self.population, None, amrscan.sha256(self.processed / "cohort_manifest.json"))
        (self.population / "evaluation_groups.tsv").write_text("isolate_id\n")
        with self.assertRaisesRegex(ValueError, "missing or changed"):
            gb.load_population(self.population, None, "x")

    def test_population_lists_must_agree(self):
        self.build_inputs()
        groups = gb.read_table(self.population / "evaluation_groups.tsv")[:-1]
        (self.population / "evaluation_groups.tsv").unlink()
        self.write_rows(self.population / "evaluation_groups.tsv", groups)
        self.sums(self.population, ("population_metadata.tsv", "evaluation_groups.tsv", "population_manifest.json"))
        with self.assertRaisesRegex(ValueError, "different isolates"):
            gb.load_population(self.population, None, amrscan.sha256(self.processed / "cohort_manifest.json"))

    def test_config_validation(self):
        self.build_inputs()
        path = self.root / "c.json"
        path.write_text(json.dumps({"mapping": {}}))
        with self.assertRaisesRegex(ValueError, "lacks"):
            gb.load_config(path)
        data = json.loads(self.config().read_text())
        data["stratum_support"]["decision"] = "accepted"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "decision record"):
            gb.load_config(path)
        data = json.loads(self.config().read_text())
        data["phenotype_corrections"].pop("decision_record")
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "phenotype_corrections needs"):
            gb.load_config(path)

    def write_corrections(self, rows):
        path = self.root / "corrections.tsv"
        self.write_rows(path, [dict(zip(gb.CORRECTION_FIELDS, r)) for r in rows])
        return path

    def test_project_phenotype_correction_swaps_binary_endpoints(self):
        self.build_inputs()
        self.run_calls()
        table = self.write_corrections([("PC1", "PRJNA2", "swap_R_S", "all_binary_endpoints", "test", "0004")])
        cfg = self.config(phenotype_corrections={"table": str(table), "sha256": amrscan.sha256(table),
                                                 "decision_record": "0004"})
        result = self.build(config=cfg)
        rows = {(r["isolate_id"], r["antimicrobial"]): r for r in self.read("genotype_phenotype_endpoints.tsv")}
        self.assertEqual(rows[("PDT000000003.1", "ampicillin")]["phenotype"], "R")
        self.assertEqual(rows[("PDT000000003.1", "ampicillin")]["phenotype_pr5"], "S")
        self.assertEqual(rows[("PDT000000003.1", "ampicillin")]["phenotype_correction"], "PC1")
        self.assertEqual(rows[("PDT000000004.1", "ciprofloxacin")]["phenotype"], "R")
        # conflicts are not binary and stay as they are; other projects are untouched
        self.assertEqual(rows[("PDT000000003.1", "ciprofloxacin")]["phenotype"], "CONFLICT")
        self.assertEqual(rows[("PDT000000003.1", "ciprofloxacin")]["phenotype_correction"], NA)
        self.assertEqual(rows[("PDT000000001.1", "ampicillin")]["phenotype"], "R")
        self.assertEqual(rows[("PDT000000001.1", "ampicillin")]["phenotype_correction"], NA)
        self.assertEqual(result["summary"]["phenotype_corrections"], {"PC1": 5})
        self.assertEqual(result["manifest"]["phenotype_corrected_endpoints"], 5)
        self.assertEqual(result["manifest"]["phenotype_corrections_sha256"], amrscan.sha256(table))

    def test_phenotype_correction_table_checks(self):
        good = ("PC1", "PRJNA2", "swap_R_S", "all_binary_endpoints", "e", "0004")
        table = self.write_corrections([good])
        with self.assertRaisesRegex(ValueError, "hash differs"):
            gb.load_corrections(table, "0" * 64)
        self.assertEqual(gb.load_corrections(table, amrscan.sha256(table)), {"PRJNA2": "PC1"})
        table = self.write_corrections([("PC1", "PRJNA2", "set_R", "all_binary_endpoints", "e", "0004")])
        with self.assertRaisesRegex(ValueError, "unsupported"):
            gb.load_corrections(table, amrscan.sha256(table))
        table = self.write_corrections([good, ("PC2", "PRJNA2", "swap_R_S", "all_binary_endpoints", "e", "0004")])
        with self.assertRaisesRegex(ValueError, "corrected twice"):
            gb.load_corrections(table, amrscan.sha256(table))

    def test_empty_cohort(self):
        self.build_inputs()
        rows = gb.read_table(self.processed / "isolate_drug_endpoints.tsv")
        (self.processed / "isolate_drug_endpoints.tsv").write_text("\t".join(rows[0]) + "\n")
        self.sums(self.processed, ("candidate_isolates.tsv", "isolate_drug_endpoints.tsv", "cohort_manifest.json"))
        self.run_calls()
        result = self.build()
        self.assertEqual(result["summary"]["phenotype_endpoints"], 0)
        self.assertEqual(self.read("genotype_phenotype_endpoints.tsv"), [])

    def test_unexpected_caller_status(self):
        self.ready()
        prov = self.calls / "PDT000000001.1" / "provenance.json"
        data = json.loads(prov.read_text())
        data["status"] = "upstream_fixture"
        prov.write_text(json.dumps(data))
        self.build()
        self.assertEqual(self.endpoints()[("PDT000000001.1", "ampicillin")]["genotype_status"], "caller_failed")

    def test_command_line(self):
        self.build_inputs()
        with patch.object(gb, "call_cohort", return_value={"isolates": 0}) as call:
            self.assertEqual(gb.main(["call", "--population", str(self.population), "--cache", str(self.cache),
                                      "--calls", str(self.calls), "--database", str(self.database)]), 0)
            call.assert_called_once()
        self.assertEqual(gb.main(["call", "--population", "x", "--cache", "x", "--calls", "x", "--database", "x",
                                  "--threads", "0"]), 2)
        self.assertEqual(gb.main(["build", "--processed", "x", "--population", "x", "--calls", "x",
                                  "--database", "x", "--outdir", "x", "--workers", "0"]), 2)
        self.run_calls()
        args = ["build", "--processed", str(self.processed), "--population", str(self.population), "--calls",
                str(self.calls), "--database", str(self.database), "--config", str(self.config()),
                "--outdir", str(self.root / "cli")]
        self.assertEqual(gb.main(args), 0)
        self.assertEqual(gb.main(args), 2)
        self.assertIsNotNone(gb.git_commit())


qc_spec = importlib.util.spec_from_file_location("baseline_qc", ROOT / "scripts/ecoli_research_v1/baseline_qc.py")
qc = importlib.util.module_from_spec(qc_spec)
qc_spec.loader.exec_module(qc)


class LabelQcTests(unittest.TestCase):
    @staticmethod
    def obs(iso, drug, phenotype, mic, n=[0]):
        n[0] += 1
        return {"record_id": f"r{n[0]}", "isolate_id": iso, "antimicrobial": drug, "phenotype": phenotype,
                "mic": mic, "mic_operator": "=", "breakpoint_standard": "EUCAST"}

    def test_inverted_project_is_flagged_and_cleared_by_correction(self):
        project_of = {"a": "P1", "b": "P2"}
        observations = []
        for mic in ("0.5", "1", "2"):
            observations += [self.obs("a", "ampicillin", "R", mic), self.obs("b", "ampicillin", "S", mic)]
        for mic in ("32", "64", "64"):
            observations += [self.obs("a", "ampicillin", "S", mic), self.obs("b", "ampicillin", "R", mic)]
        observations.append(self.obs("a", "ampicillin", "S", "not_a_number"))
        observations.append(self.obs("a", "ampicillin", "S", "0.5"))  # stays out of order after the swap
        rows, tested, inverted, flagged = qc.label_check(observations, project_of)
        self.assertEqual(flagged, ["P1"])
        fixed = qc.corrected_observations(observations, project_of, {"P1": "PC1"})
        self.assertEqual(qc.label_check(fixed, project_of)[3], [])
        order, totals = qc.out_of_order(fixed, project_of, {"P1": "PC1"})
        self.assertEqual([(r["corrected_phenotype"], r["mic"]) for r in order], [("R", "=0.5")])
        self.assertEqual(dict(totals), {"observations": 7, "out_of_order": 1})


if __name__ == "__main__":
    unittest.main()
