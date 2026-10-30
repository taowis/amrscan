"""Offline research-cohort tests using small synthetic NCBI table snapshots."""

import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
spec = importlib.util.spec_from_file_location("ecoli_cohort", ROOT / "bin/ecoli_cohort.py")
cohort = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cohort)
ast = cohort.ast_dataset


def raw_ast(pdt, drug, phenotype, suffix, **changes):
    row = {name: "" for name in ast.AST_FIELDS + sorted(ast.AST_OPTIONAL)}
    row.update({"id": f"{suffix}_{pdt}", "target_acc": pdt, "biosample_acc": "SAMN00000001",
                "scientific_name": "Escherichia coli", "antibiotic": drug, "phenotype": phenotype,
                "measurement_sign": "==", "mic": "1", "standard": "CLSI", "platform": "Test",
                "host": "Homo sapiens", "collection_date": "2020-04", "geo_loc_name": "USA: Oregon",
                "isolation_source": "urine", "epi_type": "clinical", "bioproject_acc": "PRJNA000001",
                "taxgroup_name": "E.coli and Shigella"})
    row.update(changes)
    return row


def raw_meta(pdt, number=1, **changes):
    row = {name: "" for name in ast.META_FIELDS + sorted(ast.META_OPTIONAL)}
    row.update({"target_acc": pdt, "biosample_acc": "SAMN00000001", "asm_acc": "GCA_000000001.1",
                "Run": "", "scientific_name": "Escherichia coli", "collection_date": "2020-04",
                "geo_loc_name": "USA: Oregon", "isolation_source": "urine", "host": "Homo sapiens",
                "epi_type": "clinical", "bioproject_acc": "PRJNA000001", "species_taxid": "562",
                "isolate_identifiers": f"isolate-{number}", "strain": f"strain-{number}",
                "taxid": "562", "erd_group": "PDS000001.1", "bioproject_center": "Test Center"})
    row.update(changes)
    return row


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ast_file = self.root / "ast.tsv"
        self.meta_file = self.root / "metadata.tsv"
        self.snapshot = self.root / "snapshot"
        self.config = ROOT / "benchmark/cohorts/ecoli_research_v1/cohort_config.json"

    def write_tables(self, ast_rows, meta_rows):
        for path, fields, rows in ((self.ast_file, ast.AST_FIELDS + sorted(ast.AST_OPTIONAL), ast_rows),
                                   (self.meta_file, ast.META_FIELDS + sorted(ast.META_OPTIONAL), meta_rows)):
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

    def archive(self, ast_rows, meta_rows):
        self.write_tables(ast_rows, meta_rows)
        ast.archive(self.ast_file, [self.meta_file], self.snapshot,
                    "https://fixture.invalid/ncbi-export", "2026-10-11T00:00:00Z")

    def read_output(self, out, name):
        with (out / name).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def base_population(self):
        ast_rows = [
            raw_ast("PDT000000001.1", "CIP", "resistant", "001", platform="A"),
            raw_ast("PDT000000001.1", "CIP", "resistant", "002", platform="B"),
            raw_ast("PDT000000001.1", "AMP", "resistant", "003"),
            raw_ast("PDT000000001.1", "AMP", "susceptible", "004"),
            raw_ast("PDT000000001.1", "TMP-SMX", "susceptible", "005"),
            raw_ast("PDT000000001.1", "TMP-SMX", "susceptible", "006"),
            raw_ast("PDT000000002.1", "CIP", "susceptible", "007", host="", geo_loc_name="", collection_date="", isolation_source=""),
            raw_ast("PDT000000002.1", "CIP", "intermediate", "008", host="", geo_loc_name="", collection_date="", isolation_source=""),
            raw_ast("PDT000000002.1", "nitrofurantoin", "susceptible", "015", host="", geo_loc_name="", collection_date="", isolation_source=""),
            raw_ast("PDT000000003.1", "CIP", "ND", "009"),
            raw_ast("PDT000000003.1", "moxifloxacin", "resistant", "010"),
            raw_ast("PDT000000004.1", "CIP", "not-collected-by-policy", "011"),
            raw_ast("PDT000000005.1", "CAZ", "susceptible", "012", host="", geo_loc_name="", collection_date="", isolation_source="ground beef"),
            raw_ast("PDT000000006.1", "CIP", "resistant", "013", scientific_name="Salmonella enterica"),
            raw_ast("", "CIP", "resistant", "014"),
        ]
        meta_rows = [raw_meta(f"PDT00000000{i}.1", i) for i in range(1, 6)]
        # Consecutive pipeline metadata rows are preserved but do not duplicate AST observations.
        meta_rows.append(raw_meta("PDT000000001.1", 1, erd_group="PDS000002.1"))
        # A valid-identity isolate with missing country/source but a usable collection year.
        meta_rows[1].update({"geo_loc_name": "", "host": "", "collection_date": "", "isolation_source": ""})
        # Missing assembly is analysis-specific; phenotype observations remain in the endpoint table.
        meta_rows[4]["asm_acc"] = ""
        return ast_rows, meta_rows

    def test_endpoints_duplicates_conflicts_and_nonbinary_categories(self):
        rows, metadata = self.base_population()
        self.archive(rows, metadata)
        out = self.root / "out"
        flow = cohort.build(self.snapshot, out, self.config)
        endpoints = {(r["isolate_id"], r["antimicrobial"]): r
                     for r in self.read_output(out, "isolate_drug_endpoints.tsv")}
        self.assertEqual(endpoints[("PDT000000001.1", "ciprofloxacin")]["phenotype_endpoint"], "R")
        self.assertEqual(endpoints[("PDT000000001.1", "ciprofloxacin")]["duplicate_status"], "repeated_consistent")
        self.assertEqual(endpoints[("PDT000000001.1", "ampicillin")]["phenotype_endpoint"], "CONFLICT")
        self.assertEqual(endpoints[("PDT000000001.1", "ampicillin")]["conflict_status"], "discordant_rs")
        self.assertEqual(endpoints[("PDT000000001.1", "trimethoprim-sulfamethoxazole")]["phenotype_endpoint"], "S")
        self.assertEqual(endpoints[("PDT000000002.1", "ciprofloxacin")]["phenotype_endpoint"], "NON_BINARY_OR_CONFLICT")
        self.assertEqual(endpoints[("PDT000000003.1", "ciprofloxacin")]["phenotype_endpoint"], "ND")
        self.assertEqual(endpoints[("PDT000000004.1", "ciprofloxacin")]["phenotype_endpoint"], "UNKNOWN")
        self.assertEqual(endpoints[("PDT000000003.1", "moxifloxacin")]["eligible_binary_rs"], "false")
        self.assertEqual(flow["units"]["raw_ast_observations"]["count"], 15)
        self.assertEqual(flow["units"]["ecoli_ast_observations"]["count"], 14)
        self.assertTrue(all(flow["reconciliation"].values()))

    def test_analysis_specific_missingness_and_source_preservation(self):
        rows, metadata = self.base_population()
        self.archive(rows, metadata)
        out = self.root / "out"
        cohort.build(self.snapshot, out, self.config)
        endpoints = {(r["isolate_id"], r["antimicrobial"]): r
                     for r in self.read_output(out, "isolate_drug_endpoints.tsv")}
        endpoint = endpoints[("PDT000000002.1", "ciprofloxacin")]
        self.assertEqual(endpoint["phenotype_endpoint"], "NON_BINARY_OR_CONFLICT")
        endpoint = endpoints[("PDT000000002.1", "nitrofurantoin")]
        self.assertEqual(endpoint["eligible_overall"], "true")
        self.assertEqual(endpoint["eligible_country_analysis"], "false")
        self.assertEqual(endpoint["eligible_temporal_analysis"], "false")
        self.assertEqual(endpoint["eligible_source_analysis"], "false")
        self.assertEqual(endpoints[("PDT000000005.1", "ceftazidime")]["eligible_overall"], "false")
        observations = self.read_output(out, "ast_observations.tsv")
        self.assertEqual(len(observations), 15)
        raw = json.loads(observations[0]["raw_record_json"])
        self.assertIn("antibiotic", raw)
        candidates = self.read_output(out, "candidate_isolates.tsv")
        p1 = next(r for r in candidates if r["isolate_id"] == "PDT000000001.1")
        self.assertEqual(p1["ast_record_count"], "6")
        self.assertEqual(p1["metadata_row_count"], "2")
        self.assertEqual(json.loads(p1["snp_clusters_json"]), ["PDS000001.1", "PDS000002.1"])

    def test_summaries_missingness_readiness_and_flow_reconcile(self):
        rows, metadata = self.base_population()
        self.archive(rows, metadata)
        out = self.root / "out"
        cohort.build(self.snapshot, out, self.config)
        flow = json.loads((out / "cohort_flow.json").read_text())
        summary = {r["antimicrobial"]: r for r in self.read_output(out, "antimicrobial_summary.tsv")}
        self.assertEqual(summary["ciprofloxacin"]["n_R"], "1")
        self.assertEqual(summary["ciprofloxacin"]["n_S"], "0")
        self.assertEqual(summary["ciprofloxacin"]["overall_ready"], "HUMAN_REVIEW_REQUIRED")
        self.assertEqual(summary["ampicillin"]["n_conflict"], "1")
        country = self.read_output(out, "country_summary.tsv")
        self.assertIn("NA", {r["country"] for r in country})
        year = self.read_output(out, "year_summary.tsv")
        self.assertIn("NA", {r["collection_year"] for r in year})
        source = self.read_output(out, "source_summary.tsv")
        self.assertIn("unknown", {r["source_category"] for r in source})
        missing = self.read_output(out, "missingness_summary.tsv")
        self.assertIn(("E. coli candidate isolates", "country"), {(r["unit"], r["field"]) for r in missing})
        self.assertEqual(flow["units"]["isolate_drug_endpoints"]["count"],
                         len(self.read_output(out, "isolate_drug_endpoints.tsv")))
        exclusions = self.read_output(out, "exclusions.tsv")
        self.assertTrue(any(r["analysis"] == "country" and "missing_country" in r["reason"] for r in exclusions))
        self.assertTrue((out / "SHA256SUMS").is_file())

    def test_rebuild_outputs_are_deterministic(self):
        rows, metadata = self.base_population()
        self.archive(rows, metadata)
        a, b = self.root / "a", self.root / "b"
        cohort.build(self.snapshot, a, self.config)
        cohort.build(self.snapshot, b, self.config)
        for name in sorted(p.name for p in a.iterdir()):
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)

    def test_duplicate_source_ast_ids_fail_closed(self):
        rows, metadata = self.base_population()
        rows[1]["id"] = rows[0]["id"]
        self.archive(rows, metadata)
        with self.assertRaisesRegex(ValueError, "duplicate source AST record id"):
            cohort.build(self.snapshot, self.root / "out", self.config)

    def test_conflicting_identity_is_not_overall_eligible(self):
        rows = [raw_ast("PDT000000001.1", "CIP", "resistant", "001")]
        metadata = [raw_meta("PDT000000001.1"), raw_meta("PDT000000001.1", 2, asm_acc="GCA_000000002.1")]
        self.archive(rows, metadata)
        out = self.root / "out"
        cohort.build(self.snapshot, out, self.config)
        endpoint = self.read_output(out, "isolate_drug_endpoints.tsv")[0]
        self.assertEqual(endpoint["eligible_binary_rs"], "true")
        self.assertEqual(endpoint["eligible_overall"], "false")
        self.assertIn("identity_incomplete_or_conflicting", endpoint["exclusion_reason"])

    def test_header_only_snapshot_has_zero_counts_and_valid_outputs(self):
        self.archive([], [])
        out = self.root / "out"
        flow = cohort.build(self.snapshot, out, self.config)
        self.assertEqual(flow["units"]["raw_ast_observations"]["count"], 0)
        self.assertEqual(flow["units"]["unique_ecoli_isolates"]["count"], 0)
        self.assertEqual(len(self.read_output(out, "candidate_isolates.tsv")), 0)
        self.assertEqual(len(self.read_output(out, "isolate_drug_endpoints.tsv")), 0)

    def test_bad_input_tables_fail_at_archive_boundary(self):
        self.ast_file.write_text("id\ttarget_acc\n1\tPDT000000001.1\n", encoding="utf-8")
        self.meta_file.write_text("target_acc\tbiosample_acc\nPDT000000001.1\tSAMN00000001\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "schema mismatch"):
            ast.archive(self.ast_file, [self.meta_file], self.snapshot, "fixture:bad", "2026-10-11T00:00:00Z")
        self.ast_file.write_text("id\tid\n1\t2\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate columns"):
            ast.archive(self.ast_file, [self.meta_file], self.snapshot, "fixture:bad", "2026-10-11T00:00:00Z")

    def test_empty_metadata_is_reported_without_fabricating_link(self):
        rows = [raw_ast("PDT000000001.1", "CIP", "susceptible", "001")]
        self.archive(rows, [])
        out = self.root / "out"
        cohort.build(self.snapshot, out, self.config)
        endpoint = self.read_output(out, "isolate_drug_endpoints.tsv")[0]
        self.assertEqual(endpoint["biosample_accession"], "SAMN00000001")
        self.assertEqual(endpoint["assembly_accession"], "NA")
        self.assertEqual(endpoint["assembly_available"], "false")
        self.assertEqual(endpoint["eligible_overall"], "false")

    def test_cli_entrypoint(self):
        rows, metadata = self.base_population()
        self.archive(rows, metadata)
        result = __import__("subprocess").run(
            [sys.executable, str(ROOT / "bin/ecoli_cohort.py"), "build", "--snapshot", str(self.snapshot),
             "--config", str(self.config), "--outdir", str(self.root / "cli-out")],
            capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("cohort_id", json.loads(result.stdout))

    def test_candidate_inventory_deduplicates_pipeline_rows_without_calling_them_ast(self):
        source = self.root / "candidate.tsv"
        fields = ["#target_acc", "biosample_acc", "asm_acc", "scientific_name", "strain", "erd_group", "bioproject_acc"]
        rows = [
            ["PDT000000001.1", "SAMN00000001", "GCA_000000001.1", "Escherichia coli", "A", "PDS000001.1", "PRJNA1"],
            ["PDT000000001.1", "SAMN00000001", "GCA_000000001.1", "Escherichia coli", "A", "PDS000002.1", "PRJNA1"],
            ["PDT000000002.1", "SAMN00000002", "GCA_000000002.1", "Escherichia coli", "", "", "PRJNA2"],
        ]
        with source.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(rows)
        out = self.root / "inventory"
        summary = cohort.inventory(source, out, "https://fixture.invalid", "E. coli + AST + GCA",
                                   "2026-10-11T00:00:00Z")
        self.assertEqual(summary["source_row_count"], 3)
        self.assertEqual(summary["unique_versioned_pdt_count"], 2)
        self.assertEqual(summary["candidate_rows_per_pdt"], {"1": 1, "2": 1})
        self.assertIsNone(summary["ast_observation_count"])
        self.assertIn("full AST", summary["endpoint_census_status"])
        self.assertEqual(len(self.read_output(out, "source_candidates.tsv")), 2)
        self.assertTrue((out / "SHA256SUMS").exists())

    def test_candidate_inventory_keeps_optional_taxid(self):
        source = self.root / "candidate.tsv"
        fields = ["#target_acc", "biosample_acc", "asm_acc", "scientific_name", "strain", "erd_group",
                  "bioproject_acc", "taxid"]
        rows = [["PDT000000001.1", "SAMN00000001", "GCA_000000001.1", "Escherichia coli", "A", "PDS000001.1",
                 "PRJNA1", "562"]] * 2
        with source.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(rows)
        out = self.root / "inventory-taxid"
        cohort.inventory(source, out, "https://fixture.invalid", "query", "2026-10-11T00:00:00Z")
        self.assertEqual(self.read_output(out, "source_candidates.tsv")[0]["taxid"], "562")

    def test_antimicrobial_country_count_ignores_placeholder_countries(self):
        config = cohort.validate_config(self.config)
        def endpoint(isolate, country):
            return {"antimicrobial": "ciprofloxacin", "isolate_id": isolate, "eligible_binary_rs": "true",
                    "phenotype_endpoint": "R", "country": country, "collection_year": "2020",
                    "source_category": "human", "conflict_status": "none", "assembly_available": "true"}
        rows = cohort.summarize_antimicrobials(
            [endpoint("PDT000000001.1", "USA"), endpoint("PDT000000002.1", "unknown"),
             endpoint("PDT000000003.1", "NA")], config)
        self.assertEqual(rows[0]["n_countries"], "1")
        self.assertEqual(rows[0]["country_available_n"], "1")

    def test_candidate_inventory_rejects_naive_timestamp_and_bad_rows(self):
        source = self.root / "candidate.tsv"
        fields = ["#target_acc", "biosample_acc", "asm_acc", "scientific_name", "strain", "erd_group", "bioproject_acc"]
        source.write_text("\t".join(fields) + "\nPDT000000001.1\tSAMN00000001\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "wrong column count"):
            cohort.inventory(source, self.root / "bad", "fixture:", "query", "2026-10-11T00:00:00Z")


if __name__ == "__main__":
    unittest.main()
