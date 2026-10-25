"""Offline AST parser, snapshot, cohort and transport regression tests."""

import copy
import csv
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs, urlparse, unquote

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ast_dataset", ROOT / "bin/ast_dataset.py")
ast = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ast)
FIXTURES = ROOT / "tests/fixtures"


class DatasetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.raw = self.root / "raw"
        self.out = self.root / "processed"
        self.a = ast.read_table(FIXTURES / "ast.synthetic.tsv", "ast")[0]
        self.m = ast.read_table(FIXTURES / "isolates.synthetic.tsv", "metadata")[0]
        self.policy, self.drugs = ast.load_policy()

    def snapshot(self, rows=None, metadata=None, extras=()):
        a, m = self.root / "ast.tsv", self.root / "metadata.tsv"
        ast.write_tsv(a, ast.AST_FIELDS + list(extras), [self.a] if rows is None else rows)
        ast.write_tsv(m, ast.META_FIELDS, [self.m] if metadata is None else metadata)
        ast.archive(a, [m], self.raw, "fixture:synthetic", "2026-10-10T00:00:00Z")
        return self.raw

    def normalized(self, changes=None, metadata=None):
        raw = {**self.a, **(changes or {})}
        entry = ast.source_entry(FIXTURES / "ast.synthetic.tsv", "ast", "fixture:synthetic", "2026-10-10T00:00:00Z", "ast.tsv")
        meta = [{"record": self.m}] if metadata is None else [{"record": r} for r in metadata]
        return ast.normalize(raw, entry, 1, meta, self.policy, self.drugs)

    def test_valid_record_and_full_provenance(self):
        self.snapshot()
        summary = ast.build(self.raw, self.out)
        with (self.out / "ast.normalized.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(row), ast.SCHEMA)
        self.assertEqual(row["antimicrobial"], "ciprofloxacin")
        self.assertEqual(row["phenotype"], "R")
        self.assertEqual((row["mic"], row["mic_operator"], row["mic_units"]), ("0.25", "<=", "mg/L"))
        self.assertEqual(row["collection_year"], "2020")
        self.assertEqual(json.loads(row["raw_record_json"]), self.a)
        self.assertEqual(json.loads(row["metadata_records_json"])[0]["record"], self.m)
        self.assertEqual((row["source_row"], row["source_url"]), ("1", "fixture:synthetic"))
        self.assertEqual(row["source_sha256"], ast.digest(self.raw / row["source_file"]))
        self.assertEqual(summary["eligible_isolate_drug_pairs"], 1)
        self.assertEqual(summary["ecoli_record_counts"]["phenotype"], {"R": 1})
        self.assertEqual((ROOT / "schemas/ast.tsv").read_text().rstrip().split("\t"), ast.SCHEMA)
        self.assertEqual(ast.validate(self.raw, self.out), summary)

    def test_missing_metadata_remains_missing(self):
        row = self.normalized({"target_acc": "", "phenotype": "", "host": "", "isolation_source": "",
                               "collection_date": "", "geo_loc_name": "", "mic": ""}, [])
        for field in ("isolate_id", "phenotype", "assembly_accession", "sra_accession", "country", "collection_year", "mic"):
            self.assertEqual(row[field], "NA")
        self.assertEqual(row["source_category"], "unknown")
        self.assertTrue({"missing_isolate_id", "missing_genome_link", "missing_metadata", "missing_source",
                         "missing_phenotype", "missing_country", "missing_collection_year"} <= set(row["qc_flags"].split(";")))
        row = self.normalized(metadata=[{**self.m, "asm_acc": ""}])
        self.assertIn("missing_assembly", row["qc_flags"])
        self.assertNotIn("missing_genome_link", row["qc_flags"])

    def test_all_phenotype_mappings_remain_nonbinary(self):
        expected = {"resistant": "R", "sensitive": "S", "I": "I", "intermediate": "I", "SSD": "SDD",
                    "SDD": "SDD", "susceptible-dose dependent": "SDD", "NS": "NS", "nonsusceptible": "NS",
                    "HLAR": "HLAR", "not defined": "ND", "N": "ND", "ND": "ND", "nonsense": "UNKNOWN", "": "NA"}
        for original, normalized in expected.items():
            with self.subTest(original=original):
                row = self.normalized({"phenotype": original})
                self.assertEqual((row["phenotype"], row["phenotype_original"]), (normalized, original))
        self.assertIn("unknown_phenotype", self.normalized({"phenotype": "nonsense"})["qc_flags"])

    def test_drug_mapping_is_exact_and_unknown_preserved(self):
        for original, expected in [("CIP", "ciprofloxacin"), (" CIP ", "ciprofloxacin"),
                                   ("TMP-SMX", "trimethoprim-sulfamethoxazole"), ("CRO", "ceftriaxone")]:
            row = self.normalized({"antibiotic": original})
            self.assertEqual((row["antimicrobial"], row["antimicrobial_original"]), (expected, original))
        row = self.normalized({"antibiotic": "ciprofloxacni"})
        self.assertEqual(row["antimicrobial"], "ciprofloxacni")
        self.assertIn("unknown_antimicrobial", row["qc_flags"])

    def test_numeric_mic_comparators_and_rejections(self):
        for raw, sign, value, operator in [("2", "", "2", "="), ("<=0.25", "", "0.25", "<="),
                                          (">64", ">", "64", ">"), (".5", "==", "0.5", "="),
                                          ("2.0", ">=", "2", ">="), ("==2", "=", "2", "=")]:
            self.assertEqual(ast.parse_mic(raw, sign), (value, operator, False))
        for raw, sign in [("NaN", ""), ("inf", ""), ("-2", ""), ("0", ""), ("2 mg/L", ""),
                          ("32/4", ""), ("1..2", ""), ("2", "approximately"), ("<=2", ">")]:
            self.assertEqual(ast.parse_mic(raw, sign), ("NA", "NA", True))
        self.assertEqual(ast.parse_mic("", ">"), ("NA", "NA", False))
        self.assertIn("invalid_mic", self.normalized({"mic": "oops"})["qc_flags"])

    def test_combination_mic_and_disk_diffusion_are_preserved(self):
        row = self.normalized({"mic": "32", "mic_secondary": "4", "measurement_sign": ">"})
        self.assertEqual((row["mic"], row["mic_secondary"], row["mic_operator"]), ("32", "4", ">"))
        for changes in [{"mic": "", "mic_secondary": "4"}, {"mic_secondary": "oops"}]:
            self.assertIn("invalid_mic", self.normalized(changes)["qc_flags"])
        row = self.normalized({"mic": "", "disk_diffusion": "18", "measurement_sign": "<="})
        self.assertEqual((row["mic"], row["mic_units"], row["disk_diffusion_raw"]), ("NA", "NA", "18"))
        self.assertEqual(row["measurement_sign_original"], "<=")
        self.assertEqual(row["testing_method"], "NA")

    def test_dates_are_validated_without_guessing(self):
        for value, expected in [("2020", "2020"), ("2020-03", "2020"), ("2020-02-29", "2020"),
                                ("2019-02-29", "NA"), ("2020-13", "NA"), ("2020/2021", "NA"),
                                ("03/04/2020", "NA"), ("0000", "NA"), ("", "NA")]:
            self.assertEqual(ast.year_of(value), expected)
        row = self.normalized({"collection_date": "2020/2021"}, [])
        self.assertEqual(row["collection_date"], "2020/2021")
        self.assertIn("invalid_collection_date", row["qc_flags"])

    def test_location_and_source_rules_are_conservative(self):
        row = self.normalized({"geo_loc_name": "Australia: Sydney"}, [])
        self.assertEqual((row["country"], row["location"]), ("Australia", "Australia: Sydney"))
        for value in ("Sydney", "USA;Canada", "Atlantic Ocean", ": no country"):
            row = self.normalized({"geo_loc_name": value}, [])
            self.assertEqual(row["country"], "NA")
            self.assertIn("ambiguous_country", row["qc_flags"])
        for host, source, category in [("Homo sapiens", "urine", "human"), ("pig", "", "animal"),
                                       ("", "ground beef", "food"), ("", "water", "environment"),
                                       ("", "clinical", "unknown"), ("chicken", "food", "unknown")]:
            row = self.normalized({"host": host, "isolation_source": source}, [])
            self.assertEqual(row["source_category"], category)
        self.assertIn("ambiguous_source", row["qc_flags"])

    def test_organism_filter_does_not_accept_genus_or_shigella(self):
        for value, accepted in [("Escherichia coli", True), ("Escherichia coli O157:H7", True),
                                ("Escherichia", False), ("Escherichia coli-like", False),
                                ("Escherichia fergusonii", False), ("Shigella sonnei", False)]:
            row = self.normalized({"scientific_name": value}, [])
            self.assertEqual("organism_mismatch" not in row["qc_flags"], accepted)

    def test_duplicates_are_preserved_but_not_double_counted(self):
        self.snapshot([self.a, {**self.a, "id": "another-submission"}])
        rows, endpoints, summary = ast.prepare(self.raw)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all("duplicate_record" in r["qc_flags"] for r in rows))
        self.assertEqual(endpoints[0]["distinct_observations"], "1")
        self.assertEqual(endpoints[0]["eligible"], "true")
        self.assertEqual(summary["eligible_isolate_drug_pairs"], 1)
        self.assertNotEqual(rows[0]["record_id"], rows[1]["record_id"])

    def test_conflicting_and_multiple_assays_are_not_resolved(self):
        self.snapshot([self.a, {**self.a, "id": "second", "phenotype": "susceptible", "standard": "EUCAST"}])
        rows, endpoints, summary = ast.prepare(self.raw)
        self.assertTrue(all("conflicting_ast" in r["qc_flags"] for r in rows))
        self.assertEqual(endpoints[0]["phenotype"], "NA")
        self.assertEqual(endpoints[0]["eligible"], "false")
        self.assertEqual(summary["endpoint_exclusion_counts_overlapping"]["conflicting_ast"], 1)
        second = self.normalized({"standard": "EUCAST"})
        first = self.normalized()
        endpoints = ast.annotate_groups([first, second])
        self.assertNotIn("conflicting_ast", first["qc_flags"])
        self.assertIn("multiple_ast_records", first["qc_flags"])
        self.assertEqual(endpoints[0]["eligible"], "false")

    def test_metadata_conflicts_never_choose_first(self):
        for field in ("biosample_acc", "asm_acc", "scientific_name", "geo_loc_name"):
            row = self.normalized(metadata=[self.m, {**self.m, field: "different"}])
            self.assertIn("metadata_conflict", row["qc_flags"])
        row = self.normalized(metadata=[self.m, {**self.m, "asm_acc": "GCA_900000002.1"}])
        self.assertEqual(row["assembly_accession"], "NA")
        self.assertEqual(ast.annotate_groups([row])[0]["eligible"], "false")

    def test_duplicate_metadata_does_not_expand_ast_and_blank_keys_do_not_join(self):
        self.snapshot([self.a, {**self.a, "target_acc": ""}],
                      [self.m, self.m, self.m, {**self.m, "target_acc": ""}])
        rows, _, summary = ast.prepare(self.raw)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(json.loads(rows[0]["metadata_records_json"])), 3)
        self.assertEqual(json.loads(rows[1]["metadata_records_json"]), [])
        self.assertEqual(rows[1]["assembly_accession"], "NA")
        self.assertEqual(summary["raw_metadata_records"], 4)
        self.assertLessEqual(set(summary["qc_record_counts_all"]), ast.FLAGS)

    def test_undefined_phenotypes_missing_ids_and_non_ecoli(self):
        self.snapshot([{**self.a, "phenotype": "ND"}, {**self.a, "target_acc": "", "antibiotic": "AMP"},
                       {**self.a, "target_acc": "", "antibiotic": "AMP"},
                       {**self.a, "target_acc": "OTHER", "scientific_name": "Shigella sonnei"}])
        rows, endpoints, summary = ast.prepare(self.raw)
        self.assertEqual(summary["ecoli_ast_records"], 3)
        self.assertEqual(len(endpoints), 3)  # missing identifiers never collapse together
        self.assertTrue(all(e["eligible"] == "false" for e in endpoints))
        self.assertTrue(any("undefined_phenotype" in e["exclusion_reasons"] for e in endpoints))
        self.assertEqual(summary["flow"][1]["count"], 1)
        self.assertEqual(len(rows), 4)

    def test_optional_method_and_breakpoint_version_preserved(self):
        row = self.normalized({"method": "broth microdilution", "standard_version": "2020"})
        self.assertEqual((row["testing_method"], row["breakpoint_version"]), ("broth microdilution", "2020"))
        self.assertNotIn("missing_testing_method", row["qc_flags"])
        self.snapshot([{**self.a, "method": "agar dilution", "standard_version": "2019"}], extras=ast.AST_OPTIONAL)
        self.assertEqual(ast.prepare(self.raw)[0][0]["testing_method"], "agar dilution")

    def test_header_only_build_and_rebuild_determinism(self):
        self.snapshot([], [])
        summary = ast.build(self.raw, self.out)
        self.assertEqual(summary["raw_ast_records"], 0)
        self.assertEqual(summary["eligible_isolate_drug_pairs"], 0)
        self.assertEqual((self.out / "ast.normalized.tsv").read_text(), "\t".join(ast.SCHEMA) + "\n")
        again = self.root / "again"
        ast.build(self.raw, again)
        self.assertEqual({p.name: p.read_bytes() for p in self.out.iterdir()}, {p.name: p.read_bytes() for p in again.iterdir()})

    def test_snapshot_immutable_outputs_and_no_network_on_rebuild(self):
        self.snapshot()
        before = {p.name: p.read_bytes() for p in self.raw.iterdir()}
        with patch.object(ast, "urlopen", side_effect=AssertionError("network used")):
            ast.build(self.raw, self.out)
            self.assertEqual(ast.validate(self.raw, self.out)["raw_ast_records"], 1)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.raw.iterdir()})
        with self.assertRaisesRegex(ValueError, "outside the raw snapshot"):
            ast.build(self.raw, self.raw / "processed")
        with self.assertRaisesRegex(ValueError, "already exists"):
            ast.build(self.raw, self.out)
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.snapshot()
        link = self.root / "link"
        link.symlink_to(self.root / "nonexistent")
        with self.assertRaisesRegex(ValueError, "already exists"):
            with ast.new_directory(link):
                pass
        race = self.root / "race"
        with self.assertRaisesRegex(ValueError, "appeared"):
            with ast.new_directory(race):
                race.mkdir()

    def test_schema_changes_and_malformed_tsv_fail_before_publication(self):
        path = self.root / "bad.tsv"
        for fields in [[], ast.AST_FIELDS[:-1], ast.AST_FIELDS + ["surprise"], ast.AST_FIELDS + [ast.AST_FIELDS[0]]]:
            path.write_text("\t".join(fields) + "\n")
            with self.assertRaises(ValueError):
                ast.read_table(path, "ast")
        for body in ['only-one-cell\n', '\t'.join(['x'] * (len(ast.AST_FIELDS) + 1)) + '\n', '"unclosed\n']:
            path.write_text("\t".join(ast.AST_FIELDS) + "\n" + body)
            with self.assertRaises((ValueError, csv.Error)):
                ast.read_table(path, "ast")
        path.write_text("#" + "\t".join(ast.AST_FIELDS) + "\n", encoding="utf-8-sig")
        self.assertEqual(ast.read_table(path, "ast"), [])
        with self.assertRaises(ValueError):
            ast.archive(path, [], self.raw, "x", "2026-10-10T00:00:00Z")
        with self.assertRaises(ValueError):
            ast.archive(path, [path], self.raw, "x", "2026-10-10")
        with self.assertRaises(ValueError):
            ast.archive(path, [path], self.raw, "x", "2026-10-10T00:00:00Z")
        self.assertFalse(self.raw.exists())

    def test_manifest_errors_and_checksum_tampering(self):
        self.snapshot()
        manifest_path = self.raw / "manifest.json"
        original = json.loads(manifest_path.read_text())
        for invalid in [[], None, "not-an-object"]:
            manifest_path.write_text(json.dumps(invalid))
            with self.assertRaisesRegex(ValueError, "JSON object"):
                ast.load_snapshot(self.raw)
        changes = [{"schema_version": "9"}, {"parser_version": "9"}, {"sources": []}, {"sources": "bad"},
                   {"sources": original["sources"][:1]}, {"sources": original["sources"] * 2}]
        for field, value in [("filename", "../escape.tsv"), ("filename", "."), ("role", "other"),
                             ("source_url", ""), ("source_name", 3), ("sha256", "invalid"),
                             ("schema_version", "9"), ("retrieved_at", "2026-10-10"), ("expected_rows", 9),
                             ("expected_rows", True)]:
            entry = {**original["sources"][0], field: value}
            changes.append({"sources": [entry, original["sources"][1]]})
        changes.append({"sources": [{"role": "ast"}]})
        for change in changes:
            manifest_path.write_text(json.dumps({**original, **change}))
            with self.assertRaises((ValueError, OSError)):
                ast.build(self.raw, self.out)
            self.assertFalse(self.out.exists())
        manifest_path.write_text(json.dumps(original))
        raw = self.raw / original["sources"][0]["filename"]
        raw.write_text(raw.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            ast.validate(self.raw)

    def test_processed_validation_detects_changes(self):
        self.snapshot()
        ast.build(self.raw, self.out)
        file = self.out / "cohort_summary.json"
        file.write_text("{}\n")
        with self.assertRaisesRegex(ValueError, "differs from rebuild"):
            ast.validate(self.raw, self.out)
        file.unlink()
        with self.assertRaisesRegex(ValueError, "file set"):
            ast.validate(self.raw, self.out)

    def test_antimicrobial_summary_denominators(self):
        rows = [self.normalized(), self.normalized({"target_acc": "SECOND", "phenotype": "S"}),
                self.normalized({"target_acc": "THIRD", "phenotype": "I"}),
                self.normalized({"target_acc": "FOURTH", "antibiotic": "AMP", "phenotype": "NS"})]
        endpoints = ast.annotate_groups(rows)
        summary, isolates, drugs = ast.summarize_rows(rows, endpoints)
        cip = next(r for r in drugs if r["antimicrobial"] == "ciprofloxacin")
        amp = next(r for r in drugs if r["antimicrobial"] == "ampicillin")
        self.assertEqual((cip["R"], cip["S"], cip["other"], cip["R_fraction_of_RS"]), (1, 1, 1, "0.500000"))
        self.assertEqual(amp["R_fraction_of_RS"], "NA")
        self.assertEqual((summary["ecoli_unique_isolates"], len(isolates)), (4, 4))
        rows = [self.normalized({"collection_date": "", "phenotype": "ND", "geo_loc_name": ""}, [])]
        self.assertEqual(ast.summarize_rows(rows, ast.annotate_groups(rows))[2][0]["year_min"], "NA")

    def test_ambiguous_policy_fails(self):
        policy = copy.deepcopy(self.policy)
        policy["drugs"]["new drug"] = ["CIP"]
        path = self.root / "policy.json"
        path.write_text(json.dumps(policy))
        with patch.object(ast, "POLICY", path), self.assertRaisesRegex(ValueError, "ambiguous drug alias"):
            ast.load_policy()

    def test_export_url_and_mocked_download(self):
        query = 'biosample_acc:SAMN90000001'
        params = parse_qs(urlparse(ast.export_url("ast", query, ast.AST_FIELDS)).query)
        self.assertIn(query, unquote(params["q"][0]))
        self.assertEqual(params["nolimit"], ["on"])
        self.assertIn("mic_secondary|mic_secondary", params["fields"][0])
        def response(url, timeout):
            self.assertEqual(timeout, 60)
            params = parse_qs(urlparse(url).query)
            if params["action"] == ["retrieve"]:
                return io.BytesIO(json.dumps({"success": True, "ngout": {"data": {"totalCount": 1}}}).encode())
            table = "ast.synthetic.tsv" if ".from(ast)" in params["q"][0] else "isolates.synthetic.tsv"
            return io.BytesIO((FIXTURES / table).read_bytes())
        with patch.object(ast, "urlopen", side_effect=response) as transport:
            ast.download(self.raw, query)
            self.assertEqual(transport.call_count, 4)
        manifest = json.loads((self.raw / "manifest.json").read_text())
        self.assertEqual(manifest["acquisition"], "official-table-export")
        self.assertEqual(manifest["sources"][0]["expected_rows"], 1)
        self.assertEqual(ast.validate(self.raw)["raw_ast_records"], 1)

    def test_download_limits_failures_and_empty_selection(self):
        with patch.object(ast, "urlopen", return_value=io.BytesIO(b"12345")), self.assertRaisesRegex(ValueError, "exceeds"):
            ast.retrieve("https://example.invalid", self.root / "large", 4)
        with patch.object(ast, "urlopen", side_effect=OSError("offline")), self.assertRaises(OSError):
            ast.download(self.raw)
        self.assertFalse(self.raw.exists())
        for count in [-1, "1", None]:
            payload = json.dumps({"success": True, "ngout": {"data": {"totalCount": count}}}).encode()
            with patch.object(ast, "urlopen", return_value=io.BytesIO(payload)), self.assertRaises(ValueError):
                ast.remote_count("ast", "query")
        for query, limit in [("", 100), ("x", 0)]:
            with self.assertRaises(ValueError):
                ast.download(self.raw, query, limit)
        def fake_fetch(collection, query, fields, path, role, max_bytes):
            ast.write_tsv(path, fields, [])
            return ast.source_entry(path, role, "fixture:empty", "2026-10-10T00:00:00Z", path.name), []
        with patch.object(ast, "fetch_table", side_effect=fake_fetch) as fetch:
            ast.download(self.raw)
            self.assertIn("__empty_selection__", fetch.call_args_list[1].args[1])
        self.assertEqual(ast.validate(self.raw)["raw_ast_records"], 0)

    def test_truncated_export_and_bad_accession_are_rejected(self):
        def retrieve(url, output, max_bytes):
            output.write_bytes((FIXTURES / "ast.synthetic.tsv").read_bytes())
        with patch.object(ast, "remote_count", return_value=2), patch.object(ast, "retrieve", side_effect=retrieve):
            with self.assertRaisesRegex(ValueError, "truncated"):
                ast.download(self.raw)
        self.assertFalse(self.raw.exists())
        with patch.object(ast, "fetch_table", return_value=({}, [{"target_acc": "bad query"}])):
            with self.assertRaisesRegex(ValueError, "PDT accession"):
                ast.download(self.raw)

    def test_metadata_download_batches_have_complete_accession_coverage(self):
        targets = [f"PDT{900000000 + i}.1" for i in range(101)]
        requests = []
        def fetch(collection, query, fields, path, role, max_bytes):
            if collection == "ast":
                rows = [{**self.a, "target_acc": target} for target in targets]
            else:
                requests.append(query)
                rows = []
            ast.write_tsv(path, fields, rows)
            return ast.source_entry(path, role, "fixture:batch", "2026-10-10T00:00:00Z", path.name), rows
        with patch.object(ast, "fetch_table", side_effect=fetch):
            ast.download(self.raw)
        self.assertEqual(len(requests), 2)
        self.assertEqual(sum(query.count("PDT") for query in requests), 101)
        for target in targets:
            self.assertEqual(sum(target in query for query in requests), 1)

    def test_cli_commands_and_failures(self):
        with patch("sys.stdout", new_callable=io.StringIO), patch("sys.stderr", new_callable=io.StringIO):
            self.assertEqual(ast.main(["archive", "--ast", str(FIXTURES / "ast.synthetic.tsv"), "--metadata",
                                      str(FIXTURES / "isolates.synthetic.tsv"), "--outdir", str(self.raw),
                                      "--source-url", "fixture:synthetic", "--retrieved-at", "2026-10-10T00:00:00Z"]), 0)
            self.assertEqual(ast.main(["build", "--snapshot", str(self.raw), "--outdir", str(self.out)]), 0)
            self.assertEqual(ast.main(["validate", "--snapshot", str(self.raw), "--processed", str(self.out)]), 0)
            self.assertEqual(ast.main(["summarize", "--snapshot", str(self.raw)]), 0)
            self.assertEqual(ast.main(["validate", "--snapshot", str(self.root / "absent")]), 2)
            with patch.object(ast, "download") as download:
                self.assertEqual(ast.main(["download", "--outdir", str(self.root / "other")]), 0)
                download.assert_called_once()


if __name__ == "__main__":
    unittest.main()
