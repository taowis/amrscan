"""Offline tests for the genotype/AST evidence linker and deterministic summary."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("amrfinder_validation", ROOT / "bin/amrfinder_validation.py")
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)
amr = validation.amrscan
genomes = validation.genome_dataset
FIXTURE = ROOT / "tests/fixtures/amrfinder_v4.0.23.tsv"


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.pdt = "PDT000000001.1"
        self.biosample = "SAMN00000001"
        self.assembly = "GCF_000000001.1"
        self.manifest = self.root / "cohort.tsv"
        self.genomes = self.root / "genomes"
        self.genomes.mkdir()
        self.fasta = b">contig1\nACGTACGT\n"
        fasta_path = self.genomes / f"{self.pdt}_{self.assembly}.fna"
        fasta_path.write_bytes(self.fasta)
        row = {"cohort_id": "test", "isolate_id": self.pdt, "pdt_accession": self.pdt,
               "biosample_accession": self.biosample, "assembly_accession": self.assembly,
               "organism": "Escherichia coli", "assembly_source": "NCBI Assembly",
               "ast_available": "true", "selection_reason": "synthetic fixture"}
        with self.manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=genomes.FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        self.genomes_provenance = {"schema_version": "1.0", "manifest_sha256": genomes.sha256(self.manifest),
            "assemblies": [{**row, "file": fasta_path.name, "sha256": genomes.sha256(fasta_path),
                            "source_url": "https://example.invalid/test", "retrieved_at": "2026-01-01T00:00:00+00:00"}]}
        (self.genomes / "provenance.json").write_text(json.dumps(self.genomes_provenance))

        self.ast_path = self.root / "ast.tsv"
        ast_row = {"record_id": "record-one", "isolate_id": self.pdt, "biosample": self.biosample,
                   "assembly_accession": "NA", "organism": "Escherichia coli", "antimicrobial": "ampicillin",
                   "phenotype": "R", "source_url": "https://example.invalid/ast", "source_sha256": "a" * 64,
                   "phenotype_original": "resistant", "raw_record_json": "{}"}
        fields = list(ast_row)
        with self.ast_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(ast_row)
        self.calls = self.root / "calls"
        self.calls.mkdir()
        self.out = self.root / "result"
        self.database_files = {"fixture.marker": "a" * 64}
        self.database_id = "sha256:" + hashlib.sha256(amr.json_text(self.database_files).encode()).hexdigest()
        self._write_calls([self._call_row()])

    def _call_row(self):
        return {"sample_id": self.pdt, "tool": "AMRFinderPlus", "gene": "blaTEM-1", "allele": "blaTEM-1",
                "resistance_class": "BETA-LACTAM", "drug": "NA", "identity": "99.9", "coverage": "100",
                "contig": "contig1", "start": "1", "end": "800", "strand": "+", "evidence_type": "ALLELE",
                "schema_version": "1.0", "caller_version": "4.2.7", "database_id": self.database_id,
                "thresholds_json": "{}", "element_type": "AMR", "element_subtype": "NA", "scope": "core",
                "subclass": "penicillin", "raw_file": "raw.tsv", "raw_row": "1", "raw_record_json": "{}"}

    def _write_calls(self, rows):
        raw = self.calls / f"{self.pdt}.amrfinder.tsv"
        if rows:
            with FIXTURE.open() as stream:
                source = list(csv.reader(stream, delimiter="\t"))
            with raw.open("w", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(source[0])
                writer.writerow(source[2])  # one published report row: blaTEM-156
        else:
            raw.write_text(FIXTURE.read_text().splitlines()[0] + "\n")
        provenance = {"sample_id": self.pdt, "tool": "AMRFinderPlus", "status": "success",
                      "caller_version": "4.2.7", "database_id": self.database_id,
                      "database_files_sha256": self.database_files,
                      "input_sha256": self.genomes_provenance["assemblies"][0]["sha256"],
                      "raw_sha256": amr.sha256(raw), "thresholds": {}}
        (self.calls / f"{self.pdt}.provenance.json").write_text(json.dumps(provenance))
        # Harmonize upstream-shaped test evidence, then change only sample/version/database identifiers.
        fixture_provenance = {"sample_id": "upstream", "tool": "AMRFinderPlus", "status": "upstream_fixture",
                              "caller_version": "4.2.7", "database_id": self.database_id, "thresholds": {},
                              "raw_sha256": amr.sha256(raw)}
        (self.calls / f"{self.pdt}.harmonize-provenance.json").write_text(json.dumps(fixture_provenance))
        amr.harmonize(raw, self.calls / f"{self.pdt}.harmonize-provenance.json", "upstream",
                      self.calls / f"{self.pdt}.harmonized.tsv")
        with (self.calls / f"{self.pdt}.harmonized.tsv").open() as stream:
            converted = list(csv.DictReader(stream, delimiter="\t"))
        output = self.calls / f"{self.pdt}.harmonized.tsv"
        with output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=amr.SCHEMA, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for item in converted:
                item["sample_id"] = self.pdt
                writer.writerow(item)

    def run_summary(self, exclusions=None):
        return validation.assemble(self.manifest, self.genomes, self.ast_path,
                                   self.calls, self.out, exclusions)

    def test_links_verified_assembly_to_separate_preserved_evidence_and_counts(self):
        summary = self.run_summary()
        self.assertEqual(summary["selected_cohort_isolates"], 1)
        self.assertEqual(summary["valid_assemblies"], 1)
        self.assertEqual(summary["isolates_with_determinant_calls"], 1)
        self.assertEqual(summary["isolates_linked_to_ast"], 1)
        self.assertEqual(summary["gene_family_counts"], {"blaTEM-156": 1})
        with (self.out / "genotype_evidence.tsv").open() as stream:
            genotype = list(csv.DictReader(stream, delimiter="\t"))
        with (self.out / "ast_evidence.tsv").open() as stream:
            ast = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(genotype[0]["pdt_accession"], ast[0]["isolate_id"])
        self.assertEqual(ast[0]["verified_assembly_accession"], self.assembly)
        self.assertEqual(ast[0]["assembly_accession"], "NA")  # keep source cell distinct from verified link
        self.assertEqual(ast[0]["phenotype"], "R")
        self.assertNotIn("predicted_phenotype", ast[0])
        self.assertNotIn("phenotype", genotype[0])
        self.assertEqual(json.loads((self.out / "validation_provenance.json").read_text())["inputs"]["ast_table_sha256"],
                         hashlib.sha256(self.ast_path.read_bytes()).hexdigest())

    def test_all_ast_observations_and_nonbinary_categories_are_retained(self):
        rows = validation.read_tsv(self.ast_path)
        for index, category in enumerate(("ND", "I", "SDD", "NS", "HLAR", "S"), 2):
            rows.append({**rows[0], "record_id": f"record-{index}", "phenotype": category})
        self._rewrite_ast(rows)
        summary = self.run_summary()
        self.assertEqual(summary["isolates_linked_to_ast"], 1)
        self.assertEqual(summary["ast_observation_counts_by_source_phenotype"],
                         {"HLAR": 1, "I": 1, "ND": 1, "R": 1, "S": 1, "SDD": 1, "NS": 1})
        with (self.out / "ast_evidence.tsv").open() as stream:
            linked = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(linked), 7)
        self.assertEqual({r["phenotype"] for r in linked}, {"ND", "I", "SDD", "NS", "HLAR", "S", "R"})

    def test_no_amr_records_is_reported_as_no_determinants_not_susceptible(self):
        self._write_calls([])
        summary = self.run_summary()
        self.assertEqual(summary["isolates_with_no_determinant_calls"], 1)
        self.assertEqual(summary["isolates_with_determinant_calls"], 0)
        self.assertEqual(summary["amr_determinant_records"], 0)

    def test_ast_mismatches_and_missing_exact_ast_fail_closed(self):
        rows = validation.read_tsv(self.ast_path)
        rows[0]["biosample"] = "SAMN99999999"
        self._rewrite_ast(rows)
        with self.assertRaisesRegex(ValueError, "BioSample mismatch"):
            self.run_summary()
        rows[0]["biosample"] = self.biosample
        rows[0]["organism"] = "Shigella flexneri"
        self._rewrite_ast(rows)
        with self.assertRaisesRegex(ValueError, "organism mismatch"):
            self.run_summary()
        rows[0]["organism"] = "Escherichia coli"
        rows[0]["assembly_accession"] = "GCF_999999999.1"
        self._rewrite_ast(rows)
        with self.assertRaisesRegex(ValueError, "Assembly mismatch"):
            self.run_summary()
        self._rewrite_ast([])
        with self.assertRaisesRegex(ValueError, "no exact-PDT AST"):
            self.run_summary()

    def test_unselected_ast_rows_are_ignored_and_duplicate_ids_rejected(self):
        rows = validation.read_tsv(self.ast_path)
        rows.append({**rows[0], "record_id": "external", "isolate_id": "PDT000000099.1"})
        self._rewrite_ast(rows)
        self.assertEqual(self.run_summary()["isolates_linked_to_ast"], 1)
        self.out = self.root / "duplicate-result"
        rows.append({**rows[0], "record_id": rows[0]["record_id"]})
        self._rewrite_ast(rows)
        with self.assertRaisesRegex(ValueError, "duplicate AST record id"):
            self.run_summary()

    def test_caller_genome_hash_mismatch_fails_closed(self):
        path = self.calls / f"{self.pdt}.provenance.json"
        data = json.loads(path.read_text())
        data["input_sha256"] = "wrong"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "input hash"):
            self.run_summary()

    def test_exclusion_counts_and_duplicate_selection_rejected(self):
        exclusions = self.root / "excluded.tsv"
        with exclusions.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=validation.EXCLUSION_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"pdt_accession": "PDT000000002.1", "biosample_accession": "SAMN00000002",
                             "assembly_accession": "NA", "reason": "no_assembly"})
        summary = self.run_summary(exclusions)
        self.assertEqual(summary["exclusions_by_reason"], {"no_assembly": 1})
        self.out = self.root / "result2"
        with exclusions.open("a") as stream:
            stream.write(f"{self.pdt}\t{self.biosample}\t{self.assembly}\talready_selected\n")
        with self.assertRaisesRegex(ValueError, "both selected and excluded"):
            self.run_summary(exclusions)

    def test_malformed_inputs_and_caller_provenance_are_rejected(self):
        bad_ast = self.root / "bad.tsv"
        bad_ast.write_text("a\ta\n1\n")
        with self.assertRaisesRegex(ValueError, "duplicate TSV header"):
            validation.read_tsv(bad_ast)
        bad_ast.write_text("record_id\tisolate_id\n1\tPDT\textra\n")
        with self.assertRaisesRegex(ValueError, "malformed TSV"):
            validation.read_tsv(bad_ast)
        bad_ast.write_text("record_id\n1\n")
        with self.assertRaisesRegex(ValueError, "missing columns"):
            validation.read_tsv(bad_ast, {"record_id", "phenotype"})
        excluded = self.root / "bad-exclusions.tsv"
        excluded.write_text("pdt_accession\tbiosample_accession\tassembly_accession\treason\nPDT\t\tNA\tmissing\n")
        with self.assertRaisesRegex(ValueError, "all fields are required"):
            validation.validate_exclusions(excluded)

        provenance_path = self.calls / f"{self.pdt}.provenance.json"
        original = json.loads(provenance_path.read_text())
        for change, message in [({"status": "failed"}, "caller provenance"),
                                ({"caller_version": "NA"}, "caller version"),
                                ({"raw_sha256": "wrong"}, "report hash"),
                                ({"database_id": "sha256:wrong"}, "database fingerprint")]:
            provenance_path.write_text(json.dumps({**original, **change}))
            with self.assertRaisesRegex(ValueError, message):
                self.run_summary()
        provenance_path.write_text(json.dumps(original))

        harmonized = self.calls / f"{self.pdt}.harmonized.tsv"
        with harmonized.open() as stream:
            records = list(csv.DictReader(stream, delimiter="\t"))
        records[0]["sample_id"] = "wrong"
        self._rewrite_harmonized(records)
        with self.assertRaisesRegex(ValueError, "sample/tool mismatch"):
            self.run_summary()
        records[0]["sample_id"] = self.pdt
        records[0]["database_id"] = "wrong"
        self._rewrite_harmonized(records)
        with self.assertRaisesRegex(ValueError, "provenance mismatch"):
            self.run_summary()
        records[0]["database_id"] = self.database_id
        records[0]["gene"] = "altered"
        self._rewrite_harmonized(records)
        with self.assertRaisesRegex(ValueError, "does not match raw caller output"):
            self.run_summary()

    def test_output_path_is_never_overwritten(self):
        self.run_summary()
        with self.assertRaisesRegex(ValueError, "output already exists"):
            self.run_summary()

    def test_output_race_and_write_failure_clean_staging(self):
        original_write = validation.write_tsv
        calls = 0

        def create_racing_output(path, fields, rows):
            nonlocal calls
            calls += 1
            if calls == 1:
                self.out.mkdir()
            return original_write(path, fields, rows)

        with patch.object(validation, "write_tsv", side_effect=create_racing_output):
            with self.assertRaisesRegex(ValueError, "appeared during build"):
                self.run_summary()
        self.out.rmdir()
        self.assertFalse(list(self.root.glob(".amrfinder-validation-*")))

        self.out = self.root / "failed-result"
        with patch.object(validation, "write_tsv", side_effect=OSError("fixture write failure")):
            with self.assertRaisesRegex(OSError, "fixture write failure"):
                self.run_summary()
        self.assertFalse(self.out.exists())
        self.assertFalse(list(self.root.glob(".amrfinder-validation-*")))

    def test_cli_and_module_entrypoint(self):
        args = ["summarize", "--manifest", str(self.manifest), "--genomes", str(self.genomes),
                "--ast", str(self.ast_path), "--calls", str(self.calls), "--outdir", str(self.out)]
        with patch("sys.stdout"):
            self.assertEqual(validation.main(args), 0)
        script = ROOT / "bin/amrfinder_validation.py"
        self.out = self.root / "entrypoint-result"
        argv = [str(script), *args[:-1], str(self.out)]
        with patch.object(sys, "argv", argv), patch("sys.stdout"):
            with self.assertRaises(SystemExit) as raised:
                runpy.run_path(str(script), run_name="__main__")
        self.assertEqual(raised.exception.code, 0)

    def _rewrite_ast(self, rows):
        fields = list(rows[0]) if rows else list(validation.AST_REQUIRED)
        with self.ast_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _rewrite_harmonized(self, rows):
        with (self.calls / f"{self.pdt}.harmonized.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=amr.SCHEMA, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
