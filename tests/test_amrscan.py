"""Offline regression tests; reference rows are NCBI's published expected output."""

import argparse
import csv
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("amrscan", ROOT / "bin/amrscan.py")
amrscan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(amrscan)
FIXTURE = ROOT / "tests/fixtures/amrfinder_v4.0.23.tsv"


def fixture_provenance(raw, sample="upstream"):
    # This is a source fixture, not a locally executed caller or a known database release.
    return {"sample_id": sample, "tool": "AMRFinderPlus", "status": "upstream_fixture",
            "caller_version": "NA", "database_id": "NA", "thresholds": {},
            "raw_sha256": amrscan.sha256(raw)}


class HarmonizeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.raw = self.directory / "reference.tsv"
        self.raw.write_bytes(FIXTURE.read_bytes())
        self.provenance = self.directory / "provenance.json"
        self.output = self.directory / "result.tsv"

    def convert(self, metadata=None):
        self.provenance.write_text(json.dumps(metadata or fixture_provenance(self.raw)))
        count = amrscan.harmonize(self.raw, self.provenance, "upstream", self.output)
        with self.output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return count, rows

    def mutate(self, column, value):
        with self.raw.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers, rows = reader.fieldnames, list(reader)
        rows[0][column] = value
        with self.raw.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_all_upstream_evidence_is_preserved(self):
        count, rows = self.convert()
        with FIXTURE.open() as handle:
            original = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(count, 24)
        self.assertEqual([json.loads(r["raw_record_json"]) for r in rows], original)
        self.assertEqual(list(rows[0]), amrscan.SCHEMA)
        self.assertEqual([int(r["raw_row"]) for r in rows], list(range(1, 25)))
        self.assertEqual(rows[1]["gene"], "blaTEM-156")
        self.assertEqual(rows[1]["allele"], "blaTEM-156")
        self.assertEqual(rows[2]["allele"], "NA")  # BLAST family cannot imply closest allele.
        self.assertEqual(rows[2]["identity"], "99.75")
        self.assertEqual(rows[8]["strand"], "-")
        self.assertEqual((rows[8]["start"], rows[8]["end"]), ("1", "651"))
        self.assertIn("INTERNAL_STOP", {r["evidence_type"] for r in rows})
        self.assertEqual(sum(r["gene"] == "pmrB_C84R" for r in rows), 2)
        self.assertEqual(rows[22]["coverage"], "NA")
        self.assertEqual(rows[13]["element_type"], "STRESS")
        self.assertEqual(rows[13]["resistance_class"], "NA")
        self.assertEqual({r["drug"] for r in rows}, {"NA"})
        self.assertEqual({r["caller_version"] for r in rows}, {"NA"})
        self.assertNotIn("phenotype", rows[0])

    def test_header_only_means_no_calls(self):
        self.raw.write_text(FIXTURE.read_text().splitlines()[0] + "\n")
        count, rows = self.convert()
        self.assertEqual((count, rows), (0, []))
        self.assertEqual(self.output.read_text().rstrip().split("\t"), amrscan.SCHEMA)

    def test_explicit_missing_values(self):
        for column in ("Start", "Stop", "Strand", "% Identity to reference", "Element symbol"):
            self.mutate(column, "")
        _, rows = self.convert()
        for column in ("start", "end", "strand", "identity", "gene"):
            self.assertEqual(rows[0][column], "NA")

    def test_malformed_rows_do_not_publish_partial_results(self):
        for column, value in [("Method", "NA"), ("Start", "0"), ("Start", "1.5"),
                              ("Stop", "NA"), ("Start", "1000"), ("Strand", "?"),
                              ("% Identity to reference", "NaN"),
                              ("% Identity to reference", "inf"),
                              ("% Identity to reference", "-1"),
                              ("% Coverage of reference", "101"),
                              ("% Coverage of reference", "oops")]:
            with self.subTest(column=column, value=value):
                self.raw.write_bytes(FIXTURE.read_bytes())
                self.mutate(column, value)
                with self.assertRaises(ValueError):
                    self.convert()
                self.assertFalse(self.output.exists())

    def test_invalid_headers_and_row_width(self):
        for content in ("", "bad\theader\n", FIXTURE.read_text().replace("Contig id", "Protein id", 1),
                        FIXTURE.read_text() + "short\trow\n", FIXTURE.read_text().replace("\n", "\textra\n", 2)):
            with self.subTest(content=content[:30]):
                self.raw.write_text(content)
                with self.assertRaises(ValueError):
                    self.convert()

    def test_malformed_quoting_fails_without_publishing_partial_results(self):
        lines = FIXTURE.read_text().splitlines()
        record = lines[-1].split("\t")
        for last_field in ('"unterminated', '"closed"unexpected'):
            with self.subTest(last_field=last_field):
                record[-1] = last_field
                self.raw.write_text("\n".join(lines + ["\t".join(record)]) + "\n")
                with self.assertRaises(csv.Error):
                    self.convert()
                self.assertFalse(self.output.exists())
                with patch("sys.stderr", new_callable=io.StringIO) as stderr, \
                     self.assertRaises(SystemExit) as raised:
                    amrscan.main(["harmonize", "--sample", "upstream", "--raw", str(self.raw),
                                  "--provenance", str(self.provenance), "--output", str(self.output)])
                self.assertEqual(raised.exception.code, 1)
                self.assertIn("AMRScan:", stderr.getvalue())
                self.assertNotIn("Traceback", stderr.getvalue())
                self.assertFalse(self.output.exists())

        # Strict parsing must still retain valid quoted evidence exactly.
        self.raw.write_bytes(FIXTURE.read_bytes())
        self.mutate("HMM description", 'annotation with "quotes" and\ta tab')
        count, rows = self.convert()
        self.assertEqual(count, 24)
        self.assertEqual(json.loads(rows[0]["raw_record_json"])["HMM description"],
                         'annotation with "quotes" and\ta tab')

    def test_provenance_rejects_wrong_sample_tool_hash_and_failed_run(self):
        for key, value in [("sample_id", "another"), ("tool", "BLAST"),
                           ("raw_sha256", "wrong"), ("status", "failed")]:
            with self.subTest(key=key):
                metadata = fixture_provenance(self.raw)
                metadata[key] = value
                with self.assertRaises(ValueError):
                    self.convert(metadata)

    def test_raw_evidence_cannot_be_overwritten(self):
        for output in (self.raw, self.provenance):
            with self.assertRaisesRegex(ValueError, "overwrite"):
                amrscan.harmonize(self.raw, self.provenance, "upstream", output)

    def test_unsafe_sample_ids(self):
        for sample in ("", "../x", "a b", "x;touch", "-option", "a\tb"):
            with self.subTest(sample=sample), self.assertRaises(ValueError):
                amrscan.sample_id(sample)

    def test_cli_success_and_errors(self):
        self.provenance.write_text(json.dumps(fixture_provenance(self.raw)))
        args = ["harmonize", "--sample", "upstream", "--raw", str(self.raw),
                "--provenance", str(self.provenance), "--output", str(self.output)]
        amrscan.main(args)
        self.assertTrue(self.output.exists())
        self.provenance.write_text("{}")
        with patch("sys.stderr"), self.assertRaises(SystemExit) as raised:
            amrscan.main(args)
        self.assertEqual(raised.exception.code, 1)


class CallerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.database = self.directory / "test db"
        self.database.mkdir()
        # Deliberate dummy database used only to test hashing/command construction.
        (self.database / "marker").write_text("NOT AN AMR DATABASE\n")
        self.args = argparse.Namespace(sample="upstream", input=ROOT / "tests/fixtures/tiny.fna",
            database=self.database, database_source=None, executable="amrfinder",
            threads=2, identity=-1, coverage=0.5, organism=None, plus=False,
            raw=self.directory / "raw.tsv", provenance=self.directory / "provenance.json",
            log=self.directory / "caller.log")

    def mock_caller(self, command, **kwargs):
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, stdout="TEST-DOUBLE\n")
        self.args.raw.write_bytes(FIXTURE.read_bytes())
        kwargs["stdout"].write("TEST DOUBLE replaying upstream fixture; no caller executed\n")
        return subprocess.CompletedProcess(command, 0)

    def test_command_and_content_addressed_provenance(self):
        with patch.object(amrscan.shutil, "which", return_value="/test/amrfinder"), \
             patch.object(amrscan.subprocess, "run", side_effect=self.mock_caller) as run:
            amrscan.call(self.args)
            metadata = json.loads(self.args.provenance.read_text())
            command = run.call_args.args[0]
        self.assertNotIn("shell", run.call_args.kwargs)
        self.assertIn("--report_all_equal", command)
        self.assertNotIn("--organism", command)
        self.assertNotIn("--plus", command)
        self.assertEqual(command[command.index("--database") + 1], str(self.database))
        self.assertEqual(metadata["input_sha256"], amrscan.sha256(self.args.input))
        self.assertEqual(metadata["raw_sha256"], amrscan.sha256(FIXTURE))
        self.assertEqual(metadata["caller_version"], "TEST-DOUBLE")
        self.assertEqual(metadata["status"], "success")
        self.assertEqual(metadata["database_files_sha256"]["marker"], amrscan.sha256(self.database / "marker"))
        self.assertEqual(metadata["thresholds"]["ident_min"], -1)
        first_id = metadata["database_id"]
        (self.database / "marker").write_text("changed dummy content")
        self.args.organism, self.args.plus = "Escherichia", True
        self.args.database_source = "/original/versioned-db"
        with patch.object(amrscan.shutil, "which", return_value="/test/amrfinder"), \
             patch.object(amrscan.subprocess, "run", side_effect=self.mock_caller) as run:
            amrscan.call(self.args)
        metadata = json.loads(self.args.provenance.read_text())
        self.assertNotEqual(metadata["database_id"], first_id)
        self.assertEqual(metadata["database_source"], "/original/versioned-db")
        self.assertIn("--plus", run.call_args.args[0])
        self.assertIn("Escherichia", run.call_args.args[0])

    def test_failures_propagate_and_preserve_failure_provenance(self):
        for error in (subprocess.CalledProcessError(2, ["amrfinder"]), FileNotFoundError("output")):
            with patch.object(amrscan.shutil, "which", return_value="/test/amrfinder"), \
                 patch.object(amrscan.subprocess, "run", side_effect=[
                     subprocess.CompletedProcess([], 0, stdout="TEST-DOUBLE"), error]), \
                 self.assertRaises(type(error)):
                amrscan.call(self.args)
            metadata = json.loads(self.args.provenance.read_text())
            self.assertEqual(metadata["status"], "failed")
            self.assertNotIn("raw_sha256", metadata)

    def test_invalid_settings_and_missing_tools(self):
        for key, value in [("threads", 0), ("identity", 2), ("identity", float("nan")),
                           ("coverage", -1), ("coverage", float("inf")),
                           ("database", self.directory / "absent")]:
            with self.subTest(key=key, value=value):
                args = argparse.Namespace(**vars(self.args))
                setattr(args, key, value)
                with self.assertRaises(ValueError):
                    amrscan.call(args)
        with patch.object(amrscan.shutil, "which", return_value=None), self.assertRaisesRegex(ValueError, "not found"):
            amrscan.call(self.args)
        with patch.object(amrscan.shutil, "which", return_value="/test/amrfinder"), \
             patch.object(amrscan.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, stdout="")), \
             self.assertRaisesRegex(ValueError, "version"):
            amrscan.call(self.args)
        (self.database / "marker").unlink()
        with self.assertRaisesRegex(ValueError, "empty"):
            amrscan.call(self.args)

    def test_cli_calls_runner(self):
        with patch.object(amrscan, "call") as run:
            amrscan.main(["call", "--sample", "x", "--input", "x.fa", "--database", "db",
                          "--raw", "raw.tsv", "--provenance", "p.json", "--log", "log"])
        self.assertEqual(run.call_args.args[0].sample, "x")


if __name__ == "__main__":
    unittest.main()
