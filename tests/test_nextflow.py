"""Execute real Nextflow locally without AMRFinderPlus or an AMR database."""

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not installed")
class NextflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="amrscan-nextflow-")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.env = dict(os.environ, NXF_OFFLINE="true", NXF_DISABLE_CHECK_LATEST="true",
                        NXF_ANSI_LOG="false", PATH=str(ROOT / "bin") + os.pathsep + os.environ["PATH"])
        self.db = self.directory / "structure-only-db"
        self.db.mkdir()  # Path validation only; never presented as a functional database.

    def run_nf(self, *args, success=True):
        completed = subprocess.run(["nextflow", "-log", str(self.directory / "nextflow.log"), *map(str, args)],
            cwd=self.directory, env=self.env, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=120)
        if success:
            self.assertEqual(completed.returncode, 0, completed.stdout)
        else:
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
        return completed.stdout

    def validate(self, *args, success=True):
        return self.run_nf("run", ROOT / "main.nf", "--validate_only", "--amrfinder_db", self.db,
                           *args, success=success)

    def test_config_and_module_compilation(self):
        config = self.run_nf("config", ROOT, "-flat")
        self.assertIn("params.amrfinder_identity = -1", config)
        self.assertIn("process.executor = 'local'", config)
        help_text = self.run_nf("run", ROOT / "main.nf", "--help")
        self.assertIn("Legacy: workflow/AMRScan.nf", help_text)

    def test_multiple_assemblies_and_gzip(self):
        (self.directory / "sample_a.fna").write_bytes((ROOT / "tests/fixtures/tiny.fna").read_bytes())
        with gzip.open(self.directory / "sample_b.fna.gz", "wb") as handle:
            handle.write((ROOT / "tests/fixtures/tiny.fna").read_bytes())
        output = self.validate("--input", self.directory / "sample_*.fna*")
        self.assertIn("Validated sample_a:", output)
        self.assertIn("Validated sample_b:", output)
        self.assertNotIn("process >", output)

    def test_input_and_parameter_failures(self):
        tiny = ROOT / "tests/fixtures/tiny.fna"
        invalid = self.directory / "reads.fastq"
        invalid.write_text("@read\nACGT\n+\nIIII\n")
        unsafe = self.directory / "unsafe sample.fna"
        unsafe.write_bytes(tiny.read_bytes())
        for args, message in [([], "Provide --input"),
                              (["--input", self.directory / "missing.fna"], "No such file"),
                              (["--input", invalid], "Expected assembled FASTA"),
                              (["--input", unsafe], "Unsafe sample ID"),
                              (["--input", tiny, "--threads", "0"], "positive integer"),
                              (["--input", tiny, "--amrfinder_identity", "2"], "must be -1"),
                              (["--input", tiny, "--amrfinder_coverage", "NaN"], "between 0 and 1")]:
            with self.subTest(args=args):
                self.assertIn(message, self.validate(*args, success=False))

    def test_duplicate_sample_ids_fail_before_scheduling(self):
        for name in ("same.fa", "same.fna"):
            (self.directory / name).write_bytes((ROOT / "tests/fixtures/tiny.fna").read_bytes())
        output = self.validate("--input", self.directory / "same.*", success=False)
        self.assertIn("Duplicate sample ID: same", output)

    def test_harmonizer_process_runs_on_real_upstream_fixture(self):
        raw = ROOT / "tests/fixtures/amrfinder_v4.0.23.tsv"
        provenance = self.directory / "provenance.json"
        provenance.write_text(json.dumps({"sample_id": "upstream", "tool": "AMRFinderPlus",
            "status": "upstream_fixture", "caller_version": "NA", "database_id": "NA",
            "thresholds": {}, "raw_sha256": hashlib.sha256(raw.read_bytes()).hexdigest()}))
        self.run_nf("run", ROOT / "tests/workflows/harmonize.nf", "--raw", raw,
                    "--provenance", provenance, "--outdir", self.directory / "results")
        with (self.directory / "results/harmonized/upstream.amr.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 24)
        self.assertEqual(rows[1]["allele"], "blaTEM-156")
        self.assertEqual(rows[22]["coverage"], "NA")

    def test_caller_missing_is_a_workflow_failure(self):
        (self.db / "marker").write_text("NOT AN AMR DATABASE")
        output = self.run_nf("run", ROOT / "main.nf", "--input", ROOT / "tests/fixtures/tiny.fna",
            "--amrfinder_db", self.db, "--amrfinder_bin", "missing-amrfinder-test-executable",
            "--outdir", self.directory / "results", "--threads", "1", success=False)
        self.assertIn("AMRFinderPlus executable not found", output)
        self.assertFalse((self.directory / "results/harmonized/tiny.amr.tsv").exists())

    @unittest.skipUnless(shutil.which("blastx") and shutil.which("makeblastdb") and shutil.which("Rscript"),
                         "Legacy smoke requires BLAST+ and R with dplyr/magrittr")
    def test_legacy_fasta_no_hit_smoke(self):
        deps = subprocess.run(["Rscript", "-e", "stopifnot(requireNamespace('dplyr', quietly=TRUE), requireNamespace('magrittr', quietly=TRUE))"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if deps.returncode:
            self.skipTest("Legacy R packages are unavailable")
        database = self.directory / "synthetic-protein.fa"
        database.write_text(">synthetic_not_amr\nMKWVTFISLLFLFSSAYSAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
        self.run_nf("run", ROOT / "workflow/AMRScan.nf", "--input", ROOT / "tests/fixtures/tiny.fna",
                    "--card_db", database, "--threads", "1", "--outdir", self.directory / "legacy")
        final = self.directory / "legacy/final"
        for name in ("tiny_AMR_hits_summary.csv", "tiny_AMR_hits_detailed.csv"):
            with (final / name).open() as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(list(reader), [])
                self.assertIn("Bitscore", reader.fieldnames)


if __name__ == "__main__":
    unittest.main()
