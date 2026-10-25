"""Offline tests for assembly retrieval and identity/provenance validation."""
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile
from urllib.error import URLError

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("genome_dataset", ROOT / "bin/genome_dataset.py")
genome = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genome)
ROW = {"cohort_id": "test", "isolate_id": "PDT000000001.1", "pdt_accession": "PDT000000001.1",
       "biosample_accession": "SAMN00000001", "assembly_accession": "GCF_000000001.1",
       "organism": "Escherichia coli", "assembly_source": "NCBI Assembly", "ast_available": "true",
       "selection_reason": "mocked test"}
FASTA = b">contig1\nACGTACGT\n"
REPORT = {"assemblyInfo": {"assemblyAccession": ROW["assembly_accession"]},
          "biologicalMaterial": {"biosample": {"accession": ROW["biosample_accession"]}},
          "organism": {"organismName": ROW["organism"]}}
CURRENT_REPORT = {"assemblyInfo": {"biosample": {"accession": ROW["biosample_accession"],
                   "description": {"organism": {"organismName": ROW["organism"]}}}},
                  "currentAccession": ROW["assembly_accession"]}


def package(report=REPORT, fasta=FASTA, include_fasta=True):
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as z:
        members = {}
        if include_fasta:
            members[f"ncbi_dataset/data/{ROW['assembly_accession']}/GCF_genomic.fna"] = fasta
        members["ncbi_dataset/data/assembly_data_report.jsonl"] = (json.dumps(report) + "\n").encode()
        members["ncbi_dataset/data/sequence_report.jsonl"] = b"{}\n"
        members["ncbi_dataset/data/dataset_catalog.json"] = b"{}\n"
        checksums = "".join(f"{hashlib.md5(data).hexdigest()}  {name}\n" for name, data in members.items())
        z.writestr("md5sum.txt", checksums)
        for name, data in members.items():
            z.writestr(name, data)
    return stream.getvalue()


class GenomeDatasetTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root / "cohort.tsv"
        self.write_rows([ROW])
        self.out = self.root / "genomes"

    def write_rows(self, rows):
        with self.manifest.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=genome.FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    def test_manifest_validation_and_duplicates(self):
        self.assertEqual(genome.read_manifest(self.manifest), [ROW])
        bad = dict(ROW, assembly_accession="bad")
        self.write_rows([bad])
        with self.assertRaisesRegex(ValueError, "Assembly accession"):
            genome.read_manifest(self.manifest)
        self.write_rows([ROW, ROW])
        with self.assertRaisesRegex(ValueError, "duplicate assembly"):
            genome.read_manifest(self.manifest)
        self.write_rows([ROW, dict(ROW, assembly_accession="GCF_000000002.1")])
        with self.assertRaisesRegex(ValueError, "duplicate isolate"):
            genome.read_manifest(self.manifest)

    def test_manifest_rejects_wrong_schema_missing_values_ast_and_species(self):
        self.manifest.write_text("wrong\nvalue\n")
        with self.assertRaisesRegex(ValueError, "columns"):
            genome.read_manifest(self.manifest)
        for changes, message in [({"ast_available": "false"}, "AST"),
                                 ({"organism": "Shigella flexneri"}, "organism"),
                                 ({"biosample_accession": "NA"}, "BioSample")]:
            self.write_rows([dict(ROW, **changes)])
            with self.assertRaisesRegex(ValueError, message):
                genome.read_manifest(self.manifest)
        self.write_rows([])
        with self.assertRaisesRegex(ValueError, "no validation"):
            genome.read_manifest(self.manifest)

    def test_fasta_validation(self):
        self.assertEqual(genome.parse_fasta(FASTA), (1, 8))
        for data in (b"", b">\nACG\n", b"ACG\n", b">x\n!\n", b">x\n"):
            with self.subTest(data=data), self.assertRaises(ValueError):
                genome.parse_fasta(data)
        with self.assertRaises(UnicodeDecodeError):
            genome.parse_fasta(b"\xff")

    def test_url_and_zip_member_checks(self):
        with self.assertRaisesRegex(ValueError, "HTTPS"):
            genome.request_bytes("http://bad", 10)
        with self.assertRaisesRegex(ValueError, "exactly one"):
            genome.package_fasta(package(include_fasta=False), ROW["assembly_accession"])
        with self.assertRaises(zipfile.BadZipFile):
            genome.package_fasta(b"invalid", ROW["assembly_accession"])

    def test_request_transport_success_and_byte_limit(self):
        class Response(io.BytesIO):
            def geturl(self):
                return "https://ncbi.example/final"
            def __enter__(self):
                return self
            def __exit__(self, *args):
                self.close()
        with patch.object(genome, "urlopen", return_value=Response(b"abc")) as opened:
            self.assertEqual(genome.request_bytes("https://ncbi.example/file", 3),
                             (b"abc", "https://ncbi.example/final"))
            self.assertIn("AMRScan-genome-dataset", opened.call_args.args[0].get_header("User-agent"))
        with patch.object(genome, "urlopen", return_value=Response(b"abcd")):
            with self.assertRaisesRegex(ValueError, "exceeds"):
                genome.request_bytes("https://ncbi.example/file", 3)

    def test_invalid_metadata_and_provenance_sets(self):
        for report, expected in [(dict(REPORT, assemblyInfo={"assemblyAccession": "GCF_999999999.1"}), "identity mismatch"),
                                 (dict(REPORT, organism={"organismName": "Salmonella enterica"}), "organism mismatch")]:
            with patch.object(genome, "request_bytes", return_value=(package(report), "https://ncbi.example")):
                with self.assertRaisesRegex(ValueError, expected):
                    genome.download(self.manifest, self.out)
        with patch.object(genome, "request_bytes", return_value=(package(), "https://ncbi.example")):
            genome.download(self.manifest, self.out)
        provenance = self.out / "provenance.json"
        payload = json.loads(provenance.read_text())
        payload["assemblies"].append(payload["assemblies"][0])
        provenance.write_text(json.dumps(payload))
        with self.assertRaisesRegex(ValueError, "duplicate assembly"):
            genome.validate(self.manifest, self.out)

    def test_cli_success_and_failures(self):
        with patch.object(genome, "download") as download:
            self.assertEqual(genome.main(["download", "--manifest", str(self.manifest), "--outdir", str(self.out)]), 0)
            download.assert_called_once()
        with patch("sys.stderr", new_callable=io.StringIO) as stderr:
            self.assertEqual(genome.main(["validate", "--manifest", str(self.manifest), "--genomes", str(self.out)]), 2)
            self.assertIn("error:", stderr.getvalue())

    def test_download_records_identity_hash_and_validate(self):
        with patch.object(genome, "request_bytes", return_value=(package(), "https://ncbi.example/package")):
            genome.download(self.manifest, self.out)
        item = json.loads((self.out / "provenance.json").read_text())["assemblies"][0]
        path = self.out / item["file"]
        self.assertEqual(path.read_bytes(), FASTA)
        self.assertEqual(item["sha256"], hashlib.sha256(FASTA).hexdigest())
        self.assertTrue(item["source_url"].startswith("https://"))
        self.assertEqual(genome.validate(self.manifest, self.out), 1)
        path.write_bytes(b">x\nTTTT\n")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            genome.validate(self.manifest, self.out)

    def test_current_ncbi_report_schema_and_manual_package_import(self):
        archive = package(report=CURRENT_REPORT)
        package_path = self.root / "assembly.zip"
        package_path.write_bytes(archive)
        packages = self.root / "packages.tsv"
        packages.write_text(
            "assembly_accession\tpackage_path\tsource_url\tretrieved_at\n"
            f"{ROW['assembly_accession']}\tassembly.zip\thttps://ncbi.example/archive\t2026-10-11T06:31:25Z\n",
            encoding="utf-8")
        self.assertEqual(genome.import_packages(self.manifest, packages, self.out), 1)
        self.assertEqual(genome.validate(self.manifest, self.out), 1)
        payload = json.loads((self.out / "provenance.json").read_text())
        item = payload["assemblies"][0]
        self.assertEqual(item["biosample_accession"], ROW["biosample_accession"])
        self.assertEqual(item["package_sha256"], hashlib.sha256(archive).hexdigest())
        self.assertEqual(item["retrieved_at"], "2026-10-11T06:31:25Z")

    def test_manual_import_requires_complete_packages_and_safe_metadata(self):
        archive = package(report=CURRENT_REPORT)
        package_path = self.root / "assembly.zip"
        package_path.write_bytes(archive)
        packages = self.root / "packages.tsv"
        valid = ("assembly_accession\tpackage_path\tsource_url\tretrieved_at\n"
                 f"{ROW['assembly_accession']}\tassembly.zip\thttps://ncbi.example/archive\t2026-10-11T06:31:25Z\n")
        for text, error in [(valid.replace("https://", "http://"), "HTTPS"),
                            (valid.replace("2026-10-11T06:31:25Z", "2026-10-11T06:31:25"), "timezone"),
                            (valid.replace("GCF_000000001.1", "GCF_000000002.1"), "match")]:
            packages.write_text(text, encoding="utf-8")
            with self.subTest(error=error), self.assertRaisesRegex(ValueError, error):
                genome.import_packages(self.manifest, packages, self.out)

    def test_package_import_rejects_bad_ncbi_md5_manifest(self):
        archive = package(report=CURRENT_REPORT)
        broken = io.BytesIO()
        with zipfile.ZipFile(io.BytesIO(archive)) as source, zipfile.ZipFile(broken, "w") as target:
            for info in source.infolist():
                payload = source.read(info.filename)
                if info.filename == "md5sum.txt":
                    payload = payload.replace(hashlib.md5(FASTA).hexdigest().encode(), b"0" * 32)
                target.writestr(info.filename, payload)
        package_path = self.root / "broken.zip"
        package_path.write_bytes(broken.getvalue())
        packages = self.root / "packages.tsv"
        packages.write_text(
            "assembly_accession\tpackage_path\tsource_url\tretrieved_at\n"
            f"{ROW['assembly_accession']}\tbroken.zip\thttps://ncbi.example/archive\t2026-10-11T06:31:25Z\n",
            encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
            genome.import_packages(self.manifest, packages, self.out)

    def test_manual_import_cli(self):
        with patch.object(genome, "import_packages", return_value=1) as importer:
            result = genome.main(["import-packages", "--manifest", str(self.manifest),
                                  "--packages-tsv", str(self.root / "packages.tsv"),
                                  "--outdir", str(self.out)])
        self.assertEqual(result, 0)
        importer.assert_called_once()

    def test_download_rejects_identity_fasta_existing_and_partial(self):
        mismatch = dict(REPORT)
        mismatch["biologicalMaterial"] = {"biosample": {"accession": "SAMN99999999"}}
        with patch.object(genome, "request_bytes", return_value=(package(mismatch), "https://ncbi.example")):
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                genome.download(self.manifest, self.out)
        self.assertFalse(list(self.out.glob("*.fna")))
        with patch.object(genome, "request_bytes", return_value=(package(fasta=b">x\n!\n"), "https://ncbi.example")):
            with self.assertRaisesRegex(ValueError, "invalid FASTA"):
                genome.download(self.manifest, self.out)
        with patch.object(genome, "request_bytes", side_effect=OSError("offline")):
            with self.assertRaisesRegex(OSError, "offline"):
                genome.download(self.manifest, self.out)
        self.out.mkdir(exist_ok=True)
        existing = self.out / f"{ROW['isolate_id']}_{ROW['assembly_accession']}.fna"
        existing.write_bytes(FASTA)
        with self.assertRaisesRegex(ValueError, "refusing to replace"):
            genome.download(self.manifest, self.out)

    def test_download_manifest_hash_and_biosample_validation(self):
        with patch.object(genome, "request_bytes", return_value=(package(), "https://ncbi.example")):
            genome.download(self.manifest, self.out)
        self.write_rows([dict(ROW, selection_reason="changed")])
        with self.assertRaisesRegex(ValueError, "manifest hash"):
            genome.validate(self.manifest, self.out)
        self.write_rows([ROW])
        provenance = self.out / "provenance.json"
        payload = json.loads(provenance.read_text())
        payload["assemblies"][0]["biosample_accession"] = "SAMN99999999"
        provenance.write_text(json.dumps(payload))
        with self.assertRaisesRegex(ValueError, "biosample_accession mismatch"):
            genome.validate(self.manifest, self.out)


if __name__ == "__main__":
    unittest.main()
