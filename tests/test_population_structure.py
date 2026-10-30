"""Offline tests for E. coli population structure metadata (synthetic genomes and tool output)."""

from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile
from urllib.error import URLError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
spec = importlib.util.spec_from_file_location("population_structure", ROOT / "bin/population_structure.py")
pop = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pop)
CONFIG = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/population_config.json").read_text())
LOCI = CONFIG["mlst"]["loci"]
CAND_FIELDS = ("cohort_id isolate_id pdt_accession biosample_accession assembly_accession organism strain "
               "country collection_year source_category submitter project snp_clusters_json").split()
EP_FIELDS = "isolate_id antimicrobial phenotype_endpoint eligible_overall".split()
NA = pop.NA


def contig(n, gc=True):
    return ("GC" if gc else "AT") * (n // 2)


def mlst_json(st="131", **alleles):
    values = {locus: "1" for locus in LOCI}
    values.update(alleles)
    return json.dumps([{"scheme": "ecoli_achtman_4", "sequence_type": st, "alleles": values, "filename": "g"}])


class PopulationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.cohort = self.root / "cohort"
        self.cache = self.root / "cache"
        self.config_path = self.root / "config.json"
        config = json.loads(json.dumps(CONFIG))
        config["genome_qc"].update({"fail_length_min": 100, "fail_length_max": 10000, "warn_length_min": 200,
                                   "warn_length_max": 5000, "warn_contigs_max": 3, "warn_n50_min": 50})
        config["review_flags"].update({"min_lineage_isolates": 2, "min_resistant_for_concentration": 2,
                                       "large_snp_cluster_isolates": 2})
        self.config_path.write_text(json.dumps(config))
        self.config = config

    # ------------------------------------------------------------ fixtures
    def candidate(self, n, **changes):
        row = {"cohort_id": "ecoli_research_v1", "isolate_id": f"PDT00000000{n}.1", "pdt_accession": f"PDT00000000{n}.1",
               "biosample_accession": f"SAMN0000000{n}", "assembly_accession": f"GCA_00000000{n}.1",
               "organism": "Escherichia coli", "strain": f"s{n}", "country": "USA", "collection_year": "2020",
               "source_category": "human", "submitter": "Lab A", "project": "PRJNA1",
               "snp_clusters_json": '["PDS000000001.1","PDS000000001.2"]'}
        row.update(changes)
        return row

    def write_cohort(self, candidates, endpoints, manifest_changes=None):
        self.cohort.mkdir(exist_ok=True)
        for name, fields, rows in (("candidate_isolates.tsv", CAND_FIELDS, candidates),
                                   ("isolate_drug_endpoints.tsv", EP_FIELDS, endpoints)):
            with (self.cohort / name).open("w", newline="") as h:
                w = csv.DictWriter(h, fields, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        manifest = {"cohort_id": "ecoli_research_v1", "candidate_isolate_count": len(candidates),
                    "eligible_binary_rs_endpoint_count": sum(e["eligible_overall"] == "true" for e in endpoints),
                    "source_snapshot_id": "snapshot", "code_commit": "abc"}
        manifest.update(manifest_changes or {})
        (self.cohort / "cohort_manifest.json").write_text(json.dumps(manifest))
        sums = [f"{pop.sha256(self.cohort / n)}  {n}" for n in pop.COHORT_FILES]
        (self.cohort / "SHA256SUMS").write_text("\n".join(sums) + "\n")

    def cache_genome(self, assembly, text=None):
        entry = self.cache / assembly
        entry.mkdir(parents=True, exist_ok=True)
        (entry / "genome.fna").write_text(text or f">c1\n{contig(400)}\n>c2\n{contig(200, False)}\n")
        (entry / "provenance.json").write_text(json.dumps({
            "assembly_accession": assembly, "sha256": pop.sha256(entry / "genome.fna"),
            "retrieved_at": "2026-10-11T00:00:00Z"}))

    def annotation(self, assembly, tool, raw, status="success", returncode=0):
        sha = json.loads((self.cache / assembly / "provenance.json").read_text())["sha256"]
        info = {"tool": "mlst" if tool == "mlst" else "ezclermont", "version": "fixture 1", "database_version": "db1",
                "scheme_profile_sha256": "x"}
        data = {"input_sha256": sha, "tool_info": info, "status": status, "returncode": returncode, "raw_output": raw}
        if status == "success":
            data["parsed"] = (pop.parse_mlst(raw, self.config) if tool == "mlst" else pop.parse_phylogroup(raw, self.config))
        (self.cache / assembly / f"{tool}.json").write_text(json.dumps(data))

    def standard_population(self):
        cands = [self.candidate(1), self.candidate(2, country="Canada", collection_year="2021"),
                 self.candidate(3, country="", collection_year="NA", source_category="unknown",
                                snp_clusters_json="[]"),
                 self.candidate(4, snp_clusters_json='["PDS000000002.1","PDS000000009.1"]'),
                 self.candidate(5)]
        eps = []
        for c, ph in zip(cands, "RRSSR"):
            eps.append({"isolate_id": c["isolate_id"], "antimicrobial": "ampicillin", "phenotype_endpoint": ph,
                        "eligible_overall": "true"})
            eps.append({"isolate_id": c["isolate_id"], "antimicrobial": "meropenem", "phenotype_endpoint": "ND",
                        "eligible_overall": "false"})
        self.write_cohort(cands, eps)
        for i, c in enumerate(cands[:4], 1):
            self.cache_genome(c["assembly_accession"])
        self.annotation("GCA_000000001.1", "mlst", mlst_json("131"))
        self.annotation("GCA_000000002.1", "mlst", mlst_json("131"))
        self.annotation("GCA_000000003.1", "mlst", mlst_json("-", adk="~5"))
        self.annotation("GCA_000000004.1", "mlst", mlst_json("10"))
        self.annotation("GCA_000000001.1", "phylogroup", "g\tB2\n")
        self.annotation("GCA_000000002.1", "phylogroup", "g\tU\n")
        self.annotation("GCA_000000003.1", "phylogroup", "", status="failed", returncode=1)
        return cands, eps

    def build(self, name="out", **kw):
        return pop.build(self.cohort, self.cache, self.root / name, self.config_path, **kw)

    def read(self, name, out="out"):
        with (self.root / out / name).open(newline="") as h:
            return list(csv.DictReader(h, delimiter="\t"))

    # ------------------------------------------------------------ genome QC
    def test_fasta_statistics(self):
        path = self.root / "g.fna"
        path.write_text(">a\nGGGG\nCCNN\n>b\nAATT\n>c\nAC\n")
        stats = pop.fasta_stats(path)
        self.assertEqual(stats["assembly_length"], 14)
        self.assertEqual(stats["contig_count"], 3)
        self.assertEqual(stats["n50"], 8)
        self.assertEqual(stats["longest_contig"], 8)
        self.assertEqual(stats["gc_fraction"], round(7 / 12, 6))
        self.assertEqual(stats["ambiguous_fraction"], round(2 / 14, 6))

    def test_malformed_and_empty_fasta(self):
        for text, message in (("", "no sequence"), ("ACGT\n", "invalid FASTA"), (">a\n>b\nAC\n", "empty FASTA record"),
                              (">\nAC\n", "empty FASTA header"), (">a\nAC\n>b\n", "end of file"), (">a\nAC1\n", "invalid")):
            path = self.root / "bad.fna"
            path.write_text(text)
            with self.assertRaisesRegex(ValueError, message):
                pop.fasta_stats(path)

    def test_qc_status_pass_warning_fail(self):
        q = self.config
        good = {"assembly_length": 1000, "contig_count": 2, "n50": 500, "ambiguous_fraction": 0.0}
        self.assertEqual(pop.qc_status(good, q), ("pass", ["qc_pass"]))
        status, flags = pop.qc_status(dict(good, assembly_length=6000, contig_count=9, n50=10,
                                           ambiguous_fraction=0.01), q)
        self.assertEqual(status, "warning")
        self.assertEqual(flags, ["assembly_size_warning", "excessive_ambiguous_bases", "excessive_contigs", "low_n50"])
        self.assertEqual(pop.qc_status(dict(good, assembly_length=50), q)[0], "fail")
        self.assertEqual(pop.qc_status(dict(good, ambiguous_fraction=0.5), q),
                         ("fail", ["excessive_ambiguous_bases_fail"]))
        self.assertEqual(pop.qc_status(None, q, "assembly_missing"), ("fail", ["assembly_missing"]))

    def test_genome_qc_missing_corrupt_and_malformed(self):
        self.assertEqual(pop.genome_qc(self.cache, "GCA_000000001.1", self.config)["genome_qc_flags"], "assembly_missing")
        self.cache_genome("GCA_000000001.1")
        (self.cache / "GCA_000000001.1" / "genome.fna").write_text(">x\nAC\n")
        self.assertEqual(pop.genome_qc(self.cache, "GCA_000000001.1", self.config)["genome_qc_flags"],
                         "assembly_checksum_failed")
        self.cache_genome("GCA_000000002.1", ">x\nAC1\n")
        row = pop.genome_qc(self.cache, "GCA_000000002.1", self.config)
        self.assertEqual((row["genome_qc_status"], row["genome_qc_flags"]), ("fail", "malformed_fasta"))

    # ------------------------------------------------------------ MLST / phylogroup
    def test_parse_mlst_statuses(self):
        self.assertEqual(pop.parse_mlst(mlst_json("131"), self.config)["mlst_st"], "131")
        self.assertIn("adk(1)", pop.parse_mlst(mlst_json("131"), self.config)["mlst_profile"])
        cases = {mlst_json("-", adk="~5"): "novel_or_partial_alleles", mlst_json("-", adk="5?"): "novel_or_partial_alleles",
                 mlst_json("-", adk="-"): "missing_alleles", mlst_json("-", adk="1,2"): "multiple_alleles",
                 mlst_json("-"): "novel_profile"}
        for text, status in cases.items():
            parsed = pop.parse_mlst(text, self.config)
            self.assertEqual((parsed["mlst_status"], parsed["mlst_st"]), (status, "NA"))

    def test_parse_mlst_rejects_malformed_and_other_schemes(self):
        for text in ("not json", "[]", json.dumps([{"scheme": "ecoli", "sequence_type": "1", "alleles": {"x": "1"}}])):
            with self.assertRaises(ValueError):
                pop.parse_mlst(text, self.config)

    def test_parse_phylogroup(self):
        self.assertEqual(pop.parse_phylogroup("x\tB2\n", self.config)["phylogroup_status"], "assigned")
        self.assertEqual(pop.parse_phylogroup("x\tcryptic\n", self.config)["phylogroup_status"], "cryptic_clade")
        self.assertEqual(pop.parse_phylogroup("x\tU/cryptic\n", self.config),
                         {"phylogroup": "NA", "phylogroup_status": "unresolved", "phylogroup_raw": "U/cryptic"})
        self.assertEqual(pop.parse_phylogroup("x\tEC_control_fail\n", self.config)["phylogroup_status"], "control_failed")
        with self.assertRaisesRegex(ValueError, "malformed"):
            pop.parse_phylogroup("", self.config)

    def test_annotate_runs_tool_records_provenance_and_reuses_result(self):
        self.write_cohort([self.candidate(1), self.candidate(2)], [])
        self.cache_genome("GCA_000000001.1")
        calls = []

        def runner(cmd, **kw):
            calls.append(cmd)
            Path(cmd[cmd.index("--json") + 1]).write_text(mlst_json("73"))
            return subprocess.CompletedProcess(cmd, 0, "", "")
        info = {"tool": "mlst", "version": "mlst 9", "database_version": "db", "scheme_profile_sha256": "s"}
        cohort = pop.load_cohort(self.cohort)
        result = pop.annotate("mlst", cohort, self.cache, self.root / "bin", self.config, runner=runner, info=info)
        self.assertEqual(result["status_counts"], {"missing_genome": 1, "success": 1})
        saved = json.loads((self.cache / "GCA_000000001.1" / "mlst.json").read_text())
        self.assertEqual(saved["parsed"]["mlst_st"], "73")
        self.assertEqual(saved["input_sha256"], pop.sha256(self.cache / "GCA_000000001.1" / "genome.fna"))
        self.assertIn("--scheme", saved["command"])
        self.assertEqual(pop.annotate("mlst", cohort, self.cache, self.root / "bin", self.config, runner=runner,
                                      info=info)["status_counts"]["cached"], 1)
        self.assertEqual(len(calls), 1)

    def test_annotate_failure_and_phylogroup(self):
        self.cache_genome("GCA_000000001.1")
        fail = lambda cmd, **kw: subprocess.CompletedProcess(cmd, 2, "", "boom")
        status = pop.annotate_one("phylogroup", "GCA_000000001.1", self.cache, self.root, {"tool": "ezclermont"},
                                  self.config, fail)
        self.assertEqual(status, "failed")
        ok = lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, "g\tA\n", "")
        status = pop.annotate_one("phylogroup", "GCA_000000001.1", self.cache, self.root, {"tool": "ezclermont", "v": 2},
                                  self.config, ok)
        self.assertEqual(status, "success")
        self.assertEqual(json.loads((self.cache / "GCA_000000001.1" / "phylogroup.json").read_text())["parsed"]["phylogroup"], "A")

    def test_tool_info_reads_versions_and_scheme(self):
        env_bin = self.root / "env/bin"
        scheme = self.root / "env/db/pubmlst/ecoli_achtman_4"
        scheme.mkdir(parents=True)
        env_bin.mkdir(parents=True)
        (scheme / "ecoli_achtman_4.txt").write_text("ST\tadk\n")
        (scheme / "ecoli_achtman_4_info.json").write_text('{"download_date": "2026-03-11"}')
        (self.root / "env/db/VERSION.txt").write_text("2025-12-29\n")
        fake = lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, f"{Path(cmd[0]).name} 1.0\n", "")
        with patch.object(pop.subprocess, "run", fake):
            info = pop.tool_info("mlst", env_bin, self.config)
            self.assertEqual((info["version"], info["database_version"]), ("mlst 1.0", "2025-12-29"))
            self.assertEqual(info["scheme_info"]["download_date"], "2026-03-11")
            self.assertEqual(pop.tool_info("phylogroup", env_bin, self.config)["tool"], "ezclermont")
            with self.assertRaises(ValueError):
                pop.tool_info("other", env_bin, self.config)
            (scheme / "ecoli_achtman_4.txt").unlink()
            with self.assertRaisesRegex(ValueError, "not installed"):
                pop.tool_info("mlst", env_bin, self.config)

    # ------------------------------------------------------------ SNP clusters
    def test_snp_cluster_preservation_and_missingness(self):
        self.assertEqual(pop.snp_cluster({"snp_clusters_json": '["PDS000000001.1","PDS000000001.2"]'}),
                         ("PDS000000001", "assigned", ["PDS000000001.1", "PDS000000001.2"]))
        self.assertEqual(pop.snp_cluster({"snp_clusters_json": "[]"})[:2], ("NA", "missing"))
        self.assertEqual(pop.snp_cluster({"snp_clusters_json": '["PDS1.1","PDS2.1"]'})[:2], ("NA", "conflicting"))
        with self.assertRaisesRegex(ValueError, "malformed"):
            pop.snp_cluster({"snp_clusters_json": "[", "isolate_id": "x"})

    # ------------------------------------------------------------ genome fetch
    def package(self, assembly="GCA_000000001.1", biosample="SAMN00000001"):
        stream = io.BytesIO()
        members = {f"ncbi_dataset/data/{assembly}/{assembly}_genomic.fna": f">c\n{contig(400)}\n".encode(),
                   "ncbi_dataset/data/assembly_data_report.jsonl": (json.dumps(
                       {"accession": assembly, "currentAccession": assembly, "organism": {"organismName": "Escherichia coli"},
                        "assemblyInfo": {"biosample": {"accession": biosample}}}) + "\n").encode()}
        with zipfile.ZipFile(stream, "w") as z:
            z.writestr("md5sum.txt", "".join(f"{hashlib.md5(d).hexdigest()}  {n}\n" for n, d in members.items()))
            for n, d in members.items():
                z.writestr(n, d)
        return stream.getvalue()

    def test_fetch_is_verified_resumable_and_retries(self):
        self.write_cohort([self.candidate(1), self.candidate(2)], [])
        cohort = pop.load_cohort(self.cohort)
        attempts = Counter()

        def opener(url, limit):
            assembly = url.split("/accession/")[1].split("/")[0]
            attempts[assembly] += 1
            if assembly == "GCA_000000002.1" and attempts[assembly] == 1:
                raise URLError("transient")
            return self.package(assembly, "SAMN0000000" + assembly[-3]), url
        summary = pop.fetch_genomes(cohort, self.cache, workers=2, opener=opener, pause=0)
        self.assertEqual(summary["status_counts"], {"downloaded": 2})
        self.assertIsNotNone(pop.cached_genome(self.cache, "GCA_000000002.1"))
        self.assertEqual(pop.fetch_genomes(cohort, self.cache, opener=opener, pause=0)["status_counts"], {"cached": 2})
        (self.cache / "GCA_000000001.1" / "genome.fna").write_text(">bad\nAC\n")
        self.assertEqual(pop.fetch_genomes(cohort, self.cache, opener=opener, pause=0, limit=1)["status_counts"],
                         {"downloaded": 1})
        self.assertEqual(len((self.cache / "fetch_log.jsonl").read_text().splitlines()), 3)
        with self.assertRaisesRegex(ValueError, "workers"):
            pop.fetch_genomes(cohort, self.cache, workers=9)

    def test_fetch_rejects_identity_mismatch_without_partial_files(self):
        self.write_cohort([self.candidate(1)], [])
        opener = lambda url, limit: (self.package("GCA_000000001.1", "SAMN09999999"), url)
        summary = pop.fetch_genomes(pop.load_cohort(self.cohort), self.cache, workers=1, opener=opener, pause=0)
        self.assertEqual(summary["status_counts"], {"failed": 1})
        self.assertIn("identity mismatch", summary["failures"][0]["error"])
        self.assertFalse((self.cache / "GCA_000000001.1").exists())
        self.assertEqual([p.name for p in self.cache.iterdir() if p.name != "fetch_log.jsonl"], [])

    # ------------------------------------------------------------ cohort linkage
    def test_cohort_hash_checks_and_counts(self):
        self.standard_population()
        cohort = pop.load_cohort(self.cohort, pop.sha256(self.cohort / "cohort_manifest.json"))
        self.assertEqual(len(cohort["candidates"]), 5)
        with self.assertRaisesRegex(ValueError, "expected value"):
            pop.load_cohort(self.cohort, "0" * 64)
        with (self.cohort / "isolate_drug_endpoints.tsv").open("a") as h:
            h.write("x\ty\tR\ttrue\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            pop.load_cohort(self.cohort)

    def test_cohort_rejects_bad_manifest_count_duplicates_and_malformed_rows(self):
        self.write_cohort([self.candidate(1)], [], {"candidate_isolate_count": 2})
        with self.assertRaisesRegex(ValueError, "manifest isolate count"):
            pop.load_cohort(self.cohort)
        self.write_cohort([self.candidate(1), self.candidate(1)], [])
        with self.assertRaisesRegex(ValueError, "duplicate isolate"):
            pop.load_cohort(self.cohort)
        self.write_cohort([self.candidate(1), self.candidate(2, assembly_accession="GCA_000000001.1")], [])
        self.assertEqual(pop.load_cohort(self.cohort)["duplicate_assemblies"], ["GCA_000000001.1"])
        with self.assertRaisesRegex(ValueError, "assembly shared"):
            self.build()
        (self.root / "bad.tsv").write_text("a\tb\n1\n")
        with self.assertRaisesRegex(ValueError, "malformed row width"):
            pop.read_tsv(self.root / "bad.tsv")
        (self.root / "empty.tsv").write_text("")
        with self.assertRaisesRegex(ValueError, "missing header"):
            pop.read_tsv(self.root / "empty.tsv")

    def test_header_only_and_empty_cohort(self):
        self.write_cohort([], [])
        with self.assertRaisesRegex(ValueError, "no isolates"):
            self.build()

    def test_config_must_keep_human_review(self):
        bad = dict(self.config, feasibility_rules=dict(self.config["feasibility_rules"], decision="feasible"))
        self.config_path.write_text(json.dumps(bad))
        with self.assertRaisesRegex(ValueError, "HUMAN_REVIEW_REQUIRED"):
            pop.load_config(self.config_path)
        self.config_path.write_text("{}")
        with self.assertRaisesRegex(ValueError, "lacks"):
            pop.load_config(self.config_path)

    # ------------------------------------------------------------ build
    def test_population_metadata_lineage_and_eligibility(self):
        self.standard_population()
        summary = self.build()
        rows = {r["isolate_id"]: r for r in self.read("population_metadata.tsv")}
        self.assertEqual(len(rows), 5)
        one, three, five = rows["PDT000000001.1"], rows["PDT000000003.1"], rows["PDT000000005.1"]
        self.assertEqual((one["lineage_id"], one["phylogroup"], one["ncbi_snp_cluster"]), ("ST131", "B2", "PDS000000001"))
        self.assertEqual((three["lineage_id"], three["mlst_status"], three["ncbi_snp_cluster_status"]),
                         ("UNKNOWN", "novel_or_partial_alleles", "missing"))
        self.assertEqual((three["eligible_population_analysis"], three["eligible_lineage_holdout"],
                          three["eligible_country_holdout"], three["eligible_temporal_holdout"],
                          three["eligible_source_holdout"]), ("true", "false", "false", "false", "false"))
        self.assertEqual(rows["PDT000000002.1"]["phylogroup_status"], "unresolved")
        self.assertEqual(rows["PDT000000003.1"]["phylogroup_status"], "failed")
        self.assertEqual(rows["PDT000000004.1"]["phylogroup_status"], "not_run")
        self.assertEqual(rows["PDT000000004.1"]["ncbi_snp_cluster_status"], "conflicting")
        self.assertEqual((five["genome_qc_status"], five["mlst_status"], five["eligible_population_analysis"]),
                         ("fail", "not_run", "false"))
        self.assertEqual(one["binary_rs_endpoint_count"], "1")
        self.assertEqual(one["antimicrobial_count"], "2")
        self.assertEqual(summary["unique_sts"], 2)
        self.assertEqual(summary["singleton_sts"], 1)
        self.assertEqual(summary["unresolved_st_isolates"], 2)
        self.assertEqual(summary["binary_endpoints_population"], 4)
        self.assertTrue(all(summary["reconciliation"].values()))

    def test_lineage_and_cross_tab_summaries(self):
        self.standard_population()
        self.build()
        lineages = {r["lineage_id"]: r for r in self.read("lineage_summary.tsv")}
        self.assertEqual(lineages["ST131"]["isolate_count"], "2")
        self.assertEqual(lineages["ST131"]["country_count"], "2")
        self.assertEqual(lineages["ST131"]["fraction_from_largest_country"], "0.500000")
        self.assertEqual((lineages["ST131"]["earliest_year"], lineages["ST131"]["latest_year"]), ("2020", "2021"))
        self.assertEqual(lineages["UNKNOWN"]["is_lineage"], "false")
        self.assertEqual(lineages["UNKNOWN"]["country_count"], "1")
        amr = {(r["antimicrobial"], r["lineage_id"]): r for r in self.read("lineage_antimicrobial_summary.tsv")}
        self.assertEqual((amr[("ampicillin", "ST131")]["n_R"], amr[("ampicillin", "ST131")]["resistance_fraction"]),
                         ("2", "1.000000"))
        self.assertNotIn(("meropenem", "ST131"), amr)
        country = {(r["lineage_id"], r["country"]): r["isolate_count"] for r in self.read("lineage_country_summary.tsv")}
        self.assertEqual(country[("UNKNOWN", "NA")], "1")
        years = {(r["lineage_id"], r["collection_year"]): r["isolate_count"] for r in self.read("lineage_year_summary.tsv")}
        self.assertEqual(years[("ST131", "2021")], "1")
        sources = {(r["lineage_id"], r["source_category"]): r["isolate_count"] for r in self.read("lineage_source_summary.tsv")}
        self.assertEqual(sources[("UNKNOWN", "unknown")], "1")
        clusters = {r["ncbi_snp_cluster"]: r for r in self.read("snp_cluster_summary.tsv")}
        self.assertEqual(clusters["PDS000000001"]["isolate_count"], "3")
        self.assertEqual(clusters["PDS000000001"]["binary_endpoints_R"], "2")
        projects = self.read("project_summary.tsv")
        self.assertEqual(projects[0]["isolate_count"], "5")

    def test_groups_feasibility_and_review_flags(self):
        self.standard_population()
        self.build()
        groups = {r["isolate_id"]: r for r in self.read("evaluation_groups.tsv")}
        self.assertEqual(groups["PDT000000003.1"]["lineage_group"], "UNKNOWN")
        self.assertEqual(groups["PDT000000003.1"]["country_group"], "NA")
        self.assertEqual(groups["PDT000000001.1"]["snp_cluster_group"], "PDS000000001")
        self.assertNotIn("train", (self.root / "out/evaluation_groups.tsv").read_text())
        feas = {(r["antimicrobial"], r["holdout"]): r for r in self.read("holdout_feasibility.tsv")}
        self.assertEqual(set(h for _, h in feas), {"lineage", "country", "temporal", "source"})
        self.assertTrue(all(r["status"] == "not_feasible" and r["decision"] == "HUMAN_REVIEW_REQUIRED"
                            for r in feas.values()))
        self.assertEqual(feas[("ampicillin", "lineage")]["groups"], "2")
        self.assertEqual(json.loads(feas[("ampicillin", "temporal")]["detail_json"])["cutoff_year"], 2021)
        risks = {(r["risk"], r["subject"]) for r in self.read("confounding_summary.tsv")}
        self.assertIn(("dominant_lineage", "ST131"), risks)
        self.assertIn(("resistance_lineage_concentration", "ampicillin"), risks)
        self.assertIn(("resistance_project_concentration", "ampicillin"), risks)
        self.assertIn(("lineage_project_concentration", "ST131"), risks)
        self.assertIn(("low_mlst_success", "cohort"), risks)

    def test_classify_and_cutoff_rules(self):
        rule = {"min_R": 2, "max_largest_st_fraction": 0.5}
        self.assertEqual(pop.classify({"R": 3, "largest_st_fraction": 0.4}, rule, {"min_R": 1}), "feasible")
        self.assertEqual(pop.classify({"R": 3, "largest_st_fraction": 0.6}, rule, {"min_R": 1}), "marginal")
        self.assertEqual(pop.classify({"R": 0, "largest_st_fraction": 0.6}, rule, {"min_R": 1}), "not_feasible")
        self.assertEqual(pop.recent_cutoff([2010, 2011, 2012, 2013, 2014], 0.8), 2014)
        self.assertIsNone(pop.recent_cutoff([], 0.8))

    def test_review_flags_for_country_source_time_and_clusters(self):
        rows = [{"eligible_population_analysis": "true", "lineage_id": "ST1", "country": "USA", "source_category": "food",
                 "project": f"P{i}", "collection_year": str(2010 + (i > 2) * 10), "ncbi_snp_cluster": "PDS1",
                 "genome_qc_status": "pass", "genome_sha256": "x", "mlst_status": "assigned"} for i in range(4)]
        rows += [dict(rows[0], lineage_id="ST2", collection_year="2010", ncbi_snp_cluster="NA", project="P9")
                 for _ in range(3)]
        binary = [{"antimicrobial": "x", "phenotype_endpoint": "R", "lineage_id": "ST1", "ncbi_snp_cluster": "PDS1",
                   "project": f"P{i}"} for i in range(3)]
        config = dict(self.config, review_flags=dict(self.config["review_flags"], lineage_after_cutoff_fraction=0.6))
        risks = {r["risk"] for r in pop.review_flags(rows, binary, config, 2020)}
        self.assertTrue({"lineage_country_concentration", "lineage_source_concentration", "large_snp_cluster",
                         "resistance_snp_cluster_concentration", "lineage_absent_from_recent_years"} <= risks)
        late = [dict(r, collection_year="2021") for r in rows]
        self.assertIn("lineage_year_concentration", {r["risk"] for r in pop.review_flags(late, [], config, 2020)})

    def test_lineage_components_join_sts_through_clusters(self):
        def row(iid, st, cluster, eligible="true"):
            return {"isolate_id": iid, "lineage_id": st, "ncbi_snp_cluster": cluster,
                    "eligible_population_analysis": eligible}
        rows = [row("a", "ST131", "PDS1"), row("b", "ST14207", "PDS1"), row("c", "ST131", NA),
                row("d", "UNKNOWN", "PDS1"), row("e", "UNKNOWN", "PDS7"), row("f", "UNKNOWN", "PDS7"),
                row("g", "UNKNOWN", NA), row("h", "ST10", NA), row("i", "ST131", "PDS1", eligible="false")]
        comp = pop.lineage_components(rows)
        self.assertEqual({k: comp[k] for k in "abcd"}, dict.fromkeys("abcd", "ST131+ST14207"))
        self.assertEqual((comp["e"], comp["f"], comp["g"], comp["h"]), ("PDS7", "PDS7", NA, "ST10"))
        self.assertNotIn("i", comp)

    def test_evaluation_plan_tiers_and_analyses(self):
        config = json.loads(json.dumps(self.config))
        config["evaluation_design"]["drug_tiers"] = {
            "primary": {"min_R": 2, "min_S": 2, "max_top_project_R_fraction": 0.6, "min_projects_with_both": 2,
                        "min_sts_with_both": 1},
            "secondary": {"min_R": 1, "min_S": 1}, "project_both_min": 1, "st_both_min": 1,
            "excluded": {"x": "confirmation test"}}
        config["evaluation_design"].update({"temporal_min_per_period": 1,
                                            "geography": {"reference_country": "USA", "min_pooled": 1,
                                                          "min_per_country": 1}})

        def ep(drug, ph, project="P1", st="ST1", country="USA", year="2020"):
            return {"antimicrobial": drug, "phenotype_endpoint": ph, "project": project, "lineage_id": st,
                    "country": country, "collection_year": year,
                    "eligible_temporal_holdout": "true" if year != "NA" else "false"}
        binary = [ep("a", "R"), ep("a", "S"), ep("a", "R", "P2", country="Canada", year="2022"),
                  ep("a", "S", "P2", country="Canada", year="2022"),
                  ep("b", "R"), ep("b", "R"), ep("b", "S", year="NA"), ep("b", "S", country=""),
                  ep("c", "R"), ep("c", "R"), ep("x", "R"), ep("x", "S")]
        plan = {p["antimicrobial"]: p for p in pop.evaluation_plan(binary, config, 2021)}
        a, b, c, x = plan["a"], plan["b"], plan["c"], plan["x"]
        self.assertEqual((a["tier"], a["lineage_cv"], a["project_holdout"], a["temporal"], a["geographic_transfer"]),
                         ("primary", "yes", "yes", "yes", "yes"))
        self.assertEqual((a["earlier_R"], a["recent_R"], a["top_project_R_fraction"]), ("1", "1", "0.500000"))
        self.assertEqual(json.loads(a["non_reference_countries_json"]), {"Canada": [1, 1]})
        self.assertEqual((b["tier"], b["lineage_cv"], b["project_holdout"], b["temporal"], b["geographic_transfer"]),
                         ("secondary", "yes", "no", "no", "no"))
        self.assertIn("top_project_R_fraction 1.0 > 0.6", b["tier_reason"])
        self.assertEqual((c["tier"], c["tier_reason"], c["lineage_cv"]), ("excluded", "S 0 < 1", "no"))
        self.assertEqual((x["tier"], x["tier_reason"]), ("excluded", "confirmation test"))
        self.assertTrue(all(p["source_holdout"] == "no" for p in plan.values()))
        self.assertEqual(pop.temporal_cutoff([2010, 2020], {"cutoff_year": 2021}), 2021)
        self.assertEqual(pop.temporal_cutoff([2010, 2020], {"recent_quantile": 0.5}), 2020)

    def test_build_writes_evaluation_design_outputs(self):
        self.standard_population()
        summary = self.build()
        groups = {r["isolate_id"]: r for r in self.read("evaluation_groups.tsv")}
        one, three, five = groups["PDT000000001.1"], groups["PDT000000003.1"], groups["PDT000000005.1"]
        self.assertEqual((one["lineage_component"], one["lineage_cv_test_eligible"], one["split_unit"],
                          one["geography_group"], one["period"], one["project_group"]),
                         ("ST131", "true", "PDS000000001", "reference", "earlier", "PRJNA1"))
        self.assertEqual(groups["PDT000000002.1"]["geography_group"], "non_reference")
        self.assertEqual(groups["PDT000000002.1"]["period"], "recent")
        self.assertEqual((three["lineage_component"], three["lineage_cv_test_eligible"], three["split_unit"],
                          three["geography_group"], three["period"]), (NA, "false", "PDT000000003.1", NA, NA))
        self.assertEqual(five["lineage_component"], NA)
        plan = {r["antimicrobial"]: r for r in self.read("evaluation_plan.tsv")}
        self.assertEqual(set(plan), {"ampicillin"})
        self.assertEqual(plan["ampicillin"]["tier"], "excluded")
        design = summary["evaluation_design"]
        self.assertEqual((design["temporal_cutoff_year"], design["population_without_component"],
                          design["drug_tiers"]), (2021, 1, {"excluded": 1}))
        manifest = json.loads((self.root / "out/population_manifest.json").read_text())
        self.assertEqual(manifest["evaluation_design_record"], self.config["evaluation_design"]["decision_record"])

    def test_deterministic_rebuild_and_manifest(self):
        self.standard_population()
        self.build("a")
        self.build("b")
        for path in (self.root / "a").iterdir():
            self.assertEqual(path.read_bytes(), (self.root / "b" / path.name).read_bytes(), path.name)
        manifest = json.loads((self.root / "a/population_manifest.json").read_text())
        self.assertEqual(manifest["pr5_cohort_manifest_sha256"], pop.sha256(self.cohort / "cohort_manifest.json"))
        self.assertEqual(manifest["genome_retrieval"]["genomes"], 4)
        self.assertEqual(manifest["generated_at"], "2026-10-11T00:00:00Z")
        self.assertEqual(manifest["builder_sha256"], pop.sha256(Path(pop.__file__)))
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.build("a")

    def test_phenotype_join_does_not_change_phenotype(self):
        self.standard_population()
        self.build()
        total = sum(int(r["n_total"]) for r in self.read("lineage_antimicrobial_summary.tsv"))
        self.assertEqual(total, 4)
        ep = (self.cohort / "isolate_drug_endpoints.tsv").read_text().replace("\tR\ttrue", "\tND\ttrue")
        (self.cohort / "isolate_drug_endpoints.tsv").write_text(ep)
        sums = [f"{pop.sha256(self.cohort / n)}  {n}" for n in pop.COHORT_FILES]
        (self.cohort / "SHA256SUMS").write_text("\n".join(sums) + "\n")
        with self.assertRaisesRegex(ValueError, "not binary"):
            self.build("bad")

    def test_annotation_hash_mismatch_is_rejected(self):
        self.standard_population()
        data = json.loads((self.cache / "GCA_000000001.1" / "mlst.json").read_text())
        data["input_sha256"] = "0" * 64
        (self.cache / "GCA_000000001.1" / "mlst.json").write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "does not match cached genome"):
            self.build()

    def test_cli(self):
        self.standard_population()
        code = pop.main(["build", "--cohort", str(self.cohort), "--cache", str(self.cache),
                         "--config", str(self.config_path), "--outdir", str(self.root / "cli")])
        self.assertEqual(code, 0)
        self.assertTrue((self.root / "cli/SHA256SUMS").exists())
        self.assertEqual(pop.main(["build", "--cohort", str(self.root / "none"), "--cache", str(self.cache),
                                   "--config", str(self.config_path), "--outdir", str(self.root / "x")]), 2)
        with patch.object(pop, "fetch_genomes", return_value={"ok": 1}) as fetch, \
                patch.object(pop, "annotate", return_value={"ok": 2}) as annotate:
            self.assertEqual(pop.main(["fetch-genomes", "--cohort", str(self.cohort), "--cache", str(self.cache)]), 0)
            self.assertEqual(pop.main(["annotate", "--cohort", str(self.cohort), "--cache", str(self.cache), "--tool",
                                       "mlst", "--env-bin", str(self.root), "--config", str(self.config_path)]), 0)
            fetch.assert_called_once()
            annotate.assert_called_once()


if __name__ == "__main__":
    unittest.main()
