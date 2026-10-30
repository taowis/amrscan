"""Offline tests for the error decomposition (synthetic catalogue, genomes, features and predictions)."""

from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import json
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
spec = importlib.util.spec_from_file_location("error_decomposition", ROOT / "bin/error_decomposition.py")
ed = importlib.util.module_from_spec(spec)
sys.modules["error_decomposition"] = ed
spec.loader.exec_module(ed)
gb = ed.gb
NA = ed.NA
CONFIG = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/error_decomposition_v1.json").read_text())
GCONFIG = json.loads((ROOT / "benchmark/cohorts/ecoli_research_v1/generalisation_v1.json").read_text())
STUDY = set(CONFIG["catalogue"]["antimicrobials"])


def rnd_dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def space(kept=(), sample_all=True):
    """A feature space where every k-mer is sampled (or none), and `kept` hashes are features."""
    hashes = np.array(sorted(kept), dtype=np.uint64)
    s = ed.FeatureSpace(hashes, np.arange(hashes.shape[0]), 200, 31, 42)
    s.max_hash = 2 ** 64 if sample_all else 0
    return s


def mapping_row(det, drug="ciprofloxacin", kind="point_mutation", interp="predicts_resistance", requires=NA,
                organism="Escherichia coli", family=None):
    return {"mapping_version": "ecoli_amrfinder_drug_mapping_v1", "organism": organism,
            "rule_id": f"{det}|{drug}", "determinant": det, "determinant_type": kind,
            "gene_family": family or det.split("_")[0], "allele_or_mutation": NA, "antimicrobial": drug,
            "interpretation": interp, "requires": requires, "evidence_source": "PointFinder",
            "evidence_terms": NA, "evidence_reference": "PMID 1", "evidence_level": "source_exact_term",
            "review_status": "reviewed_source", "notes": NA}


LOCI = {"gyrA_S83L": {"locus_type": "protein", "locus_reference": "WP_GYRA", "locus_position": 83,
                      "amrfinder_class": "QUINOLONE", "amrfinder_subclass": "QUINOLONE"},
        "ampC_C-42T": {"locus_type": "dna", "locus_reference": "REF_AMPC", "locus_position": 12,
                       "amrfinder_class": "BETA-LACTAM", "amrfinder_subclass": "CEPHALOSPORIN"}}


def protein(n, seed):
    r = random.Random(seed)
    return "M" + "".join(r.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 1))


def encode(prot, seed):
    r = random.Random(seed)
    back = {}
    for c, a in ed.CODON.items():
        back.setdefault(a, []).append(c)
    return "".join(r.choice(back[a]) for a in prot) + "TAA"


class CatalogueTests(unittest.TestCase):
    def build(self, rows):
        return ed.build_catalogue(rows, LOCI, {"blaTEM": "other"}, CONFIG, STUDY)

    def test_acquired_gene_and_mutation_classes(self):
        cat = self.build([mapping_row("blaTEM-1", "ampicillin", "gene", family="blaTEM"),
                          mapping_row("gyrA_S83L"), mapping_row("ampC_C-42T", "cefoxitin")])
        cls = {r["determinant"]: r["mechanism_class"] for r in cat}
        self.assertEqual(cls, {"blaTEM-1": "acquired_gene", "gyrA_S83L": "target_site_mutation",
                               "ampC_C-42T": "promoter_mutation"})
        mut = next(r for r in cat if r["determinant"] == "gyrA_S83L")
        self.assertEqual((mut["ref_allele"], mut["alt_allele"], mut["locus_position"]), ("S", "L", 83))
        self.assertEqual(next(r for r in cat if r["determinant"] == "blaTEM-1")["gene_function"], "other")

    def test_intrinsic_gene(self):
        cat = self.build([mapping_row("blaEC", "ampicillin", "gene", "does_not_predict", family="blaEC")])
        self.assertEqual(cat[0]["mechanism_class"], "intrinsic_gene")

    def test_unknown_mechanism_stops(self):
        loci = {**LOCI, "zzzX_A1T": {**LOCI["gyrA_S83L"]}}
        with self.assertRaisesRegex(ed.CatalogueError, "no mechanism class"):
            ed.build_catalogue([mapping_row("zzzX_A1T")], loci, {}, CONFIG, STUDY)

    def test_mutation_without_locus_stops(self):
        with self.assertRaisesRegex(ed.CatalogueError, "without an AMRFinderPlus locus"):
            self.build([mapping_row("gyrA_D87N")])

    def test_conflicting_definition(self):
        rows = [mapping_row("gyrA_S83L"), mapping_row("gyrA_S83L", interp="uncertain")]
        with self.assertRaisesRegex(ed.CatalogueError, "conflicting"):
            self.build(rows)

    def test_duplicate_rule(self):
        with self.assertRaisesRegex(ed.CatalogueError, "duplicate"):
            self.build([mapping_row("gyrA_S83L"), mapping_row("gyrA_S83L")])

    def test_wrong_organism(self):
        with self.assertRaisesRegex(ed.CatalogueError, "organism"):
            self.build([mapping_row("gyrA_S83L", organism="Salmonella enterica")])

    def test_wrong_antimicrobial_is_dropped_or_rejected(self):
        self.assertEqual(self.build([mapping_row("gyrA_S83L", drug="imipenem")]), [])
        row = self.build([mapping_row("gyrA_S83L")])[0]
        with self.assertRaisesRegex(ed.CatalogueError, "not a study antimicrobial"):
            ed.validate_catalogue([{**row, "antimicrobial": "imipenem"}], STUDY, CONFIG)

    def test_deterministic_ordering(self):
        rows = [mapping_row("gyrA_S83L"), mapping_row("blaTEM-1", "ampicillin", "gene", family="blaTEM"),
                mapping_row("ampC_C-42T", "cefoxitin")]
        a, b = self.build(rows), self.build(list(reversed(rows)))
        self.assertEqual(a, b)
        self.assertEqual([r["antimicrobial"] for r in a], sorted(r["antimicrobial"] for r in a))
        with self.assertRaisesRegex(ed.CatalogueError, "canonical order"):
            ed.validate_catalogue(list(reversed(a)), STUDY, CONFIG)

    def test_catalogue_version(self):
        row = self.build([mapping_row("gyrA_S83L")])[0]
        self.assertEqual(row["catalogue_version"], CONFIG["catalogue"]["version"])
        with self.assertRaisesRegex(ed.CatalogueError, "catalogue version"):
            ed.validate_catalogue([{**row, "catalogue_version": "other"}], STUDY, CONFIG)

    def test_catalogue_checksum(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "cat.tsv"
            ed.write_tsv(path, ed.CATALOGUE_FIELDS, self.build([mapping_row("gyrA_S83L")]))
            cfg = {**CONFIG, "catalogue": {**CONFIG["catalogue"], "path": str(path), "sha256": "0" * 64}}
            with self.assertRaisesRegex(ed.CatalogueError, "checksum"):
                ed.load_catalogue(cfg)
            cfg["catalogue"]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(len(ed.load_catalogue(cfg)), 1)

    def test_committed_catalogue_matches_config(self):
        self.assertTrue(ed.load_catalogue(CONFIG))

    def test_indicator_symbols(self):
        cat = self.build([mapping_row("gyrA_S83L"), mapping_row("ampC_C-42T", "cefoxitin", interp="uncertain")])
        self.assertEqual(ed.indicator_symbols(cat, ["predicts_resistance", "component"]), ["gyrA_S83L"])

    def test_parse_mutation(self):
        self.assertEqual(ed.parse_mutation("gyrA_S83L")["pos"], 83)
        self.assertEqual(ed.parse_mutation("ampC_C-42T")["pos"], -42)
        self.assertEqual(ed.parse_mutation("nfsA_Q67Ter")["alt"], "*")
        self.assertEqual(ed.parse_mutation("ampC_G-15GG")["alt"], "GG")
        self.assertIsNone(ed.parse_mutation("cirA_S90YfsTer15"))


class RepresentationTests(unittest.TestCase):
    def setUp(self):
        self.ref = protein(120, 1)
        self.ref = self.ref[:82] + "S" + self.ref[83:]
        mut = self.ref[:82] + "L" + self.ref[83:]
        self.gene = encode(mut, 2)
        left, right = rnd_dna(200, 3), rnd_dna(200, 4)
        self.contig = left + self.gene + right
        self.call = {"element_symbol": "gyrA_S83L", "element_subtype": "POINT", "contig": "c1", "method": "POINTX",
                     "start": "201", "stop": str(200 + len(self.gene)), "strand": "+"}
        self.refs = {"loci": LOCI, "proteins": {"WP_GYRA": self.ref + "*"}, "dna": {}}
        self.cfg = CONFIG["capture"]
        self.codon_at = 200 + 82 * 3

    def mutant_kmers(self, contig, at):
        codon = contig[at:at + 3]
        changed = [i for i in range(3) if min(("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
                                              key=lambda c: ed.hamming(c, codon))[i] != codon[i]]
        return ed.overlapping_kmers(contig, at + changed[0], at + changed[-1] + 1, 31)

    def hashes(self, kmers):
        s = space()
        return [s.kmer_hash(k) for k in kmers]

    def test_site_located_on_both_strands(self):
        text, off, length = ed.oriented_locus(self.contig, 201, 200 + len(self.gene), "+", 40)
        site = ed.protein_site(text, off, length, self.ref, 83, "S", "L", self.cfg)
        self.assertEqual(ed.CODON[site["codon"]], "L")
        rc = ed.revcomp(self.contig)
        start = len(self.contig) - (200 + len(self.gene)) + 1
        text, off, length = ed.oriented_locus(rc, start, start + len(self.gene) - 1, "-", 40)
        site_rc = ed.protein_site(text, off, length, self.ref, 83, "S", "L", self.cfg)
        self.assertEqual(site_rc["codon"], site["codon"])

    def test_mutation_encoded_by_sampled_feature(self):
        kept = self.hashes(self.mutant_kmers(self.contig, self.codon_at))[:1]
        r = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(kept), self.cfg)
        self.assertEqual(r["status"], "captured")
        self.assertEqual(r["mutant"]["kept"], 1)

    def test_mutation_not_encoded(self):
        r = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(sample_all=False), self.cfg)
        self.assertEqual(r["status"], "not_represented")

    def test_partial_when_sampled_but_pruned(self):
        r = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(), self.cfg)
        self.assertEqual(r["status"], "partially_represented")
        self.assertGreater(r["mutant"]["sampled"], 0)

    def test_partial_when_only_wild_type_kmer_kept(self):
        at = self.codon_at
        codon = self.contig[at:at + 3]
        wt = min((c for c in ed.CODON if ed.CODON[c] == "S"), key=lambda c: ed.hamming(c, codon))
        wt_contig = self.contig[:at] + wt + self.contig[at + 3:]
        changed = [i for i in range(3) if wt[i] != codon[i]]
        wt_kmers = ed.overlapping_kmers(wt_contig, at + changed[0], at + changed[-1] + 1, 31)
        r = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(self.hashes(wt_kmers)[:1]),
                                self.cfg)
        self.assertEqual(r["mutant"]["kept"], 0)
        self.assertGreaterEqual(r["wild_type"]["kept"], 1)
        # only the wild-type allele has a feature: the mutation shows as an absence
        unsampled = {"kmers": 31, "sampled": 0, "kept": 0}
        self.assertEqual(ed.allele_status(unsampled, {"kmers": 31, "sampled": 1, "kept": 1}),
                         "partially_represented")
        self.assertEqual(ed.allele_status(unsampled, {"kmers": 31, "sampled": 1, "kept": 0}), "not_represented")
        self.assertEqual(ed.allele_status({"kmers": 0}, {"kmers": 0}), "representation_unknown")

    def test_site_not_found_is_unknown(self):
        refs = {**self.refs, "proteins": {"WP_GYRA": protein(120, 9)}}
        r = ed.mutation_capture({"c1": self.contig}, self.call, refs, space(), self.cfg)
        self.assertEqual(r["status"], "representation_unknown")
        missing = ed.mutation_capture({}, self.call, self.refs, space(), self.cfg)
        self.assertEqual(missing["detail"], "contig_not_found")

    def test_ambiguous_bases_unknown(self):
        at = self.codon_at
        contig = self.contig[:at - 40] + "N" * 120 + self.contig[at + 80:]
        call = dict(self.call)
        r = ed.mutation_capture({"c1": contig}, call, self.refs, space(), self.cfg)
        self.assertEqual(r["status"], "representation_unknown")

    def test_dna_promoter_mutation(self):
        ref = rnd_dna(60, 5)
        ref = ref[:11] + "C" + ref[12:]
        mutant = ref[:11] + "T" + ref[12:]
        contig = rnd_dna(300, 6) + mutant + rnd_dna(300, 7)
        call = {"element_symbol": "ampC_C-42T", "element_subtype": "POINT", "contig": "c1", "method": "POINTN",
                "start": "301", "stop": "360", "strand": "+"}
        refs = {"loci": LOCI, "proteins": {}, "dna": {"REF_AMPC": ref}}
        site = 300 + 11
        kept = self.hashes(ed.overlapping_kmers(contig, site, site + 1, 31))[:1]
        r = ed.mutation_capture({"c1": contig}, call, refs, space(kept), self.cfg)
        self.assertEqual(r["status"], "captured")
        rc = ed.mutation_capture({"c1": ed.revcomp(contig)}, {**call, "start": "301", "stop": "360"}, refs,
                                 space(kept), self.cfg)
        self.assertEqual(rc["status"], "captured")

    def test_acquired_gene_represented_and_absent(self):
        contig = rnd_dna(1200, 8)
        call = {"contig": "c1", "start": "101", "stop": "900"}
        kept = self.hashes([contig[300:331]])
        self.assertEqual(ed.gene_capture({"c1": contig}, call, space(kept))["status"], "captured")
        self.assertEqual(ed.gene_capture({"c1": contig}, call, space())["status"], "partially_represented")
        self.assertEqual(ed.gene_capture({"c1": contig}, call, space(sample_all=False))["status"],
                         "not_represented")
        self.assertEqual(ed.gene_capture({}, call, space())["status"], "representation_unknown")
        outside = self.hashes([contig[1000:1031]])
        self.assertEqual(ed.gene_capture({"c1": contig}, call, space(outside))["status"], "partially_represented")

    def test_canonical_hash_and_columns(self):
        s = space()
        k = rnd_dna(31, 10)
        self.assertEqual(s.kmer_hash(k), s.kmer_hash(ed.revcomp(k)))
        self.assertIsNone(s.kmer_hash(k[:-1] + "N"))
        kept = sorted([s.kmer_hash(k), s.kmer_hash(rnd_dna(31, 11))])
        s2 = ed.FeatureSpace(np.array(kept, dtype=np.uint64), np.array([7, 7]), 200, 31, 42)
        self.assertEqual(s2.column(kept[0]), 7)
        self.assertEqual(s2.column(kept[1]), 7)  # merged identical columns share one pattern
        self.assertIsNone(s2.column(kept[0] + 1))

    def test_capture_deterministic(self):
        kept = self.hashes(self.mutant_kmers(self.contig, self.codon_at))[:2]
        a = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(kept), self.cfg)
        b = ed.mutation_capture({"c1": self.contig}, self.call, self.refs, space(kept), self.cfg)
        self.assertEqual(a, b)

    def test_genome_checks(self):
        with tempfile.TemporaryDirectory() as d:
            data = f">c1 test\n{self.contig}\n".encode()
            (Path(d) / "GCA_000000001.1").mkdir()
            (Path(d) / "GCA_000000001.1" / "genome.fna").write_bytes(data)
            row = {"isolate_id": "I1", "assembly_accession": "GCA_000000001.1",
                   "genome_sha256": hashlib.sha256(data).hexdigest()}
            args = (Path(d), self.refs, space(), self.cfg)
            ok = ed.audit_isolate(row, [self.call], *args)
            self.assertEqual(ok["genome_status"], "verified")
            self.assertEqual(ok["mutations"][0]["status"], "partially_represented")
            bad = ed.audit_isolate({**row, "genome_sha256": "0" * 64}, [self.call], *args)
            self.assertEqual((bad["genome_status"], bad["mutations"]), ("genome_hash_mismatch", []))
            gone = ed.audit_isolate({**row, "assembly_accession": "GCA_000000002.1"}, [self.call], *args)
            self.assertEqual(gone["genome_status"], "genome_missing")
            version = gb.verify_genome(Path(d), row, expected_assembly="GCA_000000001.2")
            self.assertEqual(version[0], "genome_version_mismatch")


class NoveltyTests(unittest.TestCase):
    def ep(self, n, drug="ciprofloxacin"):
        return {"isolate_id": f"I{n}", "antimicrobial": drug}

    def test_seen_unseen_and_combinations(self):
        train = [self.ep(i) for i in range(8)]
        dets = {(f"I{i}", "ciprofloxacin"): ([], ["gyrA_S83L"]) for i in range(6)}
        dets[("I6", "ciprofloxacin")] = (["qnrS1"], ["gyrA_S83L"])
        single, combo = ed.training_exposure(train, dets)
        self.assertEqual(single["gyrA_S83L"], 7)
        self.assertEqual(ed.seen_flag(["gyrA_S83L"], single, 5), "yes")
        self.assertEqual(ed.seen_flag(["qnrS1"], single, 5), "no")
        self.assertEqual(ed.seen_flag(["parC_S80I"], single, 5), "no")
        self.assertEqual(ed.seen_flag([], single, 5), NA)
        self.assertEqual(ed.mechanism_novelty("positive", [], ["gyrA_S83L"], single, 5), "seen")
        self.assertEqual(ed.mechanism_novelty("positive", ["qnrS1"], ["gyrA_S83L"], single, 5), "unseen")
        self.assertEqual(ed.mechanism_novelty("negative", [], [], single, 5), "no_established_determinant")
        self.assertEqual(ed.mechanism_novelty("uncertain", [], [], single, 5), "unresolved")
        self.assertEqual(combo[("gyrA_S83L", "qnrS1")], 1)
        self.assertEqual(combo[("gyrA_S83L", "parC_S80I")], 0)

    def test_fold_specific_and_no_test_leak(self):
        dets = {(f"I{i}", "ciprofloxacin"): ([], ["gyrA_S83L"]) for i in range(10)}
        fold1 = [self.ep(i) for i in range(5)]
        fold2 = [self.ep(i) for i in range(4)]
        self.assertEqual(ed.seen_flag(["gyrA_S83L"], ed.training_exposure(fold1, dets)[0], 5), "yes")
        self.assertEqual(ed.seen_flag(["gyrA_S83L"], ed.training_exposure(fold2, dets)[0], 5), "no")
        # endpoints of another drug or outside the training list never count
        other = fold2 + [self.ep(1, "ampicillin")]
        self.assertEqual(ed.training_exposure(other, dets)[0], ed.training_exposure(fold2, dets)[0])

    def test_split_determinants(self):
        self.assertEqual(ed.split_determinants("blaTEM-1;gyrA_S83L", "gene;point_mutation"),
                         (["blaTEM-1"], ["gyrA_S83L"]))
        self.assertEqual(ed.split_determinants(NA, NA), ([], []))

    def test_usable_columns(self):
        sub = sparse.csr_matrix(np.array([[1, 0], [1, 0], [1, 1], [0, 1], [0, 1], [0, 1]], dtype=np.int8))
        u = ed.usable_columns(sub, np.arange(6), 2, 2)
        self.assertEqual(u.tolist(), [True, True])
        self.assertEqual(ed.usable_columns(sub, np.arange(3), 2, 2).tolist(), [False, False])

    def test_shift_measures(self):
        self.assertEqual(ed.tvd(Counter(a=1), Counter(a=3)), 0)
        self.assertAlmostEqual(ed.tvd(Counter(a=1), Counter(b=1)), 1)
        self.assertAlmostEqual(ed.jsd(Counter(a=1), Counter(b=1)), 1)
        self.assertAlmostEqual(ed.jsd(Counter(a=2, b=2), Counter(a=1, b=1)), 0)
        self.assertIsNone(ed.tvd(Counter(), Counter(a=1)))


class DecompositionTests(unittest.TestCase):
    def cat(self, *a):
        return ed.error_category(*a)

    def test_resistant_categories(self):
        self.assertEqual(self.cat("R", "S", "positive", ["blaCTX-M-15"], [], "captured", NA),
                         "R_KNOWN_DET_REPRESENTED")
        self.assertEqual(self.cat("R", "S", "positive", ["blaCTX-M-15"], [], "not_represented", NA),
                         "R_KNOWN_DET_NOT_REPRESENTED")
        self.assertEqual(self.cat("R", "S", "positive", [], ["gyrA_S83L"], NA, "captured"),
                         "R_KNOWN_MUT_REPRESENTED")
        self.assertEqual(self.cat("R", "S", "positive", [], ["gyrA_S83L"], NA, "partially_represented"),
                         "R_KNOWN_MUT_NOT_REPRESENTED")
        self.assertEqual(self.cat("R", "S", "negative", [], [], NA, NA), "R_NO_ESTABLISHED_DET")
        self.assertEqual(self.cat("R", "S", "uncertain", [], [], NA, NA), "R_UNRESOLVED")
        self.assertEqual(self.cat("R", "S", NA, [], [], NA, NA), "R_GENOTYPE_UNAVAILABLE")
        self.assertEqual(self.cat("R", "S", "positive", [], ["gyrA_S83L"], NA, "representation_unknown"),
                         "R_UNRESOLVED")

    def test_multiple_mechanisms_not_forced(self):
        self.assertEqual(self.cat("R", "S", "positive", ["qnrS1"], ["gyrA_S83L"], "not_represented", "captured"),
                         "R_KNOWN_DET_AND_MUT_REPRESENTED")
        self.assertEqual(ed.mechanism_group("positive", ["qnrS1"], ["gyrA_S83L"]), "both")
        self.assertEqual(ed.mechanism_group("positive", [], ["gyrA_S83L"]), "mutation")
        self.assertEqual(ed.mechanism_group("negative", [], []), "no_established_determinant")

    def test_susceptible_categories(self):
        self.assertEqual(self.cat("S", "R", "positive", ["blaTEM-1"], [], "captured", NA), "S_KNOWN_DET_PRESENT")
        self.assertEqual(self.cat("S", "R", "negative", [], [], NA, NA), "S_NO_KNOWN_DET")
        self.assertEqual(self.cat("S", "R", "uncertain", [], [], NA, NA), "S_UNRESOLVED")
        self.assertEqual(self.cat("S", "S", "positive", [], [], NA, NA), "correct")

    def test_phenotype_uncertainty_flags(self):
        rec = {"n_ast_observations": "2", "conflict_status": "mixed_phenotypes", "methods": {NA},
               "standards": {"CLSI 2020", "EUCAST 2021"}}
        f = ed.phenotype_flags(rec, corrected=True)
        self.assertEqual(f["phenotype_flags"], "label_corrected;multiple_observations;conflict;method_missing;"
                                               "mixed_breakpoint_standards")
        self.assertEqual(ed.phenotype_flags(None, False)["phenotype_flags"], "provenance_missing")
        clean = {"n_ast_observations": "1", "conflict_status": "none", "methods": {"MIC"}, "standards": {"CLSI"}}
        self.assertEqual(ed.phenotype_flags(clean, False)["phenotype_flags"], "none")

    def test_uncertain_mapping_is_not_missing_genotype(self):
        self.assertEqual(ed.endpoint_status({"genotype_status": "determinant_uncertain",
                                             "determinant_status": "uncertain"}), "uncertain")
        self.assertEqual(ed.endpoint_status({"genotype_status": "calls_missing", "determinant_status": "NA"}), NA)

    def test_four_way_and_best_capture(self):
        self.assertEqual(ed.four_way("R", "S", "R"), "determinant_only_correct")
        self.assertEqual(ed.four_way("R", "R", "S"), "ml_only_correct")
        self.assertEqual(ed.four_way("S", "S", "S"), "both_correct")
        self.assertEqual(ed.four_way("S", "R", "R"), "both_wrong")
        self.assertEqual(ed.four_way("S", "R", NA), NA)
        self.assertEqual(ed.best_capture(["not_represented", "captured"]), "captured")
        self.assertEqual(ed.best_capture([]), NA)

    def test_deterministic_summary(self):
        rows = [{"cohort": "common", "antimicrobial": d, "split_strategy": s} for d in ("b", "a") for s in ("x", "y")]
        g1, g2 = ed._group(rows, ("cohort", "antimicrobial")), ed._group(list(reversed(rows)), ("cohort",
                                                                                               "antimicrobial"))
        self.assertEqual(sorted(g1), sorted(g2))
        self.assertEqual(list(ed._group(rows, ("cohort",), all_drugs=True)), [("common", "ALL")])


class GroupAndUpstreamTests(unittest.TestCase):
    def drop(self, drug, strat, delta, lower):
        return {"cohort": "common", "comparator": "random", "comparison": "same_test_isolates", "level": "pooled",
                "antimicrobial": drug, "holdout_strategy": strat, "delta_AUROC": str(delta),
                "delta_AUROC_ci_lower": str(lower)}

    def test_drug_groups(self):
        drops = [self.drop("a", "country", 0.2, 0.1), self.drop("a", "lineage", 0.01, 0),
                 self.drop("b", "temporal", 0.06, -0.01), self.drop("c", "lineage", 0.01, -0.01),
                 self.drop("d", "project", 0.03, 0.01), {**self.drop("e", "country", 0.3, 0.2), "cohort": "maximal"}]
        g = ed.drug_groups(drops, CONFIG["drug_groups"])
        self.assertEqual({k: v["group"] for k, v in g.items()},
                         {"a": "high_loss", "b": "intermediate", "c": "stable", "d": "intermediate"})
        self.assertEqual(g["a"]["max_delta_strategy"], "country")

    def write_upstream(self, d):
        d = Path(d)
        (d / "bench" / "matrix").mkdir(parents=True)
        files = {"benchmark_manifest": d / "bench/benchmark_manifest.json", "predictions": d / "bench/predictions.tsv.gz",
                 "feature_manifest": d / "bench/matrix/feature_manifest.json", "generalisation_config": d / "g.json",
                 "mapping_rules": d / "map.tsv", "genome_manifest": d / "genomes.tsv"}
        for name, path in files.items():
            path.write_text(json.dumps({"benchmark_version": "generalisation_v1", "genome_manifest_sha256": "g",
                                        "upstream": {}}) if name == "benchmark_manifest" else name)
        db = d / "db"
        db.mkdir()
        (db / "AMRProt-mutation.tsv").write_text("x")
        up = {**CONFIG["upstream"], "generalisation_config": str(files["generalisation_config"]),
              "mapping_rules": str(files["mapping_rules"]), "genome_manifest": str(files["genome_manifest"]),
              "amrfinder_database": str(db), "genome_manifest_sha256": "g",
              "amrfinder_database_files_sha256": {"AMRProt-mutation.tsv": hashlib.sha256(b"x").hexdigest()}}
        for name, path in files.items():
            up[f"{name}_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        up["genome_manifest_sha256"] = up["genome_manifest_sha256"]
        return {**CONFIG, "upstream": up}, files

    def test_upstream_hashes(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, files = self.write_upstream(d)
            # genome manifest hash in the benchmark manifest must match the config's
            files["benchmark_manifest"].write_text(json.dumps({"benchmark_version": "generalisation_v1",
                                                               "genome_manifest_sha256": cfg["upstream"][
                                                                   "genome_manifest_sha256"], "upstream": {}}))
            cfg["upstream"]["benchmark_manifest_sha256"] = hashlib.sha256(
                files["benchmark_manifest"].read_bytes()).hexdigest()
            got = ed.check_upstream(cfg, Path(d) / "bench")
            self.assertEqual(got["benchmark_version"], "generalisation_v1")
            for name in ("predictions", "benchmark_manifest", "generalisation_config", "genome_manifest",
                         "mapping_rules", "feature_manifest"):
                with self.subTest(name=name):
                    bad = json.loads(json.dumps(cfg))
                    bad["upstream"][f"{name}_sha256"] = "0" * 64
                    with self.assertRaisesRegex(ed.UpstreamError, name):
                        ed.check_upstream(bad, Path(d) / "bench")
            bad = json.loads(json.dumps(cfg))
            bad["upstream"]["amrfinder_database_files_sha256"]["AMRProt-mutation.tsv"] = "0" * 64
            with self.assertRaisesRegex(ed.UpstreamError, "database"):
                ed.check_upstream(bad, Path(d) / "bench")

    def task(self, drug="ciprofloxacin", strat="random", fold="fold1", test=("I1", "I2")):
        return {"cohort": "common", "antimicrobial": drug, "strategy": strat, "fold": fold,
                "task_id": f"common__{drug}__{strat}__{fold}", "test": [{"isolate_id": i} for i in test]}

    def write_preds(self, d, rows):
        path = Path(d) / "p.tsv.gz"
        fields = ["cohort", "antimicrobial", "split_strategy", "fold_or_holdout_id", "isolate_id", *ed.PRED_KEEP]
        with gzip.open(path, "wt") as h:
            w = csv.DictWriter(h, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({f: r.get(f, "NA") for f in fields})
        return path

    def pred(self, iid, drug="ciprofloxacin", fold="fold1"):
        return {"cohort": "common", "antimicrobial": drug, "split_strategy": "random", "fold_or_holdout_id": fold,
                "isolate_id": iid, "predicted_probability": "0.7", "predicted_class": "R", "correct": "yes"}

    def test_prediction_reconciliation(self):
        with tempfile.TemporaryDirectory() as d:
            preds = ed.load_predictions(self.write_preds(d, [self.pred("I1"), self.pred("I2")]))
            ed.reconcile([self.task()], preds)
            with self.assertRaisesRegex(ed.UpstreamError, "without a prediction"):
                ed.reconcile([self.task(test=("I1", "I2", "I3"))], preds)
            with self.assertRaisesRegex(ed.UpstreamError, "duplicate"):
                ed.load_predictions(self.write_preds(d, [self.pred("I1"), self.pred("I1")]))
            extra = ed.load_predictions(self.write_preds(d, [self.pred("I1"), self.pred("I2"),
                                                             self.pred("I1", "imipenem")]))
            with self.assertRaisesRegex(ed.UpstreamError, "outside the benchmark"):
                ed.reconcile([self.task()], extra)
            stray = ed.load_predictions(self.write_preds(d, [self.pred("I1"), self.pred("I2"),
                                                             self.pred("I1", fold="fold9")]))
            with self.assertRaisesRegex(ed.UpstreamError, "not frozen"):
                ed.reconcile([self.task()], stray)

    def test_config_rules(self):
        bad = {**CONFIG, "drug_groups": {**CONFIG["drug_groups"], "stable_max_delta": 0.1}}
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "c.json"
            path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, "stable_max_delta"):
                ed.load_config(path)


class CallTests(unittest.TestCase):
    def write(self, d, **prov):
        iso = Path(d) / "I1"
        iso.mkdir()
        report = "\t".join(["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name",
                             "Scope", "Type", "Subtype", "Class", "Subclass", "Method", "Target length",
                             "Reference sequence length", "% Coverage of reference", "% Identity to reference",
                             "Alignment length", "Closest reference accession", "Closest reference name",
                             "HMM accession", "HMM description", "Hierarchy node"]) + "\n"
        for sym, typ, sub in (("gyrA_S83L", "AMR", "POINT"), ("blaTEM-1", "AMR", "AMR"), ("qacEdelta1", "STRESS",
                                                                                             "BIOCIDE")):
            report += "\t".join(["NA", "c1", "10", "100", "+", sym, "x", "core", typ, sub, "C", "C", "POINTX", "30",
                                  "30", "100", "100", "30", "WP_1", "x", "NA", "NA", sym]) + "\n"
        (iso / "amrfinder.tsv").write_text(report)
        base = {"caller_version": "4.2.7", "status": "success", "input_sha256": "g" * 64,
                "raw_sha256": hashlib.sha256(report.encode()).hexdigest(),
                "database_files_sha256": dict(CONFIG["upstream"]["amrfinder_database_files_sha256"])}
        base.update(prov)
        (iso / "provenance.json").write_text(json.dumps(base))

    def test_calls_checked_against_provenance(self):
        row = {"genome_sha256": "g" * 64, "assembly_accession": "GCA_1.1"}
        cases = {"success": {}, "genome_hash_mismatch": {"input_sha256": "h" * 64},
                 "caller_version_mismatch": {"caller_version": "4.0.23"},
                 "call_not_successful": {"status": "failed"}, "report_hash_mismatch": {"raw_sha256": "0" * 64},
                 "database_mismatch": {"database_files_sha256": {}}}
        for expected, prov in cases.items():
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as d:
                self.write(d, **prov)
                status, rows = ed.read_calls(Path(d), "I1", row, CONFIG)
                self.assertEqual(status, expected)
                if expected == "success":
                    self.assertEqual([r["element_symbol"] for r in rows], ["gyrA_S83L", "blaTEM-1"])
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(ed.read_calls(Path(d), "I1", row, CONFIG)[0], "calls_missing")

    def test_known_mutation_rows(self):
        cat = {"ciprofloxacin": {"gyrA_S83L": {"determinant_type": "point_mutation", "mechanism_class":
                                               "target_site_mutation", "evidence_source": "PointFinder",
                                               "evidence_reference": "PMID 1", "interpretation":
                                               "predicts_resistance", "antimicrobial": "ciprofloxacin"}}}
        calls = [{"isolate_id": "I1", "element_symbol": "gyrA_S83L", "element_subtype": "POINT", "method": "POINTX"},
                 {"isolate_id": "I1", "element_symbol": "uhpT_E350Q", "element_subtype": "POINT",
                  "method": "POINTX"},
                 {"isolate_id": "I1", "element_symbol": "blaTEM-1", "element_subtype": "AMR", "method": "EXACTX"}]
        rows = ed.known_mutation_rows(calls, cat, CONFIG)
        self.assertEqual([(r["mutation"], r["status"]) for r in rows],
                         [("gyrA_S83L", "predicts_resistance"), ("uhpT_E350Q", "not_in_catalogue")])
        self.assertEqual(rows[0]["database_version"], "2026-08-07.1")


class AugmentationTests(unittest.TestCase):
    def calls(self):
        return {"I0": [{"element_symbol": "gyrA_S83L", "element_subtype": "POINT", "method": "POINTX"}],
                "I1": [{"element_symbol": "gyrA_S83L", "element_subtype": "POINT", "method": "PARTIALX"},
                       {"element_symbol": "blaTEM-1", "element_subtype": "AMR", "method": "EXACTX"}],
                "I2": [{"element_symbol": "parC_S80I", "element_subtype": "POINT", "method": "POINTX"},
                       {"element_symbol": "gyrA_S83L", "element_subtype": "POINT", "method": "POINTX"}]}

    def test_indicators_from_calls_only(self):
        b = ed.indicator_block(["I0", "I1", "I2", "I3"], self.calls(), ["gyrA_S83L", "parC_S80I"])
        self.assertEqual(b.toarray().tolist(), [[1, 0], [0, 0], [1, 1], [0, 0]])
        self.assertEqual(ed.indicator_block.__code__.co_varnames[:3], ("isolates", "calls_by_isolate", "symbols"))

    def test_append_columns_alignment(self):
        rng = np.random.default_rng(0)
        X = sparse.random(30, 50, density=0.2, format="csr", random_state=1)
        X.data[:] = 1
        X = X.astype(np.int8)
        B = sparse.csr_matrix((rng.random((30, 4)) < 0.3).astype(np.int8))
        got = ed.append_columns(X, B)
        self.assertEqual((got != sparse.hstack([X, B]).tocsr()).nnz, 0)
        self.assertEqual(got.shape, (30, 54))
        with self.assertRaises(ValueError):
            ed.append_columns(X, B[:10])

    def test_permuted_control_keeps_prevalence(self):
        B = ed.indicator_block([f"I{i}" for i in range(3)], self.calls(), ["gyrA_S83L", "parC_S80I"])
        P = ed.permute_rows(B, 1)
        self.assertEqual(np.asarray(P.sum(0)).tolist(), np.asarray(B.sum(0)).tolist())
        self.assertEqual(ed.permute_rows(B, 1).toarray().tolist(), P.toarray().tolist())

    def ctx(self, drug, n, captured):
        return [
            {"antimicrobial": drug, "catalogue_mutations": "gyrA_S83L" if i < n else NA,
             "catalogue_mutation_capture": ("captured" if i < captured else "not_represented") if i < n else NA}
            for i in range(n + 5)]

    def test_eligibility_and_unsupported_drug(self):
        cat = {"ciprofloxacin": {"gyrA_S83L": {"determinant_type": "point_mutation",
                                               "interpretation": "predicts_resistance"}},
               "chloramphenicol": {"catA1": {"determinant_type": "gene", "interpretation": "predicts_resistance"}},
               "levofloxacin": {"gyrA_S83L": {"determinant_type": "point_mutation",
                                              "interpretation": "predicts_resistance"}},
               "nalidixic acid": {"gyrA_S83L": {"determinant_type": "point_mutation",
                                                "interpretation": "predicts_resistance"}}}
        rows = self.ctx("ciprofloxacin", 40, 4) + self.ctx("levofloxacin", 10, 0) + \
            self.ctx("nalidixic acid", 30, 29) + self.ctx("chloramphenicol", 0, 0)
        e = ed.augmentation_eligibility(rows, CONFIG["augmentation"], cat)
        self.assertTrue(e["ciprofloxacin"]["eligible"])
        self.assertFalse(e["levofloxacin"]["eligible"])
        self.assertFalse(e["nalidixic acid"]["eligible"])
        self.assertIn("no point-mutation rule", e["chloramphenicol"]["reason"])

    def test_gap_rescue(self):
        rng = np.random.default_rng(3)
        y = np.array([1] * 30 + [0] * 30, dtype=bool)
        good = np.where(y, 0.9, 0.1) + rng.normal(0, 0.05, 60)
        weak = np.where(y, 0.6, 0.4) + rng.normal(0, 0.2, 60)
        units = [f"u{i // 2}" for i in range(60)]
        ci = ed.gap_rescue_ci(y, units, good, weak, good, good, 200, 1, 0.95)
        gap = (gb.auroc(y, good) - gb.auroc(y, weak)) - 0
        self.assertLessEqual(ci["lower"], gap + 1e-9)
        self.assertGreater(ci["upper"], 0)
        none = ed.gap_rescue_ci(y, units, good, weak, good, weak, 200, 1, 0.95)
        self.assertAlmostEqual(none["lower"], 0)
        self.assertAlmostEqual(none["upper"], 0)

    def test_augmented_model_deterministic_on_frozen_split(self):
        rng = np.random.default_rng(5)
        n = 80
        y = rng.random(n) < 0.5
        X = sparse.csr_matrix((rng.random((n, 40)) < 0.4).astype(np.int8))
        mut = sparse.csr_matrix(y.reshape(-1, 1).astype(np.int8))
        Xa = ed.append_columns(X, mut)
        eps = [{"isolate_id": f"I{i}", "phenotype": "R" if y[i] else "S", "split_unit": f"I{i}"} for i in range(n)]
        task = {"task_id": "common__ciprofloxacin__random__fold1", "cohort": "common",
                "antimicrobial": "ciprofloxacin", "strategy": "random", "fold": "fold1",
                "train": eps[:60], "test": eps[60:]}
        row_of = {e["isolate_id"]: i for i, e in enumerate(eps)}
        cfg = json.loads(json.dumps(GCONFIG))
        cfg["features"]["min_train_present"] = cfg["features"]["min_train_absent"] = 2
        s1, p1 = gb.run_task(task, Xa, row_of, cfg)
        s2, p2 = gb.run_task(task, Xa, row_of, cfg)
        self.assertEqual(p1, p2)
        self.assertEqual([p["isolate_id"] for p in p1], [e["isolate_id"] for e in task["test"]])
        self.assertEqual(s1["split_sha256"], gb.task_digest(task))
        s0, p0 = gb.run_task(task, X, row_of, cfg)
        self.assertEqual([p["isolate_id"] for p in p0], [p["isolate_id"] for p in p1])
        self.assertGreaterEqual(gb.auroc(y[60:], np.array([p["predicted_probability"] for p in p1])),
                                gb.auroc(y[60:], np.array([p["predicted_probability"] for p in p0])))


class ContextIoTests(unittest.TestCase):
    def test_context_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "ctx.tsv.gz"
            w = ed.ContextWriter(path)
            w.write({"cohort": "common", "antimicrobial": "ampicillin", "phenotype": "R", "year": 2015})
            w.write({"cohort": "maximal", "antimicrobial": "ampicillin", "phenotype": "S"})
            w.close()
            rows = ed.read_context(path)
            self.assertEqual(rows[0]["year"], "2015")
            self.assertEqual(rows[1]["four_way"], NA)
            self.assertEqual(rows[1].get("missing_field", "x"), "x")
            common = ed.read_context(path, keep=lambda r: r["cohort"] == "common")
            self.assertEqual(len(common), 1)

    def test_perf_and_support(self):
        rows = [{"phenotype": "R", "predicted_probability": "0.9", "predicted_class": "R"}] * 10 + \
               [{"phenotype": "S", "predicted_probability": "0.1", "predicted_class": "S"}] * 9
        p = ed.perf(rows, CONFIG)
        self.assertEqual((p["AUROC"], p["supported"]), (1.0, "no"))
        self.assertEqual((p["N"], p["N_R"], p["N_S"]), (19, 10, 9))
        self.assertTrue(set(ed.PERF_FIELDS) <= set(p))


# ------------------------------------------------------------ end-to-end synthetic study

def load_gb_study():
    """The 150-isolate synthetic study of the generalisation benchmark tests, run through every PR #14 command."""
    spec = importlib.util.spec_from_file_location("_gb_study", ROOT / "tests/test_generalisation_benchmark.py")
    mod = importlib.util.module_from_spec(spec)
    # the study registers its own benchmark module for pool workers; put back the one other tests use
    previous = sys.modules.get("generalisation_benchmark")
    try:
        spec.loader.exec_module(mod)

        class Study(mod.EndToEndTests):
            pass

        Study.setUpClass()
    finally:
        if previous is not None:
            sys.modules["generalisation_benchmark"] = previous
    assert Study.codes == [0] * 6, Study.codes
    return Study


def write_rows(path, fields, rows):
    with Path(path).open("w", newline="") as h:
        w = csv.DictWriter(h, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


REPORT_HEAD = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name", "Scope", "Type",
               "Subtype", "Class", "Subclass", "Method", "Target length", "Reference sequence length",
               "% Coverage of reference", "% Identity to reference", "Alignment length",
               "Closest reference accession", "Closest reference name", "HMM accession", "HMM description",
               "Hierarchy node"]


def report_line(contig, start, stop, sym, typ, sub, method):
    return "\t".join(["NA", contig, str(start), str(stop), "+", sym, "x", "core", typ, sub, "C", "C", method, "30",
                      "30", "100", "100", "30", "WP_1", "x", "NA", "NA", sym]) + "\n"


class EndToEndTests(unittest.TestCase):
    """The PR #14 synthetic study, with AMRFinderPlus calls, a small database and AST provenance, run through
    every command. Ampicillin resistance is carried by an acquired gene in the genomes; tetracycline resistance
    by a point-mutation call that no k-mer encodes, so only the augmented model can see it."""

    @classmethod
    def setUpClass(cls):
        cls.study = study = load_gb_study()
        root = cls.root = study.root
        bench, pop_dir, base_dir = study.out, root / "population", root / "baseline"
        n = 150
        gene_present = {i: ((i * 7) % 10 < 4) != (i % 23 == 0) for i in range(n)}
        mutant = {i: i % 5 < 2 and i % 17 != 0 for i in range(n)}
        # population QC fields and determinant columns that the PR #14 tests do not need
        meta = ed.pop.read_tsv(pop_dir / "population_metadata.tsv")
        for r in meta:
            i = int(r["isolate_id"][3:12])
            r.update({"genome_qc_status": "warn" if i % 29 == 0 else "pass", "genome_qc_flags": NA,
                      "assembly_length": str(2000 + i), "contig_count": "1", "n50": str(2000 + i),
                      "ambiguous_fraction": "0"})
        write_rows(pop_dir / "population_metadata.tsv", list(meta[0].keys()), meta)
        eps = ed.pop.read_tsv(base_dir / "genotype_phenotype_endpoints.tsv")
        for r in eps:
            i = int(r["isolate_id"][3:12])
            pos = r["determinant_status"] == "positive"
            det = {"ampicillin": ("blaTEM-1", "gene"), "tetracycline": ("gyrA_S83L", "point_mutation")}.get(
                r["antimicrobial"])
            r["genotype_status"] = "calls_missing" if i == 1 else "called"
            if r["antimicrobial"] == "tetracycline" and i % 13 == 0:
                r["determinant_status"] = "uncertain"
            r["determinant_ids"], r["determinant_types"] = det if det and pos else (NA, NA)
        write_rows(base_dir / "genotype_phenotype_endpoints.tsv", list(eps[0].keys()), eps)
        for d in (pop_dir, base_dir):
            names = [line.split("  ")[1].strip() for line in (d / "SHA256SUMS").read_text().splitlines()]
            (d / "SHA256SUMS").write_text("".join(f"{sha(d / x)}  {x}\n" for x in names))
        # AMRFinderPlus database, mapping and catalogue
        db = root / "db"
        db.mkdir()
        ref = protein(120, 1)
        ref = ref[:82] + "S" + ref[83:]
        (db / "AMRProt-mutation.tsv").write_text(
            "#taxgroup\taccession_version\tstandard_mutation_symbol\tmutation_position\tclass\tsubclass\n"
            "Escherichia\tWP_GYRA\tgyrA_S83L\t83\tQUINOLONE\tQUINOLONE\n"
            "Salmonella\tWP_GYRA\tgyrA_D87N\t87\tQUINOLONE\tQUINOLONE\n")
        (db / "AMR_DNA-Escherichia.tsv").write_text(
            "#accesion_version@gene_name@gene_symbol:start-stop\tstandard_mutation_symbol\tmutation_position\t"
            "class\tsubclass\nREF_AMPC\tampC_C-42T\t12\tBETA-LACTAM\tCEPHALOSPORIN\n")
        (db / "fam.tsv").write_text("#node_id\tfamily_name\nblaTEM\tclass A beta-lactamase\nacrF\tefflux pump\n")
        (db / "AMRProt.fa").write_text(f">WP_GYRA|1|gyrA\n{ref}*\n>WP_OTHER|2|x\nMAAA\n")
        (db / "AMR_DNA-Escherichia.fa").write_text(f">REF_AMPC\n{rnd_dna(40, 9)}\n")
        mapping = root / "mapping.tsv"
        rows = [mapping_row("blaTEM-1", "ampicillin", "gene", family="blaTEM"),
                mapping_row("gyrA_S83L", "tetracycline"), mapping_row("ampC_C-42T", "ampicillin", interp="uncertain")]
        write_rows(mapping, list(rows[0].keys()), rows)
        # calls: one isolate without calls, one partial hit, a stress element that is not read
        calls = root / "calls"
        manifest = {r["isolate_id"]: r for r in ed.pop.read_tsv(root / "genome_manifest.tsv")}
        files = {k: sha(db / k) for k in ("AMRProt-mutation.tsv", "AMRProt.fa", "AMR_DNA-Escherichia.tsv",
                                          "fam.tsv")}
        for i in range(1, n):
            iid, gca = f"PDT{i:09d}.1", f"GCA_{i:09d}.1"
            report = "\t".join(REPORT_HEAD) + "\n"
            if gene_present[i]:
                report += report_line(gca, 1651, 2050, "blaTEM-1", "AMR", "AMR",
                                      "PARTIALP" if i == 3 else "BLASTX")
            if mutant[i]:
                report += report_line(gca, 1, 300, "gyrA_S83L", "AMR", "POINT", "POINTX")
            report += report_line(gca, 10, 100, "qacEdelta1", "STRESS", "BIOCIDE", "BLASTX")
            (calls / iid).mkdir(parents=True)
            (calls / iid / "amrfinder.tsv").write_text(report)
            (calls / iid / "provenance.json").write_text(json.dumps({
                "caller_version": "4.2.7", "status": "success", "input_sha256": manifest[iid]["genome_sha256"],
                "raw_sha256": hashlib.sha256(report.encode()).hexdigest(), "database_files_sha256": files}))
        # AST provenance for a subset of endpoints
        ast = root / "ast"
        ast.mkdir()
        keys = [(e["isolate_id"], e["antimicrobial"]) for e in eps if e["antimicrobial"] in ("ampicillin",
                                                                                             "tetracycline")]
        keys = [k for j, k in enumerate(keys) if j % 7]
        write_rows(ast / "isolate_drug_endpoints.tsv", ["isolate_id", "antimicrobial", "n_ast_observations",
                                                        "duplicate_status", "conflict_status"],
                   [{"isolate_id": a, "antimicrobial": b, "n_ast_observations": "2" if j % 11 == 0 else "1",
                     "duplicate_status": "none", "conflict_status": "conflict" if j % 19 == 0 else "none"}
                    for j, (a, b) in enumerate(keys)])
        write_rows(ast / "ast_observations.tsv", ["isolate_id", "antimicrobial", "testing_method",
                                                  "breakpoint_standard", "breakpoint_version", "mic"],
                   [{"isolate_id": a, "antimicrobial": b, "testing_method": "" if j % 5 == 0 else "broth dilution",
                     "breakpoint_standard": ("CLSI", "EUCAST")[j % 2] if j % 9 else "",
                     "breakpoint_version": "2020" if j % 3 else "", "mic": "4"} for j, (a, b) in enumerate(keys)])
        # analysis config pointing at the synthetic inputs
        config = json.loads(json.dumps(CONFIG))
        up = config["upstream"]
        up.update({"generalisation_config": str(study.config_path), "generalisation_config_sha256":
                   sha(study.config_path), "benchmark_manifest_sha256": sha(bench / "benchmark_manifest.json"),
                   "predictions_sha256": sha(bench / "predictions.tsv.gz"),
                   "feature_manifest_sha256": sha(bench / "matrix" / "feature_manifest.json"),
                   "mapping_rules": str(mapping), "mapping_rules_sha256": sha(mapping),
                   "genome_manifest": str(root / "genome_manifest.tsv"),
                   "genome_manifest_sha256": sha(root / "genome_manifest.tsv"), "amrfinder_database": str(db),
                   "amrfinder_database_files_sha256": files})
        config["catalogue"].update({"path": str(root / "catalogue.tsv"), "sha256": "0" * 64})
        config["augmentation"]["min_carriers"] = 10
        config["bootstrap"]["replicates"] = 20
        config["support"] = {"min_R": 5, "min_S": 5}
        cls.config_path = root / "analysis.json"
        cls.config_path.write_text(json.dumps(config))
        out = cls.out = root / "decomposition"
        common = ["--config", str(cls.config_path), "--outdir", str(out)]
        inputs = ["--benchmark", str(bench), "--population", str(pop_dir), "--baseline", str(base_dir)]
        cls.args = {"common": common, "inputs": inputs, "ast": ["--ast", str(ast)], "calls": ["--calls", str(calls)],
                    "cache": ["--cache", str(study.cache)]}
        with patch("sys.stdout"):
            codes = [ed.main(["catalogue", "--config", str(cls.config_path)])]
            config["catalogue"]["sha256"] = sha(root / "catalogue.tsv")
            cls.config_path.write_text(json.dumps(config))
            codes += [ed.main(["freeze", *common, "--benchmark", str(bench)]),
                      ed.main(["calls", *common, *cls.args["calls"]]),
                      ed.main(["audit", *common, "--benchmark", str(bench), *cls.args["cache"], "--workers", "2"]),
                      ed.main(["audit", *common, "--benchmark", str(bench), *cls.args["cache"]]),
                      ed.main(["context", *common, *inputs, *cls.args["ast"]]),
                      ed.main(["augment", *common, *inputs, "--only", "__random__"]),
                      ed.main(["augment", *common, *inputs]),
                      ed.main(["report", *common, *inputs])]
        cls.codes = codes
        cls.config = config

    @classmethod
    def tearDownClass(cls):
        cls.study.tearDownClass()

    def read(self, name):
        path = self.out / name
        if name.endswith(".gz"):
            with gzip.open(path, "rt") as h:
                return list(csv.DictReader(h, delimiter="\t"))
        return ed.read_tsv(path)

    def test_commands_succeed(self):
        self.assertEqual(self.codes, [0] * 9)
        manifest = json.loads((self.out / "error_decomposition_manifest.json").read_text())
        self.assertEqual(set(manifest["drug_groups"]), {"ampicillin", "tetracycline"})
        self.assertEqual(manifest["upstream"]["predictions_sha256"], self.config["upstream"]["predictions_sha256"])

    def test_catalogue_built_from_mapping_and_database(self):
        cat = ed.read_tsv(self.root / "catalogue.tsv")
        self.assertEqual([(r["antimicrobial"], r["determinant"], r["mechanism_class"]) for r in cat],
                         [("ampicillin", "blaTEM-1", "acquired_gene"), ("ampicillin", "ampC_C-42T", "promoter_mutation"),
                          ("tetracycline", "gyrA_S83L", "target_site_mutation")])
        self.assertEqual(cat[0]["gene_function"], "other")

    def test_calls_and_capture(self):
        status = Counter(r["status"] for r in self.read("call_status.tsv"))
        self.assertEqual(status, Counter({"success": 149, "calls_missing": 1}))
        self.assertNotIn("qacEdelta1", {c["element_symbol"] for c in self.read("amr_calls.tsv.gz")})
        known = self.read("known_mutations.tsv.gz")
        self.assertTrue(known and all(r["antimicrobial"] == "tetracycline" for r in known))
        genes = Counter(r["status"] for r in self.read("determinant_capture.tsv.gz"))
        self.assertGreater(genes["captured"], 0)
        muts = self.read("mutation_capture.tsv.gz")
        self.assertTrue(muts and {r["detail"] for r in muts} == {"site_not_found"})
        self.assertTrue(all(r["detail"] for r in self.read("genome_verification.tsv")))
        verified = Counter(r["genome_status"] for r in self.read("genome_verification.tsv"))
        self.assertEqual(verified, Counter({"verified": 150}))

    def test_context_reconciles_with_predictions(self):
        ctx = self.read("prediction_error_context.tsv.gz")
        with gzip.open(self.study.out / "predictions.tsv.gz", "rt") as h:
            preds = list(csv.DictReader(h, delimiter="\t"))
        self.assertEqual(len(ctx), len(preds))
        key = lambda r: (r["cohort"], r["antimicrobial"], r["split_strategy"], r["fold_or_holdout_id"],  # noqa: E731
                         r["isolate_id"])
        self.assertEqual({key(r) for r in ctx}, {key(r) for r in preds})
        by = {key(r): r for r in preds}
        for r in ctx:
            self.assertEqual(r["predicted_probability"], by[key(r)]["predicted_probability"])
        cats = Counter(r["error_category"] for r in ctx)
        # the synthetic mutation site cannot be located, so its capture is unknown and the error unresolved
        fn_mut = [r for r in ctx if r["error_type"] == "FN" and r["mutation_determinants"] == "gyrA_S83L"]
        self.assertTrue(fn_mut and all(r["error_category"] == "R_UNRESOLVED" and
                                       r["mutation_capture"] == "representation_unknown" for r in fn_mut))
        self.assertIn("R_NO_ESTABLISHED_DET", cats)
        self.assertIn("S_KNOWN_DET_PRESENT", cats)
        self.assertTrue(any(r["established_determinant_status"] == "uncertain" for r in ctx))
        self.assertTrue(any(r["established_determinant_status"] == NA for r in ctx))
        self.assertTrue(any("conflict" in r["phenotype_flags"] for r in ctx))
        self.assertTrue(any(r["phenotype_flags"] == "provenance_missing" for r in ctx))
        self.assertTrue(any(r["partial_amrfinder_calls"] == "yes" for r in ctx))

    def test_augmentation_on_frozen_splits(self):
        plan = json.loads((self.out / "augmented" / "augmentation_plan.json").read_text())
        self.assertEqual(plan["drugs"], ["tetracycline"])
        self.assertEqual(plan["indicators"], ["gyrA_S83L"])
        self.assertIn("no point-mutation rule", plan["eligibility"]["ampicillin"]["reason"])
        manifest = json.loads((self.study.out / "benchmark_manifest.json").read_text())
        frozen = manifest["task_split_sha256"]
        with gzip.open(self.study.out / "predictions.tsv.gz", "rt") as h:
            original = Counter("__".join((r["cohort"], r["antimicrobial"], r["split_strategy"],
                                          r["fold_or_holdout_id"])) for r in csv.DictReader(h, delimiter="\t"))
        for variant in ("mutation_augmented", "permuted_control"):
            for st in sorted((self.out / "augmented" / variant / "tasks").glob("*.json")):
                rec = json.loads(st.read_text())
                self.assertEqual(rec["status"], "success")
                self.assertEqual(rec["split_sha256"], frozen[rec["task_id"]])
                rows = ed.read_tsv(st.with_suffix(".tsv"))
                self.assertEqual(len(rows), original[rec["task_id"]])

    def test_augmentation_reuses_finished_tasks(self):
        with patch("sys.stdout"):
            self.assertEqual(ed.main(["augment", *self.args["common"], *self.args["inputs"],
                                      "--variants", "mutation_augmented"]), 0)
        for st in (self.out / "augmented" / "mutation_augmented" / "tasks").glob("*.json"):
            self.assertEqual(json.loads(st.read_text())["status"], "success")

    def test_augmented_metrics_and_rescue(self):
        metrics = self.read("mutation_augmented_metrics.tsv")
        self.assertTrue(metrics and {r["variant"] for r in metrics} == {"mutation_augmented", "permuted_control"})
        ctx = [r for r in self.read("prediction_error_context.tsv.gz") if r["cohort"] == "common"
               and r["antimicrobial"] == "tetracycline" and r["split_strategy"] == "random"]
        y = np.array([r["phenotype"] == "R" for r in ctx])
        p = np.array([float(r["predicted_probability"]) for r in ctx])
        row = next(r for r in metrics if r["variant"] == "mutation_augmented" and r["split_strategy"] == "random")
        self.assertAlmostEqual(float(row["original_AUROC"]), gb.auroc(y, p), places=5)
        self.assertGreater(float(row["augmented_AUROC"]), float(row["original_AUROC"]))
        control = next(r for r in metrics if r["variant"] == "permuted_control" and r["split_strategy"] == "random")
        self.assertLess(float(control["augmented_AUROC"]), float(row["augmented_AUROC"]))
        for r in self.read("generalisation_gap_rescue.tsv"):
            gap = float(r["original_delta_AUROC"]) - float(r["augmented_delta_AUROC"])
            self.assertAlmostEqual(float(r["gap_rescue"]), gap, places=5)

    def test_report_tables(self):
        for name in ("representation_audit.tsv", "mutation_capture.tsv", "determinant_capture_summary.tsv",
                     "error_decomposition.tsv.gz", "error_decomposition_summary.tsv", "four_way_comparison.tsv",
                     "determinant_correct_ml_wrong.tsv.gz", "ml_correct_determinant_wrong.tsv.gz",
                     "determinant_correct_ml_wrong_summary.tsv", "population_vs_mechanism_novelty.tsv",
                     "performance_by_mechanism.tsv", "mechanism_train_test_overlap.tsv",
                     "mechanism_prevalence_by_task.tsv", "mechanism_shift_summary.tsv", "phenotype_quality_audit.tsv",
                     "genome_quality_audit.tsv", "explain_generalisation_drop.tsv", "manuscript_table.tsv",
                     "mechanism_catalogue_snapshot.tsv"):
            with self.subTest(name=name):
                self.assertTrue(self.read(name), name)
        summary = json.loads((self.out / "error_decomposition_summary.json").read_text())
        ctx = self.read("prediction_error_context.tsv.gz")
        self.assertEqual(summary["prediction_rows"], len(ctx))
        self.assertEqual(summary["errors"], sum(r["error_type"] in ("FN", "FP") for r in ctx))
        summ = [r for r in self.read("error_decomposition_summary.tsv") if r["antimicrobial"] != "ALL"]
        self.assertEqual(sum(int(r["errors"]) for r in summ), summary["errors"])
        fw = [r for r in self.read("four_way_comparison.tsv") if r["antimicrobial"] != "ALL"]
        self.assertEqual(sum(int(r["determinant_only_correct"]) for r in fw),
                         len(self.read("determinant_correct_ml_wrong.tsv.gz")))
        table = self.read("manuscript_table.tsv")
        amp = [r for r in table if r["antimicrobial"] == "ampicillin"]
        self.assertTrue(amp and all(r["augmented_AUROC"] == NA for r in amp))

    def test_report_is_deterministic(self):
        before = (self.out / "error_decomposition_summary.tsv").read_bytes()
        with patch("sys.stdout"):
            ed.main(["report", *self.args["common"], *self.args["inputs"]])
        self.assertEqual((self.out / "error_decomposition_summary.tsv").read_bytes(), before)

    def test_changed_config_or_inputs_stop(self):
        changed = json.loads(self.config_path.read_text())
        changed["novelty"]["min_training_carriers"] = 3
        path = self.root / "changed.json"
        path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError, "changed after freezing"):
            ed.main(["report", "--config", str(path), "--outdir", str(self.out), *self.args["inputs"]])
        bad = json.loads(self.config_path.read_text())
        bad["upstream"]["predictions_sha256"] = "0" * 64
        path.write_text(json.dumps(bad))
        with self.assertRaises(ed.UpstreamError), patch("sys.stdout"):
            ed.main(["freeze", "--config", str(path), "--outdir", str(self.root / "other"),
                     "--benchmark", str(self.study.out)])

    def test_cli_requires_inputs(self):
        with patch("sys.stderr"):
            with self.assertRaises(SystemExit):
                ed.main(["context", *self.args["common"], *self.args["inputs"]])
            with self.assertRaises(SystemExit):
                ed.main(["report", "--config", str(self.config_path)])


if __name__ == "__main__":
    unittest.main()
