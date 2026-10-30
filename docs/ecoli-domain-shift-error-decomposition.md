# E. coli domain-shift error decomposition (error_decomposition_v1)

This analysis takes the frozen predictions of the generalisation benchmark and
asks why some antimicrobials lose performance under domain holdouts and others
do not. It does not change a split, label, endpoint or model of that benchmark.
Design: [decision 0008](decisions/0008-ecoli-error-decomposition.md).

PR #15 tests whether antimicrobial-specific domain-shift failures can be
explained by resistance-mechanism distribution and genomic feature
representation. The analysis does not establish causal resistance mechanisms.

## 1. What the generalisation benchmark found

The benchmark ([guide](ecoli-generalisation-benchmark.md),
[decision 0007](decisions/0007-ecoli-generalisation-benchmark.md)) compared a
random split with lineage, country, project, USA-to-other and temporal
holdouts for an L1 logistic regression on FracMinHash k-mers. Random splitting
was mildly optimistic overall (median common-cohort ΔAUROC 0.005 to 0.025 per
holdout). Clonal/lineage leakage was measurable in PR #14 but explained little
of the observed performance difference. The larger drops were confined to a
few drugs: nalidixic acid, cefazolin, cefepime, ceftazidime, ciprofloxacin,
cefoxitin and piperacillin-tazobactam. The curated determinant calls of the
PR #7 baseline tied or beat the model under every holdout.

## 2. Why this analysis

Several of the drugs with large drops have resistance that comes from point
mutations (gyrA, parC and parE for the quinolones; the ampC promoter for some
cephalosporin resistance), while the features are k-mers sampled at 1 in 200,
which see about one single-base change in seven. Three explanations were open
and are not exclusive:

- feature-representation failure: the model cannot see a determinant that is
  in the genome;
- resistance-mechanism distribution shift: the held-out population carries
  mechanisms that are rare or absent in training;
- domain generalisation failure: the mechanisms are shared and visible, but
  the model still does worse on the new population (background, sampling,
  phenotype practice).

## 3. Frozen upstream

| Input | SHA-256 |
| --- | --- |
| generalisation config (`generalisation_v1.json`) | `00fa371a…c65f` |
| benchmark manifest (file) | `f0b51c76…f303` |
| predictions (`predictions.tsv.gz`, 1,030,746 rows) | `63b1a2c0…0ce0` |
| feature manifest | `ec85fc5c…b15f` |
| determinant mapping v1 | `67f0eb67…795a` |
| genome manifest (PR #8, 10,521 rows) | `e85a45fc…f2d6` |
| PR #7 baseline manifest | `eccf2a71…c0d9` (through the benchmark manifest) |
| PR #5 cohort manifest | `21348c7d…7caf` (through the benchmark manifest) |
| AMRFinderPlus database 2026-08-07.1 | `AMRProt-mutation.tsv` `b6f44197…`, `AMRProt.fa` `41d5ebf4…`, `AMR_DNA-Escherichia.tsv` `495c2c35…`, `fam.tsv` `00323cb1…` |

Every stage checks these hashes and stops on a difference. The frozen PR #14
tasks are rebuilt from the same inputs and compared with the benchmark
manifest's split hashes; every prediction row is matched to exactly one test
endpoint of one frozen task, and a duplicate, missing or extra row stops the
analysis. The analysis config (`error_decomposition_v1.json`, SHA-256
`e97bd57d…af75`) and the catalogue were committed before any capture status
or error category was computed (commit `6024a4f`).

## 4. Mechanism catalogue

`benchmark/mappings/ecoli_amr_mechanisms_v1.tsv` (SHA-256 `c41ae6fd…eb34`)
takes every rule of determinant mapping v1 for the 24 study antimicrobials,
unchanged, and adds a mechanism class, the AMRFinderPlus locus of each point
mutation and the gene function. It has 43,890 rows: 5,428 genes (3,563 with a
`predicts_resistance` or `component` rule) and 63 point mutations (62 with
such a rule; 151 target-site and 53 promoter rows across drugs). Classes come
from fixed rules in the config: genes are `acquired_gene` except the
intrinsic blaEC family; mutations take the class of their gene; a mutation in
an unlisted gene stops the build. Nothing is inferred from model features or
phenotypes.

## 5. Mutation extraction

The AMRFinderPlus 4.2.7 calls behind PR #7 (`--organism Escherichia`,
database 2026-08-07.1) were read again. Each call set was checked against its
provenance record: caller version, database file hashes, the input genome's
SHA-256 against the PR #8 manifest, call status, and the report hash. All
10,521 passed. They hold 99,319 AMR calls, of which 16,680 are point
mutations (16,023 protein, 657 nucleotide; 105 symbols).
`known_mutations.tsv.gz` lists each called mutation with its catalogue status
per drug (36,524 rows).

## 6. Feature representation audit

The PR #14 features are sourmash 4.9.4 FracMinHash 31-mers (`scaled=200`,
`seed=42`), presence only. Of 921,736 observed hashes, 419,872 passed the
global 5/5 presence filter and were merged into 271,797 identical-presence
patterns. Each task then keeps only features present in at least 5 and absent
from at least 5 of its training isolates.

At scaled 200 a k-mer is kept with probability 1/200, so a single changed base
covered by 31 k-mers has about a 14% chance that at least one of them is
sampled, and less that the sampled one also passes the presence filter. Capture
is computed in each genome after checking its SHA-256 against the PR #8
manifest; phenotypes and model outputs are not read.

## 7. Determinant capture

A called gene is `captured` when at least one 31-mer inside it is a kept PR #14
feature, `partially_represented` when some of its k-mers are sampled but all
were pruned, and `not_represented` when none is sampled. Of 80,486 gene calls,
79,882 were captured, 88 partially represented and 516 not represented. The
acquired genes of every study drug were captured in more than 99% of
common-cohort carriers (`representation_audit.tsv`). A few short genes were
not represented: aac(6')-Ib-cr in 33 of 79 calls, arr-3 in 45 of 65.

## 8. Mutation capture

For a point mutation the site is located from the AMRFinderPlus gene
coordinates and the reference protein (or reference DNA for promoter and rRNA
mutations), the changed bases are found against the nearest wild-type codon,
and the 31-mers overlapping a changed base are hashed for the mutant allele
and for the wild-type allele. `captured` means a mutant-allele k-mer is a kept
feature, so the model can see the mutation directly. `partially_represented`
means no such k-mer, but a mutant k-mer was sampled and pruned or a wild-type
k-mer is kept, so the mutation shows only as a missing wild-type feature.
`not_represented` means no sampled k-mer at the site for either allele.
`representation_unknown` covers sites that could not be located.

Of 16,680 mutation calls, 257 were captured, 1,988 partially represented,
14,310 not represented and 125 unknown (non-catalogue insertions and
frameshifts, and five sites not found). The main quinolone mutations were
almost never captured: gyrA S83L in 0 of 2,298 calls (5 partial), parC S80I in
0 of 1,685 (24 partial), gyrA D87N 1,278 partial and 349 not represented,
parC E84V and parE I529L never sampled. Per drug, over common-cohort carriers:

| Antimicrobial | Mutation carriers | Captured | Partial | Not represented |
| --- | --- | --- | --- | --- |
| nalidixic acid | 1,230 | 1.1% | 14.0% | 84.9% |
| ciprofloxacin | 10,662 | 0.3% | 16.2% | 83.5% |
| levofloxacin | 8,294 | 0.2% | 15.6% | 84.2% |
| cefoxitin (ampC promoter) | 178 | 0 | 88.8% | 11.2% |
| ceftazidime | 196 | 0 | 52.0% | 48.0% |
| cefazolin | 106 | 0 | 22.6% | 77.4% |

(Counts are endpoint-mutation pairs, so an isolate with three mutations counts
three times.) For the same drugs the acquired genes were captured in 99 to
100% of carriers.

## 9. Mechanism novelty

For each task a determinant is seen in training when at least 5 training
endpoints carry it, whatever their phenotype. A test endpoint's mechanism is
seen when all its established determinants are, and unseen otherwise. Pooled
over common-cohort drugs, the share of resistant determinant-positive test
endpoints with an unseen mechanism was 2.1% under random splitting, 2.7%
lineage, 3.1% project, 7.0% geography, 8.2% temporal and 8.8% country. The
point mutations themselves were nearly always seen: 98% of resistant
mutation carriers under random, 95% lineage, 94% project, 87% geography, 99%
temporal and 82% country.

## 10. Population novelty against mechanism novelty

Sensitivity on resistant determinant-positive endpoints, pooled over drugs and
holdouts (common cohort, `population_vs_mechanism_novelty.tsv`):

| Population axis | Population seen | Mechanism seen | Mechanism unseen |
| --- | --- | --- | --- |
| lineage (ST) | yes | 0.950 (n = 54,931) | 0.636 (n = 2,236) |
| lineage (ST) | no | 0.938 (n = 29,328) | 0.564 (n = 1,764) |
| country | yes | 0.957 (n = 58,355) | 0.715 (n = 1,621) |
| country | no | 0.921 (n = 26,247) | 0.529 (n = 2,380) |
| project | yes | 0.958 (n = 42,582) | 0.821 (n = 1,030) |
| project | no | 0.934 (n = 42,020) | 0.529 (n = 2,971) |

An unseen lineage, country or project costs one to four points of
sensitivity when the mechanism is seen; an unseen mechanism costs 14 to 40
points whether or not the population is new. Across the 84 drug and holdout
pairs of the common cohort, the size of the AUROC drop correlated with the
unseen-mechanism fraction (Spearman ρ = 0.66) more than with mechanism
composition shift (0.36), mutation-associated fraction (0.36), prevalence
shift (0.33) or lineage shift (−0.18).

## 11. Country, project, temporal and lineage mechanism shift

Mechanism composition of resistant test endpoints differed from training
under every holdout (median total variation distance 0.15 random, 0.31
lineage, 0.36 temporal, 0.38 geography, 0.43 country, 0.45 project), and the
share of R endpoints moved much more under the domain holdouts (median
absolute shift 0.22 to 0.32) than under lineage (0.001).

- Country. The cefazolin drop (0.26) comes almost entirely from the USA
  holdout, whose training set is the 161 non-USA isolates: 552 of 725
  resistant determinant-positive test endpoints (76%) carry a mechanism with
  fewer than 5 training carriers. Canada as the held-out country has none.
- Project. Several projects enrich for resistance: in the cefepime holdout of
  PRJNA885502 (the Netherlands) 84% of test endpoints are R against 12% in
  training; in PRJNA809394 4% against 38%.
- Temporal. Later isolates (2021 on) are much more often resistant for
  cefepime (43% against 8%), ceftazidime (70% against 14%) and ciprofloxacin
  (57% against 10%); unseen mechanisms are 5% (ciprofloxacin) to 18%
  (cefepime) of the resistant test endpoints.
- Lineage. Mechanisms are mostly shared (2.7% unseen), and the drops are
  small except for nalidixic acid (0.10) and piperacillin-tazobactam (0.11;
  70 resistant test endpoints, 20% unseen).

## 12. Error decomposition

Each false negative and false positive of the frozen predictions is given one
category from the frozen rules, and the evidence columns (capture, novelty,
QC) are kept beside it, so an endpoint with both a gene and a mutation is not
forced into one explanation. Over all cohorts and tasks there are 50,506
errors in 1,030,746 predictions:

| Category | Errors |
| --- | --- |
| R_KNOWN_DET_REPRESENTED | 8,726 |
| R_NO_ESTABLISHED_DET | 7,183 |
| R_KNOWN_MUT_NOT_REPRESENTED | 4,796 |
| R_UNRESOLVED | 805 |
| R_KNOWN_DET_AND_MUT_REPRESENTED | 290 |
| R_KNOWN_DET_NOT_REPRESENTED | 41 |
| R_KNOWN_MUT_REPRESENTED | 27 |
| S_NO_KNOWN_DET | 13,280 |
| S_KNOWN_DET_PRESENT | 13,139 |
| S_UNRESOLVED | 2,219 |

Share of false negatives by category, common cohort, pooled over drugs:

| Strategy | FN | known determinant, captured | known mutation, not represented | no established determinant | unresolved |
| --- | --- | --- | --- | --- | --- |
| random | 1,301 | 0.31 | 0.16 | 0.48 | 0.04 |
| lineage | 1,857 | 0.36 | 0.24 | 0.35 | 0.04 |
| project | 1,743 | 0.37 | 0.23 | 0.34 | 0.04 |
| country | 2,497 | 0.50 | 0.20 | 0.25 | 0.03 |
| geography | 1,140 | 0.44 | 0.27 | 0.24 | 0.04 |
| temporal | 676 | 0.48 | 0.28 | 0.15 | 0.06 |

The holdouts add false negatives in two places: resistant endpoints whose
mutation the features do not encode, and resistant endpoints whose
determinant is captured but which the model still calls S. Errors without an
established determinant are about as many in number under random and holdout
splits, so their share falls. Per drug, pooled over holdouts, the mutation
category dominates the false negatives of nalidixic acid (77%) and
ciprofloxacin (47%); for cefepime (84%), ceftazidime (72%) and cefazolin (68%)
the determinant was present and captured; for chloramphenicol, gentamicin and
tetracycline most false negatives (59 to 64%) have no established
determinant.

Endpoints the determinant rule got right and the model got wrong
(`determinant_correct_ml_wrong.tsv.gz`) grow from 1,745 under random to 3,268
under country holdout, while the converse (`ml_correct_determinant_wrong.tsv.gz`)
falls from 2,311 to 1,342 (common cohort, `four_way_comparison.tsv`). Both
were wrong on 1,605 to 1,835 endpoints under the random, lineage, country and
project strategies.

Resistant endpoints without an established determinant are reported as
unexplained resistance: not explained by the current established-determinant
mapping. They are not evidence of a new mechanism.

## 13. Stable-drug controls

The same audit was run for every drug. For chloramphenicol, gentamicin and
tetracycline (largest drop 0.002, 0.013 and 0.008) resistance is carried by
acquired genes, which the features capture in more than 99% of carriers; no
catalogue mutation is involved; at most 5% of resistant determinant-positive
test endpoints have an unseen mechanism under any holdout (tetracycline at
most 1.1%); and the determinant baseline AUROC is close to the model's (0.95
to 0.99). This is case D of the frozen decision tree.

Levofloxacin is the counter-example. It is in the stable group (largest drop
0.007) although its mutations are as poorly represented as ciprofloxacin's (84%
not represented). Its country holdout has only 86 test endpoints, its random
AUROC is already 0.987, and augmentation changes little. Poor representation of a mutation therefore does
not by itself produce a drop; it does when the population the mutation sits in
changes between training and test.

## 14. Mutation augmentation

This is a sensitivity analysis labelled `mutation_augmented`; the PR #14
benchmark results are unchanged and remain the primary results. The PR #14
feature matrix was extended by 62 columns, one 0/1 indicator per catalogue
point mutation with a `predicts_resistance` or `component` rule, set from the
AMRFinderPlus calls (partial hits excluded). The indicators come from
established external knowledge (mapping v1 and the AMRFinderPlus database),
not from phenotypes. The same columns are used for every drug.

A drug was refitted if its catalogue has such a rule, at least 20 common-cohort
endpoints carry the mutation, and fewer than 90% of carriers have it captured.
Ten drugs qualified: amoxicillin-clavulanic acid, ampicillin, cefazolin,
cefoxitin, ceftazidime, ceftriaxone, ciprofloxacin, levofloxacin, nalidixic
acid and trimethoprim-sulfamethoxazole. Sulfisoxazole (19 carriers) and
ceftiofur (14) did not. The 262 common-cohort tasks of these drugs were refitted
with the exact frozen splits (checked against the PR #14 split hashes), the
same model, C grid, inner CV, seeds, feature filter and threshold, with no new
tuning. A permuted control repeats the refit with the indicator rows shuffled
across isolates (seed 20261029), keeping prevalence but breaking the link to
each genome. All 524 refits succeeded.

## 15. Generalisation-gap rescue

Gap rescue is the original (random − holdout) ΔAUROC minus the augmented one,
on the same test isolates, with a 95% interval from 1,000 paired bootstrap
resamples of split units (seed 20261028). A rescue counts only where its
interval is above zero and the permuted control is near zero.

| Antimicrobial | Holdout | Test N | Original Δ | Augmented Δ | Rescue (95% CI) | Control rescue |
| --- | --- | --- | --- | --- | --- | --- |
| nalidixic acid | country | 4,509 | 0.110 | 0.053 | 0.057 (0.018 to 0.096) | 0.001 |
| nalidixic acid | geography | 1,705 | 0.175 | 0.060 | 0.115 (0.053 to 0.186) | 0.004 |
| nalidixic acid | lineage | 4,522 | 0.104 | 0.006 | 0.098 (0.054 to 0.151) | 0.001 |
| nalidixic acid | project | 4,252 | 0.177 | 0.058 | 0.119 (0.074 to 0.167) | 0.001 |
| nalidixic acid | temporal | 343 | 0.131 | 0.008 | 0.123 (0.051 to 0.197) | 0.015 |
| ciprofloxacin | country | 6,173 | 0.067 | 0.021 | 0.046 (0.035 to 0.059) | 0.002 |
| ciprofloxacin | geography | 1,784 | 0.079 | 0.023 | 0.056 (0.023 to 0.087) | −0.010 |
| ciprofloxacin | lineage | 7,529 | 0.006 | 0.001 | 0.004 (−0.000 to 0.009) | 0.000 |
| ciprofloxacin | project | 6,106 | 0.027 | 0.014 | 0.013 (0.005 to 0.022) | 0.001 |
| ciprofloxacin | temporal | 1,064 | 0.068 | 0.038 | 0.030 (0.016 to 0.046) | 0.000 |
| cefoxitin | temporal | 751 | 0.075 | 0.031 | 0.044 (0.021 to 0.066) | −0.007 |
| cefoxitin | project | 4,811 | 0.040 | 0.017 | 0.023 (0.014 to 0.033) | 0.017 |
| amoxicillin-clavulanic acid | lineage | 6,356 | 0.027 | 0.007 | 0.020 (0.014 to 0.027) | −0.001 |
| amoxicillin-clavulanic acid | project | 6,239 | 0.040 | 0.013 | 0.026 (0.017 to 0.037) | 0.018 |
| cefazolin | country | 2,519 | 0.260 | 0.263 | −0.004 (−0.007 to −0.001) | 0.022 |
| ceftazidime | geography | 733 | 0.166 | 0.159 | 0.007 (0.002 to 0.014) | −0.011 |
| ceftazidime | temporal | 783 | 0.222 | 0.218 | 0.005 (−0.002 to 0.011) | −0.000 |

The full table, with every drug and holdout, is `generalisation_gap_rescue.tsv`.
Holdout performance before and after augmentation for the quinolones
(`mutation_augmented_metrics.tsv`):

| Antimicrobial | Holdout | AUROC | AUPRC | Balanced accuracy | MCC |
| --- | --- | --- | --- | --- | --- |
| nalidixic acid | project | 0.698 → 0.915 | 0.326 → 0.863 | 0.654 → 0.887 | 0.287 → 0.833 |
| nalidixic acid | geography | 0.693 → 0.911 | 0.361 → 0.835 | 0.616 → 0.876 | 0.273 → 0.834 |
| nalidixic acid | lineage | 0.776 → 0.968 | 0.402 → 0.898 | 0.699 → 0.913 | 0.324 → 0.864 |
| nalidixic acid | random | 0.880 → 0.974 | 0.631 → 0.918 | 0.820 → 0.932 | 0.500 → 0.867 |
| ciprofloxacin | country | 0.906 → 0.964 | 0.797 → 0.905 | 0.856 → 0.924 | 0.658 → 0.807 |
| ciprofloxacin | geography | 0.892 → 0.962 | 0.816 → 0.935 | 0.790 → 0.820 | 0.673 → 0.756 |
| ciprofloxacin | temporal | 0.886 → 0.936 | 0.914 → 0.950 | 0.826 → 0.866 | 0.648 → 0.729 |

For nalidixic acid the mutation indicators remove 50 to 95% of the holdout gap
and bring the holdout AUROC (0.91 to 0.97) to about the level of the
determinant baseline (0.89 to 0.92) or above it. They also raise the random-split AUROC (0.88 to
0.97): the representation limit costs performance everywhere and more under
holdouts. For ciprofloxacin they remove 44 to 71% of the country, geography
and temporal gaps. For cefoxitin (ampC promoter) they remove 59% of the
temporal gap. Two rescues (cefoxitin and amoxicillin-clavulanic acid, project)
are matched by the permuted control and are read as an effect of adding
columns, not of the mutations. For cefazolin, cefepime (not eligible: no
catalogue mutation) and ceftazidime there is no meaningful rescue.

## 16. Phenotype QC

Over the common-cohort holdouts, errors are more often from a project whose
labels were corrected (decision 0004; 4.4% of errors against 1.7% of correct
predictions) and more often lack a breakpoint standard (3.4% against 1.4%).
Conflicting or repeated observations are rare (under 0.1%). Testing method is
not recorded in the processed AST table for any endpoint, so it cannot be
audited. For ceftazidime the country-holdout errors were more often tested
under EUCAST than the correct predictions (35% against 12%); for cefepime and
ceftazidime the temporal-holdout errors and correct predictions both mix CLSI
and EUCAST (33 to 63% EUCAST). Breakpoints were not reinterpreted. These
differences are consistent with, but do not show, a phenotype-definition
shift for the extended-spectrum cephalosporins.

## 17. Genome QC

Errors and correct predictions have the same median assembly length (5.14 and
5.11 Mb), contig count (147 and 146) and N50 (131 and 129 kb). Errors are
slightly more often QC-flagged (0.44% against 0.27%) and more often have a
partial AMRFinderPlus hit for a catalogue gene (5.1% against 1.9%). Genome
quality does not account for a material share of errors. No genome was
removed.

## 18. Limitations

- The representation result holds for this feature set: FracMinHash at scaled
  200 with a 5/5 presence filter. A denser k-mer, unitig or gene-plus-variant
  representation would see the mutations, and the gap it leaves has not been
  measured here.
- One model, one cohort and one species. The cohort is dominated by USA
  surveillance and food-animal projects.
- Capture describes the features, not the fitted model. A captured k-mer may
  still have zero weight, and the model may use linked k-mers as proxies.
- Mechanism novelty depends on mapping v1. Unmapped or uncertain determinants
  count as no established determinant or unresolved; unseen mechanisms are
  unseen within these frozen tasks, not new to biology.
- Country, project, source and year are confounded; the holdout effects are
  for these test populations.
- The augmentation adds known mutations only. It cannot help where the
  missing signal is an unmapped mechanism or a phenotype difference.
- Testing method is missing for all endpoints, and the breakpoint audit is
  descriptive.
- The decomposition is descriptive. Associations across 84 drug and holdout
  pairs share isolates and are not independent.

## 19. Interpretation

Against the frozen decision tree:

- Nalidixic acid and ciprofloxacin are case A. The large drops come with
  mutations that the features almost never encode, which are nevertheless
  seen in training, and adding them as indicators rescues most of the
  nalidixic-acid gap and roughly half of the ciprofloxacin gap while the
  permuted control does not. A meaningful part of the apparent domain-shift
  failure for these drugs is a feature-representation failure.
- Cefazolin is case B. Its country drop is driven by a holdout whose training
  set is small and does not contain most of the test population's mechanisms.
- Cefepime and ceftazidime are closest to case C, with a phenotype component
  that this analysis cannot resolve. Their determinants are captured and
  mostly seen; augmentation does not help; the share of resistant isolates
  moves strongly between domains, and many errors are susceptible isolates
  carrying a known determinant (S_KNOWN_DET_PRESENT). The determinant baseline
  still beats the model under most of their holdouts.
- Piperacillin-tazobactam is unresolved: 70 resistant test endpoints in its
  lineage holdout, 20% with unseen mechanisms, no augmentation.
- Chloramphenicol, gentamicin and tetracycline are case D.

So the three explanations each hold for different drugs. Lineage novelty on
its own explains little; mechanism novelty and representation explain more.

## 20. Implications

- For mutation-driven resistance, a sampled k-mer representation is not a
  fair test of genomic prediction. A next model should carry known variants
  explicitly or use a representation that sees every base, and the frozen
  benchmark should be repeated with it as a confirmatory analysis.
- For acquired-gene resistance the k-mer features already see the
  determinants; there, better models are unlikely to beat curated rules
  without new biological information.
- The R/determinant-negative endpoints (7,183 false negatives) and the
  S/determinant-positive ones (13,139 false positives) are the material for
  discovery work, which needs independent replication and is not started
  here.

## 21. Reproducibility

With the PR #14 output in `datasets/ecoli_research_v1/generalisation_v1`, the
PR #6, PR #7 and PR #5 builds, the genome cache and the AMRFinderPlus calls in
`datasets/ecoli_research_v1/`, run the commands in
`benchmark/cohorts/ecoli_research_v1/error_decomposition_v1/run_record.json`
with the environment in `environment/generalisation-benchmark.txt`. The audit
and augmentation resume from their checkpoints; the augmentation reuses a
task only if its split hash still matches. The run took about 4 hours on a
16 GB laptop (augmentation 3.6 hours, peak 6.9 GB). Deviations and the code
corrections after freezing are in the run record and in decision 0008.

## 22. Output files

Committed in `benchmark/cohorts/ecoli_research_v1/error_decomposition_v1/`
(with `SHA256SUMS`):

| File | Content |
| --- | --- |
| `error_decomposition_manifest.json` | upstream hashes, drug groups and case studies, frozen before any mechanism table |
| `call_status.tsv`, `known_mutations.tsv.gz` | provenance check per isolate; called mutations with catalogue status |
| `determinant_capture.tsv.gz`, `mutation_capture.tsv`, `determinant_capture_summary.tsv` | k-mer level capture per call |
| `representation_audit.tsv` | capture per drug and determinant over common-cohort endpoints |
| `genome_verification.tsv` | genome hash check |
| `mechanism_train_test_overlap.tsv`, `mechanism_prevalence_by_task.tsv`, `mechanism_shift_summary.tsv` | mechanism overlap and shift per task |
| `error_decomposition.tsv.gz`, `error_decomposition_summary.tsv` | every error with its context; counts per drug and strategy |
| `four_way_comparison.tsv`, `determinant_correct_ml_wrong.tsv.gz`, `ml_correct_determinant_wrong.tsv.gz`, `determinant_correct_ml_wrong_summary.tsv` | determinant rule against model |
| `population_vs_mechanism_novelty.tsv`, `performance_by_mechanism.tsv` | performance by population and mechanism novelty, mechanism group and capture |
| `phenotype_quality_audit.tsv`, `genome_quality_audit.tsv` | QC of errors against correct predictions |
| `augmentation_plan.json`, `mutation_augmented_metrics.tsv`, `generalisation_gap_rescue.tsv` | sensitivity analysis |
| `explain_generalisation_drop.tsv`, `manuscript_table.tsv`, `error_decomposition_summary.json` | drug-level summary |
| `run_record.json` | commands, timings, deviations |

The full context table (`prediction_error_context.tsv.gz`, about 76 MB,
1,030,746 rows), the call table, the audit checkpoint and the augmented
predictions stay in `datasets/ecoli_research_v1/error_decomposition_v1/`.
