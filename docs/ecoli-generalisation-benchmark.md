# E. coli AMR generalisation benchmark (generalisation_v1)

## 1. Question and hypothesis

How much does a conventional random train/test split overstate genomic AMR
prediction, compared with test sets made of unseen lineages, unseen
countries, unseen submitting projects and later isolates?

The hypothesis under test is that a model evaluated on a random split scores
better than the same model evaluated on biologically or epidemiologically
separated test populations. The main quantity is the drop

    Δmetric = metric(random) − metric(holdout)

for each drug and holdout, computed on the same test isolates. The benchmark
can come out against the hypothesis: no drug, threshold or rule was chosen
from model results.

Random-split and holdout performance are compared on a common eligible
cohort as the primary analysis wherever feasible. Split definitions and
eligibility rules were frozen independently of observed model performance.

## 2. Version

Benchmark `generalisation_v1`, config
`benchmark/cohorts/ecoli_research_v1/generalisation_v1.json` 1.0.0, code
`bin/generalisation_benchmark.py` 1.0.0, design in
[decision 0007](decisions/0007-ecoli-generalisation-benchmark.md). Any change of
rule, drug list, feature or model is a new benchmark version.

## 3. Upstream frozen data

| Input | Identity |
| --- | --- |
| PR #5 cohort | `cohort_manifest.json` SHA-256 `21348c7d…7caf` |
| PR #6 population build | `population_manifest.json` SHA-256 `ca4c0f65…3468`; `evaluation_groups.tsv` and `population_metadata.tsv` checked against its `SHA256SUMS` |
| PR #7 baseline (primary, corrected labels) | `baseline_manifest.json` SHA-256 `eccf2a71…c0d9`, mapping v1, correction table `a8a92d9c…eb71`; `genotype_phenotype_endpoints.tsv` checked against its `SHA256SUMS` |
| PR #8 genome manifest | `genome_manifest.tsv` SHA-256 `e85a45fc…f2d6`, 10,521 rows |
| PR #6 population config | `population_config.json` SHA-256 `df55c31d…9c0d`, for the decision 0001 tier rules |

The builder checks each hash and stops on any difference. It also checks that
the PR #6 and PR #7 builds name the same PR #5 manifest, that the PR #7 build
used the configured correction table and mapping, and that every endpoint's
assembly accession is the same in PR #6, PR #7 and PR #8.

Phenotypes are the `phenotype` column of the PR #7 endpoint table: the PR #5
label, with R and S swapped for project PRJNA1322038 (decision 0004). An
endpoint enters if it is binary R/S and its isolate is population-eligible
(PR #5 eligible, genome QC not failed).

## 4. Genome verification

`sketch` reads every genome in the PR #8 manifest once, computes its SHA-256 and
compares it with the manifest. Failures are classed as `genome_missing`,
`genome_hash_mismatch`, `genome_unreadable` or `genome_version_mismatch` (the
cohort's assembly differs from the manifest's versioned GCA) and are written to
`genome_verification.tsv`; they are never silently dropped. The exact GCA
version of PR #8 is used; newer NCBI versions are not looked up. Identical
genome bytes under two isolates are listed as duplicates.

## 5. Feature representation

Each verified genome is sketched with sourmash 4.9.4 (`MinHash(n=0, ksize=31,
scaled=200, seed=42)`, all records, non-ACGT k-mers skipped). The feature for a
k-mer is its presence. Which k-mers are features is fixed by the hash function,
so no vocabulary is learned from genomes or labels.

`matrix` builds an isolate × feature matrix over all verified genomes. It drops
features present in fewer than 5 or absent from fewer than 5 genomes, and
merges columns with identical presence in every genome into one pattern. Both
steps are exact with respect to every task: the per-task filter (below)
would drop the same features, and an L1-penalised model treats identical
columns as one. The matrix manifest records the sketch parameters, the
genome-verification hash, the row list, and hashes of the k-mer list, the
k-mer-to-pattern map and the matrix arrays.

## 6. Leakage prevention

- Features for a task are those present in at least 5 and absent from at least
  5 of its training isolates; test genomes are projected onto that set. A
  k-mer seen only in test genomes is never a feature. The task's feature list is
  recorded by hash (`vocabulary_sha256`).
- No label is used to choose features.
- C is chosen by inner cross-validation on the training rows only, with the
  feature filter recomputed in each inner fold.
- The threshold is fixed at 0.5.
- Every task's training and test isolate lists are hashed into the benchmark
  manifest at freeze time; `run` and `report` rebuild them and stop on any
  difference.

## 7. Model

L1-penalised logistic regression (scikit-learn 1.7.2, liblinear solver,
`class_weight=balanced`, tol 1e-4, up to 1,000 iterations, intercept scaling 1,
fixed seed). One model, no second model, no neural networks.

## 8. Hyperparameter policy

C ∈ {0.03, 0.1, 0.3}, chosen per task by 3-fold inner cross-validation on the
training rows, folds grouped by split unit (SNP cluster or isolate) and
balanced for R and S, by mean inner AUROC; ties go to the smaller C. The same
policy is used for every strategy. The grid was set from a timing run on the
largest training set (ampicillin, 6,229 isolates, 5–6 s per fit, 109 to 885
non-zero coefficients across the grid) before any prediction was scored.

## 9. Common cohort

For each drug: endpoints with a verified genome, a lineage component, a known
country and a valid collection year. All seven strategies run on this set. The
excluded endpoints and reasons are in `cohort_eligibility.tsv`.

## 10. Maximal cohorts

Each strategy also runs on every endpoint with a verified genome and the
metadata it needs: all for the random splits and project, a known country for
country and geography, a valid year for temporal. For lineage, isolates
without a component train but are never tested. Maximal results are
secondary, are labelled `cohort = maximal`, and are never mixed with
common-cohort results in one comparison.

## 11. Random splits

`random`: stratified 5-fold cross-validation over isolates. This is the
conventional split, and SNP clusters can fall on both sides. `random_cluster`:
the same over split units, so a cluster is never divided (the decision 0001
reference). Every common-cohort isolate gets exactly one out-of-fold prediction
from each. Units are assigned in order of size, ties in order of a SHA-256 key
of seed, drug and unit, each to the fold that keeps R and S loads most even.

## 12. Lineage holdout

5-fold cross-validation over decision 0001 lineage components (connected
components of shared ST or shared NCBI SNP cluster), assigned the same way.
Isolates without a component (no ST and no cluster) are training-only in the
maximal cohort and absent from the common cohort. Per-lineage descriptive
results come from the out-of-fold predictions.

## 13. Country holdout

Leave one country out for each country with at least 10 R and 10 S for the drug
(the decision 0001 per-country rule), training on all other countries.
Isolates without a country are not used. Primary drugs only. A separate
`geography` task trains on USA isolates and tests on all other countries
pooled, for primary drugs with `geographic_transfer = yes`. Because country,
project, source and year move together in this cohort, results are
geographic/domain-shift results, not country effects.

A leave-one-project-out task (`project`) is also run for each project with at
least 10 R and 10 S, primary drugs with `project_holdout = yes`. Decision 0001
names this as its main external-validity axis.

## 14. Temporal holdout

Train on isolates collected before 2021, test on 2021 on (decision 0001
cutoff), for primary drugs with `temporal = yes` (at least 20 R and 20 S in
each period). One cutoff, fixed before any model was run. Training years
always precede test years.

## 15. Split validation

`freeze` fails if any task breaks its rule:
- no isolate in both partitions, no duplicate isolate or assembly;
- every isolate has a verified genome;
- random_cluster: no split unit on both sides;
- lineage: no component, ST or SNP cluster on both sides, and no test isolate
  without a component;
- country: no test country or unknown country in training;
- project: no test project in training;
- geography: training all USA, test none;
- temporal: max(training year) < min(test year).

For CV strategies the folds must cover each test-eligible endpoint exactly
once and training plus test must equal the cohort. Feasibility minimums: 20 R
and 20 S in training for every holdout task; 10 R and 10 S in a held-out
country or project; 20 R and 20 S in the geography and temporal test sets.
Tasks that fail are in `not_feasible.tsv` with the reason.

## 16. Leakage audit

`split_leakage_audit.tsv` gives, for every task, the share of test isolates
whose ST, SNP cluster, lineage component, project, country and period also occur
in training, plus training and test year ranges. Country, project, geography
and temporal holdouts do not remove training isolates that share a SNP cluster
with the test set (that matches deployment); the overlap is reported.

## 17. Metrics

For every fold and for pooled predictions: N, R, S, AUROC, AUPRC (average
precision), sensitivity, specificity, PPV, NPV, balanced accuracy, F1, MCC,
Brier score, and calibration intercept and slope (logistic fit of the
phenotype on logit p, an evaluation measure only; no recalibration). Metrics
undefined for one-class sets are NA. Rows with fewer than 10 R or 10 S are kept
with `supported = no`. Balanced class weights centre the probabilities on 50%
prevalence, so calibration is reported but not interpreted as a main result.
PPV and NPV depend on the sampled prevalence and do not transfer.

## 18. Uncertainty

95% percentile intervals from 1,000 bootstrap replicates resampling split units
(SNP clusters, or isolates when unclustered), seed 20261027. Drops are
bootstrapped in pairs: the same resample for the random and the holdout
predictions of the same isolates. Fold-level results are kept with median, IQR
and range. The unit for comparing strategies is the drug or fold, not the
isolate; no significance tests are pooled over isolates.

## 19. Known-determinant comparison

The PR #7 mapping v1 calls, frozen, give a determinant prediction (positive =
R) for endpoints with an interpretable determinant status. `determinant_vs_ml.tsv`
compares it with the model on the same test endpoints for each strategy.
The baseline does not train, so on a fixed set of isolates it scores the same
under every strategy. Differences between strategies in its row come only from
the test population.

## 20. R/determinant-negative subgroup

`unexplained_resistance_ml.tsv` gives, per strategy, the number of R/DET- test
endpoints and how many the model calls R, overall and split by whether the
isolate's ST was seen in training; and the same for S/DET+ endpoints called S.
A correct call on an R/DET- endpoint is not evidence of a new mechanism: the
model may use lineage-linked k-mers.

## 21. Domain novelty

Each prediction row records whether the isolate's ST, SNP cluster and country
occur in that task's training set and whether its year is after the last
training year. `unseen_lineage_summary.tsv` gives metrics by ST seen or unseen
and by cluster seen or unseen. No composite novelty score is made.

## 22. Results

Run on 2026-10-24 and 2026-10-25 (`run_record.json`). All 10,521 genomes
matched their PR #8 SHA-256, with no duplicates. 8,194 isolates and 109,865
endpoints entered; the common cohorts hold 107,205 of them (excluded: 121
endpoints without a country, 10 without a lineage component, 3 without a
year). 1,069 tasks were frozen and all succeeded; 72 strategy and drug pairs
were not feasible (`not_feasible.tsv`: the 64 secondary-drug group holdouts,
no temporal task for cefazolin, streptomycin and sulfisoxazole, and no USA
country holdout or USA-to-other task for levofloxacin, at 19 R / 79 S). No
pooled result is unsupported. Three random, three random_cluster and four
lineage folds have fewer than 10 R or 10 S and are marked in
`fold_metrics.tsv`.

The tables below are for the common cohort. Maximal-cohort results are in the
same files with `cohort = maximal`; the medians and the drugs with the
largest drops are much the same.

### 22.1 Drop against the random split

Δ is pooled AUROC of the random-split predictions minus that of the holdout
predictions, on the same test isolates, with 95% paired bootstrap intervals
(`generalisation_drop.tsv`). Positive values mean the random split scored
higher. (s) marks secondary drugs, which have only the lineage holdout.

| Drug | N (R) | Random AUROC | Δ lineage | Δ country | Δ project | Δ geography | Δ temporal |
| --- | --- | --- | --- | --- | --- | --- | --- |
| amikacin (s) | 3215 (47) | 0.891 | 0.009 (−0.064 to 0.075) |  |  |  |  |
| amoxicillin-clavulanic acid | 6356 (870) | 0.969 | 0.027 (0.017 to 0.037) | 0.046 (0.032 to 0.061) | 0.040 (0.028 to 0.053) | 0.002 (−0.004 to 0.008) | 0.041 (0.014 to 0.071) |
| ampicillin | 7566 (3572) | 0.986 | 0.004 (0.002 to 0.007) | 0.038 (0.031 to 0.045) | 0.022 (0.016 to 0.027) | 0.011 (0.005 to 0.017) | 0.018 (0.010 to 0.025) |
| aztreonam (s) | 2262 (259) | 0.990 | 0.004 (0.001 to 0.007) |  |  |  |  |
| cefazolin | 2520 (870) | 0.960 | 0.022 (0.012 to 0.031) | 0.260 (0.223 to 0.298) | 0.080 (0.060 to 0.099) | 0.006 (−0.023 to 0.035) |  |
| cefepime | 3532 (656) | 0.965 | 0.017 (0.009 to 0.024) | 0.187 (0.154 to 0.224) | 0.115 (0.096 to 0.135) | 0.112 (0.065 to 0.166) | 0.185 (0.152 to 0.220) |
| cefoxitin | 5027 (875) | 0.975 | 0.002 (−0.005 to 0.009) | 0.036 (0.026 to 0.047) | 0.040 (0.026 to 0.056) | 0.014 (0.005 to 0.024) | 0.075 (0.040 to 0.109) |
| ceftazidime | 3245 (893) | 0.979 | 0.007 (0.004 to 0.010) | 0.071 (0.057 to 0.085) | 0.028 (0.021 to 0.036) | 0.166 (0.129 to 0.202) | 0.222 (0.188 to 0.260) |
| ceftiofur (s) | 1300 (515) | 0.979 | 0.000 (−0.003 to 0.004) |  |  |  |  |
| ceftriaxone | 6917 (1475) | 0.987 | 0.003 (−0.001 to 0.006) | 0.011 (0.006 to 0.016) | 0.016 (0.010 to 0.021) | 0.002 (−0.001 to 0.005) | 0.010 (−0.035 to 0.049) |
| chloramphenicol | 4386 (353) | 0.990 | 0.002 (−0.001 to 0.007) | −0.003 (−0.009 to 0.002) | −0.002 (−0.009 to 0.002) | −0.004 (−0.012 to 0.000) | −0.001 (−0.004 to 0.001) |
| ciprofloxacin | 7529 (1275) | 0.975 | 0.006 (0.002 to 0.010) | 0.067 (0.053 to 0.082) | 0.027 (0.019 to 0.036) | 0.079 (0.054 to 0.106) | 0.068 (0.048 to 0.091) |
| ertapenem (s) | 3465 (65) | 0.818 | 0.078 (−0.012 to 0.195) |  |  |  |  |
| gentamicin | 7216 (1176) | 0.985 | 0.000 (−0.002 to 0.003) | 0.010 (0.005 to 0.014) | 0.001 (−0.001 to 0.003) | 0.010 (0.002 to 0.019) | 0.013 (0.000 to 0.029) |
| kanamycin (s) | 1183 (161) | 0.975 | 0.015 (0.001 to 0.032) |  |  |  |  |
| levofloxacin | 2941 (840) | 0.987 | 0.007 (0.004 to 0.010) | −0.033 (−0.097 to 0.010) | 0.005 (0.002 to 0.011) |  | 0.004 (−0.003 to 0.012) |
| meropenem (s) | 6139 (57) | 0.922 | −0.032 (−0.070 to −0.000) |  |  |  |  |
| nalidixic acid | 4522 (289) | 0.880 | 0.104 (0.066 to 0.149) | 0.110 (0.075 to 0.143) | 0.177 (0.134 to 0.224) | 0.175 (0.115 to 0.240) | 0.131 (0.063 to 0.201) |
| piperacillin-tazobactam (s) | 2817 (107) | 0.831 | 0.112 (0.046 to 0.184) |  |  |  |  |
| streptomycin | 4096 (1474) | 0.969 | 0.004 (0.001 to 0.007) | 0.015 (0.011 to 0.019) | 0.017 (0.012 to 0.022) | 0.024 (0.014 to 0.033) |  |
| sulfisoxazole | 4035 (1445) | 0.992 | 0.003 (0.001 to 0.007) | 0.014 (0.008 to 0.020) | 0.013 (0.007 to 0.019) | 0.008 (0.001 to 0.016) |  |
| tetracycline | 6940 (2909) | 0.991 | 0.001 (−0.001 to 0.004) | 0.007 (0.004 to 0.010) | 0.008 (0.005 to 0.013) | 0.001 (−0.001 to 0.004) | −0.001 (−0.005 to 0.001) |
| tobramycin (s) | 2850 (261) | 0.980 | 0.002 (−0.003 to 0.009) |  |  |  |  |
| trimethoprim-sulfamethoxazole | 7146 (1366) | 0.985 | 0.006 (0.003 to 0.009) | 0.012 (0.008 to 0.015) | 0.011 (0.008 to 0.015) | 0.015 (0.007 to 0.027) | 0.001 (−0.004 to 0.005) |

Across drugs (median, interquartile range, and the number of drugs whose
interval lies above or below zero):

| Holdout | Drugs | Median ΔAUROC | IQR | Above 0 | Below 0 | Median Δ balanced accuracy |
| --- | --- | --- | --- | --- | --- | --- |
| lineage | 24 | 0.005 | 0.002 to 0.016 | 14 | 1 | 0.013 |
| country | 16 | 0.025 | 0.010 to 0.070 | 14 | 0 | 0.039 |
| project | 16 | 0.019 | 0.009 to 0.040 | 14 | 0 | 0.024 |
| geography | 15 | 0.011 | 0.002 to 0.079 | 10 | 0 | 0.029 |
| temporal | 13 | 0.018 | 0.002 to 0.103 | 8 | 0 | 0.021 |

The random split scored higher than every holdout for most drugs, but by a
small margin for many of them. The lineage drop was small (median 0.005):
with a lineage component held out the model still ranked isolates almost as
well as with it in training. The larger drops came with country, project,
geography and temporal holdouts, and they were concentrated in a few drugs:
nalidixic acid under every holdout, cefepime, ceftazidime and cefazolin under
the group and temporal holdouts, and ciprofloxacin. Chloramphenicol,
tetracycline and gentamicin showed almost no drop. The one interval below
zero was meropenem under lineage (57 R).

Against `random_cluster` the drops are almost the same (lineage median
0.003, country 0.025, project 0.018, geography 0.009, temporal 0.015). Over
the whole common cohort, random minus random_cluster pooled AUROC had median
0.001 (range −0.006 to 0.055); dividing SNP clusters between training and test
added little except for the drugs with few R (ertapenem 0.055, amikacin 0.045,
kanamycin 0.027).

### 22.2 Leakage in the random split

In the random split, 90.9% of test isolates had their ST and 67.0% of
clustered test isolates had their SNP cluster in training; random_cluster
brought the cluster share to 0, and lineage brought both to 0
(`split_leakage_audit.tsv`). The country, project and geography holdouts still
had 3 to 9% of clustered test isolates sharing a cluster with training, and the
temporal holdout 22%; these were not removed (decision 0007 §5).

### 22.3 Heterogeneity

Fold-level results vary much more than the pooled drops suggest
(`fold_metrics.tsv`, `country_shift_summary.tsv`). The lowest held-out AUROCs
were nalidixic acid in Canada (0.64, 80 R) and in PRJNA966974 (0.65), cefazolin
in the USA (0.66, 772 R) and in PRJNA278886 (0.68), and cefoxitin in Zambia
(0.67, 15 R). Country, project, source and year move together in this cohort,
so these are differences between test populations, not effects of country.

### 22.4 Seen and unseen lineages

Within each strategy, predictions for test isolates whose ST was in training
were a little better than for those whose ST was not: the median difference in
balanced accuracy across drugs was 0.014 for random, 0.001 for country, −0.002
for project, 0.006 for geography and 0.016 for temporal, with wide ranges
(`unseen_lineage_summary.tsv`). Being in an unseen ST did not by itself mark
an isolate as hard to predict. Under the group holdouts few test isolates
share a cluster with training, so 20 of the 82 cluster-seen rows are
unsupported.

### 22.5 Determinant baseline against the model

On the endpoints with an interpretable determinant status
(`determinant_vs_ml.tsv`), the median balanced accuracy over drugs was:

| Strategy | Drugs | Determinant | Model | Model better | Determinant better |
| --- | --- | --- | --- | --- | --- |
| random | 24 | 0.953 | 0.948 | 6 | 5 |
| random_cluster | 24 | 0.953 | 0.944 | 4 | 8 |
| lineage | 24 | 0.953 | 0.941 | 3 | 12 |
| country | 16 | 0.957 | 0.911 | 0 | 13 |
| project | 16 | 0.957 | 0.933 | 1 | 12 |
| geography | 15 | 0.967 | 0.930 | 2 | 9 |
| temporal | 13 | 0.914 | 0.848 | 0 | 8 |

"Better" counts drugs whose paired interval for the difference excludes zero.
Under a random split the model and the curated determinant calls were close;
under every holdout the determinant calls did as well or better for most
drugs. The determinant column differs between strategies only because the
test populations differ.

### 22.6 R/determinant-negative endpoints

The model called R for 309 of 928 R/DET- endpoints under the random split, 276
of 928 under lineage, 76 of 708 under country, 102 of 702 under project, 21 of
298 under geography and 39 of 141 under temporal
(`unexplained_resistance_ml.tsv`). Under the group holdouts most of these calls
were for isolates whose ST was in training (country 66 of 76, project 93 of 102,
temporal 36 of 39). The model also called S for 2,002 of 3,167 S/DET+ endpoints
under random. Neither count is evidence of a new mechanism (section 24).

### 22.7 Model and compute

C was chosen as 0.03 in 168 tasks, 0.1 in 444 and 0.3 in 457. Models kept a
median of 224 non-zero coefficients (5 to 1,121) out of a median of 219,491
features per task. The run took 15.6 hours of task time on one worker (median
27 s per task, peak resident memory 7.0 GB); the report took 12 minutes.

## 23. Limitations

- One feature representation and one model. The drops describe this model;
  another model or representation could lose more or less.
- FracMinHash at scaled 200 sees about 14% of single-nucleotide changes, so
  point-mutation resistance (quinolones, part of the cephalosporin
  resistance) is only partly visible. The drugs with the largest drops include
  nalidixic acid and ciprofloxacin, and this may be part of the reason; the
  benchmark cannot separate it from population shift.
- C = 0.3, the top of the grid, was chosen in 43% of tasks. A wider grid could
  change absolute scores; it was fixed before any result and is not changed
  in this version.
- Country, project, source and year are confounded. Holdout drops are drops
  for these particular test populations; they do not isolate one factor.
- The cohort is the public NCBI Pathogen Detection AST collection, dominated by
  USA surveillance projects and food-animal sources; other settings may differ.
- Group holdouts keep training isolates that share a SNP cluster with the test
  set (3 to 22% of clustered test isolates), so they are less strict than a
  purge would be.
- Secondary drugs have few R (47 to 515 in the common cohort) and only the
  lineage holdout; their intervals are wide.
- Balanced class weights make the probabilities uncalibrated for the cohort
  prevalence; Brier score and calibration are not interpreted.
- One run. The model is deterministic given the seeds, so there is no
  run-to-run noise to report; the bootstrap intervals cover resampling of
  split units, not choice of seed or fold assignment.

## 24. Interpretation guardrails

- A drop is the difference between two evaluations of one model on the same
  isolates. It is not a property of the drug, and it says nothing about other
  models.
- A small or absent drop does not show that a model generalises to every new
  population; it shows that this held-out group was not much harder than a
  random test set.
- A correct model call on an R/DET- endpoint is not evidence of a new
  resistance mechanism. The model may use lineage-linked k-mers, and most such
  calls under group holdouts were for STs seen in training. No mechanism claim
  is made from these results.
- Country results are geographic or domain-shift results.
- Maximal-cohort and common-cohort results are not compared with each other.
- Results with `supported = no` are reported for completeness and not
  interpreted.

## 25. Reproducibility

Install the environment in `environment/generalisation-benchmark.txt`
(Python 3.12), then run, with the PR #6 population build, the PR #7 corrected
baseline and the PR #8 genome cache in `datasets/ecoli_research_v1/`:

    for step in sketch matrix freeze run report; do
      .venv/bin/python bin/generalisation_benchmark.py $step \
        --population datasets/ecoli_research_v1/population_20261015 \
        --baseline datasets/ecoli_research_v1/genotype_baseline_20261024 \
        --outdir datasets/ecoli_research_v1/generalisation_v1 \
        $( [ $step = sketch ] && echo --cache datasets/ecoli_research_v1/genome_cache )
    done

`sketch` and `run` resume from their checkpoints. `freeze` on the same inputs
gives the same `task_split_sha256` values as the committed
`benchmark_manifest.json`; `run` and `report` stop if any task differs. The
committed copy of the summary tables is in
`benchmark/cohorts/ecoli_research_v1/generalisation_v1/` with `SHA256SUMS`
and `run_record.json` (commands, timings, deviations). The prediction table
(`predictions.tsv.gz`, about 13 MB, 1,030,746 rows), the split lists, the task
checkpoints, the sketches and the feature matrix (about 1 GB) stay in
`datasets/ecoli_research_v1/generalisation_v1/` and are rebuilt by the commands
above.

## 26. Output files

| File | Content |
| --- | --- |
| `benchmark_metrics.tsv` | pooled metrics with intervals, and fold median, IQR and range, per cohort, drug and strategy |
| `fold_metrics.tsv` | one row per task: training size, chosen C, features, non-zero coefficients, metrics |
| `generalisation_drop.tsv` | Δ for every metric against random and random_cluster, pooled and per fold |
| `random_split_optimism.tsv` | one row per drug and metric, all drops side by side |
| `manuscript_table.tsv` | common cohort, pooled, AUROC and Δ against random |
| `determinant_vs_ml.tsv` | determinant calls against the model on the same endpoints |
| `unexplained_resistance_ml.tsv` | R/DET- and S/DET+ endpoints by ST seen in training |
| `unseen_lineage_summary.tsv` | metrics by ST and SNP cluster seen or unseen |
| `country_shift_summary.tsv`, `temporal_shift_summary.tsv` | composition of each held-out group |
| `split_summary.tsv`, `split_leakage_audit.tsv`, `not_feasible.tsv` | the frozen tasks and their audit |
| `cohort_eligibility.tsv`, `drug_plan_corrected.tsv` | cohorts and drug tiers |
| `benchmark_manifest.json`, `benchmark_summary.json`, `feature_manifest.json`, `genome_verification.json` | hashes and counts |
| `figure_data/` | data behind the six summary figures |
