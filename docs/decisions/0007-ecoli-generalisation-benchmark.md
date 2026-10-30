# 0007. Design of the E. coli generalisation benchmark (generalisation_v1)

- Status: accepted
- Date: 2026-10-24
- Applies to: `benchmark/cohorts/ecoli_research_v1/generalisation_v1.json`
  (config 1.0.0), read by `bin/generalisation_benchmark.py` 1.0.0
- Written before: any model was fitted on cohort data. The only model fits
  made before this record was frozen were timing runs, whose predictions
  were not scored.
- Numbering: written as 0005 and frozen under that name. Numbers 0005 and
  0006 went to the container and caller-wrapper records on main while the
  benchmark was running, so the file was renamed to 0007 after the run. The
  frozen config (whose SHA-256 is in the benchmark manifest) and one docstring
  in the builder still say 0005; they were left unchanged so the recorded
  hashes stay valid.

## Context

The benchmark asks how much a random train/test split overstates genomic AMR
prediction compared with test sets drawn from unseen lineages, countries,
submitting projects and later years. Decision 0001 already fixed how the
cohort is grouped (lineage components, SNP clusters as the split unit, a 2021
temporal cutoff, leave-one-project-out as the main external axis, USA against
the rest for geography) and which drugs are studied. Decision 0004 then
swapped R and S for PRJNA1322038 and asked for the 0001 drug tiers to be
looked at again before any model comparison.

## Decisions

### 1. Drug tiers on corrected labels

The decision 0001 tier rules are run unchanged
(`population_structure.evaluation_plan`) on the binary endpoints of
population-eligible isolates with the corrected labels from the PR #7 primary
build. The result is written to `drug_plan_corrected.tsv` with every run.

The 16 primary drugs and their temporal and geographic flags do not change.
Imipenem falls from secondary to excluded (46 R, below the 50 R floor);
amikacin (59 R), ertapenem (90 R) and meropenem (78 R) stay secondary. The
correction moved the carbapenem R counts down because most of the former
PRJNA1322038 "resistant" carbapenem records were susceptible. The benchmark
covers 16 primary and 8 secondary drugs.

### 2. Endpoints

Every binary R/S endpoint of a population-eligible isolate is used,
including the 2,826 whose PR #7 determinant mapping is uncertain; those have
no determinant comparator (`established_determinant_status = NA`). An endpoint
enters only if its genome matches the PR #8 manifest (SHA-256 and versioned
GCA, which must also equal the PR #6 and PR #7 assembly).

### 3. Feature representation

Canonical 31-mer presence, sampled by FracMinHash (sourmash 4.9.4, scaled 200,
seed 42): a k-mer is a feature if its 64-bit hash is below 2^64/200. The kept
set is fixed by the hash function, so no vocabulary is learned from genomes or
labels, and test genomes are projected onto the same coordinates.

For every task the model uses only features present in at least 5 and absent
from at least 5 training isolates. Before that, the matrix build drops
features present in fewer than 5 or absent from fewer than 5 of all genomes,
and merges columns with identical presence over all genomes. Neither step
changes any task: a feature dropped globally would fail every per-task filter,
and two columns identical in every genome are identical in every training set,
where an L1 penalty treats them as one.

Scaled 200 was chosen from a 200-genome sizing run before any model fit: about
25,000 sampled k-mers per genome, about 260 million non-zeros before pruning,
which fits a 16 GB machine. A 1 kb acquired gene contributes at least one
sampled k-mer with probability above 0.99; a single SNP is represented only
when one of its 31 overlapping k-mers is sampled (about 14%). Point-mutation
resistance (for example gyrA and parC for the fluoroquinolones) is therefore
only partly visible, and the model may substitute lineage-linked k-mers. This
is a limitation of the representation, stated with the results.

Alternatives: unitigs or full k-mer sets built from all genomes (rejected for
this version: memory on the available machine, and a vocabulary built from all
genomes would include test genomes); per-fold vocabularies built from training
genomes (equivalent here, since the hash-defined universe needs no
vocabulary); AMRFinderPlus calls as features (rejected: the comparison with the
determinant baseline would then compare a curated database with itself).

### 4. Model

One model: L1-penalised logistic regression (scikit-learn 1.7.2, liblinear),
balanced class weights, threshold 0.5 for every task. C is chosen from
{0.03, 0.1, 0.3} by 3-fold inner cross-validation on the training rows only,
grouped by split unit, by mean AUROC (ties to the smaller C). No second model
is fitted. The same feature rule, model, C policy and threshold are used for
every split strategy.

Balanced weights make the predicted probabilities reflect a 50% prevalence, so
the Brier score and calibration are reported but are not a main result.

### 5. Split strategies

The comparison has seven strategies. All of them run on the same
per-drug common cohort, and each also runs on its own maximal cohort.

| Strategy | Rule | Drugs |
| --- | --- | --- |
| `random` | stratified 5-fold over isolates; SNP clusters can be divided | primary and secondary |
| `random_cluster` | stratified 5-fold over split units (decision 0001 reference) | primary and secondary |
| `lineage` | 5-fold over lineage components (decision 0001); isolates without a component train only | primary and secondary |
| `country` | leave one country out, countries with at least 10 R and 10 S | primary |
| `project` | leave one project out, projects with at least 10 R and 10 S | primary with `project_holdout = yes` |
| `geography` | train on USA, test on all other countries | primary with `geographic_transfer = yes` |
| `temporal` | train before 2021, test 2021 on | primary with `temporal = yes` |

Countries and projects get a task only if they meet their test minimum. Every
task also needs at least 20 R and 20 S in training, and the geography and
temporal tasks need at least 20 R and 20 S in test. A task that fails is
written to `not_feasible.tsv` with its reason.

Fold assignment is deterministic: units in order of size, then of a SHA-256
key built from the seed, the drug and the unit, each placed where the R and S
load stays most even.

The primary comparator for a holdout is the `random` prediction for the
same test isolates, since every common-cohort isolate gets one out-of-fold
random prediction. The difference then comes only from what the model was
trained on. `random_cluster` is the second comparator, and comparing against
the pooled random metric over all isolates is reported separately with its own
label.

Departures from decision 0001, for review:

- Decision 0001 §2 keeps SNP clusters whole in every split and allows a fully
  isolate-level random split only as an illustration of leakage. This benchmark's
  question is about the conventional random split, so `random` is the
  comparator and `random_cluster` is reported next to it. Their difference
  measures how much cluster division adds on its own.
- For the country, project, geography and temporal holdouts, the group alone
  defines the split. Training isolates that share a SNP cluster with the test
  set are not removed, because that matches deployment (a clone seen
  in the USA before it reaches Canada). The overlap is reported per task
  (`split_leakage_audit.tsv`), and predictions are split by cluster seen or
  unseen.
- Decision 0001 does not make country a primary axis. Leave-one-country-out is
  added because the benchmark question asks for it, and it is reported as geographic/domain
  shift, since in this cohort country is largely confounded with project,
  source and year.

### 6. Metrics and uncertainty

AUROC, AUPRC, sensitivity, specificity, PPV, NPV, balanced accuracy, F1, MCC,
Brier score and calibration intercept and slope are computed for every fold
and for the pooled predictions of each drug and strategy. 95% intervals come
from 1,000 bootstrap replicates that resample split units, not single
isolates (seed fixed). Drops are bootstrapped in pairs, with the same resample
for the random and the holdout predictions. Rows with fewer than 10 R or 10 S
are kept and marked `supported = no`. Fold-level results are kept, with
median, IQR and range per drug and strategy. No significance tests across
drugs are made; the unit of replication for strategy comparisons is the drug
or fold, not the isolate.

### 7. Seeds

Split 20261024, inner CV 20261025, model 20261026, bootstrap 20261027.

## Freezing

`freeze` writes `benchmark_manifest.json` with the hash of this config, the
upstream hashes and a SHA-256 of every task's training and test isolate list.
`run` and `report` rebuild the tasks and stop if any differs. Any later
change of rule, drug list, feature or model is a new benchmark version.

## Consequences

- Results are reported for all eligible drugs and tasks, including undefined
  or unsupported ones.
- The determinant baseline does not train, so on the same isolates it scores
  the same under every strategy. Its comparison with the model is reported
  per strategy, on the endpoints with an interpretable determinant status.
- Model predictions on R/DET- endpoints are reported overall and by whether
  the isolate's ST was seen in training. Correct predictions there are not
  evidence of new mechanisms.
