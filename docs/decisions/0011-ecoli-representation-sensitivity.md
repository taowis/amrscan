# 0011. Design of the E. coli representation sensitivity analysis (representation_sensitivity_v1)

- Status: accepted
- Date: 2026-10-29
- Applies to: `benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1.json`
  (config 1.0.0), read by `bin/representation_sensitivity.py` 1.0.0
- Written before: any denser sketch, feature count, capture status or model
  fit of this analysis. The only results read were the generalisation_v1 and
  error_decomposition_v1 outputs already reported in decisions 0007 and 0008.

## Context

The error decomposition (decision 0008) found that the scaled 200 FracMinHash
features almost never contain a k-mer at the main quinolone resistance
mutations (gyrA S83L in 0 of 2,298 calls), and that adding 62 known-mutation
indicators removed most of the nalidixic-acid holdout gap and about half of
the ciprofloxacin gap while a permuted control did not. Its limitations note
that a denser representation would see the mutations and that the gap it
leaves was not measured. A reviewer will ask whether denser k-mers alone would
have fixed the quinolone result, and whether the main conclusion (mechanism
novelty explains more of the loss than lineage novelty) depends on the feature
set.

This analysis answers that by changing one thing, the feature density, and
keeping the cohort, genomes, endpoints, splits, model, tuning, threshold,
metrics, determinant baseline and mutation panel of the frozen benchmark.

## Hypotheses

Stated before any result and open to being falsified:

1. A denser FracMinHash representation captures a larger share of known
   point-mutation sites than scaled 200.
2. More mutation capture improves quinolone holdout performance.
3. Explicit mutation indicators add robustness beyond denser k-mers.
4. Drugs whose resistance is carried by acquired genes that scaled 200 already
   captures (chloramphenicol, tetracycline, gentamicin) stay stable.
5. Denser features do not rescue the cefazolin country drop, which decision
   0008 attributed to mechanisms absent from training.
6. Cefepime and ceftazidime stay unresolved if density is not what limits
   them.
7. Levofloxacin stays a counter-example: poor mutation representation without
   a performance drop.

## Decisions

### 1. Representations

- R0 is the PR #14 representation, reused as it is: the generalisation_v1
  matrix and predictions are read, not rebuilt or refitted.
- R1 is the same FracMinHash family (sourmash 4.9.4, canonical 31-mers, seed
  42, presence, global 5/5 prune, identical columns merged) at one smaller
  scaled value. Nothing else changes.
- R2 and its permuted control are the error_decomposition_v1 runs, reused:
  they already use the frozen splits, model and panel.
- R3 is R1 plus the same 62 indicators, and R3_permuted is R1 plus the same
  permuted block, on the R2 task set. They show whether the indicators still
  help once the k-mers are denser.

### 2. Sketching

Each genome is verified against the PR #8 manifest and sketched once at scaled
20, the densest candidate. A FracMinHash sketch at a coarser scaled value is
exactly the subset of hashes under its threshold, so all candidates and R0 come
from the same sketch. Downsampling to 200 must give the PR #14 sketch of every
genome hash for hash; a difference stops the analysis.

### 3. Choosing R1 without looking at performance

Candidates are scaled 100, 50 and 20. A pilot over the full cohort measures
only resources and capture: observed and kept features, merged columns,
non-zeros, matrix size, build time and memory, the time and memory of model
fits on a few tasks (no prediction is scored), and the capture of called genes
and mutations. No metric is computed and no test phenotype is read.

R1 is the densest candidate for which every genome is sketched, the largest
task's projected peak memory fits one available machine with 25% left free,
and the projected compute for R1, R3 and R3_permuted fits a budget that the
project owner sets after seeing the pilot resource table. The available
machines are the 16 GB laptop used for PR #14 and PR #15, and NCI Gadi nodes.
The selection, the budget and the projections are committed in
`selection.json` before the first R1 fit. A cheaper candidate is not chosen
because it might do better, and no candidate other than R1 is fitted. Scaled
values below 20 are not tried.

### 4. Benchmark

All 1,069 frozen tasks are refitted for R1 (common and maximal cohorts; the
analysis tables use the common cohort, as in decision 0008). Each task is
rebuilt from the frozen inputs and must match the manifest's split hash, and
its train and test lists must equal those of the R0 task and of the R2 task
where one exists. Task statuses are completed, not_applicable, failed or
resource_failed, with a reason; no task is dropped.

### 5. Comparisons

Metrics, pooling, gap, gap rescue with a paired bootstrap, the 84-pair
Spearman correlations and the four-way comparison with the determinant rule
are those of decision 0008, recomputed for each representation on identical
endpoints. The mechanism-novelty and lineage-shift values do not depend on the
features and are read from error_decomposition_v1.

### 6. Reading the result

- R1 captures the mutations and rescues the quinolone holdouts: sparse
  sampling contributed to the original failure.
- R1 captures them but rescues only part, and R3 or R2 does better: denser
  k-mers help, explicit mutations are more reliable.
- R1 still misses them and R2 rescues: the density tested is not enough.
- R1 captures them and nothing improves: missing mutations did not explain the
  loss.

Cefazolin staying poor supports the mechanism-distribution reading; cefepime
and ceftazidime staying poor keeps them unresolved; the stable drugs staying
stable supports their use as controls. Results that fit none of these are
reported as they are. The comparisons are interventions on the representation
but not a causal mediation analysis.

## Consequences

- PR #14 and PR #15 results are not changed or overwritten; this analysis
  writes to its own directory.
- The denser matrix may not fit the laptop; R1 can then run on Gadi, with the
  project and paths of decision 0008 (Gadi runs).
- If R1 rescues nothing, that is the result; density is not tuned further.

## Alternatives considered

- Choosing the scaled value by benchmark performance: rejected, it would turn
  a confirmatory analysis into a search.
- A unitig or full k-mer representation: a different feature family, so a
  change would not be attributable to density alone.
- Refitting R2 instead of reusing it: the reused runs already satisfy the
  comparison and refitting cannot change them.

## Outcome (2026-10-30)

R1 is scaled 50, chosen by the rule from the resource pilot (commit a00cd62,
before any R1 fit). All 1,069 R1 tasks and all 262 R3 and R3_permuted tasks
completed on the frozen splits; the R0 platform check reproduced the PR #14
predictions exactly. Results are in
[the analysis guide](../ecoli-representation-sensitivity.md).

- Hypotheses 1 and 2 are supported: gyrA S83L capture went from 0 to 99.9%,
  and the nalidixic-acid holdout gaps fell by 55 to 97%, about as much as with
  the explicit indicators. Ciprofloxacin was rescued in part (21 to 42% of the
  gaps).
- Hypothesis 3 is supported for ciprofloxacin (R2 beats R1 under country,
  geography and project holdouts) and not for nalidixic acid. The indicators
  add at most 0.016 AUROC on top of the denser k-mers.
- Hypotheses 4 and 7 are supported: the stable drugs and levofloxacin did not
  change.
- Hypothesis 5 is not supported as stated: the denser k-mers alone rescued
  0.106 (0.075 to 0.136) of the cefazolin country gap of 0.260, and the
  catalogue mutations did not.
- Hypothesis 6 is supported in part: cefepime country improved (rescue
  0.073); cefepime and ceftazidime still have the largest gaps left.
- The mechanism-novelty correlation fell from 0.66 to 0.48 and stayed above
  the lineage-shift correlation (−0.14).

Deviations from this record: the fit pilot was resubmitted after a failed
first run; the pilot read the feature spaces from the built matrices; capture
takes the best copy of a determinant; R1 ran with memory-packed tasks on
256 GB nodes; the budget was raised from 1,300 to 1,600 SU without changing
the selection; and direct comparisons between R1, R2, R3 and R3_permuted were
added after the R0-referenced tables were read. Details are in the guide.
