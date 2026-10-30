# E. coli representation sensitivity (representation_sensitivity_v1)

This analysis refits the frozen generalisation benchmark with denser k-mer
features, and with known-mutation indicators on top of them, to see whether
the conclusions of the error decomposition depend on the feature set. It is a
confirmatory sensitivity analysis. It does not replace the PR #14 benchmark or
the PR #15 decomposition, whose results stay as they are. Design: [decision
0011](decisions/0011-ecoli-representation-sensitivity.md).

Representation density was selected without examining downstream holdout
predictive performance. All representation comparisons used frozen benchmark
splits and identical test isolates.

## 1. Why this analysis

The PR #14 features are FracMinHash k-mers sampled at 1 in 200 (scaled 200).
The error decomposition ([guide](ecoli-domain-shift-error-decomposition.md),
[decision 0008](decisions/0008-ecoli-error-decomposition.md)) found that these
features almost never contain a k-mer spanning the main quinolone resistance
mutations: gyrA S83L was captured in 0 of 2,298 calls. Adding 62 known-mutation
indicators to the same features removed most of the nalidixic-acid holdout gap
and about half of the ciprofloxacin gap, while a permuted control did not. Two
questions were left open. Would denser k-mers alone have fixed the quinolone
result? And does the main conclusion, that mechanism novelty explains more of
the domain-shift loss than lineage novelty, depend on the feature set?

## 2. What PR #15 found

- Drops under domain holdouts were drug-specific and large for a few drugs
  (nalidixic acid, cefazolin, cefepime, ceftazidime, ciprofloxacin, cefoxitin,
  piperacillin-tazobactam).
- Across the 84 drug by holdout pairs of the common cohort, the AUROC drop
  correlated with the fraction of test resistance carried by mechanisms unseen
  in training (Spearman ρ = 0.66) and not with lineage shift (ρ = −0.18).
- Nalidixic acid and ciprofloxacin were read as feature-representation
  failures, cefazolin (country) as a mechanism-distribution shift, cefepime and
  ceftazidime as unresolved, and chloramphenicol, gentamicin and tetracycline
  as stable controls. Levofloxacin was a counter-example: its mutations were
  just as poorly captured but it lost little performance.
- The curated determinant rule tied or beat the model under most holdouts;
  rule correct and model wrong grew from 1,745 endpoints under random splits
  to 3,268 under country holdouts.

## 3. Hypotheses

Stated in decision 0011 before any denser sketch was made:

1. Denser FracMinHash features capture a larger share of known point-mutation
   sites than scaled 200.
2. More mutation capture improves quinolone holdout performance.
3. Explicit mutation indicators add robustness beyond denser k-mers.
4. Chloramphenicol, tetracycline and gentamicin, whose acquired genes scaled
   200 already captures, stay stable.
5. Denser features do not rescue the cefazolin country drop.
6. Cefepime and ceftazidime stay unresolved if density is not what limits
   them.
7. Levofloxacin stays a counter-example.

## 4. Frozen inputs

| Input | SHA-256 |
| --- | --- |
| analysis config (`representation_sensitivity_v1.json`) | `f20203e3…1566` |
| generalisation config | `00fa371a…c65f` |
| benchmark manifest (split hashes of all 1,069 tasks) | `f0b51c76…f303` |
| PR #14 predictions | `63b1a2c0…0ce0` |
| PR #14 feature manifest | `ec85fc5c…b15f` |
| genome manifest (PR #8, 10,521 genomes) | `e85a45fc…f2d6` |
| error decomposition config | `e97bd57d…af75` |
| error decomposition `SHA256SUMS` | `f9bceaca…7b6b` |
| augmentation plan | `39e5a31f…0f72` |
| mutation panel (62 indicators) | `60aa963e…5f` |
| mechanism catalogue | `c41ae6fd…4eb34` |

`check` verifies these and every later stage stops if the analysis config
changes. The cohort, genomes, phenotype endpoints, splits, model (L1 logistic
regression, liblinear), C grid, inner cross-validation, seeds, feature filter,
threshold, metrics, determinant baseline and mutation panel are those of PR #14
and PR #15. Only the k-mer density changes.

All 10,521 genomes were verified against the PR #8 manifest (no hash failure,
none missing). Each was sketched once at scaled 20, and downsampling that
sketch to scaled 200 reproduced the PR #14 sketch hash for hash for every
genome.

Every refitted task was rebuilt from the frozen inputs and had to match the
manifest's split hash, and its train and test lists had to equal those of the
R0 task and, where one exists, the R2 task. All 1,069 R1 tasks and all 262 R3
and 262 R3_permuted tasks passed. The report also checks that each
representation predicts exactly the R0 test endpoints of every drug and
strategy; none differed.

## 5. Representations

| Name | Features | Tasks |
| --- | --- | --- |
| R0 | PR #14 FracMinHash, scaled 200, reused as it is | 1,069 (PR #14) |
| R1 | the same family at scaled 50 | 1,069 refitted |
| R2 | R0 plus the 62 indicators (PR #15 `mutation_augmented`), reused | 262 |
| R2_permuted | R0 plus the indicator rows shuffled across isolates, reused | 262 |
| R3 | R1 plus the 62 indicators | 262 refitted |
| R3_permuted | R1 plus the same permuted block | 262 refitted |

R1 keeps sourmash 4.9.4, canonical 31-mers, seed 42, presence coding, the
global 5/5 prune and the merging of identical columns. The indicators are 0/1
calls from AMRFinderPlus (database 2026-08-07.1) for catalogue point mutations
with a `predicts_resistance` or `component` rule; they come from external
knowledge, not from phenotypes. The permuted block keeps each column's
prevalence and breaks its link to the genome (seed 20261029). The analysis
tables use the common cohort, as PR #15 did (532 R0/R1 tasks, 132 drug and
strategy pairs; 262 tasks and 68 pairs for the R2 drug set).

## 6. Choosing the density

The candidates were scaled 100, 50 and 20. A pilot over the full cohort
measured resources and capture only; no prediction was scored and no test
phenotype was read.

| | R0 (200) | 100 | 50 | 20 |
| --- | --- | --- | --- | --- |
| kept features | 419,872 | 839,758 | 1,679,419 | 4,197,009 |
| non-zeros | 253.0 M | 501.5 M | 977.2 M | 2,283.8 M |
| peak memory of the largest fit | 11 GB | 21 GB | 40 GB | 103 GB |
| projected SU (R1, R3, R3_permuted) | | 312 | 1,239 | 11,197 |
| catalogue mutation calls captured | 1.5% | 21% | 53% | 70% |
| gyrA S83L captured | 0% | 55% | 99.9% | 99.9% |
| acquired gene calls captured | 99.25% | 99.82% | 99.91% | 99.91% |

The frozen rule takes the densest candidate whose largest task fits one
machine with 25% memory left and whose projected compute fits the budget the
project owner set after seeing this table. With 1,300 SU the rule gave scaled
50; the budget was later raised to 1,600 SU after a schedule simulation, which
does not change the choice (any budget from 1,239 to 11,197 SU gives scaled
50). `selection.json` was committed (a00cd62) before the first R1 fit. No other
density was fitted.

## 7. Capture

At scaled 50 the share of carriers whose mutation is spanned by a kept k-mer
rose for every quinolone:

| Antimicrobial | Mutation carrier endpoints | R0 | R1 |
| --- | --- | --- | --- |
| nalidixic acid | 615 | 1.1% | 85.0% |
| ciprofloxacin | 5,331 | 0.3% | 76.0% |
| levofloxacin | 4,147 | 0.2% | 74.8% |

By mutation, gyrA S83L went from 0 to 2,295 of 2,298 calls, gyrA D87N from 0
to 1,610 of 1,627, and parC S80I from 0 to 865 of 1,685. Cefazolin (53
endpoints with a catalogue mutation, 21% captured) and ceftazidime (98, 51%)
gained less. Acquired-gene capture was already complete or nearly so and
stayed that way: chloramphenicol 397/397 under both, gentamicin 1,139/1,139,
tetracycline 3,170 and 3,172 of 3,175, cefazolin 1,450/1,450, cefepime 1,351
of 1,352, ceftazidime 1,184/1,184.

## 8. Quinolones

Holdout AUROC (common cohort; gap = random minus holdout AUROC on the
holdout's test isolates):

| Antimicrobial | Holdout | R0 | R1 | R2 | R3 | R2_perm | R0 gap | R1 gap | R2 gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nalidixic acid | lineage | 0.776 | 0.969 | 0.968 | 0.973 | 0.777 | 0.104 | 0.003 | 0.006 |
| nalidixic acid | country | 0.765 | 0.954 | 0.920 | 0.963 | 0.765 | 0.110 | 0.017 | 0.053 |
| nalidixic acid | project | 0.698 | 0.911 | 0.915 | 0.926 | 0.699 | 0.177 | 0.062 | 0.058 |
| nalidixic acid | geography | 0.693 | 0.890 | 0.911 | 0.906 | 0.693 | 0.175 | 0.080 | 0.060 |
| nalidixic acid | temporal | 0.714 | 0.936 | 0.946 | 0.939 | 0.714 | 0.131 | 0.016 | 0.008 |
| ciprofloxacin | lineage | 0.969 | 0.985 | 0.984 | 0.987 | 0.969 | 0.006 | 0.001 | 0.001 |
| ciprofloxacin | country | 0.906 | 0.947 | 0.964 | 0.948 | 0.908 | 0.067 | 0.038 | 0.021 |
| ciprofloxacin | project | 0.947 | 0.965 | 0.972 | 0.970 | 0.949 | 0.027 | 0.021 | 0.014 |
| ciprofloxacin | geography | 0.892 | 0.930 | 0.962 | 0.932 | 0.882 | 0.079 | 0.056 | 0.023 |
| ciprofloxacin | temporal | 0.886 | 0.935 | 0.936 | 0.948 | 0.887 | 0.068 | 0.040 | 0.038 |

Gap rescue against R0 with 95% intervals from 1,000 paired bootstrap resamples
of split units (`generalisation_gap_rescue.tsv`):

| Antimicrobial | Holdout | R1 rescue | R2 rescue | R3 rescue |
| --- | --- | --- | --- | --- |
| nalidixic acid | lineage | 0.101 (0.059 to 0.153) | 0.098 (0.054 to 0.151) | 0.098 (0.056 to 0.152) |
| nalidixic acid | country | 0.093 (0.058 to 0.129) | 0.057 (0.018 to 0.096) | 0.096 (0.062 to 0.133) |
| nalidixic acid | project | 0.114 (0.069 to 0.163) | 0.119 (0.074 to 0.167) | 0.127 (0.081 to 0.172) |
| nalidixic acid | geography | 0.095 (0.025 to 0.178) | 0.115 (0.053 to 0.186) | 0.106 (0.040 to 0.184) |
| nalidixic acid | temporal | 0.115 (0.045 to 0.187) | 0.123 (0.051 to 0.197) | 0.107 (0.030 to 0.183) |
| ciprofloxacin | country | 0.028 (0.019 to 0.039) | 0.046 (0.035 to 0.059) | 0.029 (0.019 to 0.041) |
| ciprofloxacin | project | 0.006 (−0.001 to 0.012) | 0.013 (0.005 to 0.022) | 0.010 (0.003 to 0.018) |
| ciprofloxacin | geography | 0.023 (0.004 to 0.043) | 0.056 (0.023 to 0.087) | 0.026 (0.004 to 0.048) |
| ciprofloxacin | temporal | 0.029 (0.015 to 0.045) | 0.030 (0.016 to 0.046) | 0.041 (0.025 to 0.062) |

The permuted control R2_permuted rescued nothing for either drug (−0.010 to
0.015). R3_permuted, which has the R1 k-mers, rescued as much as R1 did, as
expected if the rescue comes from the k-mers.

For nalidixic acid, denser k-mers on their own removed 55 to 97% of every
holdout gap, about as much as the explicit indicators; holdout sensitivity
rose from 0.29 to 0.47 under R0 to 0.73 to 0.85 under R1 (0.74 to 0.83 under
R2). The random-split AUROC rose too (0.880 to 0.972), so the sparse
representation cost performance everywhere and more under holdouts. For
ciprofloxacin R1 removed 21 to 42% of the country, project, geography and
temporal gaps, against 44 to 71% for R2, and R1 did not raise the geography
holdout sensitivity (0.60 to 0.58; R2 0.64).

The direct comparisons between the new representations
(`representation_pairwise.tsv`, holdout AUROC difference with a paired
bootstrap interval) separate the two kinds of information:

| Antimicrobial | Holdout | R2 − R1 | R3 − R1 | R3 − R3_permuted | R3 − R2 |
| --- | --- | --- | --- | --- | --- |
| nalidixic acid | country | −0.034 (−0.053 to −0.015) | 0.009 (0.000 to 0.019) | 0.009 (−0.000 to 0.019) | 0.043 (0.025 to 0.061) |
| nalidixic acid | geography | 0.021 (0.004 to 0.040) | 0.016 (0.002 to 0.034) | 0.019 (0.003 to 0.039) | −0.005 (−0.017 to 0.008) |
| nalidixic acid | project | 0.004 (−0.016 to 0.026) | 0.016 (0.001 to 0.030) | 0.014 (0.001 to 0.029) | 0.012 (−0.002 to 0.025) |
| ciprofloxacin | country | 0.017 (0.012 to 0.022) | 0.001 (−0.001 to 0.003) | 0.003 (0.000 to 0.005) | −0.015 (−0.021 to −0.010) |
| ciprofloxacin | geography | 0.032 (0.009 to 0.053) | 0.002 (−0.006 to 0.009) | −0.002 (−0.008 to 0.004) | −0.030 (−0.051 to −0.011) |
| ciprofloxacin | project | 0.007 (0.002 to 0.012) | 0.005 (0.002 to 0.008) | 0.005 (0.001 to 0.008) | −0.002 (−0.004 to 0.001) |
| ciprofloxacin | temporal | 0.001 (−0.008 to 0.009) | 0.013 (0.008 to 0.020) | 0.013 (0.008 to 0.020) | 0.012 (0.003 to 0.024) |

For nalidixic acid the explicit indicators (R2) and the denser k-mers (R1) do
about equally well; R1 is better under the country holdout and R2 under
geography. For ciprofloxacin R2 is better than R1 under the country,
geography and project holdouts. Adding the indicators to the denser k-mers
(R3 against R1, and against its permuted control) adds a little in places,
up to 0.016 AUROC for nalidixic acid and 0.013 for the ciprofloxacin temporal
holdout, and nothing elsewhere. For ciprofloxacin R3 is worse than R2 under
the country and geography holdouts (−0.015 and −0.030). One possible reason is
that with four times as many k-mer columns beside the 62 indicators, the L1
model shares the signal between them and the k-mers that track the same
mutations; this was not examined. The explicit mutations remain the more reliable
representation for ciprofloxacin, and they help most when the k-mer block
beside them is small.

## 9. Levofloxacin

Levofloxacin's mutations went from 0.2% to 74.8% captured, yet nothing
changed: holdout AUROC 0.878 to 0.980 under R0 and 0.871 to 0.984 under R1,
gaps of −0.033 to 0.007 under R0 and 0.002 to 0.017 under R1, and every R1
rescue interval includes zero or is below 0.007. Its country holdout has 86
test isolates and 10 resistant, too few to read. Levofloxacin remains a
counter-example: poor mutation representation was not enough on its own to
cause a generalisation loss. The acquired genes (344 endpoints, captured
under both) and k-mers linked to the mutations may already have carried
enough information under R0.

## 10. Stable-drug controls

Chloramphenicol, tetracycline and gentamicin stayed where they were. Their
holdout gaps under R1 ranged from −0.006 to 0.014 (R0 −0.004 to 0.013), every
R1 rescue was within ±0.005, and AUROC changed by less than 0.004 everywhere
(`stable_drug_controls.tsv`). The controls hold.

## 11. Cefazolin

Hypothesis 5 is not supported as stated. R1 rescued a substantial part of the
cefazolin country drop: holdout AUROC 0.700 to 0.806, gap 0.260 to 0.154,
rescue 0.106 (0.075 to 0.136), MCC 0.43 to 0.66, sensitivity 0.46 to 0.63. R3
and R3_permuted give the same rescue (0.103 and 0.112), and R2 gives none
(−0.004), so it comes from the denser k-mers, not from the catalogue mutations
(only 53 cefazolin endpoints carry one). The other cefazolin holdouts did not
improve (project 0.003, lineage −0.007, geography −0.035 with an interval
including zero).

What the denser k-mers see is not known from this analysis. Gene-level capture
was already complete under R0, so a likely candidate is finer resolution
within genes (alleles of blaCTX-M, blaTEM or blaCMY, or the ampC promoter
region) or of neighbouring sequence; this needs a separate look and is listed
for review. The mechanism-distribution reading of PR #15 still stands in
part: the country gap left under R1 (0.154) is the largest of any drug and
holdout, and the training set of that holdout is small. But part of what
PR #15 attributed to unseen mechanisms is recovered by a denser
representation of the same genes.

## 12. Cefepime and ceftazidime

Both remain largely unresolved. R1 reduced the cefepime country gap from 0.187
to 0.114 (rescue 0.073, 0.058 to 0.089) and smaller amounts for project (0.026)
and temporal (0.028), but not geography (−0.015, interval including zero).
Cefepime has no catalogue mutation, so it has no R2 or R3. For ceftazidime the
rescues were small (0.003 to 0.029) and the gaps stayed large: geography 0.137,
temporal 0.208 (R0 0.166 and 0.222). Under the temporal holdout the R1
balanced accuracy fell (0.58 to 0.54) while AUROC rose slightly, a threshold
effect of a shift in resistance prevalence. The determinant baseline still
beats the model on ceftazidime geography and on cefepime country, project and
geography.

## 13. All drugs

Across all 132 common-cohort drug and strategy pairs, R1 raised AUROC with an
interval above zero in 50 and lowered it with an interval below zero in 4
(largest fall 0.007, cefazolin lineage). Increases above 0.02 were confined to
nalidixic acid, ciprofloxacin, cefazolin, cefepime, ceftazidime, cefoxitin,
amoxicillin-clavulanic acid, ertapenem and piperacillin-tazobactam.

| | R0 | R1 |
| --- | --- | --- |
| median holdout gap (84 pairs) | 0.013 | 0.012 |
| gaps above zero | 77 | 77 |
| gap intervals above zero | 60 | 53 |
| largest gap | 0.260 (cefazolin country) | 0.208 (ceftazidime temporal) |
| ρ(gap, unseen-mechanism fraction) | 0.66 | 0.48 |
| ρ(gap, lineage shift) | −0.18 | −0.14 |

On the 48 pairs of the R2 drug set the unseen-mechanism correlation was 0.68
under R0, 0.45 under R1, 0.58 under R2 and 0.54 under R3, and the lineage-shift
correlation stayed between −0.18 and −0.27. Mechanism novelty remains more
strongly associated with the domain-shift loss than lineage novelty under
every representation, but the association is weaker once the k-mers are
denser, probably in part because the cefazolin and nalidixic-acid drops,
which combined large gaps with high unseen fractions, shrink.

## 14. Determinant baseline

The determinant baseline (PR #7 rules, unchanged) was compared with each
representation on the same endpoints (`determinant_vs_representation.tsv`).
Under the 84 holdouts the model's balanced accuracy beat the rule in 16 pairs
under R0 and 16 under R1, with a median difference of −0.021 and −0.016. Under
random splits R1 reached parity (median +0.001, 13 of 24 drugs). Denser k-mers
narrow the gap to the rules under domain shift a little but do not close it.

Rule correct and model wrong, pooled over drugs (`rule_ml_discordance_pooled.tsv`):

| Split | R0 | R1 | R2 (10 drugs) | R0 on the same 10 drugs | R3 |
| --- | --- | --- | --- | --- | --- |
| random | 1,745 | 1,346 | 506 | 994 | 531 |
| country | 3,268 | 1,833 | 1,612 | 2,432 | 1,260 |

Denser features cut the country excess of rule-correct, model-wrong endpoints
by 44% (3,268 to 1,833); model correct, rule wrong fell slightly (1,342 to
1,077).

## 15. Cost

| | R0 | R1 |
| --- | --- | --- |
| kept features | 419,872 | 1,679,419 |
| non-zeros (density) | 253.0 M (8.8%) | 977.2 M (10.7%) |
| matrix on disk | 0.97 GB | 3.7 GB |
| matrix build | 28 min | 31 min |
| model fits, 532 common tasks on the same Gadi nodes | 16.9 task-hours | 65.2 task-hours (median 4.1 times per task) |
| peak memory of the largest fit | 11 GB | 40 GB |

R1 does not fit the 16 GB laptop used for PR #14 and PR #15 and was run on
Gadi (normalbw, 28 cores, 256 GB). The Gadi jobs used 1,296 SU in all: 79.5
for the failed first fit pilot, 62.0 for the fit pilot, 2.4 for the R1
matrix, 28.5 for the R0 platform check, 699.6 for R1 (338.3 and 361.3), 214.6
for R3 and 209.7 for R3_permuted. The 62 indicators of R2 cost almost nothing:
R2 tasks took a median of 26 s on the laptop against 559 s for R3 on Gadi.

For nalidixic acid and ciprofloxacin, R2 reaches the same or better holdout
performance as R1 at a small fraction of the cost. Where the information sits
in a few known variants, encoding them is cheaper and at least as reliable as
sampling the genome more densely; where it does not (cefazolin, cefepime), only
the denser k-mers helped.

## 16. Noise baseline

The R0 tasks of the common cohort (532) were refitted on Gadi with the PR #14
matrix (`R0_gadi`). Every task gave byte-identical predictions to the
workstation run: same split, isolates, C and features used, no class
disagreement, maximum probability difference 0 (`r0_platform_check.tsv`,
`scripts/ecoli_research_v1/compare_platform_runs.py`). The fits are
deterministic across the two platforms, so the run-to-run noise of a
representation is zero and the R1 against R0 differences are read against the
paired bootstrap intervals, which measure sampling uncertainty over split
units.

## 17. Limitations

- One denser setting was tested. The choice was made from resources, and
  denser settings were not tried afterwards; how far capture and rescue would
  continue below scaled 50 is not known.
- Scaled 50 still misses about a quarter of the ciprofloxacin and
  levofloxacin mutation calls (parC S80I is captured in about half), so R1 is
  not a full-genome representation.
- The comparisons change the representation and keep everything else, which
  makes them interventions on the features, but they are not a causal
  mediation analysis; that more capture comes with rescue is consistent with,
  not proof of, representation failure contributing to the loss.
- Which sequence gives the cefazolin and cefepime rescue was not examined.
- The R1-against-R2 and R3-against-R1 comparisons in section 8 were added
  after the R0-referenced tables had been read (see Deviations).
- The bootstrap intervals are per pair and are not corrected for the number
  of pairs.
- One model family, as in PR #14.

## 18. Reading the result

Against the decision matrix of decision 0011:

- Nalidixic acid is case A: R1 captures the mutations and rescues the holdout
  gaps about as fully as the explicit indicators. Sparse FracMinHash sampling
  made a substantial contribution to the original nalidixic-acid failure.
- Ciprofloxacin is case B: R1 captures most mutation calls and rescues part of
  the gaps, and the explicit indicators (R2) rescue more.
- Levofloxacin remains a counter-example (hypothesis 7 supported).
- The stable drugs remain stable (hypothesis 4 supported).
- Cefazolin does not fit case E: its country drop is partly rescued by the
  denser k-mers alone (hypothesis 5 not supported as stated).
- Cefepime and ceftazidime are case F in the main: partial rescue for the
  cefepime country holdout, large gaps left for both (hypothesis 6 supported
  in part).
- Hypotheses 1 and 2 are supported; hypothesis 3 is supported for
  ciprofloxacin and not for nalidixic acid.

The overall PR #15 conclusion survives with a change of weight. Domain-shift
losses remain drug-specific, lineage novelty still explains little, and
mechanism novelty is still the stronger correlate of the loss; but more of
the original loss than PR #15 could show was a matter of representation,
including part of the cefazolin country drop.

## 19. Manuscript implications

- The PR #14 benchmark stays the primary result. This analysis is reported as
  a confirmatory sensitivity analysis that answers the reviewer question
  "would denser k-mers fix the quinolone result?": for nalidixic acid largely
  yes, for ciprofloxacin in part, and explicit mutation features are cheaper.
- The cefazolin country drop should be described as partly representation
  and partly mechanism distribution, not as a mechanism shift alone.
- The mechanism-novelty correlation should be reported under both
  representations (0.66 and 0.48).
- No further primary experiment is needed for the manuscript; the cefazolin
  sequence question can be described as a limitation or examined as a
  secondary analysis.

## 20. Reproducibility

The commands, job records and SU are in
`benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1/run_record.json`.
The genome cache and the PR #14, PR #15, PR #6 and PR #7 outputs in
`datasets/ecoli_research_v1/` are needed. Sketching, the pilot and the report
ran on the workstation; the fit pilot, the R1 matrix and the R1, R3 and
R3_permuted refits ran on Gadi with the scripts in
`scripts/gadi/representation_sensitivity/` (submitted with `qsub -P <project>
-l storage=scratch/<project>`). The report took 32 minutes and 1.5 GB on the
laptop. The R1 matrix (3.7 GB) and the predictions are not committed.

## 21. Deviations

- Fit pilot, first submission: all four jobs measured the fits and then failed
  writing their JSON because `pilot/` did not exist on the cluster copy; 79.5
  SU were spent and the numbers lost. Fixed (the directory is created and
  the record is also printed to the job log) and the jobs resubmitted; the
  matrices from the first submission were kept, so the second ran fits only.
- Resource pilot: the feature spaces were read from the matrices the fit pilot
  had built (`pilot_scaled_*/matrix`) instead of the one-pass bucket scan,
  which needs one random read per sketch and bucket and would have taken about
  three days on the workstation's external disk. Both routes give the same
  spaces (tests); the R0 space was checked hash for hash against PR #14.
- Several copies of one determinant in an isolate: capture takes the best
  status of the copies, as the PR #15 error decomposition does. The PR #15
  `representation_audit.tsv` took the last copy instead; this changes one
  endpoint (trimethoprim-sulfamethoxazole, dfrA12, an isolate with 22 copies:
  210 vs 211 of 211 captured). Per-call R0 statuses match PR #15 exactly.
- Task memory records: `max_rss_bytes` in the R1, R3 and R3_permuted task
  checkpoints is the child process's `ru_maxrss`, which on Linux starts from
  the parent's resident size at fork (about 3.6 GB in the R1 jobs and 17 GB in
  the R3 jobs, whose parent builds the indicator matrix first). Those pages
  are shared copy-on-write, so small tasks read high; the values are upper
  bounds. Large tasks matched the fit-pilot estimate (median 0.9 of it).
- Budget: the owner set 1,300 SU after the resource pilot, then raised it to
  1,600 SU after a schedule simulation projected about 1,490 SU; the selection
  is the same for any budget from 1,239 to 11,197 SU
  (`budget_revision.json`). The runs used 1,154.7 SU.
- R1 ran on 256 GB normalbw nodes with tasks packed by estimated memory
  (`run --memory-gb`), not one task per core; predictions are byte-identical
  to the unpacked runner (tests).
- Pairwise comparisons: decision 0011 says R3 shows whether the indicators
  still help once the k-mers are denser, but its tables compare every
  representation with R0 only. The direct comparisons (R2 and R3 against R1,
  R3 against R3_permuted and R2; `representation_pairwise.tsv`, command
  `pairwise`) were added after the report had been read. They use the same
  paired bootstrap and change no other table.

## 22. Output files

Committed in `benchmark/cohorts/ecoli_research_v1/representation_sensitivity_v1/`
(with `SHA256SUMS`):

| File | Content |
| --- | --- |
| `representation_sensitivity_manifest.json`, `selection.json`, `budget.json`, `budget_revision.json` | upstream hashes; density selection and budget, committed before any R1 fit |
| `genome_verification.tsv.gz` | genome check and scaled 200 reproduction per genome |
| `pilot/` | resource pilot, capture by candidate, fit-pilot records |
| `representation_metrics.tsv`, `representation_metric_deltas.tsv` | metrics per representation and paired deltas against R0 |
| `generalisation_gap_by_representation.tsv`, `generalisation_gap_rescue.tsv` | gaps and rescue against R0 |
| `representation_pairwise.tsv` | direct comparisons between the new representations |
| `quinolone_representation_analysis.tsv`, `levofloxacin_counterexample.tsv`, `stable_drug_controls.tsv`, `cefazolin_control.tsv`, `unresolved_drug_controls.tsv` | drug views |
| `determinant_vs_representation.tsv`, `rule_ml_discordance.tsv`, `rule_ml_discordance_pooled.tsv` | determinant rule against each representation |
| `all_drug_robustness.tsv`, `correlation_robustness.tsv` | gap summary and correlations per representation |
| `manuscript_primary_table.tsv`, `manuscript_quinolone_table.tsv`, `figure_data/` | manuscript tables and plotting data |
| `task_status.tsv` | status of every task of every representation |
| `representation_sensitivity_summary.json` | summary |
| `R1_benchmark/` | the generalisation_v1 report tables for R1 (all 1,069 tasks) |
| `r0_platform_check.tsv` | R0 workstation against Gadi, per task |
| `run_record.json`, `logs/` | commands, PBS job records and SU |
