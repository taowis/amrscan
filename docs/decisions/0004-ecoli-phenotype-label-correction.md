# 0004. Swap R and S for project PRJNA1322038

- Status: accepted
- Date: 2026-10-24
- Applies to: `benchmark/cohorts/ecoli_research_v1/phenotype_corrections_v1.tsv`
  (correction `PC0001`), read by `bin/genotype_baseline.py` 1.1.0 through the
  `phenotype_corrections` entry of `baseline_config.json` and
  `baseline_config_mapping_v2.json`
- Written after: the first full baseline (`genotype_baseline_20261018`) had
  been built and its discordant endpoints reviewed.

## Context

Reviewing the carbapenem and amikacin R_DET_NEG endpoints of the first build
showed that most came from one submitting project, PRJNA1322038 (137
isolates, all submitted by the University of Queensland; mostly Pakistan,
Zambia and Thailand; 107 collected in 2022; all tested on Vitek). Its R/S
labels run against its own MICs: ampicillin MIC >32 mg/L is labelled S in 118
records and MIC <2 mg/L R in 16, and meropenem MIC <0.25 mg/L is labelled R in
120 records.

`scripts/ecoli_research_v1/baseline_qc.py` compares, for each project and
antimicrobial with at least 3 R and 3 S MIC observations, the median MIC of R
and S records. The check uses the submitted MICs and labels only, not the
genotype calls.

- PRJNA1322038: the R median is below the S median for 26 of 26
  antimicrobials.
- The other 21 projects that report MICs: 0 inverted.
- In PRJNA1322038 the I labels sit between the two groups (for example
  meropenem I at 4 and 8 mg/L, R at 0.5 mg/L or below, S at 16 mg/L or
  above), which fits an R/S swap and not a general scrambling of records.
- After swapping R and S, no antimicrobial in the project is inverted, and 3 of
  its 3,281 R/S MIC observations remain on the wrong side of the other
  group's median (one each for amoxicillin-clavulanic acid, cefoxitin and
  ticarcillin-clavulanic acid; `qc/corrected_mic_out_of_order.tsv`).

## Decision

For every binary R/S endpoint of an isolate in PRJNA1322038, use S where PR #5
has R and R where PR #5 has S. I, other categories and conflicts are not
touched. The rule is applied to the whole project, not record by record, and
is stated in a versioned table with its hash in the baseline config.

The PR #5 files stay frozen and unchanged; the swap happens when the baseline
reads the endpoints. Every endpoint row keeps the PR #5 label
(`phenotype_pr5`) next to the one used (`phenotype`) and the correction ID
(`phenotype_correction`). The corrected build, `genotype_baseline_20261024`,
replaces `genotype_baseline_20261018` as the primary baseline. The earlier
build stays in the repository as the record of the uncorrected labels.

## Effect

3,281 endpoints change label, 2,041 of them metric-eligible. Over all study
drugs (primary variant, mapping v1), R_DET_NEG goes from 1,876 to 1,013 and
S_DET_POS from 4,204 to 3,344. The largest changes are in the secondary drugs
that the project tested heavily:

| Antimicrobial | Sensitivity before → after | Specificity before → after |
| --- | --- | --- |
| ertapenem | 0.158 → 0.506 | 0.996 → 1.000 |
| imipenem | 0.214 → 0.886 | 0.997 → 0.998 |
| meropenem | 0.306 → 0.857 | 0.996 → 0.998 |
| amikacin | 0.423 → 0.797 | 0.934 → 0.927 |
| piperacillin-tazobactam | 0.443 → 0.731 | 0.900 → 0.906 |
| ampicillin | 0.940 → 0.945 | 0.959 → 0.986 |
| gentamicin | 0.869 → 0.921 | 0.986 → 0.995 |

The `blaNDM-5` rules, flagged in the first audit only because of this
project's 14 susceptible carriers, are no longer flagged; the audit now flags
22 rules instead of 31.

## Why this rule and not another

The correction rests on the MICs and the order of R, I and S within the
project, not on agreement with the genotype. It was found while looking at
genotype-discordant endpoints, so the corrected baseline is not blind to
that review. It is still not tuned to the baseline: it is one rule, fixed
before the corrected build was run, and it would be the same whatever the
genotype calls were.

## Alternatives considered

- Keep the labels and report a post hoc table without the project (the
  previous state): rejected. The primary metrics for the carbapenems and
  amikacin then mostly describe a labelling error.
- Drop the project: rejected. Its MICs are usable, and dropping it removes
  all Pakistan (54) and Zambia (36) isolates and half of the Thai ones (25
  of 50) from the cohort.
- Re-derive every label from MIC and breakpoints: rejected for now. The
  project does not record breakpoint versions, and re-interpreting one
  project under different rules from the rest would make its labels
  inconsistent with the cohort. This would be a cohort-wide change for a
  later PR #5 revision.
- Edit the PR #5 build: rejected. It is frozen and its hash is checked by
  PR #6 and PR #7. A future PR #5 revision can take up the same correction.

## Consequences for decision 0001

Decision 0001 kept the carbapenems secondary partly because 62–76% of their
resistant isolates came from PRJNA1322038, and excluded polymyxin B because
its R was almost entirely from that project. Those counts used the
uncorrected labels. The study tiers are not changed here, but they should be
looked at again before the machine-learning comparison.

The error should also be reported to the submitter and to NCBI Pathogen
Detection, so that the source table can be fixed.
