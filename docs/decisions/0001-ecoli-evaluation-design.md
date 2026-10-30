# 0001. Evaluation design for the E. coli research cohort

- Status: accepted
- Date: 2026-10-15
- Applies to: `benchmark/cohorts/ecoli_research_v1/population_config.json`
  (version 1.1.0, section `evaluation_design` and `feasibility_rules.temporal.cutoff_year`),
  `bin/population_structure.py` 1.1.0 outputs `evaluation_groups.tsv` and
  `evaluation_plan.tsv`
- Evidence: population build `datasets/ecoli_research_v1/population_20261015`
  (PR #5 cohort manifest SHA-256 `21348c7d…7caf`); all counts below are
  population-eligible isolates (8,202) or their binary R/S endpoints unless
  stated otherwise

## Context

PR #6 described the population structure of the research cohort and screened
which holdout designs each antimicrobial could support, leaving six questions
for review: how to group lineages, whether SNP clusters may be split, which
drugs to evaluate, how to use country, whether to keep a source holdout, and
where to put the temporal cutoff.

The census shows that country, collection year, isolation source and project
are largely the same variable. The study that submitted an isolate decides
where and when it was collected, what it was isolated from, which drugs were
tested and how resistant the sample is:

- 1,515 of the 1,570 isolates collected in 2020 come from one project
  (PRJNA809394, Walter Reed Army Institute of Research).
- Before 2020 the cohort is dominated by PRJNA292663 (FDA, 2,648 of 5,231) and
  PRJNA966974 (Canadian Food Inspection Agency, 1,128). From 2021 on it comes
  mainly from projects that do not appear earlier: PRJNA885502 (338),
  PRJNA1271593 (330), PRJNA1322038 (107), plus PRJNA278886 (338, with 170
  earlier).
- The 2021-on isolates are resistance-enriched: ampicillin is 909 R to 229 S
  (80% R) against 3,621 R of 7,787 (46%) over the whole population.
- 1,128 of 1,292 Canadian isolates are from PRJNA966974, and 338 of 505 Dutch
  isolates from PRJNA885502.
- 2,076 of 3,373 human isolates are from PRJNA809394, and all 528 food isolates
  from PRJNA292663; 4,171 isolates have no source category.

These decisions were made before any genotype–phenotype comparison or model
was run, and do not use model performance.

## Decisions

### 1. Lineage grouping

Lineage holdouts group isolates by connected components of shared assigned ST
or shared NCBI SNP cluster (`lineage_component` in `evaluation_groups.tsv`).
Components are formed over population-eligible isolates only. Primary analysis
is grouped cross-validation over components. A secondary analysis splits the
large STs internally by SNP cluster to estimate within-lineage performance;
both are reported.

UNKNOWN isolates (no assigned ST) join the component of their SNP cluster.
Those without a cluster have no component (`lineage_cv_test_eligible =
false`): they may be used for training but are never test data. They are not
assigned to a nearest ST. A sensitivity analysis drops all UNKNOWN isolates.

Evidence: seven components join more than one ST, each a large ST absorbing a
rare single-locus variant through a shared cluster (ST131 + ST14207, ST10 +
ST16873, ST101 + ST8604, ST69 + ST10065 + ST13763, ST3580 + ST9584, ST115 +
ST16941, ST1594 + ST10007). No two common STs are chained. Of 103 UNKNOWN
population isolates, 30 have a cluster (19 in a component with an ST, 11 in
cluster-only components) and 73 have none.

Alternatives: ST alone (rejected: the seven cross-ST clusters would put
near-identical genomes on both sides); SNP cluster alone (rejected: 56% of
population isolates have no cluster, and clusters within one ST still share
lineage-level resistance); clonal complexes (not available from the 7-gene
call without an extra, versioned grouping step).

### 2. SNP clusters are never split

In every split, including random and cross-validation splits, all isolates of
an NCBI SNP cluster stay together (`split_unit`: the cluster, or the isolate
when it has no cluster). A cluster-grouped random split is the optimistic
reference. A fully isolate-level random split may be shown only as an
illustration of leakage, labelled as such.

Evidence: 941 clusters contain two or more cohort isolates (2,969 isolates);
the largest has 61, so grouping costs little. Unclustered isolates are not
assumed to be independent; NCBI clusters only within its distance threshold.

### 3. Antimicrobials

Tiers are computed from population endpoints by `evaluation_plan`
(`drug_tiers` in the config) and written to `evaluation_plan.tsv`:

- primary: at least 100 R and 100 S; one project contributes at most 50% of R;
  at least three projects with 10 R and 10 S; at least ten STs with 5 R and
  5 S;
- secondary: not primary, at least 50 R and 50 S;
- excluded: below secondary, or named in `drug_tiers.excluded` with a reason.

Result: 16 primary (ampicillin, amoxicillin-clavulanic acid, cefoxitin,
cefazolin, ceftriaxone, ceftazidime, cefepime, ciprofloxacin, levofloxacin,
nalidixic acid, gentamicin, trimethoprim-sulfamethoxazole, tetracycline,
chloramphenicol, streptomycin, sulfisoxazole); 9 secondary (meropenem,
ertapenem, imipenem, amikacin, tobramycin, piperacillin-tazobactam, aztreonam,
ceftiofur, kanamycin); the rest excluded.

Named exclusions: cefotaxime (118 S to 641 R, resistance-enriched sampling);
cefotaxime-clavulanic acid and ceftazidime-clavulanic acid (ESBL confirmation
tests, nearly all R from PRJNA885502); colistin and polymyxin B (unreliable
outside reference broth microdilution; polymyxin B R almost entirely from
PRJNA1322038); trimethoprim (74% of R from one project); ampicillin-sulbactam
(one project with both phenotypes).

The carbapenems are clinically the most important drugs here but stay
secondary: 62–76% of their resistant isolates come from PRJNA1322038, so
project and phenotype cannot be separated. Secondary drugs enter lineage
cross-validation only, are reported with confidence intervals and are not used
for model selection.

Ciprofloxacin, levofloxacin and nalidixic acid share gyrA/parC mechanisms, and
ceftriaxone, ceftazidime and cefepime are also strongly correlated; results
within these groups are not independent evidence.

### 4. Geography

Country is not a primary holdout axis and is not restricted to three
countries. The primary external-validity analysis is leave-one-project-out,
which covers country, source and time shifts together. Geographic transfer is
one analysis: train on USA isolates (the reference country), test on all other
countries pooled (`geography_group`), for primary drugs with at least 20 R and
20 S outside the USA. Per-country results are descriptive and only for
countries with at least 10 R and 10 S for the drug
(`non_reference_countries_json`). Reports state that country is confounded with
project and source.

### 5. Source

There is no source holdout. Source is used only to stratify results where it
is known.

Evidence: half the population has no source category, known sources follow
single projects, and only tetracycline passed the source screen in PR #6.

### 6. Temporal cutoff

The cutoff is fixed at 2021 (`feasibility_rules.temporal.cutoff_year`):
`period` is `earlier` before 2021 and `recent` from 2021 on (1,197 isolates,
15% of those with a year). The analysis is reported as a combined time and
study shift, not as pure temporal extrapolation.

A primary drug enters the temporal analysis only with at least 20 R and 20 S
in both periods. Drugs that fail this (streptomycin 3 R / 13 S recent,
sulfisoxazole 0 / 1, cefazolin 257 / 12) are left out of the temporal
analysis. They are not imputed, and related drugs are not substituted for
them.

Interpretive breakpoints changed over 2009–2025 (for example CLSI lowered the
Enterobacterales fluoroquinolone breakpoints in 2019). Where MICs are
available, the genotype–phenotype work should add a sensitivity analysis that
reinterprets all MICs with one breakpoint version.

Alternatives: the 80th-percentile year used in the PR #6 screen gives 2020,
which makes the test period 55% one project from one year and 35% of dated
isolates rather than about 20%.

## Consequences

- PR #7 and any later model use `evaluation_groups.tsv` and
  `evaluation_plan.tsv` as built, without redefining groups or tiers.
- Performance is reported by project as well as overall. Because many projects
  sequenced resistant isolates selectively, PPV and NPV do not transfer to other
  populations; sensitivity, specificity and balanced accuracy are the main
  measures.
- `feasibility_rules.decision` stays `HUMAN_REVIEW_REQUIRED`: the screen in
  `holdout_feasibility.tsv` is kept as a description, and this record is the
  decision.
- Changing any rule here needs a new decision record and a new config version.
