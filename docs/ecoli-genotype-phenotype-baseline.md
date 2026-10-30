# E. coli established-determinant genotype–phenotype baseline

## Purpose

This baseline asks how well established AMR determinants, detected by
AMRFinderPlus and interpreted with an explicit mapping, explain the R/S
phenotypes of the frozen research cohort. Any later model is compared against
it. Nothing here is a model; nothing is trained or tuned.

Three statements hold throughout:

- Absence of a mapped AMR determinant is not equivalent to biological
  susceptibility.
- R/determinant-negative isolates are candidates for unexplained resistance,
  not evidence of novel resistance mechanisms.
- The determinant-to-antimicrobial mapping is defined independently of
  downstream machine-learning performance.

## Inputs and linkage

| Input | Identity |
| --- | --- |
| PR #5 cohort (`ecoli_research_v1`) | `cohort_manifest.json` SHA-256 `21348c7d…7caf`; the builder checks the PR #5 `SHA256SUMS` and this hash |
| PR #6 population build (`ecoli_research_v1_population`) | `population_manifest.json` SHA-256 `ca4c0f65…3ba3`, code `397329b`; the builder checks its `SHA256SUMS`, this hash, and that it was built from the same PR #5 manifest |
| Genomes | PR #6 genome cache; each call's input hash must equal the PR #6 `genome_sha256` |
| AMRFinderPlus | 4.2.7, database 2026-08-07.1, fingerprint `sha256:748425a0…978b` (SHA-256 over the hashes of all 162 database files) |
| Mapping | `ecoli_amrfinder_drug_mapping_v1`, see decision 0002 |
| Settings | `benchmark/cohorts/ecoli_research_v1/baseline_config.json` 1.0.0 |

Endpoints are joined on the PR #5 isolate ID (versioned PDT accession), and the
PDT, BioSample and Assembly accessions must agree between the PR #5 endpoint and
the PR #6 row; a mismatch, a duplicate endpoint, or isolate lists that differ
stop the build. No name-based join is used. PR #5 endpoint definitions,
eligibility and conflict handling, and PR #6 lineage, geography, period and
source fields are used as built. The country, year and source on each PR #5
endpoint agree with PR #6 for all 162,629 endpoints.

## Calling

`bin/genotype_baseline.py call` runs the PR #4 wrapper (`bin/amrscan.py call`,
then `harmonize`) on each cached genome with the PR #4 options:
`--organism Escherichia`, `--ident_min -1`, `--coverage_min 0.5`,
`--translation_table 11`, `--report_all_equal`, `--print_node`, nucleotide
input, no `--plus`. Each isolate gets its own directory with the raw report,
harmonized report, log and provenance. A call is reused only if its
provenance names the isolate, reports success, matches the genome hash, the
database fingerprint, the caller version and the option set, the raw report
hash matches, and re-harmonizing the raw report reproduces the stored
harmonized file byte for byte. Anything else is moved to `.quarantine/` and
called again, so a failed or interrupted isolate never forces a restart and
validated calls are never overwritten. Calls are logged to `call_log.jsonl`.

The AMRFinderPlus environment was recreated from the PR #4 explicit package
list (`datasets/amrfinderplus-db/env-4.2.7-osx-arm64.explicit.txt`) after the
earlier copy under `/private/tmp` had been partly removed by the system; the
database files were rechecked against their recorded SHA-256 list (162 of 162
match). Re-running isolate `PDT000041778.1` reproduced the PR #4 raw report
byte for byte.

## Determinant evidence

Only `element_type == AMR` rows count: acquired genes, alleles, and point
mutations (`subtype POINT`). STRESS and VIRULENCE rows are kept in the
harmonized files but never make a drug positive. Hits by `INTERNAL_STOP`,
`PARTIALP` or `PARTIALX` turn a rule into uncertain; a full-length hit for the
same symbol takes precedence. AMR elements without any rule for a study drug
are listed in `unmapped_determinants.tsv` and do not change any status.

## Mapping

Detection (AMRFinderPlus) and antimicrobial interpretation (the mapping) are
separate steps. The mapping is a versioned file of explicit rules, one per
AMRFinderPlus symbol and study antimicrobial, built from pinned ResFinder and
PointFinder phenotype tables; the method, sources, inferences and flagged
rules are in [decision 0002](decisions/0002-ecoli-determinant-mapping-v1.md).
It was committed (`cfe27b0`) before any cohort call was read. No rule was
changed after looking at the data; rules whose carriers disagree with them are
flagged in `mapping_audit.tsv` for review.

The study antimicrobials are the 16 primary and 9 secondary drugs of
[decision 0001](decisions/0001-ecoli-evaluation-design.md). PR #5 did not
freeze a drug panel; decision 0001 did, from phenotype counts and project
structure only. The other 19 drugs have no rules and their endpoints carry
`mapping_unavailable`.

## Definitions

Genotype states: `determinant_positive` (a predicting rule applies, with any
required co-mutation present, or both parts of the trimethoprim-sulfamethoxazole
combination are present), `determinant_uncertain` (only uncertain, partial or
incomplete evidence), `determinant_negative` (a successful, provenance-checked
call with no determinant satisfying the mapping). Genome unavailable, genome QC
failure, missing or failed calls, provenance mismatch and harmonization failure
are separate states and are never counted as negative.

Phenotype states are PR #5 endpoints: R and S are binary; conflicting,
nonbinary and missing endpoints are excluded and stay in the table.

Concordance categories are `R_DET_POS`, `R_DET_NEG`, `S_DET_POS`, `S_DET_NEG`.
For metrics, phenotype R is the positive reference state and
determinant-positive the baseline positive call (TP = R_DET_POS, FN =
R_DET_NEG, FP = S_DET_POS, TN = S_DET_NEG).

## Metrics and intervals

For each study drug: N, R, S, TP, FN, FP, TN, sensitivity, specificity, PPV,
NPV, accuracy, balanced accuracy, F1 and MCC. Undefined values are `NA`.
Sensitivity, specificity, PPV and NPV have Wilson 95% intervals. Each isolate
contributes at most one endpoint per drug, but isolates are not independent:
clonal lineages and SNP clusters share determinants and phenotypes, so the
intervals are too narrow and are descriptive. No AUROC is given; the rule makes
a binary call and no score is invented for it. Metrics are per drug; no pooled
cross-drug accuracy is reported.

Metric-eligible endpoints are PR #5 binary R/S endpoints of population-eligible
isolates (PR #6) for study drugs with a positive or negative genotype state.
The primary variant leaves uncertain endpoints out; the prespecified variants
count them as positive or negative, and the two-mechanism variant applies a
stricter fluoroquinolone rule (decision 0002).

## Stratification

Counts and metrics by lineage (Achtman ST) and phylogroup, country and
reference/non-reference geography, period (before 2021 / from 2021, decision
0001) and calendar year, and source category, plus project and SNP cluster for
error analysis. Metrics are given for a stratum only with at least 10 R and
10 S endpoints (`baseline_config.json`, marked for review); smaller strata keep
their counts with `NA` metrics, and missing values form an `unknown` stratum.
This is description, not a holdout experiment, and differences between strata
are not causal: country, period, source and project are largely the same
variable in this cohort (decision 0001).

## Discordance tables

`unexplained_resistance.tsv` holds every metric-eligible R endpoint without a
mapped determinant, with lineage, cluster, place, time, source, QC, caller,
database and mapping provenance, and the other AMR determinants the isolate
carries (which, under the mapping, do not explain the target drug). Before any
such isolate could support a novel-mechanism claim one would need to examine
mapping completeness, known but unmapped mutations, regulatory and promoter
variation, gene copy number, structural variation, assembly errors, AST quality
and breakpoint context, rare determinants, database coverage and population
structure.

`susceptible_with_determinant.tsv` holds every metric-eligible S endpoint with
a mapped determinant. Possible reasons include an inactive or disrupted gene,
allele-specific substrate range, low expression, breakpoint effects, phenotype
error, an over-broad mapping rule, and assembly or caller artefacts.

## Real-cohort execution

All 10,521 PR #5 isolates had a PR #6 genome (10,500 QC pass, 21 warning, none
failed) and a successful, provenance-checked AMRFinderPlus call; no call
failed, mismatched provenance or failed harmonization. The calls ran on one
workstation (Apple M4, 10 cores, 16 GB) with one thread per call and up to
eight calls at a time, in about 11.4 hours of active run time across two
sessions; 22 incomplete outputs from interrupted runs were quarantined and
called again. The run record is
`benchmark/cohorts/ecoli_research_v1/genotype_baseline_20261018/call_run.json`.

There are two builds from the same calls. `genotype_baseline_20261018` used
the PR #5 labels as they are. `genotype_baseline_20261024` swaps R and S for
project PRJNA1322038 ([decision 0004](decisions/0004-ecoli-phenotype-label-correction.md))
and is the primary baseline. The build reads 42,000 call files; it took between 17 seconds and 9 minutes
depending on how much of the external drive was in the file cache (13 to 26
seconds of CPU), under 1 GB of memory. The corrected build's record is
`genotype_baseline_20261024/build_run.json`. A second build of each
into a new directory gave an identical `SHA256SUMS`.

| Count | Value |
| --- | ---: |
| Cohort isolates / unique assemblies | 10,521 / 10,521 |
| PR #5 phenotype endpoints | 162,629 |
| Binary R/S endpoints | 121,689 |
| Endpoints of study antimicrobials with a binary phenotype | 113,143 |
| Determinant-positive / negative / uncertain (study endpoints) | 25,306 / 85,011 / 2,826 |
| Genotype-unavailable endpoints | 0 |
| Metric-eligible endpoints (primary) | 110,317 |
| Endpoints with a corrected label (all / metric-eligible) | 3,281 / 2,041 |
| R_DET_POS / R_DET_NEG / S_DET_POS / S_DET_NEG | 21,962 / 1,013 / 3,344 / 83,998 |

Of the 10,521 isolates, 8,202 are population-eligible in PR #6. The other
2,319 contribute no binary endpoint, so the population-ineligible exclusion is
empty for every drug. `baseline_denominators.tsv` reconciles every drug from
the PR #5 total to the metric-eligible set: non-binary or conflicting 40,940,
not a study drug 8,546, mapping uncertain 2,826, genotype unavailable 0. The
label swap changes no denominator. Metric-eligible endpoints are all single
binary observations or repeated consistent ones (40); no conflicting endpoint
reaches a denominator.

## Phenotype labels in project PRJNA1322038

One submitting project, PRJNA1322038 (137 isolates, mostly from Pakistan,
Zambia and Thailand, 107 collected in 2022, all tested on Vitek), has R/S
labels that run opposite to its own MIC values. Ampicillin MIC >32 mg/L is
labelled susceptible in 118 records and MIC <2 mg/L resistant in 16;
meropenem MIC <0.25 mg/L is labelled resistant in 120 records. The labels
are inverted relative to the MICs for all 26 antimicrobials tested in this
project, and for none in the other 21 projects that report MICs
(`qc/phenotype_label_check.tsv`, made by
`scripts/ecoli_research_v1/baseline_qc.py`). Projects that report only disk
zones or categories cannot be checked this way.

Decision 0004 swaps R and S for every binary endpoint of the project, as one
rule in `benchmark/cohorts/ecoli_research_v1/phenotype_corrections_v1.tsv`.
The PR #5 files are not changed; each endpoint row keeps `phenotype_pr5` and
`phenotype_correction` next to the label used. After the swap no antimicrobial
in the project is inverted (`qc/phenotype_label_check_corrected.tsv`), and 3
of its 3,281 R/S MIC observations still sit on the wrong side of the other
group's median (`qc/corrected_mic_out_of_order.tsv`).

Before the correction this project supplied 118 of 160 ertapenem, 120 of 125
imipenem and 119 of 129 meropenem R_DET_NEG endpoints, and 97 of 109 amikacin
ones. After it, 2 ertapenem, 1 meropenem and no imipenem or amikacin R_DET_NEG
endpoints come from the project. The uncorrected build and its QC tables stay
in `genotype_baseline_20261018/` for the record. Sensitivity before → after
for the most affected drugs: ertapenem 0.158 → 0.506, imipenem 0.214 → 0.886,
meropenem 0.306 → 0.857, amikacin 0.423 → 0.797, piperacillin-tazobactam
0.443 → 0.731, gentamicin 0.869 → 0.921. The correction was found while
reviewing discordant endpoints, but the rule depends only on the MICs.

## Per-antimicrobial baseline (primary variant, mapping v1, corrected labels)

Phenotype R is the positive reference state. FN is R_DET_NEG and FP is
S_DET_POS. Intervals are Wilson 95% and too narrow for clustered isolates.

| Antimicrobial | Tier | N | R | S | TP | FN | FP | TN | Sensitivity (95% CI) | Specificity (95% CI) | PPV | NPV | Bal. acc. | F1 | MCC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| amoxicillin-clavulanic acid | primary | 6433 | 911 | 5522 | 847 | 64 | 42 | 5480 | 0.930 (0.911–0.945) | 0.992 (0.990–0.994) | 0.953 | 0.989 | 0.961 | 0.941 | 0.932 |
| ampicillin | primary | 7773 | 3714 | 4059 | 3511 | 203 | 58 | 4001 | 0.945 (0.938–0.952) | 0.986 (0.982–0.989) | 0.984 | 0.952 | 0.966 | 0.964 | 0.933 |
| cefazolin | primary | 2593 | 958 | 1635 | 913 | 45 | 438 | 1197 | 0.953 (0.938–0.965) | 0.732 (0.710–0.753) | 0.676 | 0.964 | 0.843 | 0.791 | 0.662 |
| cefepime | primary | 3588 | 711 | 2877 | 690 | 21 | 483 | 2394 | 0.971 (0.955–0.981) | 0.832 (0.818–0.845) | 0.588 | 0.991 | 0.901 | 0.733 | 0.682 |
| cefoxitin | primary | 5097 | 967 | 4130 | 854 | 113 | 24 | 4106 | 0.883 (0.861–0.902) | 0.994 (0.991–0.996) | 0.973 | 0.973 | 0.939 | 0.926 | 0.911 |
| ceftazidime | primary | 3339 | 967 | 2372 | 956 | 11 | 371 | 2001 | 0.989 (0.980–0.994) | 0.844 (0.828–0.858) | 0.720 | 0.995 | 0.916 | 0.834 | 0.771 |
| ceftriaxone | primary | 6918 | 1535 | 5383 | 1505 | 30 | 126 | 5257 | 0.981 (0.972–0.986) | 0.977 (0.972–0.980) | 0.923 | 0.994 | 0.979 | 0.951 | 0.937 |
| chloramphenicol | primary | 4456 | 369 | 4087 | 362 | 7 | 16 | 4071 | 0.981 (0.961–0.991) | 0.996 (0.994–0.998) | 0.958 | 0.998 | 0.989 | 0.969 | 0.967 |
| ciprofloxacin | primary | 7424 | 1372 | 6052 | 1339 | 33 | 385 | 5667 | 0.976 (0.966–0.983) | 0.936 (0.930–0.942) | 0.777 | 0.994 | 0.956 | 0.865 | 0.839 |
| gentamicin | primary | 7358 | 1226 | 6132 | 1129 | 97 | 30 | 6102 | 0.921 (0.904–0.935) | 0.995 (0.993–0.997) | 0.974 | 0.984 | 0.958 | 0.947 | 0.937 |
| levofloxacin | primary | 2893 | 902 | 1991 | 890 | 12 | 305 | 1686 | 0.987 (0.977–0.992) | 0.847 (0.830–0.862) | 0.745 | 0.993 | 0.917 | 0.849 | 0.784 |
| nalidixic acid | primary | 4472 | 331 | 4141 | 287 | 44 | 20 | 4121 | 0.867 (0.826–0.899) | 0.995 (0.993–0.997) | 0.935 | 0.989 | 0.931 | 0.900 | 0.893 |
| streptomycin | primary | 4127 | 1486 | 2641 | 1451 | 35 | 304 | 2337 | 0.976 (0.967–0.983) | 0.885 (0.872–0.896) | 0.827 | 0.985 | 0.931 | 0.895 | 0.836 |
| sulfisoxazole | primary | 4066 | 1458 | 2608 | 1426 | 32 | 16 | 2592 | 0.978 (0.969–0.984) | 0.994 (0.990–0.996) | 0.989 | 0.988 | 0.986 | 0.983 | 0.974 |
| tetracycline | primary | 7078 | 3005 | 4073 | 2952 | 53 | 45 | 4028 | 0.982 (0.977–0.987) | 0.989 (0.985–0.992) | 0.985 | 0.987 | 0.986 | 0.984 | 0.972 |
| trimethoprim-sulfamethoxazole | primary | 5853 | 1372 | 4481 | 1343 | 29 | 50 | 4431 | 0.979 (0.970–0.985) | 0.989 (0.985–0.992) | 0.964 | 0.994 | 0.984 | 0.971 | 0.963 |
| amikacin | secondary | 3347 | 59 | 3288 | 47 | 12 | 240 | 3048 | 0.797 (0.677–0.880) | 0.927 (0.918–0.935) | 0.164 | 0.996 | 0.862 | 0.272 | 0.340 |
| aztreonam | secondary | 2284 | 306 | 1978 | 273 | 33 | 47 | 1931 | 0.892 (0.852–0.922) | 0.976 (0.969–0.982) | 0.853 | 0.983 | 0.934 | 0.872 | 0.852 |
| ceftiofur | secondary | 1297 | 517 | 780 | 502 | 15 | 29 | 751 | 0.971 (0.953–0.982) | 0.963 (0.947–0.974) | 0.945 | 0.980 | 0.967 | 0.958 | 0.930 |
| ertapenem | secondary | 3506 | 89 | 3417 | 45 | 44 | 0 | 3417 | 0.506 (0.404–0.607) | 1.000 (0.999–1.000) | 1.000 | 0.987 | 0.753 | 0.672 | 0.707 |
| imipenem | secondary | 3250 | 44 | 3206 | 39 | 5 | 5 | 3201 | 0.886 (0.760–0.951) | 0.998 (0.996–0.999) | 0.886 | 0.998 | 0.942 | 0.886 | 0.885 |
| kanamycin | secondary | 1232 | 178 | 1054 | 167 | 11 | 5 | 1049 | 0.938 (0.893–0.965) | 0.995 (0.989–0.998) | 0.971 | 0.990 | 0.967 | 0.954 | 0.947 |
| meropenem | secondary | 6196 | 77 | 6119 | 66 | 11 | 15 | 6104 | 0.857 (0.762–0.918) | 0.998 (0.996–0.999) | 0.815 | 0.998 | 0.927 | 0.835 | 0.834 |
| piperacillin-tazobactam | secondary | 2814 | 130 | 2684 | 95 | 35 | 251 | 2433 | 0.731 (0.649–0.800) | 0.906 (0.895–0.917) | 0.275 | 0.986 | 0.819 | 0.399 | 0.407 |
| tobramycin | secondary | 2923 | 291 | 2632 | 273 | 18 | 39 | 2593 | 0.938 (0.904–0.961) | 0.985 (0.980–0.989) | 0.875 | 0.993 | 0.962 | 0.905 | 0.895 |

Sensitivity is at least 0.93 for amoxicillin-clavulanic acid, ampicillin, all
cephalosporins except cefoxitin, chloramphenicol, the fluoroquinolones,
streptomycin, sulfisoxazole, tetracycline, kanamycin, tobramycin and
trimethoprim-sulfamethoxazole. Specificity is lowest for cefazolin (0.73),
cefepime (0.83), ceftazidime (0.84), levofloxacin (0.85) and streptomycin
(0.89). The weakest drugs are ertapenem (sensitivity 0.51, 89 R),
piperacillin-tazobactam (0.73, specificity 0.91) and amikacin (0.80, PPV
0.16); cefoxitin (0.88) and nalidixic acid (0.87) are the lowest primary
drugs.

The prespecified variants move the numbers where they should.
Trimethoprim-sulfamethoxazole has 1,451 uncertain endpoints (one component
only); counting them as positive drops specificity from 0.989 to 0.757, and
counting them as negative drops sensitivity from 0.979 to 0.926. The
two-mechanism fluoroquinolone variant raises ciprofloxacin specificity from
0.936 to 0.978 and levofloxacin from 0.847 to 0.937, at almost no cost in
sensitivity, with 722 and 340 single-mechanism endpoints set aside as
uncertain.

| Antimicrobial | uncertain_as_positive sens/spec | uncertain_as_negative sens/spec | fq_two_mechanisms sens/spec (N) |
| --- | --- | --- | --- |
| amoxicillin-clavulanic acid | 0.931 / 0.971 | 0.916 / 0.993 |  |
| cefepime | 0.971 / 0.818 | 0.957 / 0.835 |  |
| cefoxitin | 0.885 / 0.967 | 0.866 / 0.994 |  |
| ceftazidime | 0.989 / 0.830 | 0.977 / 0.846 |  |
| ceftriaxone | 0.981 / 0.959 | 0.976 / 0.977 |  |
| ciprofloxacin | 0.976 / 0.895 | 0.972 / 0.939 | 0.972 / 0.978 (6985) |
| levofloxacin | 0.987 / 0.793 | 0.987 / 0.857 | 0.987 / 0.937 (2689) |
| nalidixic acid | 0.868 / 0.964 | 0.865 / 0.995 |  |
| trimethoprim-sulfamethoxazole | 0.980 / 0.757 | 0.926 / 0.992 |  |
| piperacillin-tazobactam | 0.739 / 0.884 | 0.709 / 0.909 |  |

## Discordance and error analysis

`qc/discordance_concentration.tsv` gives, for each drug and concordance
category, the largest lineage, country, project, source, period and (for
determinant-positive endpoints) determinant set.

- Determinant-positive but susceptible, where one determinant dominates:
  `aac(6')-Ib-cr5` for amikacin (97% of 240 S_DET_POS, mostly USA projects
  PRJNA278886 and PRJNA809394); `blaTEM-1` for cefazolin (91% of 438, almost
  all USA human isolates from PRJNA809394); CTX-M-27 and CTX-M-15 for cefepime
  and ceftazidime (mainly PRJNA278886); `blaOXA-1` (62% of 251) and `blaCMY-2`
  for piperacillin-tazobactam; aadA genes for streptomycin (PRJNA966974,
  Canada); `gyrA_S83L` alone or with `parE_I529L` for the fluoroquinolones.
  For cefepime, ceftazidime and piperacillin-tazobactam the pattern is the
  expected gap between enzyme-level resistance in the sources and clinical
  breakpoints; for cefazolin it concerns a rule inferred from cephalothin.
- Resistant without a mapped determinant: ampicillin (203, mostly older USA
  isolates from PRJNA292663), cefoxitin (113), gentamicin (97). For
  gentamicin, 22 animal isolates from PRJNA1142533 (USA) carry `rmtE1`, which
  v1 leaves unmapped (see below). Of the 44 ertapenem R_DET_NEG endpoints, the
  largest group (19) is from the Netherlands (PRJNA885502).
- Lineage: ciprofloxacin specificity in ST131 is 0.34 (130 S), and 70 of its
  86 S_DET_POS carry `gyrA_S83L;parE_I529L`; cefepime
  specificity is 0.57 in ST131 and 0.56 in ST38.

Many error groups sit mostly in one country and one project: ampicillin
R_DET_NEG (USA 90%, PRJNA292663 69%), cefazolin S_DET_POS (USA 96%,
PRJNA809394 92%), sulfisoxazole R_DET_NEG (Canada 91%). Project, country and
period overlap heavily in this cohort (decision 0001), so these groups cannot
be attributed to any one of them.

## Stratified results

Strata need at least 10 R and 10 S endpoints for metrics; 464 of 16,527
lineage and phylogroup rows, 146 of 445 country rows, 59 of 72 period rows and
58 of 87 source rows reach that. The others keep counts with `NA` metrics.

- Lineage: across reported lineages, balanced accuracy for tetracycline stays
  above 0.91; for ciprofloxacin it ranges from 0.67 (ST131) to 0.997, and for
  cefepime from 0.75 (ST648).
- Geography: reference and non-reference countries differ by at most about
  six points in sensitivity or specificity for ampicillin, ceftriaxone and
  cefepime; gentamicin sensitivity is 0.90 in non-reference and 0.93 in
  reference countries. Before the correction Pakistan had balanced accuracy
  near zero for several drugs; it is now 0.87 to 1.00.
- Time: recent (from 2021) specificity is low for ceftriaxone (0.67, 55 S) and
  cefepime (0.53, 622 S), driven by ESBL carriers in PRJNA278886.
- Source: gentamicin sensitivity in animal isolates is 0.52 (56 R), mostly the
  `rmtE1` isolates; human isolates have the highest balanced accuracy for most
  drugs.

These are descriptive. Strata overlap, and differences between them are not
causal.

## Mapping audit and exploratory mapping v2

`mapping_audit.tsv` flags 22 rules with at least 20 carriers whose R fraction
is below 0.5, for review only (31 before the label correction; the `blaNDM-5`
flags were caused by PRJNA1322038 and are gone). They include
`aac(6')-Ib-cr5` for amikacin (24 R of 257 carriers), `blaTEM-1` and
`blaSHV-1` for cefazolin, CTX-M-14 and CTX-M-27 for cefepime, ceftazidime and
aztreonam, `blaOXA-1` and `blaCMY-2` for piperacillin-tazobactam,
`ampC_T-32A` and `ampC_C-42T` for the cephalosporins, `parC_S57T`,
`parE_I355T` and `qnrB19` for the quinolones, and `sul2` and `dfrA16` as
components of trimethoprim-sulfamethoxazole. No rule was changed because of
the audit.

The review of R_DET_NEG endpoints did find one mapping error: alleles such as
`rmtE1`, whose AMRFinderPlus parent node is a ResFinder gene, got no rule in
v1. [Decision 0003](decisions/0003-ecoli-determinant-mapping-v2.md) defines
mapping v2, which fixes only this, and the corrected baseline was rebuilt with
it (`genotype_baseline_20261024/mapping_v2/`). Gentamicin goes from 1,129 TP / 97 FN to 1,151 / 75 (sensitivity
0.921 to 0.939) and kanamycin gains one FP; no other drug changes. v2 is exploratory and
the primary baseline stays on v1.

## Known caller limitation: blaEC control difference

The PR #4 upstream software control matched 25 of 26 rows; the remaining
difference is the closest reference reported for `blaEC` (EC-15 expected,
EC-28 observed, same symbol and class). It is still unresolved. It does not
affect this baseline: `blaEC` and its alleles are `does_not_predict` for every
beta-lactam in mapping v1, and the reported symbol, not the closest reference,
is what the mapping reads.

## Limitations

- Assemblies, not reads. Copy number, low-frequency variants and some
  structural changes cannot be seen, and assembly gaps can hide a gene.
- The mapping follows ResFinder and PointFinder, which describe departure
  from the wild type, while PR #5 phenotypes are clinical categories under
  mixed CLSI and EUCAST versions. Breakpoint differences account for much of
  the S_DET_POS for cefepime, ceftazidime, piperacillin-tazobactam and the
  fluoroquinolones.
- Five drugs are inferred from related drugs and all their rules are flagged
  for review; cefazolin in particular has low specificity.
- Porin loss, efflux and regulatory changes outside the PointFinder table are
  not mapped. This is the likely reason for many carbapenem, cefoxitin and
  amoxicillin-clavulanic acid R_DET_NEG endpoints, but it has not been
  checked.
- Phenotype quality varies by submitter, as PRJNA1322038 shows. Its labels
  were corrected by a whole-project swap, not re-interpreted from MICs and
  breakpoints. Projects with only categorical data could not be checked
  against MICs at all.
- Isolates are clustered by lineage, SNP cluster and project, so the
  confidence intervals are too narrow. Metrics are per drug; nothing is
  pooled across drugs and no significance tests are run.
- One AMRFinderPlus database (2026-08-07.1) was used. Database-version
  sensitivity was not run; the fingerprint is recorded so it can be.

## Future ML comparison

A later model will be compared against this baseline drug by drug, on the
same PR #5 endpoints with the decision 0004 correction and the PR #6 groups
(lineage, country and period holdouts from decision 0001). The comparison
should report, at least, the primary v1 baseline and the same metrics for the
model on identical endpoints, and say separately how the model does on
R_DET_NEG and S_DET_POS endpoints. Mapping v2 was informed by this cohort and
should not be used as the comparison baseline for it.

## Outputs and reproduction

The full builds (`genotype_phenotype_endpoints.tsv` is 63 MB) are in
`datasets/ecoli_research_v1/genotype_baseline_20261018/` and
`genotype_baseline_20261024/` and are not committed. The small tables
(metrics, denominators, summary, manifest, checksums, mapping audit, unmapped
determinants, country, temporal and source strata, the QC tables and the v2
metrics) are in `benchmark/cohorts/ecoli_research_v1/genotype_baseline_20261024/`
(primary) and `genotype_baseline_20261018/` (uncorrected labels). The build,
QC and v2 commands are in `call_run.json` and in the script headers; the
calls can be redone on NCI Gadi with `scripts/gadi/`.
