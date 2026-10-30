# E. coli genotype–phenotype baseline output schema

Outputs of `bin/genotype_baseline.py build`. All tables are UTF-8 TSV with a
header and `NA` for unavailable values. A build is written to a new directory
and never overwrites an earlier one; `SHA256SUMS` covers every other file.

| File | Row unit | Meaning |
| --- | --- | --- |
| `genotype_phenotype_endpoints.tsv` | PR #5 isolate × antimicrobial endpoint | Every PR #5 endpoint with its genotype status, determinant evidence, concordance category, PR #6 grouping fields and exclusion reason. |
| `baseline_metrics.tsv` | Study antimicrobial × rule variant | Confusion counts and metrics with Wilson 95% intervals. |
| `baseline_denominators.tsv` | Antimicrobial | Endpoint counts from the PR #5 total to the metric-eligible set; the exclusions add up to the total. |
| `lineage_baseline.tsv` | Dimension (`lineage`, `phylogroup`) × stratum × antimicrobial | Stratified counts and, with enough support, metrics. |
| `country_baseline.tsv` | (`country`, `geography_group`) × stratum × antimicrobial | As above. |
| `temporal_baseline.tsv` | (`period`, `collection_year`) × stratum × antimicrobial | As above. |
| `source_baseline.tsv` | `source` × stratum × antimicrobial | As above. |
| `error_distribution.tsv` | (`project`, `snp_cluster`) × stratum × antimicrobial | The same counts by submitting project and NCBI SNP cluster. |
| `unexplained_resistance.tsv` | Metric-eligible endpoint with phenotype R and no mapped determinant | Candidates for unexplained resistance with lineage, place, time, source, QC and provenance. |
| `susceptible_with_determinant.tsv` | Metric-eligible endpoint with phenotype S and a mapped determinant | The determinants that made the call, with the same context. |
| `mapping_audit.tsv` | Antimicrobial × determinant with a rule | Carriers among population-eligible binary endpoints, their R/S split, sole-support counts, top lineages and countries, and review flags. |
| `unmapped_determinants.tsv` | Detected AMR element with no rule for any study drug | Isolate counts; these never change a determinant status. |
| `baseline_summary.json` | Build | Counts by caller status, genotype status and concordance category, and the definitions used. |
| `baseline_manifest.json` | Build | Upstream PR #5 and PR #6 hashes, caller and database versions and fingerprint, mapping version and hashes, configuration hash, code commit and builder hash. |

Post-build checks from `scripts/ecoli_research_v1/baseline_qc.py` (not part
of the build, and not covered by its `SHA256SUMS`):

| File | Row unit | Meaning |
| --- | --- | --- |
| `qc/phenotype_label_check.tsv` | Project × antimicrobial with at least 3 R and 3 S MICs | Median MIC of R and S observations under the PR #5 labels; `inverted` when R is lower. |
| `qc/phenotype_label_check_corrected.tsv` | As above | The same check after the corrections in `--corrections` are applied. |
| `qc/corrected_mic_out_of_order.tsv` | MIC observation of a corrected project | Corrected R at or below the median MIC of corrected S, or corrected S at or above the median of corrected R, within project, antimicrobial and breakpoint standard. |
| `qc/discordance_concentration.tsv` | Study antimicrobial × concordance category | Largest lineage, country, project, source, period and determinant set, with shares. |
| `qc/metrics_without_flagged_projects.tsv` | Study antimicrobial × subset | Primary metrics for all projects and without projects flagged by the label check; post hoc. |
| `qc/qc_summary.json` | Check run | Input hashes, the label-check rule, the flagged projects and, with corrections, the counts after correction. |

## Endpoint fields

`phenotype` is the label used for metrics. `phenotype_pr5` is the PR #5
endpoint as frozen, and `phenotype_correction` is the correction ID from the
configured correction table when the two differ, otherwise `NA`. The only
correction type is `swap_R_S` for every binary endpoint of a submitting
project (decision 0004).

`genotype_status` is one of:

- `determinant_positive`, `determinant_negative`, `determinant_uncertain`:
  a successful, provenance-checked call interpreted with the mapping;
- `genome_unavailable`, `genome_qc_failure`, `caller_missing`, `caller_failed`,
  `provenance_mismatch`, `harmonization_failed`: no interpretable call. None of
  these is ever counted as negative;
- `mapping_unavailable`: the antimicrobial is not a study drug in decision 0001,
  so the mapping has no rules for it.

`determinant_status` repeats the interpreted state (`positive`, `negative`,
`uncertain`) or is `NA`. `determinant_ids` lists the determinants that make the
call positive, sorted, and `determinant_types` gives `gene` or
`point_mutation` for each in the same order. `uncertain_determinant_ids`
lists determinants whose rules were uncertain, partial, or missing a required
co-mutation.

`concordance_category` is `R_DET_POS`, `R_DET_NEG`, `S_DET_POS` or `S_DET_NEG`
for binary phenotypes with a positive or negative state, otherwise `NA`. These
are neutral labels. In `baseline_metrics.tsv` phenotype R is the positive
reference state and determinant-positive the baseline positive call, so TP =
R_DET_POS, FN = R_DET_NEG, FP = S_DET_POS and TN = S_DET_NEG.

`exclusion_reason` gives the first reason, in this order, an endpoint is not
metric-eligible: `non_binary_or_conflict` (PR #5), `population_ineligible`
(PR #6), `not_study_antimicrobial`, `genotype_unavailable:<status>`,
`mapping_uncertain`.

`fq_two_mechanisms_status` is the ciprofloxacin/levofloxacin state under the
prespecified two-mechanism variant; `NA` for other drugs.

## Rule variants

| `rule_variant` | Uncertain endpoints |
| --- | --- |
| `primary` | excluded from N |
| `uncertain_as_positive` | counted as determinant-positive |
| `uncertain_as_negative` | counted as determinant-negative |
| `fq_two_mechanisms` | ciprofloxacin and levofloxacin only; positive needs two predicting determinants |

## Stratum fields

`support_status` is `reported` when the stratum has at least `min_R` R and
`min_S` S endpoints (`baseline_config.json`); otherwise metrics are `NA` and
the reason is stated. `r_det_neg_fraction` is R_DET_NEG / R. `flags` may
contain `only_R`, `only_S`, `extreme_R_S_imbalance` (minority class below 5%),
`one_lineage_majority` and `one_project_majority` (one lineage or project holds
at least half the stratum). Missing values are kept as the stratum `unknown`.

## Mapping files

`benchmark/mappings/ecoli_amrfinder_drug_mapping_v1.tsv` has one row per
determinant × antimicrobial with columns `mapping_version`, `organism`,
`rule_id` (`determinant|antimicrobial`), `determinant` (exact AMRFinderPlus
symbol), `determinant_type`, `gene_family`, `allele_or_mutation`,
`antimicrobial`, `interpretation` (`predicts_resistance`, `uncertain`,
`does_not_predict`, `component`, `HUMAN_REVIEW_REQUIRED`), `requires`
(co-mutations, or the component name for a combination rule),
`evidence_source`, `evidence_terms`, `evidence_reference`, `evidence_level`,
`review_status` and `notes`. The JSON beside it records the SHA-256 of the TSV,
the study drugs, the combination rule, the rule variants, the source commits
and hashes, and the references. See decision 0002.
