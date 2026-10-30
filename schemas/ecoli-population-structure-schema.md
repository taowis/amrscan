# E. coli population structure output schema

All tables are UTF-8 TSV with a header and `NA` for unavailable scalar values.
JSON-valued cells are compact, sorted-key JSON. Files are immutable outputs of
`bin/population_structure.py build`; rebuild into a new directory.

| File | Row unit | Meaning |
| --- | --- | --- |
| `population_metadata.tsv` | PR #5 isolate | Identity, genome hash and QC, MLST, phylogroup, NCBI SNP cluster, lineage, metadata, phenotype counts and genomic eligibility flags. Every PR #5 isolate appears once. |
| `genome_qc.tsv` | PR #5 isolate | Assembly statistics, `pass`/`warning`/`fail` and the flags behind them. |
| `mlst.tsv` | PR #5 isolate | Scheme, ST (only when assigned), status, allele profile, tool and database versions, input SHA-256. |
| `phylogroup.tsv` | PR #5 isolate | Phylogroup, status, raw tool call, method, version, input SHA-256. |
| `evaluation_groups.tsv` | PR #5 isolate | Grouping variables and holdout eligibility, plus `lineage_component`, `lineage_cv_test_eligible`, `split_unit` (SNP cluster or isolate, never divided), `project_group`, `geography_group` (`reference`/`non_reference`) and `period` (`earlier`/`recent`). No train/validation/test assignment. |
| `lineage_summary.tsv` | Lineage (ST or UNKNOWN) | Isolate count, cohort fraction, country/year/source spread, binary endpoints, phylogroups. |
| `lineage_antimicrobial_summary.tsv` | Antimicrobial × lineage | Binary R/S counts of population-eligible endpoints. |
| `lineage_country_summary.tsv` | Lineage × country | Isolate count. |
| `lineage_year_summary.tsv` | Lineage × collection year | Isolate count. |
| `lineage_source_summary.tsv` | Lineage × source category | Isolate count. |
| `snp_cluster_summary.tsv` | NCBI SNP cluster | Isolates, countries, years, lineages, projects, binary R/S endpoints. |
| `project_summary.tsv` | BioProject × submitter | Isolates, countries, lineages, binary R/S endpoints. |
| `holdout_feasibility.tsv` | Antimicrobial × holdout axis | Group counts, R/S balance, screening status from explicit config rules, and `HUMAN_REVIEW_REQUIRED`. |
| `evaluation_plan.tsv` | Antimicrobial | Tier (`primary`/`secondary`/`excluded`) and reason, the numbers behind it, period and non-reference R/S counts, and `yes`/`no` for each analysis, following `evaluation_design` in the config and decision 0001. |
| `confounding_summary.tsv` | Review flag | Risk type, subject, metric, value, threshold and detail. |
| `population_summary.json` | Population | Counts and reconciliation checks. |
| `population_manifest.json` | Build | PR #5 cohort hashes, genome retrieval, tool versions, config hash, code commit, annotation counts. |
| `SHA256SUMS` | File | SHA-256 for every other output file. |

`lineage_id` is `ST<n>` for an assigned Achtman 7-gene ST and `UNKNOWN` for any
other MLST outcome; `UNKNOWN` is a missing label, not a lineage
(`is_lineage = false`). `ncbi_snp_cluster` is the unversioned PDS accession;
the original versions are in `ncbi_snp_cluster_versions_json`. MLST,
phylogroup and SNP cluster are separate fields and are never substituted for
one another. Phenotype columns are counts joined from PR #5 and never alter
PR #5 endpoints.
