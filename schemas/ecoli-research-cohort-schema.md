# E. coli research cohort output schema

All tables are UTF-8 TSV with a header and `NA` for unavailable scalar values.
JSON-valued cells are compact, sorted-key JSON. Files are immutable outputs of
`bin/ecoli_cohort.py`; rebuild into a new directory.

| File | Row unit | Meaning |
| --- | --- | --- |
| `ast_observations.tsv` | Source AST observation | Normalized source rows, unchanged phenotype evidence and source provenance. |
| `candidate_isolates.tsv` | Versioned PDT accession | Consensus-linked isolate metadata, source hashes, and identity/metadata conflicts. |
| `isolate_drug_endpoints.tsv` | Isolate × antimicrobial | Derived phenotype category, observation counts, duplicate/conflict status, and analysis eligibility flags. |
| `eligible_isolates.tsv` | Overall-eligible isolate | Distinct isolates with at least one overall-eligible endpoint. |
| `exclusions.tsv` | Isolate × antimicrobial × analysis | Explicit exclusion reason for each failed eligibility dimension. |
| `antimicrobial_summary.tsv` | Antimicrobial | Isolate-level R/S counts, conflicts, metadata availability, and readiness status. |
| `country_summary.tsv` | Country | Isolate and endpoint counts with binary R/S counts. |
| `year_summary.tsv` | Collection year | Isolate and endpoint counts with binary R/S counts. |
| `source_summary.tsv` | Source category × antimicrobial | Isolate and endpoint counts with binary R/S counts. |
| `missingness_summary.tsv` | Unit × field | Explicit missing count and fraction for observation and isolate metadata. |
| `conflict_summary.tsv` | Conflicted isolate × antimicrobial | Conflicted endpoint categories and contributing source AST identifiers. |
| `cohort_flow.json` | Cohort | Reconciled counts with the unit named for every count. |
| `cohort_manifest.json` | Build | Source snapshot, file hashes, parser/config/builder versions and code commit. |
| `SHA256SUMS` | File | SHA-256 checksums for every other output file. |

Phenotype categories R and S are the only binary values. Repeated same-category
observations produce one endpoint; discordant R/S observations produce
`CONFLICT`; mixed binary/nonbinary observations produce
`NON_BINARY_OR_CONFLICT`. Nonbinary, missing, and unknown categories are retained
as observed. No AMR genotype is used to define or resolve an AST endpoint.
