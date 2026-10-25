# AST cohort schema, version 1.0.0

`schemas/ast.tsv` is the machine-readable header contract for both normalized
and E. coli record tables. It replaces the unimplemented foundation placeholder.
All tables use UTF-8, tab delimiters, standard CSV quoting and LF line endings.

Normalized missing values are the literal `NA`. Recognized input missing tokens
are empty/whitespace, NA, N/A, null, missing, not collected, not provided and
not applicable (case-insensitive). Original fields and JSON retain exact source
cells, including empty strings; JSON arrays may be empty. `unknown` is an explicit
source category and `UNKNOWN` an unmapped phenotype. Neither implies a phenotype.

## Normalized fields, in output order

| Field | Meaning |
| --- | --- |
| `record_id` | SHA-256 of AST source hash plus one-based data-row number; stable within the snapshot. |
| `isolate_id` | Versioned PDT accession, the exact metadata join key. |
| `biosample` | Consensus BioSample accession. |
| `assembly_accession` | Consensus Assembly accession; unresolved conflicts remain missing. |
| `sra_accession` | Consensus source Run value, preserved as a string, including supplied multi-accession lists. |
| `organism` | Consensus scientific name; the taxgroup is not a species label. |
| `antimicrobial` | Explicit canonical drug name, or original trimmed unknown name. |
| `antimicrobial_original` | Exact original antibiotic cell, including an empty cell. |
| `phenotype` | R, S, I, SDD, NS, HLAR, ND, UNKNOWN or NA; never inferred from MIC. |
| `phenotype_original` | Exact original phenotype cell. |
| `mic` | Positive decimal first-component concentration, if safely parsed. |
| `mic_operator` | <, <=, =, >=, >, or NA. |
| `mic_units` | mg/L from the source-column definition when a MIC was supplied, otherwise NA. |
| `mic_raw` | Exact source mic cell; comparator may be embedded in archived input. |
| `mic_secondary` | Positive decimal second-component concentration, if safely parsed. |
| `mic_secondary_raw` | Exact source mic_secondary cell. |
| `measurement_sign_original` | Exact measurement_sign cell; also retains the sign for disk diffusion. |
| `disk_diffusion_raw` | Exact source disk-diffusion measurement in mm, separate from MIC. |
| `testing_method` | Source method only if supplied in an explicitly enriched archive; current download leaves NA. |
| `breakpoint_standard` | Submitted standard, retained without reconciling standards. |
| `breakpoint_version` | Explicit standard_version in an enriched archive; never derived from reagent or date. |
| `platform` | Submitted laboratory typing platform. |
| `vendor` | Submitted test vendor. |
| `reagent` | Submitted kit/method version or reagent; not the breakpoint version. |
| `collection_date` | Consensus source date, precision preserved; exact cells also remain in raw JSON. |
| `collection_year` | Year from a valid YYYY, YYYY-MM or YYYY-MM-DD only. |
| `country` | Allowlisted explicit country prefix of geo_loc_name; otherwise NA. |
| `location` | Consensus source geo_loc_name, retaining locality. |
| `isolation_source` | Consensus source specimen/environment description. |
| `host` | Consensus source host value. |
| `source_category` | human, animal, food, environment or unknown under explicit policy rules. |
| `source_file` | AST snapshot basename. |
| `source_row` | One-based data-row ordinal after the header, not a physical text line number. |
| `source_url` | Exact AST source endpoint or archived-source description from the manifest. |
| `source_sha256` | Hash of the complete original AST file. |
| `retrieved_at` | AST source retrieval timestamp with timezone. |
| `source_record_id` | Source AST id, distinct from the snapshot record_id. |
| `metadata_records_json` | All matching metadata rows with source filename/row/hash/URL/retrieval time. |
| `raw_record_json` | Complete original AST row keyed by source field names; no original cell is normalized. |
| `schema_version` | 1.0.0, the normalized contract version. |
| `policy_version` | Version of the explicit mapping and cohort policy. |
| `qc_flags` | Sorted semicolon-separated flags; empty means no flags. |

## QC vocabulary

`qc_flags.tsv` has `record_id` and `flag`. Every flag below is also present in
the corresponding normalized row. Flags report evidence quality; only the
explicit endpoint policy determines structural eligibility.

| Flag | Meaning |
| --- | --- |
| `ambiguous_country` | Location prefix is absent from the country allowlist; includes valid unlisted countries. |
| `ambiguous_source` | Explicit host and source rules imply different source categories. |
| `conflicting_ast` | More than one distinct recognized defined phenotype in the isolate–drug group. |
| `duplicate_record` | Raw AST observation excluding id repeats within the isolate–drug group; every copy flagged. |
| `invalid_collection_date` | Nonmissing date is ambiguous, non-ISO, or invalid; raw value retained. |
| `invalid_mic` | Invalid concentration/comparator or second component without a valid first component. |
| `metadata_conflict` | Multiple nonmissing values disagree in a normalized metadata field. |
| `missing_assembly` | No resolved Assembly accession; blocks assembly-based endpoint eligibility. |
| `missing_breakpoint_standard` | Submitted testing standard absent. |
| `missing_breakpoint_version` | Explicit breakpoint standard version absent. |
| `missing_collection_year` | No year can be derived under the date policy. |
| `missing_country` | No safely normalized country. |
| `missing_genome_link` | Neither an Assembly nor an SRA link is available. |
| `missing_isolate_id` | No PDT identifier; cannot link metadata or collapse with other missing identifiers. |
| `missing_metadata` | No exact versioned PDT match in the archived metadata. |
| `missing_phenotype` | Phenotype is a missing token. |
| `missing_source` | Both host and isolation source are missing. |
| `missing_testing_method` | Explicit testing method absent. |
| `multiple_ast_records` | The group has more than one distinct raw AST observation. |
| `organism_mismatch` | Resolved scientific name does not pass the E. coli filter. |
| `unknown_antimicrobial` | The supplied drug is outside the versioned mapping, including missing drugs. |
| `unknown_phenotype` | A nonmissing phenotype has no explicit mapping. |

## Derived tables and summary

`ecoli_isolates.tsv` fields, in order:

```text
isolate_id
ast_records
biosample_values_json
assembly_accession_values_json
sra_accession_values_json
collection_year_values_json
country_values_json
source_category_values_json
```

`ecoli_endpoints.tsv` fields, in order:

```text
isolate_id
antimicrobial
phenotype
eligible
exclusion_reasons
record_ids_json
record_count
distinct_observations
```

`antimicrobial_summary.tsv` fields, in order:

```text
antimicrobial
ast_records
isolates
eligible_pairs
R
S
other
R_fraction_of_RS
countries
year_min
year_max
assembly_linked_isolates
country_known_isolates
year_known_isolates
source_classified_isolates
```

Isolate metadata sets are sorted JSON arrays containing all distinct nonmissing
values among retained E. coli rows. They do not assert a single agreed country
or date across different AST records. Missing identifiers are counted in the
record layer but do not become a fictitious isolate row.

Endpoint `eligible` is the string `true` or `false`. `record_count` includes
all retained copies; `distinct_observations` excludes repeats identical apart
from source id. `record_ids_json` links all evidence. `exclusion_reasons` is a
sorted semicolon-separated set of blocking flags, plus `undefined_phenotype`
for NA/UNKNOWN/ND. The latter is an endpoint reason, not a row QC flag.

`cohort_summary.json` records all-source QC counts, E. coli record distributions,
unique isolates, assembly linkage, defined phenotypes, source-file/metadata
counts, policy version/hash and the manifest hash. `flow` entries state stage,
unit and count. Pair exclusion reasons overlap; never add them to infer the
number of excluded groups. The full eligibility rules and drug-summary
denominators are in the [builder guide](../docs/ast-dataset-builder.md).
