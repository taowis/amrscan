# NCBI Pathogen Detection AST adapter

`bin/ast_dataset.py` is the executable adapter. The caller and Nextflow workflow
remain separate. Python 3.10+ and its standard library suffice at runtime.

`policy.json` version 1.0.0 contains explicit drug aliases, phenotype categories,
source-category rules and a conservative country allowlist. The policy hash is
written to every cohort summary. Review mappings as scientific inputs when
changing this file; increment its version when semantics change.

- [User guide](../../docs/ast-dataset-builder.md)
- [Raw source contract](../../schemas/ast_raw_schema.md)
- [Normalized/cohort contract](../../schemas/ast_cohort_schema.md)
- [Validation evidence](../../docs/ast-validation.md)

The 44 drug names and their aliases use the
[NCBI BioSample vocabulary](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/)
inspected on 2026-10-10. This is a selected vocabulary, not an exhaustive list.
`TMP-SMX` is an additional explicit project alias for
trimethoprim-sulfamethoxazole. NCBI's `SSD` and the commonly used `SDD` are
explicit aliases for the separate SDD category. No fuzzy matching is used.
Unlisted drugs are retained and flagged, including valid drugs outside this
initial vocabulary.

The 67 country strings are a selected subset of the
[INSDC geographic vocabulary](https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/)
inspected on 2026-10-10. They retain source spelling; this is not a complete
country registry. Unlisted names, historical names, oceans and free-text places
keep their raw location but have missing normalized country. Expand this list
only after reviewing the source vocabulary, without assigning locations to
modern countries by inference.

Host and isolation-source mappings are deliberately narrow project rules. They
use exact trimmed, case-insensitive matches. If host and source imply different
categories, the result is `unknown` with `ambiguous_source`; neither wins by
precedence. `clinical`, `environmental/other`, specimen type and taxgroup labels
alone do not imply human or environmental origin.
