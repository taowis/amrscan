# AMRScan result TSV schema 1.0

One row represents one caller-reported element/reference match. Rows are not
deduplicated or ranked by AMRScan. File encoding is UTF-8, with a tab delimiter,
a header, and standard CSV quoting (including doubled quotes in JSON cells).
Use a TSV/CSV parser, not string splitting, to recover `raw_record_json`.
The literal `NA` represents missing or not applicable values. A header-only
file means the caller reported no elements, not phenotypic susceptibility.

| Column | Type / interpretation |
| --- | --- |
| sample_id | Validated filename-derived sample identifier |
| tool | `AMRFinderPlus`; reserved for future caller adapters |
| gene | Caller element symbol, verbatim; may be a family, allele or mutation label |
| allele | Element symbol only for `ALLELE`, `ALLELEX`, `ALLELEP`; otherwise `NA` |
| resistance_class | Caller `Class` for `Type=AMR`; otherwise `NA` |
| drug | `NA` in this adapter; no individual drug or AST endpoint is inferred |
| identity | Caller percent identity, 0–100, or `NA`; reference sequence type matters |
| coverage | Caller percent of reference covered, 0–100, or `NA`; not read depth |
| contig | Caller contig ID, or `NA` |
| start / end | Caller Start / Stop: 1-based inclusive nucleotide coordinates, start <= end, or both `NA` |
| strand | `+`, `-`, or `NA`; independent of coordinate order |
| evidence_type | Original caller `Method`, including partial, stop and mutation methods |
| schema_version | `1.0` |
| caller_version | Actual `amrfinder --version` output; fixture-only conversions use `NA` |
| database_id | `sha256:` digest of canonical JSON of database relative paths/file hashes |
| thresholds_json | Requested identity/coverage, translation table, organism, plus and equal-hit policy |
| element_type / element_subtype | Caller `Type` / `Subtype`, retaining AMR/STRESS/VIRULENCE distinctions |
| scope | Caller core/plus assignment |
| subclass | Caller `Subclass`, without asserting it is an individual drug |
| raw_file | Original report basename, found under `raw/amrfinderplus/<sample_id>/` |
| raw_row | One-based data-record index, excluding the header |
| raw_record_json | Every original column/value, including future unknown columns, preserved as JSON |

The adapter targets AMRFinderPlus v4 nucleotide reports. Missing/duplicate
required headers, malformed rows, invalid numbers/coordinates, wrong provenance
or changed report checksums fail explicitly. It does not silently drop rows,
infer gene families from allele names, assign the nearest reference allele to a
BLAST hit, or translate subclass strings into susceptibility interpretations.
For mutations, `gene` retains the full mutation label rather than guessing a
gene/variant split. `raw_record_json` retains the closest reference, hierarchy,
HMM and other annotations even when they cannot fit the common columns.

The adjacent provenance JSON also stores the exact argument vector, input
SHA-256, database source path and per-file hashes, raw report hash and completion
status. The untouched caller log retains runtime database/software information.
The database digest is a content identifier, not a fabricated release number.
Keep a release-named database snapshot and the original log for publication.

Future RGI/CARD and ResFinder adapters must populate these semantics explicitly,
use `NA` for unavailable values, and preserve their own original records and
provenance. Cross-caller reconciliation is a separate,
versioned layer. No resistant/susceptible classification is produced here.

Sources: [NCBI report fields](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format)
and [interpretation](https://github.com/ncbi/amr/wiki/Interpreting-results).
