# AST dataset builder

The builder preserves public AST observations and genome-linking metadata for
*Escherichia coli* isolates. It does not apply breakpoints or infer phenotype
from genotype.

AMRScan does not infer phenotypic resistance from the presence of an AMR gene.
AST phenotype and genomic AMR determinants remain separate evidence layers.

## Download and rebuild

Run from the repository root with Python 3.10+:

```bash
# Small public selection for checking the source interface first.
python3 bin/ast_dataset.py download \
  --query 'biosample_acc:SAMN05170351' --outdir datasets/check/raw

# Default selection: NCBI scientific_name:"Escherichia coli".
python3 bin/ast_dataset.py download --outdir datasets/ecoli/raw

# These operations are offline.
python3 bin/ast_dataset.py build \
  --snapshot datasets/ecoli/raw --outdir datasets/ecoli/processed
python3 bin/ast_dataset.py validate \
  --snapshot datasets/ecoli/raw --processed datasets/ecoli/processed
python3 bin/ast_dataset.py summarize --snapshot datasets/ecoli/raw
```

`download` retrieves measurement-level AST and associated isolate tables from
NCBI's official table exporter. The [raw schema](../schemas/ast_raw_schema.md)
documents URLs, parameters, field names and schema-change behavior. Requests
are sequential, use a 60-second socket timeout, and have a 128 MiB per-file
limit (`--max-bytes`). Narrow the query when this is insufficient. Failed
requests do not publish a snapshot. There is no automatic retry or resume;
retry explicitly into a new path. No genome sequence is downloaded.

`validate --snapshot` checks raw provenance/schema and computes the cohort.
With `--processed`, it rebuilds offline and compares the complete expected
output set byte for byte, detecting stale or modified results. `summarize`
recomputes a summary from the snapshot, not from a potentially stale report.
Processing loads the selected tables into memory; full-scale streaming and
cross-release reconciliation are outside this version.

## Offline demonstration

The two fixtures are invented test data and are not breakpoint examples:

```bash
python3 bin/ast_dataset.py archive \
  --ast tests/fixtures/ast.synthetic.tsv \
  --metadata tests/fixtures/isolates.synthetic.tsv \
  --source-url fixture:synthetic \
  --retrieved-at 2026-10-10T00:00:00Z \
  --outdir datasets/demo/raw
python3 bin/ast_dataset.py build \
  --snapshot datasets/demo/raw --outdir datasets/demo/processed
python3 bin/ast_dataset.py validate \
  --snapshot datasets/demo/raw --processed datasets/demo/processed
```

For previously downloaded source tables use `archive` with their true source
and retrieval time. It requires the documented projection, not an arbitrary
CSV downloaded with different browser display columns. Keep all manifest
entries when transporting the raw directory. `datasets/`, root `raw/` and
root `processed/` are ignored by Git; use these locations for generated data.

## Interpretation and cohort policy

Every source row appears in `ast.normalized.tsv`. The *E. coli* subset retains
rows whose resolved scientific name is exactly `Escherichia coli` or starts
with that binomial followed by a space and a strain/subspecies descriptor.
Generic `Escherichia`, Shigella, the joint taxgroup label and unresolved species
conflicts are excluded from this subset with `organism_mismatch`. This is a
source-name filter, not a new taxonomic classification.

Metadata joins use the entire versioned PDT accession. There is no fuzzy
identifier match, cross-version collapse, BioSample-only genome assignment or
Cartesian expansion. All matching metadata rows and their provenance are
retained. A single nonmissing consensus value fills a normalized field;
disagreeing nonmissing values make that field missing and add
`metadata_conflict`. This also applies to disagreement between AST and isolate
metadata. Multiple identical metadata rows do not duplicate AST observations.

The [versioned policy](../ingestion/ncbi_pathogen_detection/policy.json)
contains every drug and phenotype alias. Only surrounding whitespace and case
are normalized before exact lookup. Unmapped drug names remain visible with
`unknown_antimicrobial`. The original drug and phenotype are always retained.

| Source phenotype, case-insensitive | Normalized |
| --- | --- |
| R, resistant | R |
| S, susceptible, sensitive | S |
| I, intermediate | I |
| SDD, SSD, susceptible-dose dependent | SDD |
| NS, nonsusceptible | NS |
| HLAR, high level aminoglycoside resistance | HLAR |
| N, ND, not defined | ND |
| Missing token | NA |
| Any other nonmissing value | UNKNOWN |

I, SDD, NS and HLAR remain separate nonbinary categories. ND is a recognized
undefined category, not missing data and not a binary label. Submitter
phenotypes are never recalculated from MIC or resistance genes.

MIC parsing accepts positive decimal values and `<`, `<=`, `=`, `>=`, `>`;
`==` becomes `=`. An embedded comparator is accepted if it agrees with a
separate supplied comparator. A bare number denotes equality. Invalid values,
negative/zero concentrations, nonfinite values, conflicting signs or a second
component without a valid first component receive `invalid_mic`; raw values
are retained. Combined concentrations use the source's separate `mic` and
`mic_secondary` columns. A single `32/4` string is preserved but flagged rather
than split without a verified source contract. No clinical range or breakpoint
is applied. Disk diffusion remains separate, with its raw measurement and sign;
it is not converted to MIC. Detailed disk-diffusion QC is future work.

Collection dates retain their source precision. Only valid `YYYY`, `YYYY-MM`
and `YYYY-MM-DD` yield a year. Date intervals, non-ISO dates and invalid calendar
dates remain available but are flagged with missing normalized year. No month
or day is imputed. Country is the explicit prefix of NCBI `geo_loc_name` before
`:`, accepted only when present in the policy's country allowlist. Locality is
never geocoded. The unmodified source location remains available. Source
categories use the narrow exact host/source mappings described in the
[adapter notes](../ingestion/ncbi_pathogen_detection/README.md).

## Repeated measurements and endpoint eligibility

The observation identity is the complete raw AST row excluding only the
source's `id`. Exact repeats within a versioned isolate–normalized-drug group
are flagged on **all** repeated rows, retained, and counted once when assessing
endpoint eligibility. Distinct rows, even with agreeing phenotypes, receive
`multiple_ast_records` and are not collapsed to a preferred assay. Differences
in MIC, standards, methods, platform, reagent or supplied metadata therefore
require review. Alias-equivalent but otherwise nonidentical raw rows also stay
separate under this conservative policy.

More than one distinct recognized phenotype (R/S/I/SDD/NS/HLAR) within the group
adds `conflicting_ast`. ND, UNKNOWN and missing phenotypes are not themselves
assertions of susceptibility or resistance. They still block endpoint
eligibility as undefined evidence. A group containing any mismatched organism
row is also ineligible even if another row passes the species filter.

`ecoli_endpoints.tsv` records one group per isolate–drug combination with all
record IDs and reasons. Missing isolate identifiers each form a separate
unresolved group and are ineligible. An eligible group requires:

- explicit *E. coli* identity, an isolate ID and an Assembly accession;
- a recognized drug and a recognized, defined phenotype;
- one distinct raw observation (exact repeats are allowed);
- no invalid MIC or unresolved AST/metadata conflict.

An SRA link alone is retained but is insufficient for this assembly-based
endpoint. Missing MIC is allowed when the submitter provided an explicit
phenotype. Missing country, year, source category, method or breakpoint version
is reported but does not erase the observation or automatically block this
structural eligibility. **Eligible does not mean assay comparability has been
established.** Method/standard/version adjudication and any later binary-label
policy require a separately reviewed analysis plan. No rule treats CLSI and
EUCAST as interchangeable.

## Outputs and accounting

| Output | Unit/purpose |
| --- | --- |
| `ast.normalized.tsv` | Every AST source record, with provenance and QC |
| `ecoli_ast.tsv` | All records passing the species filter, including unresolved AST |
| `ecoli_isolates.tsv` | Distinct nonmissing PDT IDs, metadata value sets, record counts |
| `ecoli_endpoints.tsv` | Isolate–drug groups, eligibility and all exclusion reasons |
| `qc_flags.tsv` | One row per record–flag pair |
| `cohort_summary.json` | Record counts, pair counts, distributions, flow, policy/source hashes |
| `antimicrobial_summary.tsv` | Drug-specific record/isolate counts and structural readiness |

The [cohort schema](../schemas/ast_cohort_schema.md) defines columns and flags.
Summary distributions count *E. coli* AST records, not independent isolates.
The flow lists record exclusions and then switches explicitly to isolate–drug
group units. Exclusion-reason counts overlap; the disjoint total is
`ineligible_pairs`. Original rows are never deleted merely for missing study
covariates.

In the drug summary, R/S/other count eligible isolate–drug endpoints. Other
means I/SDD/NS/HLAR; it is not merged with R or S. `R_fraction_of_RS` uses only
R+S as denominator and is missing when that denominator is zero. Country,
year, assembly and source-completeness statistics describe all retained *E.
coli* records for that drug, counting distinct nonmissing isolate IDs. These
columns provide selection evidence, not population resistance estimates or an
automatically selected antibiotic panel.

## Limits and next work

The [validation record](ast-validation.md) distinguishes offline tests, actual
NCBI table retrieval, and the downloader's unvalidated live Python transport.
NCBI indexes change; the interface is not transactionally frozen. Source data
may lack method, standard version, Assembly accession or complete metadata.
The initial country/drug/source dictionaries are intentionally incomplete and
must be reviewed before use on a full cohort. Source-name taxonomy filtering
and versioned PDT identity do not resolve duplicate biological isolates across
releases. The current tables are held in memory and the full public cohort has
not been downloaded or performance-tested.

Human review must decide assay comparability, handling of nonbinary categories,
conflict adjudication, additional vocabulary entries, antibiotic inclusion and
minimum metadata requirements. This layer supplies no lineage labels. Existing
`isolate_identifiers` are preserved as source metadata, not parsed into a
validated lineage assignment.

The next step is a pinned AMRFinderPlus environment with known positive/negative
controls and a small, accession-verified assembly-linked cohort. Preserve the
separate AST and genotype evidence while checking sample/assembly consistency.
