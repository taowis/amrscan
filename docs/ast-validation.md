# AST builder validation — 2026-10-10

Baseline: `main` at `4cd894afef126b6d3dade2b381cf976b3a0fd86b`.
Implementation branch: `feature/ast-dataset-builder`.
Runtime for local validation: Python 3.14.6, coverage.py 7.10.6,
Nextflow 25.04.6 and installed Java, BLAST+ and R dependencies.
CI retains Python 3.12, Java 17 and Nextflow 25.04.6.

## Automated checks

| Group | Command | Fixtures and assertions | Result |
| --- | --- | --- | --- |
| AST unit/regression | `python3 -m unittest discover -s tests -p test_ast_dataset.py -v` | Two one-row synthetic tables; derived malformed, missing, duplicate/conflicting and header-only cases; mocked HTTP/counts; 101-accession metadata batching; provenance and deterministic offline rebuild | PASS, 27 tests |
| Existing caller/harmonizer | `python3 -m unittest discover -s tests -p test_amrscan.py -v` | 24 upstream expected AMRFinderPlus records, tiny FASTA, mocked caller, invalid inputs and failure propagation | PASS, 14 tests |
| Combined coverage | `python3 -m coverage run --branch --source=bin -m unittest discover -s tests -p 'test_a*.py' -q` followed by `python3 -m coverage report --show-missing --fail-under=90` | Both Python unit suites; line and branch measurement of production code | PASS, 41 tests; aggregate 99.49% |
| Complete lightweight suite | `python3 -m unittest discover -s tests -p 'test_*.py' -v` | All unit, regression and workflow scenarios below | PASS, 48 tests, no skips |
| Nextflow/legacy regression | `python3 -m unittest discover -s tests -p test_nextflow.py -v` | Real module/config compilation, invalid inputs, duplicate IDs, gzip/multiple assemblies, real fixture harmonization, missing caller failure and BLAST+/R no-hit smoke | PASS, 7 tests, no skips |

Locally coverage.py was installed in a separate tools directory added through
`PYTHONPATH`. It is a test dependency, not a runtime dependency. Python coverage details:

| Module | Statements executed | Branches executed | Combined coverage |
| --- | --- | --- | --- |
| `bin/ast_dataset.py` | 416/417 (99.76%) | 179/180 (99.44%) | 99.66% |
| `bin/amrscan.py` | 136/137 (99.27%) | 57/58 (98.28%) | 98.97% |
| Total | 552/554 (99.64%) | 236/238 (99.16%) | 99.49% |

The unmeasured statements are subprocess entry guards; tests call `main`
directly. Real CLI build/validate smoke checks were run separately. No line
coverage percentage is claimed for Nextflow or R; those have scenario-based
integration checks. No unrelated flaky test failure was observed.

CI extends the existing coverage invocation to `test_a*.py` so it includes both
`test_amrscan.py` and `test_ast_dataset.py`. Workflow validation and the 90%
coverage floor remain. CI makes no NCBI requests.

CI runs `git diff --check` on the committed checkout.

## Current source and live small-selection check

Official source documentation, the browser export configuration and actual
export headers were inspected before implementation. See
[raw schema](../schemas/ast_raw_schema.md) for the exact endpoints/projections.
The test selection was the public BioSample `SAMN05170351`, associated with
`PDT000145880.5`, not the full *E. coli* cohort.

`curl -L --fail --get` retrieved both `action=solr2txt` tables and
`action=retrieve` counts. The queries were `biosample_acc:SAMN05170351` for AST
and `target_acc:(PDT000145880.5)` for metadata, with the documented fields,
`type=tsv` and `nolimit=on`. Actual counts matched the official count responses:

| Source | Rows | SHA-256 |
| --- | --- | --- |
| AST | 26 | `2ddb8c1931ac1c21260088a61e1932291c802eb880049eaead66ab63019c99ce` |
| Isolate metadata | 3 | `df699fce9707e801060699d699f8f39bd9c0450e6a0b7d6e0c6e7cef7fd1f67b` |

The three metadata rows shared a PDT accession and identical selected values.
The join retained all three source rows without expanding 26 AST observations
into 78. Combined MIC components survived normalization; the amoxicillin–
clavulanic acid source measurement was `>32.0` with a second component `16.0`.
There were 16 R, 6 S and 4 ND records. All 26 were *E. coli*. The metadata
supplied SRA `SRR4065706`, but no Assembly accession. Method and breakpoint
version were also absent. No Assembly accession or breakpoint version was
invented: assembly-linked isolates and eligible assembly-based endpoints were
both zero. All original rows remained available.

The archived snapshot was built and checked offline with:

```bash
python3 bin/ast_dataset.py build \
  --snapshot snapshot/raw \
  --outdir snapshot/processed
python3 bin/ast_dataset.py validate \
  --snapshot snapshot/raw \
  --processed snapshot/processed
```

Both returned success; validation rebuilt and compared every output byte. Raw
and processed real data stayed outside Git. The tiny committed fixtures are
synthetic and are not copies of this live sample.

## Validation limits and review decisions

The production Python `download` command was attempted with the same narrow
BioSample query. The validation machine returned a DNS resolution error
from `urllib`; it published no partial snapshot. Thus **official live table
retrieval and offline integration passed, but a full live invocation of the
Python downloader has not passed**. Its HTTP handling, byte limits, counts,
metadata batching and error cleanup are tested offline with mocks. Re-run the
narrow download command on a host where Python can resolve NCBI before using
it for a study. Full-cohort scale and future NCBI releases remain unvalidated.

I also checked that source rows and nonbinary categories survive, metadata
joins do not multiply observations, conflicts never choose a preferred row,
MIC parsing never assigns phenotype, and all exclusions have explicit units
and reasons.

Scientific review is still required for assay comparability, standards/versions,
nonbinary endpoints, conflict adjudication, the initial drug/country/source
vocabularies, species-name filtering, biological duplicates across PDT versions,
assembly eligibility and antibiotic inclusion. These are explicit research
policy decisions, not silently implemented statistical assumptions.
