# Lightweight tests

Run from the repository root:

```bash
python3 -m unittest discover -s tests -p test_amrscan.py -v
python3 -m unittest discover -s tests -p test_ast_dataset.py -v
python3 -m unittest discover -s tests -p test_amrfinder_validation.py -v
python3 -m unittest discover -s tests -p test_nextflow.py -v
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The unit suite validates 24 authentic upstream records, including exact alleles,
BLAST family assignments, partials, internal stops, negative strands, multiple
mutations and non-AMR elements. Tests deliberately corrupt copies to exercise
failure paths and derive a header-only no-hit case. Caller-wrapper tests mock
subprocess execution and replay that reference file; their version marker is
`TEST-DOUBLE`, never a claimed installed caller version. They test argument
construction, database/input/report checksums and failed-command propagation.

Nextflow tests invoke the actual installed runtime in isolated temporary
directories. They compile DSL2 modules, validate config and assembly inputs,
check gzip/sample IDs and errors, run the real harmonizer process, and verify
that a missing AMRFinderPlus executable fails the workflow. The optional legacy
test runs actual BLAST+ and R against a tiny synthetic no-hit input/reference.
Missing Nextflow or legacy dependencies are reported as skips.

`test_amrfinder_validation.py` builds a synthetic one-isolate package, mocked
caller provenance and AST rows. It checks exact-PDT linkage, Assembly/BioSample
and caller hashes, conflict/nonbinary observation retention, explicit no-call
status, deterministic counts and safe output publication. The genome retrieval
tests also exercise the current NCBI assembly-report schema and offline package
import, including package MD5 validation. These tests are separate from the
real-data run documented in `docs/amrfinder-validation.md`.

## nf-test

`tests/nf-test/` holds nf-test cases that run in the pinned container
(`nf-test test --profile docker`). `modules/harmonize_amrfinder.nf.test`
converts the upstream 4.0.23 report and checks a snapshot of the result.
`modules/amrfinderplus.nf.test` and the end-to-end case in `main.nf.test` run
AMRFinderPlus 4.2.7 against the pinned database on
`fixtures/blaTEM-1_NG_050145.1.fna`, the NCBI reference sequence for
*blaTEM-1* (public domain, see `fixtures/NCBI-LICENSE.txt`), and expect exactly
one `blaTEM-1` row. These cases are tagged `database`. The other pipeline
cases check `--validate_only` and a missing database argument.

`legacy/tests/test-run.R` belongs to version 1 and is not part of the
lightweight suite; see [legacy/README.md](../legacy/README.md).

## Coverage

CI uses coverage.py 7.10.6 (development dependency only):

```bash
python3 -m pip install coverage==7.10.6
python3 -m coverage run --branch --source=bin -m unittest discover -s tests -p 'test_*.py' -v
python3 -m coverage report --show-missing --fail-under=90
```

When that package is unavailable, Python's standard `trace` tool measures line
coverage (not branch coverage):

```bash
python3 -m trace --count --summary --missing --coverdir .test-results/coverage \
  --module unittest discover -s tests -p 'test_a*.py'
```

Inspect `bin.amrscan` in the output and `bin.amrscan.cover`; imported standard
library and test-module coverage is not production-code coverage. Nextflow and
legacy R are covered by executed scenarios, not by this Python percentage.

## Fixture provenance and integration limitations

`fixtures/amrfinder_v4.0.23.tsv` is an **unmodified** copy of NCBI's
[`test_dna.expected` at tag `amrfinder_v4.0.23`](https://github.com/ncbi/amr/blob/amrfinder_v4.0.23/test_dna.expected).
Upstream Git blob: `37ddb555d52dacba5ffaff07f04a9a5f30541f8e`.
The test checks conversion of this upstream expectation, not AMRFinderPlus
execution against today's database. Neither the runtime database release nor
local caller version is asserted; fixture provenance explicitly uses
`upstream_fixture` and `NA` for those fields. NCBI's public-domain notice is
retained in `fixtures/NCBI-LICENSE.txt`.

`fixtures/tiny.fna` is a short synthetic sequence for path and no-hit plumbing
tests, not a positive AMR control. The dummy database marker used by unit tests
is not a biological database. No mock calls or fixtures are accepted by the
production Nextflow caller path as successful real AMRFinderPlus executions.

A real integration run requires AMRFinderPlus, its dependencies, an immutable
indexed database release, and caller controls. Real calls, harmonization and
the AST-linked summary for a five-isolate cohort are documented in
`docs/amrfinder-validation.md`; they run outside this test suite.
The validation CLI audits evidence and does not substitute for the scientific
study design.




## AST fixtures and scenarios

`fixtures/ast.synthetic.tsv` and `fixtures/isolates.synthetic.tsv` each contain
one invented record. Accessions and phenotype/MIC combinations are test values,
not clinical observations or breakpoint examples. Their headers match the
explicit NCBI export projections inspected on 2026-10-10. Tests derive tiny
mutations and combinations for missing values, S/I/R/SDD/NS/HLAR/ND categories,
unknown names, MIC syntax, partial dates, country/source rules, duplicates,
conflicts, strict schema failures, provenance, transport failure and rebuilding.
Transport tests mock HTTP bytes and reported counts; no unit test uses a network.

The real, narrowly selected source check is recorded in
[AST validation](../docs/ast-validation.md). Downloaded tables remain outside Git.
