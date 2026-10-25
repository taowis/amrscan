# 0006. AMRFinderPlus runs through a tested Python wrapper

- Status: accepted
- Date: 2026-10-25
- Applies to: `modules/local/amrfinderplus.nf`,
  `modules/local/harmonize_amrfinder.nf`, `bin/amrscan.py`

## Context

The Nextflow modules are short: `AMRFINDERPLUS` runs `amrscan.py call` and
`HARMONIZE_AMRFINDER` runs `amrscan.py harmonize`. The nf-core convention is
different: a module calls the tool directly in its script block and writes a
`versions.yml`.

AMRScan needs more than a version string for each sample. A result is only
useful if it can be traced to the exact assembly, database files, thresholds
and caller version, and a failed call must leave a record rather than nothing.
The same calls also need to run outside Nextflow, in batch jobs on a
workstation, and give the same output.

## Decision

- `amrscan.py call` checks the arguments, hashes the input assembly and every
  database file, reads the AMRFinderPlus version, runs AMRFinderPlus with a
  fixed argument list, hashes the report and writes `<sample>.provenance.json`.
  The provenance file is written in a `finally` block, so a failed call is
  recorded with `"status": "failed"`.
- `amrscan.py harmonize` converts the raw report to the documented schema
  ([result schema](../result-schema.md)). It refuses failed runs and reports
  whose hash does not match the provenance record.
- The Nextflow modules stay thin: they stage files, set CPUs, memory and time
  per process, and publish outputs. Scheduling, containers and executors
  stay in Nextflow.

## Evidence

- `tests/test_amrscan.py` covers the wrapper and harmonizer, including invalid
  input, a missing caller and failure propagation; CI enforces a 90% branch
  coverage floor on `bin/`.
- nf-test runs both modules and the whole pipeline in the pinned container
  against the real database, and the `test` profile finds `blaTEM-1` in
  NCBI reference sequence NG_050145.1.
- The same `call` and `harmonize` commands have been run outside Nextflow,
  with a process pool, on cohorts of thousands of assemblies.

## Consequences

- Provenance logic is written once and tested in Python, and both Nextflow
  and batch runs produce the same files.
- The AMRFinderPlus command line is inside `bin/amrscan.py`, not in the
  module, so readers who know nf-core have to look one level further.
- The modules cannot be dropped into nf-core as they are.

## Alternatives considered

- **Call `amrfinder` directly in the module and add a separate provenance
  process:** closer to nf-core and the likely next step. It splits one
  record across two tasks, so the provenance step has to trust that the files
  it hashes are the ones the caller used.
- **Use the nf-core `amrfinderplus/run` module:** records tool and database
  version strings, not input and database file hashes, and does not cover
  runs outside Nextflow.
- **Provenance in Groovy inside the workflow:** harder to test and not
  reusable outside Nextflow.
