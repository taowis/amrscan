# 0007. nf-core conventions for parameters, metadata, versions and resources

- Status: accepted
- Date: 2026-10-28
- Applies to: `main.nf`, `nextflow_schema.json`, `nextflow.config`,
  `conf/base.config`, `workflows/amrscan.nf`, `modules/local/*.nf`

## Context

The workflow checked its parameters with hand-written Groovy, passed a bare
sample ID string between processes, set resources per process name and wrote
no execution reports. Anyone used to nf-core pipelines has to learn a
different layout, and a large run on a cluster had no per-task record of
memory and time to size the next one from. Every task error stopped the run,
including a task the scheduler killed for exceeding its memory or walltime,
which on a shared cluster is an expected event, not an error in the data.

## Decision

- **Parameters.** `nextflow_schema.json` (JSON Schema draft 2020-12) describes
  every parameter. `main.nf` calls `validateParameters()` from the nf-schema
  plugin, pinned at 2.5.1, before building any channel, and `--help` is
  generated from the schema. Checks the schema cannot express stay in
  `main.nf`: identity must be -1 or in [0, 1], and the FASTA name, sample ID
  and duplicate checks.
- **Metadata.** Samples travel as `[meta, files]` with `meta.id` the sample
  ID, so later inputs (a samplesheet, per-sample organism) add keys instead
  of changing every process signature.
- **Versions.** Each process writes `versions.yml`; the workflow collects them
  into `pipeline_info/software_versions.yml`. The per-sample provenance JSON
  from decision 0006 is unchanged and remains the record to cite; the
  versions file is a run-level summary.
- **Resources.** Processes carry the labels `process_single` and
  `process_medium` instead of per-name settings. Memory and time scale with
  `task.attempt`.
- **Errors.** A task that exits with 130 to 145, 104 or 175 (killed by a
  signal, the kernel or the scheduler, the nf-core set) is retried up to
  twice with more memory and time. Every other failure, including a failed
  AMRFinderPlus call, still stops the run, so a bad input or database is never
  skipped.
- **Reports.** Every run writes the Nextflow report, timeline, trace and DAG
  to `pipeline_info/` with a timestamp.

## Evidence

- `tests/test_nextflow.py` runs the installed Nextflow offline and checks
  the schema messages for a missing input or database, a database path that
  is a file, a non-integer or zero thread count and out-of-range thresholds,
  and that the Groovy identity check still catches values between -1 and 0.
- nf-test checks the meta map through both modules, `versions.yml` with
  AMRFinderPlus 4.2.7, and that the end-to-end run writes
  `software_versions.yml` (nf-test replaces the trace with its own).
- The `test` profile run with Docker on 2026-10-28 produced the same
  harmonized result as before and the four reports.

## Consequences

- The first run downloads the nf-schema plugin. Offline machines and cluster
  compute nodes need it fetched beforehand
  (`nextflow plugin install nf-schema@2.5.1`); CI does this before the
  offline tests.
- Error messages for parameters now come from nf-schema
  (`Missing required parameter(s): input`) instead of the earlier wording.
- `--threads` still sets the CPUs of AMRFinderPlus tasks, through the
  `process_medium` label.

## Alternatives considered

- **Keep the Groovy checks:** no dependency, but the checks, the help text
  and the README table would keep drifting apart.
- **Start from the nf-core pipeline template:** brings MultiQC, iGenomes and
  a large utility layer the workflow does not use, and nf-core lint rules
  that assume a pipeline in the nf-core organisation.
- **Retry every failure:** hides genuine caller and input errors behind
  repeated attempts.
