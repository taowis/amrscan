# 0008. Running on NCI Gadi: samplesheet input, one PBS job per caller task

- Status: accepted
- Date: 2026-10-28
- Applies to: `conf/gadi.config`, `nextflow.config` (Singularity profiles),
  `main.nf`, `assets/schema_input.json`, `bin/amrscan.py call
  --expected-sha256`, `bin/genome_dataset.py samplesheet`,
  `examples/gadi-scale/`

## Context

The `gadi` profile had been written from the nf-core/configs `nci_gadi`
profile but never run. The first runs on Gadi, on 2026-10-28, found three
faults:

- it enabled Apptainer, but Gadi provides SingularityCE 3.11.3 (`module load
  singularity`) and no `apptainer` command;
- the `oras://` SIF tag could not be pulled (decision 0005, amended);
- the project and storage were read from `params.gadi_project` inside the
  profile. Command-line params are not visible while a config file is read,
  so `--gadi_project` never reached the paths built from it, and the run used
  the login default project (here one with no allocation).

Large runs also needed an input that can be checked before anything is
scheduled and tied to the exact files downloaded.

## Decision

- **Engine and image.** The profile enables Singularity and uses the Docker
  image of decision 0005, converted once into
  `/scratch/$PROJECT/$USER/nxf_singularity_cache`.
- **Project and storage** come from the environment, as in `nci_gadi`:
  `PROJECT`, and `GADI_STORAGE` if set (default
  `gdata/$PROJECT+scratch/$PROJECT`). The `gadi_project` and `gadi_storage`
  parameters are removed. `examples/gadi-scale/run_amrscan.pbs` sets
  `PROJECT` to the project the head job was submitted under.
- **Executors.** Each AMRFinderPlus task is a PBS job on `normalbw`.
  Harmonizing a report takes about a second, so it runs with the `local`
  executor inside the head job instead of as a PBS job per sample.
- **Head job.** Nextflow runs as a one-CPU PBS job on `copyq`, which has
  internet access for the first image and plugin download, a 10-hour
  walltime, and can submit jobs.
- **Samplesheet.** `--input` also takes a CSV with `sample`, `fasta` and an
  optional `sha256`, validated with nf-schema against
  `assets/schema_input.json`. The hash travels in `meta.sha256` and
  `amrscan.py call --expected-sha256` stops the task if the file differs.
  `genome_dataset.py samplesheet` writes the sheet from downloaded NCBI
  packages.

## Evidence

- Pilot on 2026-10-28, 20 *E. coli* assemblies, `--threads 4`: 20 of 20
  succeeded; AMRFinderPlus tasks used 14 to 21 s of real time, at most about
  465 MB of memory and about 200 to 230% CPU, and 0.65 SU in total. The head
  job used 0.74 SU over 11 minutes, most of it converting the image on the
  first run.
- Full run on 2026-10-28, 1,000 assemblies, `--threads 2`: 2,000 tasks
  succeeded without a retry in 1 h 10 min; AMRFinderPlus tasks took a median
  24 s and at most 304 MB; 20.69 SU for the task jobs and 4.68 SU for the head
  job. All 1,000 reports are byte-identical to an earlier run of the same
  files on a laptop. Details in `docs/gadi-scale-run.md`, records in
  `examples/gadi-scale/`.
- `tests/test_nextflow.py` checks the samplesheet schema messages and that a
  wrong `sha256` stops the task; `tests/test_amrscan.py` and
  `tests/test_genome_dataset.py` cover the wrapper check and the samplesheet
  writer.

## Consequences

- A run on another project needs `PROJECT` set (it is at login) or the head
  job submitted with `-P`; nothing in the repository names a project.
- One PBS job per assembly costs queueing time rather than SU: the pilot
  tasks ran for about 20 s each.
- The first run on a new account spends a few minutes converting the image.

## Alternatives considered

- **One large PBS job with the local executor:** fewer jobs and no queueing
  per task, but the node is charged for its full width while tasks start and
  finish, and the run is limited to one node.
- **PBS job arrays:** not available on Gadi.
- **Keep `--gadi_project` and resolve paths lazily:** the cache directory and
  directives are read when the config is parsed, so a parameter cannot set
  them without a wrapper script.
