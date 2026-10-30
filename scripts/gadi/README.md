# Running the PR #7 AMRFinderPlus calls on NCI Gadi

These scripts run the same calls as `bin/genotype_baseline.py call` on a
workstation, but on one Gadi `normal` node. Genomes are downloaded again from
NCBI on Gadi; the small inputs are copied from the workstation. Every genome
must match the SHA-256 in the PR #6 build before it is called, and a control
isolate must reproduce the PR #4 report byte for byte before the cohort run
starts.

Replace `ab12` with your NCI project code throughout.

## 1. Workstation: bundle the inputs

```bash
scripts/gadi/pack_inputs.sh            # writes datasets/gadi_inputs_<date>.tar and .sha256
rsync -av datasets/gadi_inputs_*.tar* $USER@gadi-dm.nci.org.au:/scratch/ab12/$USER/amrscan/
```

The bundle (about 390 MB) holds the three PR #5 cohort files and their
`SHA256SUMS`, the PR #6 build `population_20261015`, the AMRFinderPlus database
2026-08-07.1 with its checksum list, and the PR #4 report for the control
isolate. It is written without macOS metadata files, which would otherwise
change the database fingerprint.

## 2. Gadi login node: code, environment, checks

```bash
mkdir -p /scratch/ab12/$USER/amrscan && cd /scratch/ab12/$USER/amrscan
sha256sum -c gadi_inputs_*.tar.sha256
git clone --branch feature/ecoli-genotype-phenotype-baseline https://github.com/taowis/amrscan-dev.git
cd amrscan-dev
tar -xf ../gadi_inputs_*.tar
export PROJECT=ab12
scripts/gadi/setup_login.sh
```

`setup_login.sh` needs `conda`. If you do not have it, install Miniforge
(conda-forge only, no extra channel terms) under `/g/data`:

```bash
curl -fsSLO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /g/data/ab12/$USER/miniforge3
source /g/data/ab12/$USER/miniforge3/etc/profile.d/conda.sh
```

The environment goes to `/g/data/ab12/$USER/envs/amrscan-amrfinder-4.2.7`
(from `environment/amrfinderplus.yml`, AMRFinderPlus 4.2.7). Its explicit
package list is written to `/scratch/ab12/$USER/amrscan/env-4.2.7-linux-64.explicit.txt`.
The script stops if the AMRFinderPlus version is not 4.2.7 or any transferred
file fails its checksum.

## 3. Download genomes (copyq, about 2-3 hours)

```bash
cd /scratch/ab12/$USER/amrscan/amrscan-dev
qsub -P ab12 -l storage=scratch/ab12+gdata/ab12 scripts/gadi/fetch_genomes.pbs
```

At most three concurrent NCBI requests, as in PR #6. The cache is resumable:
if the job stops, submit it again. The job log ends with the number of cached
assemblies (10,521 expected); failures are listed in `genome_cache/fetch_log.jsonl`.

## 4. AMRFinderPlus calls (normal, 48 cores, about 4-6 hours)

```bash
qsub -P ab12 -l storage=scratch/ab12+gdata/ab12 scripts/gadi/amrfinder_calls.pbs
```

The job copies the database to node-local `$PBS_JOBFS`, rechecks it, calls the
control isolate and compares its raw report with the PR #4 one. A mismatch is
recorded in `platform_control.jsonl` and stops the job (exit code 3) before any
cohort call. Otherwise it runs the cohort with 48 workers and prints the status
counts from `call_log.jsonl`. If the walltime runs out, submit again: valid
calls are reused, anything incomplete is moved to `.quarantine/` and redone.

## 5. Baseline (normal, a few minutes)

```bash
qsub -P ab12 -l storage=scratch/ab12+gdata/ab12 -v BASELINE_TAG=gadi scripts/gadi/build_baseline.pbs
```

## 6. Copy results back

```bash
rsync -av $USER@gadi-dm.nci.org.au:/scratch/ab12/$USER/amrscan/amrscan-dev/datasets/ecoli_research_v1/amrfinder_calls_4.2.7_db2026-08-07.1 \
  datasets/ecoli_research_v1/amrfinder_calls_4.2.7_db2026-08-07.1_gadi
rsync -av $USER@gadi-dm.nci.org.au:/scratch/ab12/$USER/amrscan/amrscan-dev/datasets/ecoli_research_v1/genotype_baseline_gadi \
  datasets/ecoli_research_v1/
rsync -av $USER@gadi-dm.nci.org.au:/scratch/ab12/$USER/amrscan/env-4.2.7-linux-64.explicit.txt datasets/amrfinderplus-db/
```

The calls take a few hundred MB and about 42,000 files; check your `/scratch`
inode quota with `lquota` first. Files on `/scratch` that are not accessed for
100 days are removed, so copy the results back soon after the run.
