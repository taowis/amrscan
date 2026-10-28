# 5. Running on Gadi

The `gadi` profile submits each AMRFinderPlus task as a PBS job on
`normalbw`, charges it to `$NCI_PROJECT` and runs it in the pinned image under
SingularityCE. It follows the nf-core/configs `nci_gadi` profile and has been
run on 1,000 *E. coli* assemblies ([record](../gadi-scale-run.md)).

Work in `/scratch/<project>/$USER` throughout; the commands below use qn73
as the project. The home quota is 10 GB, so
keep Nextflow's home and the image cache on scratch too (the scripts below do
this).

## One-time setup on a login node

```bash
module load singularity
W=/scratch/qn73/$USER/amrscan-run
mkdir -p $W && cd $W
git clone https://github.com/taowis/amrscan.git
export SINGULARITY_CACHEDIR=$W/singularity_tmp SINGULARITY_TMPDIR=$W/singularity_tmp
mkdir -p $SINGULARITY_TMPDIR
amrscan/scripts/fetch_amrfinder_db.sh --engine singularity 2026-08-07.1 $W/db
```

This downloads and indexes the database release (about 125 MB) with the
pinned image.

## Inputs

For more than a few files, use a samplesheet (`sample,fasta,sha256`; see
[your own assemblies](your-data.md)). To download public assemblies by
accession, `bin/genome_dataset.py download` takes a manifest and
`bin/genome_dataset.py samplesheet` writes the sheet with each file's
SHA-256, which every task then checks. `examples/gadi-scale/fetch_genomes.pbs`
does both in a `copyq` job, the queue with internet access.

## The head job

Run Nextflow itself as a small PBS job on `copyq`: it has internet access for
the first image and plugin download, a 10-hour walltime, and can submit the
task jobs. `examples/gadi-scale/run_amrscan.pbs` is a working example:

```bash
cd $W/amrscan
qsub -P qn73 -l storage=scratch/qn73 \
  -v NCI_PROJECT=qn73,WORK=$W,THREADS=2 -o $W/head.log examples/gadi-scale/run_amrscan.pbs
```

It expects `samplesheet.csv` and `db/2026-08-07.1` in `$WORK`. Give the same
project to `-P` and `NCI_PROJECT`: the first pays for the head job, the
second for the task jobs. It passes `-resume`, so
resubmitting after a failure only reruns what did not finish.

## What the profile sets

| Setting | Value |
| --- | --- |
| Executor | `pbspro` for AMRFinderPlus, up to 300 queued jobs; harmonizing runs inside the head job |
| Project | `$NCI_PROJECT`, default `qn73` |
| Storage | `$GADI_STORAGE` if set, else `scratch/<project>` |
| Queue | `normalbw` below 128 GB, `normal` up to 190 GB, `hugemembw` above |
| Image | the Docker image, converted once into `$NXF_SINGULARITY_CACHEDIR`, else `singularity/` in the launch directory |

Command-line parameters cannot set the project or storage: the profile is
read before they are applied. Export `NCI_PROJECT` or `GADI_STORAGE` instead.

## Sizing a run

Each run writes a trace to `pipeline_info/`. On the pilot of 20 *E. coli*
assemblies, an AMRFinderPlus task took about 20 s and peaked below 500 MB
while using about two of the four CPUs it asked for, so the 1,000-assembly
run used `--threads 2`. Try 20 assemblies, read the trace, then scale up.
