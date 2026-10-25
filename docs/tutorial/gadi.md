# 5. Running on Gadi

!!! warning
    The `gadi` profile follows the nf-core/configs `nci_gadi` profile but has
    not yet been run on Gadi as part of this project. Check it on a small
    input first and report problems as GitHub issues.

The profile submits each task as a PBS job, charges it to your project, and
runs it in the Apptainer/Singularity image.

## One-time setup on a login node

```bash
module load nextflow singularity
cd /scratch/$PROJECT/$USER
git clone https://github.com/taowis/amrscan.git
cd amrscan
scripts/fetch_amrfinder_db.sh --engine apptainer 2026-08-07.1 /g/data/$PROJECT/$USER/amrfinderplus
```

If `apptainer` is not on the path after loading the module, use
`singularity`; edit the last line of the `apptainer)` case in the script or
run `amrfinder_index` through `singularity exec` yourself.

Check that the Nextflow module is 25.04.6 or later (`nextflow -version`). The
profile uses the `project` and `storage` directives, which work with the
Nextflow module NCI provides.

## Pull the image before running

Compute nodes cannot reach the internet. Run the head job where there is
network access (a login node for a short run, or a `copyq` job), so the image
is pulled into the cache once:

```bash
nextflow run main.nf -profile gadi \
  --input '/g/data/'$PROJECT'/assemblies/*.fna' \
  --amrfinder_db /g/data/$PROJECT/$USER/amrfinderplus/2026-08-07.1 \
  --outdir /scratch/$PROJECT/$USER/amrscan_results
```

The image is stored in `/scratch/$PROJECT/$USER/nxf_apptainer_cache`. Later
runs reuse it.

## What the profile sets

| Setting | Value |
| --- | --- |
| Executor | `pbspro`, up to 300 queued jobs |
| Project | `$PROJECT` (override with `--gadi_project`) |
| Storage | `gdata/<project>+scratch/<project>` (override with `--gadi_storage`) |
| Queue | `normalbw` below 128 GB, `normal` up to 190 GB, `hugemembw` above |
| Image | the SIF from `oras://community.wave.seqera.io/...` |

AMRFinderPlus tasks ask for 4 CPUs, 4 GB and 2 hours by default, which fits a
single assembly with room to spare. For thousands of assemblies the per-job
overhead dominates; see the [exercises](exercises.md) for a discussion.
