# AMRScan

AMRScan runs AMRFinderPlus on assembled bacterial genomes with Nextflow and
keeps a provenance record for every call: the software version, the database
release and a hash of every database file, the input hash and the exact
command. The raw caller report is kept untouched next to a harmonized table
with fixed columns.

AMRScan reports resistance determinants. It does not predict phenotype: a gene
hit is not proof of resistance, and no hit is not proof of susceptibility.

## Where to start

- New to the pipeline: work through the [workshop](tutorial/index.md). It takes
  about two hours and needs only a laptop with Docker.
- Running on NCI Gadi: [Running on Gadi](tutorial/gadi.md).
- Reading the output: [result schema](result-schema.md).
- Why the software environment is set up this way:
  [decision 0005](decisions/0005-containers.md).

## Quick start

```bash
git clone https://github.com/taowis/amrscan.git
cd amrscan
scripts/fetch_amrfinder_db.sh --engine docker
nextflow run main.nf -profile test,docker
```

Version 1, the BLAST and CARD toolkit described in the
[preprint](https://arxiv.org/abs/2507.08062), is kept under
[`legacy/`](https://github.com/taowis/amrscan/tree/main/legacy) in the
repository.
