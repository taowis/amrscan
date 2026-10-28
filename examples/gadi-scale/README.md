# 1,000 E. coli assemblies on Gadi

Scripts and records for the run described in
[docs/gadi-scale-run.md](../../docs/gadi-scale-run.md).

| File | What it is |
| --- | --- |
| `ecoli_1000.tsv` | `bin/genome_dataset.py` manifest: 1,000 NCBI assemblies with PDT and BioSample accessions |
| `ecoli_1000.sha256` | expected SHA-256 of each downloaded FASTA |
| `fetch_genomes.pbs` | `copyq` job: download in chunks of 100, write and check `samplesheet.csv` |
| `run_amrscan.pbs` | `copyq` head job: `main.nf -profile gadi` on the samplesheet |
| `pilot/` | 20-assembly pilot: trace, versions, PBS usage |
| `full/` | 1,000-assembly run: trace, gzipped report and timeline, DAG, versions, PBS usage, `retrieval.tsv`, database hashes |

To repeat it (qn73 here; use your own project):

```bash
W=/scratch/qn73/$USER/amrscan-scale
git clone https://github.com/taowis/amrscan.git $W/amrscan && cd $W/amrscan
module load singularity
export SINGULARITY_CACHEDIR=$W/singularity_tmp SINGULARITY_TMPDIR=$W/singularity_tmp
mkdir -p $SINGULARITY_TMPDIR $W/full
scripts/fetch_amrfinder_db.sh --engine singularity 2026-08-07.1 $W/db
ln -s ../db $W/full/db
qsub -P qn73 -l storage=scratch/qn73 -v WORK=$W/full -o $W/full/fetch.log examples/gadi-scale/fetch_genomes.pbs
# when it has finished:
qsub -P qn73 -l storage=scratch/qn73 -v NCI_PROJECT=qn73,WORK=$W/full,THREADS=2 -o $W/full/head.log examples/gadi-scale/run_amrscan.pbs
```

Expect about 30 SU and a little over two hours, most of it queueing.
