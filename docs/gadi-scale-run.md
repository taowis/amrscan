# 1,000 E. coli assemblies on NCI Gadi

A record of running the workflow with the `gadi` profile on 1,000 public
*Escherichia coli* assemblies, on 2026-10-28. It shows the profile working at
scale, what it cost, and that the results match an independent run of the same
files. It is a software run; no phenotype is used. Design: [decision
0008](decisions/0008-gadi-runs.md). Scripts and raw records:
[`examples/gadi-scale/`](https://github.com/taowis/amrscan/tree/main/examples/gadi-scale).

## Inputs

- **Assemblies.** `ecoli_1000.tsv`: a random 1,000 (seed 20261028) of 10,521
  NCBI Pathogen Detection *E. coli* isolates that have susceptibility data and
  a GenBank assembly. `fetch_genomes.pbs` downloaded them from the NCBI
  Datasets v2 API in a `copyq` job (48 min, 1.60 SU), checked each package's
  MD5 list, accession, BioSample and organism, and wrote `samplesheet.csv`.
  All 1,000 FASTA files have the SHA-256 listed in `ecoli_1000.sha256`, which
  was recorded when the same assemblies were first downloaded on 2026-10-11.
  `full/retrieval.tsv` lists the URL, time and hashes of every package. 5.1
  Gb of sequence in total.
- **Database.** AMRFinderPlus database 2026-08-07.1, downloaded from NCBI FTP
  and indexed on a login node with `scripts/fetch_amrfinder_db.sh --engine
  singularity` (`full/db_2026-08-07.1.SHA256SUMS`).
- **Software.** AMRFinderPlus 4.2.7 and Python 3.12.15 in
  `community.wave.seqera.io/library/ncbi-amrfinderplus_python:d86c65c0b0a702da`
  (digest `sha256:32374e0d1cdd…c38c`), converted by SingularityCE 3.11.3 to a
  SIF of 491 MB; Nextflow 25.04.6 (Gadi module), nf-schema 2.5.1.

## Command

```bash
qsub -P qn73 -l storage=scratch/qn73 \
  -v NCI_PROJECT=qn73,WORK=/scratch/qn73/$USER/amrscan-scale/full,THREADS=2 \
  examples/gadi-scale/run_amrscan.pbs
```

which runs `main.nf -profile gadi --input samplesheet.csv --amrfinder_organism
Escherichia --threads 2`. Two CPUs per task came from the pilot below.

## Pilot

20 assemblies with `--threads 4`: 20 of 20 succeeded, tasks took 14 to 21 s
and peaked at about 465 MB, and used 200 to 230% CPU, so half the CPUs asked
for were idle. Task jobs 0.65 SU, head job 0.74 SU (11 min, most of it the
first image conversion). The pilot also found the three profile faults listed
in decision 0008. Its `software_versions.yml` predates the fix that writes one
entry per process.

## Result

| | |
| --- | --- |
| Tasks | 2,000 succeeded: 1,000 AMRFinderPlus PBS jobs on `normalbw`, 1,000 harmonizations in the head job; no retries, no failures |
| Wall time | 1 h 10 min (head job, `copyq`) |
| AMRFinderPlus task time | median 24 s, 95th percentile 37 s, maximum 5 min 35 s |
| Peak memory per task | median 268 MB, maximum 304 MB (4 GB requested) |
| CPU use per task | median 154% of 2 CPUs |
| Queue wait per task | median 4 min 10 s, maximum 7 min |
| Service units | task jobs 20.69 SU (998 of 1,000 job logs report it); head job 4.68 SU; download 1.60 SU |
| Output | 1,000 harmonized tables, 9,441 rows; 187 assemblies with no reported element |

Each task spent longer queueing than running, which is the cost of one PBS
job per assembly; it is time, not SU. A smaller memory request (1 GB would
cover the peak) does not change the SU on `normalbw`, where CPUs are charged
first.

## Agreement with an earlier run

The same 1,000 assemblies had been run on 2026-10-17 on an Apple silicon
laptop with a native osx-arm64 Conda build of AMRFinderPlus 4.2.7 and the same
database release, outside Nextflow, through the same `amrscan.py call`
wrapper and settings. For all 1,000 the input SHA-256 is the same and the
AMRFinderPlus reports are byte-identical.

The database fingerprints differ (`sha256:017e00ac…` on Gadi,
`sha256:748425a0…` on the laptop). The release directory on NCBI FTP now
holds five files that were not there on the first download
(`ReferenceGeneCatalog.txt`, `ReferenceGeneHierarchy.txt`, `amr_targets.fa`,
`changelog.txt`, `mapgenelist.txt`), and the BLAST and HMMER index files
differ between the two indexing runs. Every file present in both downloads has
the same hash. The fingerprint does what it is meant to do: it tells the two
directories apart, and the reports show the difference did not change a call.

## Files

`examples/gadi-scale/`: `ecoli_1000.tsv`, `ecoli_1000.sha256`,
`fetch_genomes.pbs`, `run_amrscan.pbs`; `pilot/` and `full/` with the
Nextflow trace, versions and PBS usage, and for the full run the gzipped
report and timeline, the DAG, `retrieval.tsv` and the database hashes. The
assemblies, work directories and per-sample results were not kept in the
repository; the manifest and hashes are enough to fetch and check the same
files again.
