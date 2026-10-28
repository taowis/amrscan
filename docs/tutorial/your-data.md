# 4. Your own assemblies

## Input

AMRScan takes assembled nucleotide FASTA: `.fa`, `.fna` or `.fasta`,
optionally gzipped. It does not take reads. The sample name is the file name
without those suffixes, so `isolate_07.fna.gz` becomes `isolate_07`. Names
must be unique and use only letters, digits, `_`, `.` and `-`.

```bash
nextflow run main.nf -profile docker \
  --input 'assemblies/*.fna.gz' \
  --amrfinder_db "$PWD/db/amrfinderplus/2026-08-07.1" \
  --outdir results/my_run
```

Keep the quotes around the glob so Nextflow, not the shell, expands it.

For more than a handful of files, a CSV samplesheet with columns `sample`,
`fasta` and an optional `sha256` is easier to check and to keep with the
results; see the README. `--input samples.csv` takes it in place of the glob.

Before a large run, check the inputs without running anything:

```bash
nextflow run main.nf --validate_only \
  --input 'assemblies/*.fna.gz' \
  --amrfinder_db "$PWD/db/amrfinderplus/2026-08-07.1"
```

## Organism-specific screening

If all assemblies in a run are the same supported taxon, add
`--amrfinder_organism`, for example `Escherichia` or `Klebsiella_pneumoniae`.
This turns on curated point-mutation screening and taxon-specific filtering.
Without it you get gene screening only. List supported names with:

```bash
docker run --rm community.wave.seqera.io/library/ncbi-amrfinderplus_python:d86c65c0b0a702da amrfinder -l
```

Run each taxon separately; the option applies to every sample in a run.

## Reading the results

Each row of the harmonized table is one element reported by AMRFinderPlus.
Columns are described in the [result schema](../result-schema.md). A few
things to keep in mind:

- `evidence_type` matters. `EXACTX` and `ALLELEX` are exact matches;
  `PARTIALX` and `PARTIAL_CONTIG_END` mean the hit is incomplete, often
  because the gene runs off a contig edge.
- `element_type` separates AMR from STRESS and VIRULENCE elements (the latter
  only with `--amrfinder_plus`).
- An empty table with only a header means AMRFinderPlus found nothing. It
  does not mean the isolate is susceptible.
- `drug` is always `NA`. Turning genes into drug calls needs a separate,
  explicit mapping and phenotype data.

## Resources

Each AMRFinderPlus task asks for `--threads` CPUs (default 4) and 4 GB of
memory. On a laptop, lower the threads:

```bash
nextflow run main.nf -profile docker --threads 2 ...
```
