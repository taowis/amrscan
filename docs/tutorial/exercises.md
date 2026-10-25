# 6. Exercises

## 1. A negative result

Run the pipeline on `tests/fixtures/tiny.fna`, a short synthetic sequence.

```bash
nextflow run main.nf -profile docker \
  --input tests/fixtures/tiny.fna \
  --amrfinder_db "$PWD/db/amrfinderplus/2026-08-07.1" --outdir results/ex1
```

What does the harmonized file contain? What can you say about this sample?

??? note "Answer"
    Only the header. AMRFinderPlus found no elements. That says nothing about
    susceptibility: the sequence is 49 bases of synthetic DNA. A real "no hit"
    also only means nothing in the database matched under the thresholds used.

## 2. Same input, different database

Fetch a second release, run the positive control against it and compare the
`database_id` and the reported rows.

```bash
scripts/fetch_amrfinder_db.sh --engine docker 2026-05-15.1
```

??? note "Answer"
    The `database_id` changes even if the call does not. For a paper, report
    the software version and database release, and keep the provenance JSON
    with your results.

## 3. Organism option

Download an *E. coli* assembly from NCBI (for example GCF_037966445.1), run it
with and without `--amrfinder_organism Escherichia`, and compare
`evidence_type` across the two runs.

??? note "Answer"
    With the organism set, AMRFinderPlus also screens curated point
    mutations, which appear with `POINTX` or `POINTN` in `evidence_type` when
    present. Without it the run is a gene screen only, so a resistance
    mutation in a chromosomal gene would be missed.

## 4. Planning a large run

You have 5,000 assemblies to screen on Gadi. Each AMRFinderPlus task takes
about a minute. What would you change?

??? note "Answer"
    One PBS job per assembly means 5,000 jobs that each queue for longer than
    they run. Options: lower `--threads` to 1 or 2 so more jobs fit, raise the
    Nextflow `executor.queueSize` only if your project allows it, or run the
    head job inside a single large PBS job with the `local` executor so tasks
    share one node. Always try 20 assemblies first and look at the Nextflow
    report (`-with-report`) before scaling up.
