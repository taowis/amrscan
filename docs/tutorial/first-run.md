# 3. First run

The `test` profile runs one known positive: the *blaTEM-1* gene from NCBI
reference sequence NG_050145.1. If the setup is right, AMRScan finds it.

```bash
nextflow run main.nf -profile test,docker
```

Replace `docker` with your engine from part 1. It should finish in about a
minute after the image is downloaded.

## What came out

```text
results/test/
├── raw/amrfinderplus/blaTEM-1_NG_050145.1/
│   ├── blaTEM-1_NG_050145.1.amrfinder.tsv      the AMRFinderPlus report, unchanged
│   ├── blaTEM-1_NG_050145.1.amrfinder.log      caller messages and versions
│   └── blaTEM-1_NG_050145.1.provenance.json    command, versions, hashes
└── harmonized/blaTEM-1_NG_050145.1.amr.tsv     one row per reported element
```

Look at the harmonized table:

```bash
cut -f1-5,9 results/test/harmonized/blaTEM-1_NG_050145.1.amr.tsv
```

You should see one row: gene `blaTEM-1`, class `BETA-LACTAM`, coverage
`100.00`.

Then look at the provenance file. This is what you cite in a methods section:

```bash
python3 -m json.tool results/test/raw/amrfinderplus/*/*.provenance.json | head -20
```

`caller_version` is the AMRFinderPlus version, `database_id` the database
fingerprint, and `input_sha256` identifies the exact FASTA that was analysed.

## Rerunning

Run the same command again with `-resume`. Nextflow reuses the finished tasks
and nothing is recomputed. Change a parameter and only the affected tasks run
again.
