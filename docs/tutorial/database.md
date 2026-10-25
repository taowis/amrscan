# 2. Get the database

AMRFinderPlus needs its own database, and results depend on which release you
use. AMRScan never downloads or updates a database during a run. You fetch one
release once, keep it unchanged, and point every run at it.

```bash
scripts/fetch_amrfinder_db.sh --engine docker
```

This downloads release `2026-08-07.1` (about 260 MB once indexed) into
`db/amrfinderplus/2026-08-07.1`. Use `--engine apptainer` on HPC or
`--engine local` if `amrfinder_index` is already on your `PATH` (for example
inside the Conda environment). To use a different release, give it as an
argument:

```bash
scripts/fetch_amrfinder_db.sh --engine docker 2026-05-15.1
```

Next to the release directory the script writes:

- `2026-08-07.1.SHA256SUMS`: hashes of the downloaded files;
- `2026-08-07.1.PROVENANCE.txt`: source URL and download time.

Running the script again checks the hashes instead of downloading.

!!! note
    Do not edit files inside the release directory. AMRScan fingerprints every
    file in it, and the fingerprint (`database_id`) appears in each result
    row. A changed file means a different `database_id`, which is the point.
