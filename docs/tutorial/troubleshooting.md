# Troubleshooting

**`Provide --amrfinder_db`**: the database path is missing. Use the release
directory itself, for example `db/amrfinderplus/2026-08-07.1`, not its parent.

**`Expected assembled FASTA`**: the input is not `.fa`, `.fna` or `.fasta`
(optionally `.gz`). AMRScan does not accept reads.

**`Duplicate sample ID`**: two files give the same name after the suffix is
removed, for example `a.fa` and `a.fna`. Rename one.

**Docker: `permission denied` on the socket**: Docker is not running, or your
user is not in the `docker` group.

**Apple silicon: a warning about the platform**: expected; the image is
`linux/amd64` and runs under emulation.

**The database fingerprint changed between runs**: something in the release
directory changed. Check with
`shasum -a 256 -c db/amrfinderplus/2026-08-07.1.SHA256SUMS` from inside the
release directory. The fingerprint also covers the index files written by
`amrfinder_index`, so two separate downloads of the same release may not give
the same fingerprint; keep one copy and reuse it.

**Gadi: image pull fails inside a job**: compute nodes have no internet. Pull
the image first as described in [Running on Gadi](gadi.md).
