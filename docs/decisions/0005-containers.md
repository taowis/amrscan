# 0005. Software environments: one pinned image, three engines

- Status: accepted
- Date: 2026-10-25
- Applies to: `conf/base.config`, the `docker`, `apptainer`, `singularity`,
  `conda` and `gadi` profiles in `nextflow.config`, `conf/gadi.config`,
  `environment/amrfinderplus.yml`, `scripts/fetch_amrfinder_db.sh`

## Context

Until now the v2 workflow ran with whatever `amrfinder` was on `PATH`, so a
result depended on the user's machine. Both processes call `bin/amrscan.py`,
which needs Python 3.10+ in the same environment as AMRFinderPlus. The
BioContainers image `quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0`
has AMRFinderPlus but no Python (checked 2026-10-25), so it cannot run the
wrapper.

People will run AMRScan on laptops (Docker), on HPC where only Apptainer or
Singularity is allowed (NCI Gadi), and on hosts where containers are not
available at all (Conda).

## Decision

- Build one image from the existing Conda specification
  (`ncbi-amrfinderplus=4.2.7`, `python=3.12`, conda-forge and bioconda) with
  Seqera Containers and pin it by tag:
  - Docker/OCI: `community.wave.seqera.io/library/ncbi-amrfinderplus_python:d86c65c0b0a702da`
  - Apptainer/Singularity SIF: `oras://community.wave.seqera.io/library/ncbi-amrfinderplus_python:03c3120e1c088608`
    (replaced on 2026-10-28, see below)
- The `conda` profile builds from `environment/amrfinderplus.yml`, the same
  specification, so the three engines carry the same tool versions.
- The database is not baked into the image. `scripts/fetch_amrfinder_db.sh`
  downloads one named release from NCBI FTP (default `2026-08-07.1`), writes
  SHA-256 sums and provenance next to the release directory, and indexes it
  with the pinned image. The workflow still fingerprints every database file.
- The Gadi profile follows the nf-core/configs `nci_gadi` profile (PBS Pro,
  project and storage directives from Gadi's Nextflow module, Singularity
  module, lenient cache) and uses the SIF.

Amended 2026-10-28: the first run on Gadi (SingularityCE 3.11.3) could not pull
the `oras://` tag; the registry returns a Docker manifest for it
(`could not get image manifest, received mediaType:
application/vnd.docker.distribution.manifest.v2+json`). The `apptainer`,
`singularity` and `gadi` profiles now use the Docker image above, which
Singularity and Apptainer convert to a SIF once and keep in their cache. All
engines therefore run the same image digest. See decision 0008.

## Evidence

- The Docker image reports AMRFinderPlus 4.2.7 and runs `amrfinder_index`;
  the pinned database indexed in 32 s on an Apple silicon laptop under
  emulation.
- `nextflow run main.nf -profile test,docker` detects `blaTEM-1` (100%
  coverage) in the NCBI reference sequence NG_050145.1; the same case is part
  of the nf-test suite run in CI.
- The `conda` profile was run end to end on macOS ARM64 on 2026-10-25.
- The Gadi profile has not yet been run on Gadi.

## Consequences

- Results no longer depend on what is installed on the host, and the image tag
  is recorded in the Nextflow trace.
- Updating AMRFinderPlus means a new image request, a new tag in two places in
  the configuration, and rerunning the database-tagged tests.
- Seqera Containers is a free public registry. If it became unavailable the
  image can be rebuilt from `environment/amrfinderplus.yml`.

## Alternatives considered

- **BioContainers AMRFinderPlus image:** no Python, so the wrapper cannot run.
- **Own Dockerfile pushed to GHCR:** full control, but adds a build workflow
  and a registry to maintain; GHCR packages start private and need manual
  visibility changes.
- **Database inside the image:** simpler first run, but couples the software
  and database release cycles and makes the image several hundred MB larger.
- **Wave on the fly (`wave.enabled`):** needs network access from the run
  host, which Gadi compute nodes do not have.
