#!/bin/bash
# Run locally from anywhere: bundle the small inputs Gadi needs (PR #5 cohort
# tables, PR #6 build, AMRFinderPlus database, PR #4 control report) into one
# tar file with the repository-relative paths kept. Genomes are not included.
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT="${1:-datasets/gadi_inputs_$(date -u +%Y%m%d).tar}"
P5=datasets/ecoli_research_v1/ast_census_20261012/processed
# No macOS metadata files: they would change the database fingerprint on Gadi.
COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs -cf "$OUT" \
  "$P5/SHA256SUMS" "$P5/candidate_isolates.tsv" "$P5/isolate_drug_endpoints.tsv" "$P5/cohort_manifest.json" \
  datasets/ecoli_research_v1/population_20261015 \
  datasets/amrfinderplus-db/2026-08-07.1 datasets/amrfinderplus-db/2026-08-07.1.SHA256SUMS \
  datasets/amrfinderplus-db/PROVENANCE.txt \
  datasets/ecoli_validation/calls/PDT000041778.1.amrfinder.tsv
shasum -a 256 "$OUT" | tee "$OUT.sha256"
ls -lh "$OUT"
