#!/bin/bash
# Run once on a Gadi login node (needs internet), from the repository root,
# after the inputs have been copied in: creates the pinned AMRFinderPlus
# environment and checks every transferred file against its checksum.
set -euo pipefail
cd "$(dirname "$0")/../.."
export AMRSCAN_REPO="$PWD"
source scripts/gadi/env.sh

if [ ! -x "$AMRSCAN_ENV/bin/amrfinder" ]; then
  command -v conda >/dev/null || { echo "conda not found: install Miniforge first (see README)"; exit 1; }
  conda env create -p "$AMRSCAN_ENV" -f environment/amrfinderplus.yml
fi
conda list --explicit -p "$AMRSCAN_ENV" > "$AMRSCAN_ROOT/env-4.2.7-linux-64.explicit.txt"
test "$(amrfinder --version)" = "4.2.7" || { echo "unexpected AMRFinderPlus version: $(amrfinder --version)"; exit 1; }
python3 --version

echo "checking transferred inputs"
(cd "$PR5_DIR" && grep -E ' (candidate_isolates.tsv|isolate_drug_endpoints.tsv|cohort_manifest.json)$' SHA256SUMS \
  | sha256sum -c -)
test "$(sha256sum "$PR5_DIR/cohort_manifest.json" | cut -d' ' -f1)" = "$PR5_MANIFEST_SHA256"
(cd "$PR6_DIR" && sha256sum -c --quiet SHA256SUMS) && echo "PR #6 build: all files match"
(cd "$DB_DIR" && sha256sum -c --quiet "../$DB_NAME.SHA256SUMS") && echo "database: all files match"
test "$(sha256sum "$CONTROL_EXPECTED" | cut -d' ' -f1)" = "$CONTROL_RAW_SHA256"
amrfinder --database "$DB_DIR" --database_version
echo "setup complete"
