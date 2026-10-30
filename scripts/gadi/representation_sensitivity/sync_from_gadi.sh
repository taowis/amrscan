#!/bin/bash
# Run on the workstation: copy the Gadi run outputs back. Matrices stay on Gadi.
set -euo pipefail
: "${GADI_DIR:?set GADI_DIR to the repository checkout on Gadi}"
GADI_HOST="${GADI_HOST:-gadi-dm.nci.org.au}"
cd "$(dirname "$0")/../../.."
R=datasets/ecoli_research_v1/representation_sensitivity_v1
rsync -a --relative --exclude matrix --exclude sketches \
  "$GADI_HOST:$GADI_DIR/./$R/" ./
rsync -a --relative "$GADI_HOST:$GADI_DIR/./$R/logs" ./ 2>/dev/null || true
echo "copied back into $R"
