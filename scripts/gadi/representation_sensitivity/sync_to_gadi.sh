#!/bin/bash
# Run on the workstation from anywhere: copy the code and the inputs the Gadi
# stages read to the repository checkout on Gadi, keeping repository-relative
# paths. GADI_DIR is the checkout on Gadi, e.g.
#   GADI_DIR=/scratch/qn73/$USER/Projects/amrscan-dev scripts/gadi/representation_sensitivity/sync_to_gadi.sh
set -euo pipefail
: "${GADI_DIR:?set GADI_DIR to the repository checkout on Gadi}"
GADI_HOST="${GADI_HOST:-gadi-dm.nci.org.au}"
cd "$(dirname "$0")/../../.."
D=datasets/ecoli_research_v1
G=$D/generalisation_v1
E=$D/error_decomposition_v1
R=$D/representation_sensitivity_v1
# code and committed records
rsync -a --relative --exclude __pycache__ bin benchmark environment scripts "$GADI_HOST:$GADI_DIR/"
# frozen inputs: PR #6 and PR #7 builds, PR #14 manifest, predictions and task records, PR #15 calls and R2 runs
rsync -a --relative \
  $D/population_20261015 $D/genotype_baseline_20261024 \
  $G/benchmark_manifest.json $G/predictions.tsv.gz $G/genome_verification.tsv $G/matrix/feature_manifest.json \
  $G/tasks $E/amr_calls.tsv.gz $E/augmented \
  $R/representation_sensitivity_manifest.json $R/genome_verification.tsv $R/sketches \
  "$GADI_HOST:$GADI_DIR/"
# the PR #14 matrix, for the R0 platform check
rsync -a --relative $G/matrix "$GADI_HOST:$GADI_DIR/"
[ -f $R/selection.json ] && rsync -a --relative $R/selection.json $R/pilot "$GADI_HOST:$GADI_DIR/"
echo "synced to $GADI_HOST:$GADI_DIR"
