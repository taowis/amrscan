# Shared settings for the Gadi scripts. Source it; override any value by
# exporting it first.
: "${PROJECT:?PROJECT (NCI project code) is not set}"
export AMRSCAN_ROOT="${AMRSCAN_ROOT:-/scratch/$PROJECT/$USER/amrscan}"
export AMRSCAN_REPO="${AMRSCAN_REPO:-$AMRSCAN_ROOT/amrscan-dev}"
export AMRSCAN_ENV="${AMRSCAN_ENV:-/g/data/$PROJECT/$USER/envs/amrscan-amrfinder-4.2.7}"
export AMRSCAN_BRANCH="${AMRSCAN_BRANCH:-feature/ecoli-genotype-phenotype-baseline}"

# Same relative layout as the local runs, so recorded commands match.
export PR5_DIR=datasets/ecoli_research_v1/ast_census_20261012/processed
export PR6_DIR=datasets/ecoli_research_v1/population_20261015
export GENOME_CACHE=datasets/ecoli_research_v1/genome_cache
export DB_NAME=2026-08-07.1
export DB_DIR=datasets/amrfinderplus-db/$DB_NAME
export DB_SUMS=datasets/amrfinderplus-db/$DB_NAME.SHA256SUMS
export CALLS_DIR=datasets/ecoli_research_v1/amrfinder_calls_4.2.7_db2026-08-07.1

export PR5_MANIFEST_SHA256=21348c7da115dc0306db07d32839bd91263328a2409c4a1fb65cce8b327d3caf
export DB_ID=sha256:748425a026512fda44008b9d106a4546518cc57dac59eb6432232eb3ac2c978b
# PR #4 report for PDT000041778.1 (GCA_000770275.1), AMRFinderPlus 4.2.7 on macOS ARM64.
export CONTROL_ISOLATE=PDT000041778.1
export CONTROL_ASSEMBLY=GCA_000770275.1
export CONTROL_EXPECTED=datasets/ecoli_validation/calls/PDT000041778.1.amrfinder.tsv
export CONTROL_RAW_SHA256=6be131e7b9ca35b78732e69b01e8153c62bbc0a3226f53ae77cf4a5dbe7569ba

export PATH="$AMRSCAN_ENV/bin:$PATH"
