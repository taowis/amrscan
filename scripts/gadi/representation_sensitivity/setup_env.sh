#!/bin/bash
# Run once on a Gadi login node (internet access) from the repository root:
# Python 3.12 virtual environment with the pinned packages of
# environment/generalisation-benchmark.txt, and a record of what was installed.
set -euo pipefail
module load python3/3.12.1
# keep pip's cache and temporary files off the home quota
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$PWD/.pip-cache}" TMPDIR="${TMPDIR_SETUP:-$PWD/.tmp}"
mkdir -p "$TMPDIR"
[ -x .venv-gadi/bin/python ] || python3 -m venv .venv-gadi
.venv-gadi/bin/pip install --quiet --upgrade pip
.venv-gadi/bin/pip install --quiet -r environment/generalisation-benchmark.txt
mkdir -p datasets/ecoli_research_v1/representation_sensitivity_v1/logs
{ .venv-gadi/bin/python --version; .venv-gadi/bin/pip freeze; } \
  > datasets/ecoli_research_v1/representation_sensitivity_v1/logs/gadi_environment.txt
.venv-gadi/bin/python -c "import sourmash, sklearn; assert sourmash.VERSION == '4.9.4' and sklearn.__version__ == '1.7.2'"
echo "environment ready"
