#!/usr/bin/env bash
# Download a fixed AMRFinderPlus database release from NCBI FTP, record file
# hashes and index it, so every run points at the same immutable directory.
# Hashes and provenance are written next to the release directory, not inside
# it, because AMRScan fingerprints every file in the database directory.
#
# usage: scripts/fetch_amrfinder_db.sh [--engine local|docker|apptainer] [VERSION] [OUTROOT]
#   VERSION  database release, default 2026-08-07.1 (format 4.2, AMRFinderPlus 4.2.x)
#   OUTROOT  parent directory, default db/amrfinderplus; the release lands in OUTROOT/VERSION
#   --engine how to run amrfinder_index: on PATH (local, e.g. the Conda environment),
#            or inside the pipeline container with Docker or Apptainer
set -euo pipefail

engine=local
if [[ ${1:-} == --engine ]]; then engine=$2; shift 2; fi
version=${1:-2026-08-07.1}
outroot=${2:-db/amrfinderplus}
format=4.2
image=community.wave.seqera.io/library/ncbi-amrfinderplus_python:d86c65c0b0a702da
base="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/${format}/${version}"

dest="${outroot}/${version}"
sums="${outroot}/${version}.SHA256SUMS"
if [[ -f $sums && -f "${outroot}/${version}.indexed" ]]; then
    echo "Database ${version} already present in ${dest}; checking hashes"
    (cd "$dest" && shasum -a 256 -c --quiet "$(cd .. && pwd)/${version}.SHA256SUMS")
    exit 0
fi

mkdir -p "$dest"
files=$(curl -fsSL "${base}/" | grep -o 'href="[^"/?]*"' | sed 's/href="//; s/"$//' | grep -v '^index')
[[ -n $files ]] || { echo "No files listed at ${base}/" >&2; exit 1; }
for f in $files; do
    curl -fsSL --retry 3 -o "${dest}/${f}" "${base}/${f}"
done
(cd "$dest" && shasum -a 256 $files) > "$sums"

abs=$(cd "$dest" && pwd)
case $engine in
    local)     amrfinder_index "$abs" ;;
    docker)    docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$abs:$abs" "$image" amrfinder_index "$abs" ;;
    apptainer) apptainer exec -B "$abs" "docker://${image}" amrfinder_index "$abs" ;;
    *)         echo "Unknown engine: $engine" >&2; exit 2 ;;
esac
touch "${outroot}/${version}.indexed"
{
    echo "source: ${base}/"
    echo "downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "indexed_with: ${engine}"
} > "${outroot}/${version}.PROVENANCE.txt"
echo "AMRFinderPlus database ${version} ready in ${dest}"
