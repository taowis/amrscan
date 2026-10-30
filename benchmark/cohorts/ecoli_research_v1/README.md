# E. coli research cohort configuration

`cohort_config.json` is version 1.0.0 of the phenotype-only E. coli cohort
rules. It requires a versioned Pathogen Detection accession, identity-consistent
BioSample and GenBank assembly links for overall eligibility, and one
unambiguous R or S endpoint per recognized antimicrobial. Repeated consistent
AST observations are collapsed only for the derived endpoint; every observation
is retained. Discordance and nonbinary values remain explicit and are excluded
from the binary denominator.

Study-drug panel thresholds are intentionally unset pending scientific review of
the population census. See [the cohort methods and current census status](../../../docs/ecoli-research-cohort.md).

`population_config.json` is version 1.1.0 of the population-structure rules used
by `bin/population_structure.py`: genome QC bounds, the MLST scheme, the
phylogroup method, the lineage definition, holdout eligibility, the holdout
screening rules, the review-flag thresholds and, in `evaluation_design`, the
evaluation design decided in
[decision 0001](../../../docs/decisions/0001-ecoli-evaluation-design.md). None of
these values were chosen with reference to genotype–phenotype or model
results. See [the population structure methods and census](../../../docs/ecoli-population-structure.md).

`baseline_config.json` is version 1.1.0 of the genotype–phenotype baseline
settings read by `bin/genotype_baseline.py`: the mapping version and files, the
expected PR #5 and PR #6 manifest hashes, the AMRFinderPlus version and database
fingerprint, the minimum support for reporting stratum metrics (10 R and 10 S,
marked for review), the temporal periods from decision 0001 and the thresholds
for flagging mapping rules in the audit. The mapping itself is in
`benchmark/mappings/` and is described in
[decision 0002](../../../docs/decisions/0002-ecoli-determinant-mapping-v1.md).
`baseline_config_mapping_v2.json` is the same configuration pointing at the
exploratory mapping v2 of
[decision 0003](../../../docs/decisions/0003-ecoli-determinant-mapping-v2.md).

`phenotype_corrections_v1.tsv` lists the phenotype corrections read by both
configurations: R and S are swapped for every binary endpoint of PRJNA1322038
([decision 0004](../../../docs/decisions/0004-ecoli-phenotype-label-correction.md)).
The PR #5 build itself is not changed.

`genotype_baseline_20261024/` is the primary baseline, built with the
correction. `genotype_baseline_20261018/` is the first build, with the PR #5
labels as they are, kept as the record of the uncorrected labels. Each holds
the small tables of the full-cohort baseline (metrics, denominators, summary,
manifest, checksums, mapping audit, unmapped determinants, and country,
temporal and source strata), a run record (`call_run.json` for the calls and
the first build, `build_run.json` for the corrected build), the post-build
checks from
`scripts/ecoli_research_v1/baseline_qc.py` in `qc/`, and the v2 rebuild in
`mapping_v2/`. The endpoint, lineage and error-distribution tables and the
discordance lists are larger and stay in `datasets/`; they are rebuilt from
the calls with the command in `call_run.json`.

## Assembly list

`genome_manifest.tsv` lists the 10,521 GenBank assemblies used by the cohort, one
row per isolate: Pathogen Detection isolate ID, versioned GCA accession,
BioSample, SHA-256 of the extracted `genome.fna` and of the downloaded NCBI
Datasets package, the package member the FASTA was taken from, record and base
counts, and the UTC retrieval time. It was built from the `provenance.json`
files in the local genome cache and matches the `genome_sha256` column of the
2026-10-15 population build. `genome_fetch_log.jsonl` is the cache's fetch log
(pilot of 20, then the main run of 10,501, 2026-10-11 19:07 to 21:14 UTC).

The genome cache itself (about 51 GB, `datasets/ecoli_research_v1/genome_cache/`)
is not in git. Each assembly can be fetched again from

    https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<GCA>/download?include_annotation_type=GENOME_FASTA,SEQUENCE_REPORT&filename=assembly.zip

with `bin/population_structure.py fetch-genomes` (see
[the population structure guide](../../../docs/ecoli-population-structure.md)).
After a re-download, compare each `genome.fna` against `genome_sha256` here; the
accessions are versioned, so the bytes should be identical.

## Generalisation benchmark

`generalisation_v1.json` is the frozen configuration of the generalisation
benchmark
([decision 0007](../../../docs/decisions/0007-ecoli-generalisation-benchmark.md),
[guide](../../../docs/ecoli-generalisation-benchmark.md)). Its `decision_record`
field still names the record by its first number, 0005; the file is hashed into
the benchmark manifest and was not edited after freezing.

`generalisation_v1/` holds the summary tables of the 2026-10-24 run: metrics,
drops, leakage audit, determinant comparison, figure data, the benchmark and
feature manifests, `run_record.json` (commands, timings, deviations) and
`SHA256SUMS`. The predictions (about 13 MB), split lists, task checkpoints,
sketches and feature matrix (about 1 GB) stay in
`datasets/ecoli_research_v1/generalisation_v1/` and are rebuilt from the genome
cache with the commands in `run_record.json`.

## Error decomposition

`error_decomposition_v1.json` is the frozen configuration of the error
decomposition
([decision 0008](../../../docs/decisions/0008-ecoli-error-decomposition.md),
[guide](../../../docs/ecoli-domain-shift-error-decomposition.md)). It pins the
generalisation_v1 outputs, the AMRFinderPlus 2026-08-07.1 database files and
the mechanism catalogue by hash.

`error_decomposition_v1/` holds the tables of the 2026-10-25 run:
`error_decomposition_manifest.json` (drug groups and case studies, written
before any mechanism table), representation audit and capture tables, error
decomposition, mechanism overlap and shift, QC audits, the mutation
augmentation and gap-rescue tables, `run_record.json` and `SHA256SUMS`. The
full context table (about 76 MB), the call table, the audit checkpoint and the
augmented predictions stay in `datasets/ecoli_research_v1/error_decomposition_v1/`
and are rebuilt with the commands in `run_record.json` from the PR #14 output,
the genome cache and the AMRFinderPlus calls in
`datasets/ecoli_research_v1/amrfinder_calls_4.2.7_db2026-08-07.1/` (written by
the `call` command of `bin/genotype_baseline.py`, PR #7). The AMRFinderPlus database is fetched with
`scripts/fetch_amrfinder_db.sh` (version 2026-08-07.1) and checked against the
file hashes in the config.
