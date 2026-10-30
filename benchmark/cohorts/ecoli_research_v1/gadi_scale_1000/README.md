# Gadi run of 1,000 cohort assemblies (2026-10-28)

A software run, not an analysis: 1,000 of the 10,521 `ecoli_research_v1`
assemblies called with the workflow's `gadi` profile on NCI Gadi (project
qn73), to check that profile at scale and that it reproduces the cohort calls.
No phenotype is used. Scripts, inputs and run records are in
[`examples/gadi-scale/`](../../../../examples/gadi-scale/); the write-up is
[docs/gadi-scale-run.md](../../../../docs/gadi-scale-run.md); design in
[decision 0008](../../../../docs/decisions/0008-gadi-runs.md).

| File | What it is |
| --- | --- |
| `comparison_with_cohort_calls.tsv` | per assembly: input SHA-256, SHA-256 of the Gadi report and of the PR #7 cohort report, whether they are identical, harmonized row count, both database fingerprints |
| `MANIFEST.tsv` | where the uncommitted outputs, database and image are and how to get them again |
| `head_job_qstat_xf.txt` | `qstat -xf` of the Nextflow head job |

All 1,000 reports are byte-identical to the cohort calls in
`datasets/ecoli_research_v1/amrfinder_calls_4.2.7_db2026-08-07.1` (made on
2026-10-17 with the osx-arm64 Conda build), so the Gadi copy adds no new
calls to the study. The database fingerprints differ for the reasons given in
the write-up.

Not committed, kept in `datasets/gadi_scale_1000/` (55 MB): the full results
(raw reports, provenance, harmonized tables, `pipeline_info`), the pilot
results, `samplesheet.csv`, the chunked manifests, the per-chunk retrieval
`provenance.json`, the download and head-job logs and the per-task SU lines.
The 1,000 FASTA files and the Nextflow work directories were deleted from
Gadi after the run; `examples/gadi-scale/fetch_genomes.pbs` downloads the
same files again and stops unless every one matches `ecoli_1000.sha256`.
