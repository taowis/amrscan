# AMRScan version 1

This directory holds version 1 of AMRScan, as described in the 2025 preprint
([arXiv:2507.08062](https://arxiv.org/abs/2507.08062)): `workflow/AMRScan.nf`
and its config, `scripts/AMRScan.R`, the R Markdown demonstrations and their
rendered reports in `docs/`, the R package `DESCRIPTION`, and the preprint
source in `paper/`. The code is kept as it was so the preprint can be
reproduced. New work should use version 2, described in the
[main README](../README.md).

The example genome and the CARD protein homolog model are not stored in the
repository. Download the assembly
[GCF_037966445.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_037966445.1/)
into `data/` at the repository root and
`protein_fasta_protein_homolog_model.fasta` from the
[CARD download page](https://card.mcmaster.ca/download) into `db/`, noting the
CARD release you used.

```bash
nextflow run legacy/workflow/AMRScan.nf \
  --input data/GCF_037966445.1_ASM3796644v1_genomic.fna \
  --card_db db/protein_fasta_protein_homolog_model.fasta \
  --outdir results/legacy --threads 4
```

Version 1 needs BLAST+, R with `dplyr`/`magrittr`, and for FASTQ input
`Biostrings`/`ShortRead`. Its config does not give a complete container
environment. The highest-bitscore summary is a display ranking, not a phenotype
call; inspect the detailed output. The HTML report process is disabled, as in
the original workflow.

The Nextflow FASTA no-hit path has an optional smoke test in
`tests/test_nextflow.py` (skipped unless BLAST+ and R are installed). The
standalone R script and `tests/test-run.R` are demonstrations: the script
changes to a fixed local directory and expects the historical example inputs,
so it is not used as an automated test.
