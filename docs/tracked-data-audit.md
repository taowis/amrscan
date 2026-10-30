# Historical tracked artifacts retained in this PR

The following files were already tracked before the v2 foundation work. They
remain unchanged; `.gitignore` prevents new runtime files from being added but
does not remove existing Git history. No new genome or AMR database is added.

| File/group | Approximate size | Later cleanup proposal |
| --- | --- | --- |
| `data/GCF_037966445.1_ASM3796644v1_genomic.fna` | 5.48 MB | Versioned download manifest/checksum |
| `results/preprocessed/GCF_037966445.1_ASM3796644v1_genomic.fasta` | 5.48 MB | Remove duplicate generated sequence in a separate PR |
| `db/protein_fasta_protein_homolog_model.fasta` | 2.33 MB | External database installation with version/license metadata |
| `results/blast/*.tsv`, `results/final/*.csv` | 1.44 MB total | Archived reproducible example output |
| `docs/AMRScan_R.html`, `docs/AMRScan_Nextflow.html` | 1.25 MB total | Separate rendered historical documentation |
| `paper/paper_for_arXiv.pdf` | 0.30 MB | Retain manuscript; decide archival policy separately |
| `workflow/.nextflow.log`, `results/preprocessed/*log` | Small | Stop tracking runtime logs in cleanup PR |

Do not rewrite history as part of this refactoring. `DESCRIPTION`, `LICENSE` and
manuscript metadata also require a separate review: local user edits in those
files were explicitly excluded from the foundation commit. The tracked baseline
LICENSE contains unresolved merge markers; this pre-existing issue was not
silently combined with the user's local license repair.
