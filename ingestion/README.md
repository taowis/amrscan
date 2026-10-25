# AST ingestion

The [NCBI Pathogen Detection adapter](ncbi_pathogen_detection/README.md) uses
`bin/ast_dataset.py` to archive public AST and isolate tables, verify raw hashes,
and build the normalized contract in `schemas/ast.tsv`. Its versioned policy
defines explicit name mappings and conservative metadata rules.

Downloads and offline transformations are separate operations. AST observations
remain independent of the genomic determinant caller. Future adapters must
preserve original source records, assay metadata and snapshot provenance, and
must document their own source contracts before producing this interface.
