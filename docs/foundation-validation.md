# Foundation validation — 2026-10-10

## Scope and review

The v2 entry point adds modular AMRFinderPlus calling, raw evidence/provenance,
and schema 1.0 harmonization. Existing legacy algorithms and historical data are
retained. Legacy source changes are explanatory comments and final newlines only.
Review covered input collisions, shell quoting, missing and malformed results,
provenance mismatch, partial/mutation/non-AMR evidence, failure propagation and
separation of genotypes from phenotypes. No further blocking issue was found
within this foundation scope. Biological integration validation remains pending.

The release-candidate review found that unterminated quotes and unexpected text
after a closing quote could be accepted silently. Strict TSV parsing now rejects
both cases, returns a concise CLI error and writes no partial result. The new
regression test failed for both cases before the fix and passes afterward; valid
quoted evidence still round-trips without changes.

## Tests executed locally

Environment: macOS; Nextflow 25.04.6; Python 3.14.6; R 4.5.0; BLAST+ 2.16.0.

| Group | Command | Fixtures and assertions | Coverage | Result |
| --- | --- | --- | --- | --- |
| Unit/regression | `python3 -m unittest discover -s tests -p test_amrscan.py -v` | 24 NCBI upstream records; corrupted/header-only variants; subprocess mocks; complete evidence retention, no allele/phenotype inference, checksums, invalid input and failure propagation | See line coverage below | PASS: 14 tests |
| Workflow | `python3 -m unittest discover -s tests -p test_nextflow.py -v` | Tiny synthetic FASTA/FASTQ, gzip and duplicate-ID files; upstream report; dummy path-only database; actual DSL2 compilation, config, validation, harmonization process and missing-tool failure | Scenario coverage; no numeric Nextflow metric | PASS: 7 tests, including actual legacy BLAST/R no-hit smoke |
| Full lightweight suite | `python3 -m unittest discover -s tests -p 'test_*.py' -v` | Both groups above | 21/21 scenarios; not a line metric | PASS: 21 tests in 32.422 s; no skips or flaky failures |
| Python line coverage | `python3 -m trace --count --summary --missing --coverdir .test-results/coverage --module unittest discover -s tests -p test_amrscan.py` | Unit/regression suite above | `bin/amrscan.py`: 162/163 executable lines, 99.4%; entry-point invocation line unmeasured by in-process tests | PASS |
| R syntax | `Rscript -e 'invisible(parse(file="scripts/AMRScan.R")); invisible(parse(file="tests/test-run.R"))'` | Historical R source and test script | Syntax only, not functional coverage | PASS |
| Upstream fixture identity | `git hash-object tests/fixtures/amrfinder_v4.0.23.tsv` | Exact upstream report bytes | Git blob `37ddb555d52dacba5ffaff07f04a9a5f30541f8e` | PASS |

Nextflow config flattening and `--help` compilation also passed directly and are
repeated inside the workflow suite. Coverage.py is not installed locally, so
branch coverage is not claimed here. CI is configured to run coverage.py 7.10.6
with a 90% combined statement/branch threshold; its outcome is recorded on the PR.

## Limitations and next PR

- AMRFinderPlus and an indexed AMRFinderPlus database are absent locally. No real
  caller integration or biological accuracy benchmark was executed.
- The upstream TSV is a published expectation, not output from this environment.
  Mock subprocess checks validate wrapper behavior only.
- The legacy standalone R demonstration was parsed, not executed: its hard-coded
  working directory/output defects remain. Legacy FASTQ/report paths require
  broader compatibility testing.
- Nucleotide-only calling, one optional taxon per run, local execution, and
  per-sample database hashing define this initial interface. Containers, combined
  protein/GFF input and scalable manifest caching remain future work.
- Historical large tracked artifacts and baseline license merge markers are
  listed in `tracked-data-audit.md`; no history cleanup was performed.
- PR #2: pin a reproducible caller/database environment, validate known positive
  and negative controls and multi-sample execution, then add AST ingestion with
  explicit assay/breakpoint provenance. ML remains a separate stage.

## Files added

- `.github/workflows/ci.yml`
- `benchmark/README.md`
- `bin/amrscan.py`
- `conf/base.config`
- `docs/foundation-validation.md`
- `docs/result-schema.md`
- `docs/tracked-data-audit.md`
- `main.nf`
- `modules/local/amrfinderplus.nf`
- `modules/local/harmonize_amrfinder.nf`
- `nextflow.config`
- `schemas/ast.tsv`
- `tests/README.md`
- `tests/fixtures/NCBI-LICENSE.txt`
- `tests/fixtures/amrfinder_v4.0.23.tsv`
- `tests/fixtures/tiny.fna`
- `tests/test_amrscan.py`
- `tests/test_nextflow.py`
- `tests/workflows/harmonize.nf`
- `training/README.md`
- `training/ingestion/README.md`
- `workflows/amrscan.nf`

## Files modified

- `.gitignore`
- `README.md`
- `scripts/AMRScan.R`
- `workflow/AMRScan.nf`

Pre-existing local changes in `DESCRIPTION`, `LICENSE` and `paper/paper.md` are
excluded from this commit. No new genome or AMR database is included.

## Final repository tree

```text
amrscan-dev/
├── .github/
│   └── workflows/
│       └── ci.yml
├── .gitignore
├── AMRScan_Nextflow.Rmd
├── AMRScan_R.Rmd
├── DESCRIPTION
├── LICENSE
├── README.md
├── benchmark/
│   └── README.md
├── bin/
│   └── amrscan.py
├── conf/
│   └── base.config
├── data/
│   └── GCF_037966445.1_ASM3796644v1_genomic.fna
├── db/
│   └── protein_fasta_protein_homolog_model.fasta
├── docs/
│   ├── AMRScan_Nextflow.html
│   ├── AMRScan_R.html
│   ├── foundation-validation.md
│   ├── index.html
│   ├── result-schema.md
│   └── tracked-data-audit.md
├── main.nf
├── modules/
│   └── local/
│       ├── amrfinderplus.nf
│       └── harmonize_amrfinder.nf
├── nextflow.config
├── paper/
│   ├── paper.bib
│   ├── paper.md
│   ├── paper_for_arXiv.Rmd
│   ├── paper_for_arXiv.pdf
│   └── paper_for_arXiv.tex
├── results/
│   ├── blast/
│   │   └── GCF_037966445.1_ASM3796644v1_genomic_blast_results.tsv
│   ├── final/
│   │   ├── GCF_037966445.1_ASM3796644v1_genomic_AMR_hits_detailed.csv
│   │   └── GCF_037966445.1_ASM3796644v1_genomic_AMR_hits_summary.csv
│   └── preprocessed/
│       ├── GCF_037966445.1_ASM3796644v1_genomic.fasta
│       └── GCF_037966445.1_ASM3796644v1_genomic_preprocess.log
├── schemas/
│   └── ast.tsv
├── scripts/
│   └── AMRScan.R
├── tests/
│   ├── README.md
│   ├── fixtures/
│   │   ├── NCBI-LICENSE.txt
│   │   ├── amrfinder_v4.0.23.tsv
│   │   └── tiny.fna
│   ├── test-run.R
│   ├── test_amrscan.py
│   ├── test_nextflow.py
│   └── workflows/
│       └── harmonize.nf
├── training/
│   ├── README.md
│   └── ingestion/
│       └── README.md
├── workflow/
│   ├── .nextflow.log
│   ├── AMRScan.nf
│   └── nextflow.config
└── workflows/
    └── amrscan.nf
```
