# AMRFinderPlus assembly validation

This validation layer keeps source AST observations separate from genomic
caller evidence. It does not infer phenotype from genes, predict R/S, apply
CLSI/EUCAST breakpoints or calculate genotype–phenotype performance. A caller
no-hit result does not establish susceptibility.

## Cohort selection and retrieval

`benchmark/cohorts/ecoli_amrfinder_validation.tsv` lists the five verified
isolates; `benchmark/cohorts/ecoli_amrfinder_exclusions.tsv` lists the two
screened candidates that were not admitted. An isolate is admitted only after
its exact versioned PDT accession, BioSample and Assembly agree across the
Pathogen Detection AST and isolate records and the NCBI Datasets assembly
report, including *Escherichia coli* and strain. AST rows are retained at
observation level, including duplicates, conflicts and nonbinary categories.
Inclusion is based on source availability and identity, never on AMRFinderPlus
calls or agreement with AST.

The BioSample used in the AST builder examples, `SAMN05170351`, had AST linked to SRA
`SRR4065706` but no Assembly, so it was not used.

### Selection (2026-10-11)

`PDT000041778.1` was carried forward from an earlier lookup. To add isolates
reproducibly, the isolates collection was exported for
`scientific_name:"Escherichia coli" AND AST_phenotypes:* AND asm_acc:GCA_*`
(21,042 rows) and archived, because the service ignores sort parameters. PDTs
were then examined in ascending numeric `target_acc` order. A candidate was
accepted when its metadata rows agreed; it had exact *E. coli*, a BioSample and
a GenBank Assembly; its BioSample, strain and NCBI SNP cluster were not already
in the cohort; and its AST rows, retrieved by exact PDT, had one row per drug,
phenotypes in R/S/I/SDD/NS and drug names known to the ingestion policy.
Screening stopped after four acceptances. No genotype field or caller output
was consulted. Six PDTs were examined:

| PDT | BioSample | Assembly | Decision |
| --- | --- | --- | --- |
| `PDT000026734.2` | `SAMN02581401` | `GCA_000633675.2` | Excluded: AST includes ND phenotypes. |
| `PDT000041780.1` | `SAMN03075589` | `GCA_000770285.1` | Excluded: same SNP cluster (`PDS000117091`) as `PDT000041778.1`. |
| `PDT000045580.1` | `SAMN03177674` | `GCA_000797605.1` | Admitted. |
| `PDT000045581.1` | `SAMN03177676` | `GCA_000797575.1` | Admitted. |
| `PDT000045582.1` | `SAMN03177659` | `GCA_000797595.1` | Admitted. |
| `PDT000045583.1` | `SAMN03177675` | `GCA_000797585.1` | Admitted. |

The AST rule prefers isolates whose observations can be processed, so the four
added isolates have cleaner AST than a random draw would. They come from one
FDA BioProject (`PRJNA266657`) and were submitted together; the cohort is not
a population sample.

### Verified cohort

| PDT | BioSample | Assembly | Strain | Submitter | Sequences | Bases | AST rows |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `PDT000041778.1` | `SAMN03075588` | `GCA_000770275.1` | MRSN17749 | WRAIR (`PRJNA261723`) | 92 | 5,046,460 | 48 |
| `PDT000045580.1` | `SAMN03177674` | `GCA_000797605.1` | CVM N36963PS | FDA (`PRJNA266657`) | 85 | 4,836,410 | 15 |
| `PDT000045581.1` | `SAMN03177676` | `GCA_000797575.1` | CVM N37067PS | FDA (`PRJNA266657`) | 104 | 4,776,649 | 15 |
| `PDT000045582.1` | `SAMN03177659` | `GCA_000797595.1` | CVM N36225PS | FDA (`PRJNA266657`) | 154 | 5,080,401 | 15 |
| `PDT000045583.1` | `SAMN03177675` | `GCA_000797585.1` | CVM N37059PS | FDA (`PRJNA266657`) | 80 | 4,920,814 | 15 |

AST and isolate metadata were exported on 2026-10-11 at 11:29 UTC with
`bin/ast_dataset.py download` and query
`target_acc:(PDT000041778.1 OR PDT000045580.1 OR PDT000045581.1 OR PDT000045582.1 OR PDT000045583.1)`:
108 AST rows (manifest SHA-256 of `asts.tsv`
`4d1401c23adabf8c93e52f87e097a32bd54a68a0988b6ff181d987b6371fa6bb`) and 10
isolate rows, two per PDT from consecutive NCBI pipeline runs that differ only
in run and cluster fields. A wide-field TSV and JSON export, count responses
and HTTP headers were saved alongside. The 48 `PDT000041778.1` rows are
identical by id and value to a separate export taken at 11:03 UTC. Assembly
packages were downloaded at 11:30 UTC with `bin/genome_dataset.py download`,
and NCBI Datasets `dataset_report` responses were saved to check strain and
taxid. All five FASTA SHA-256 values are in the genome provenance and in the
link table below.

Raw files stay outside Git under `datasets/ncbi_ast_manual/` (exports,
candidate list, screening log, README and `SHA256SUMS` per snapshot) and
`datasets/ecoli_validation/cohort5_20261011/` (genomes, assembly reports,
processed AST, summary).

The downloader retrieves only explicitly listed accessions from NCBI Datasets,
then checks the returned Assembly accession, BioSample, exact organism, FASTA
structure and SHA-256. Its provenance records the manifest hash, source URL,
retrieval time, package member and FASTA hash. It refuses to overwrite prior
outputs. Genome files and provenance are ignored by Git.

If the command-line environment cannot reach NCBI, retrieve the complete
official NCBI Datasets ZIP in a browser or another network-enabled environment.
Keep `md5sum.txt`, the Assembly report and FASTA in the package. Create a
tab-separated package manifest with columns `assembly_accession`,
`package_path`, `source_url` and timezone-qualified `retrieved_at`; paths may
be absolute or relative to that table. The offline importer checks every
member listed in `md5sum.txt`, current and legacy NCBI report identity fields,
FASTA structure, package hash and cohort mapping before publishing files.

```sh
python3 bin/genome_dataset.py import-packages \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --packages-tsv /path/to/ncbi-packages.tsv \
  --outdir datasets/ecoli_validation/genomes
python3 bin/genome_dataset.py validate \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes datasets/ecoli_validation/genomes
```

```sh
python3 bin/genome_dataset.py download \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --outdir datasets/ecoli_validation/genomes
python3 bin/genome_dataset.py validate \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes datasets/ecoli_validation/genomes
```

Once candidate accessions are selected, archive official AST evidence through
the existing builder. The query is the live NCBI browser export interface and
may return a large, changing table; inspect the resulting PDT/BioSample/
Assembly mappings before editing the cohort manifest.

```sh
python3 bin/ast_dataset.py download \
  --outdir datasets/ecoli_validation/ast_raw \
  --query 'scientific_name:"Escherichia coli"' --max-bytes 536870912
python3 bin/ast_dataset.py build \
  --snapshot datasets/ecoli_validation/ast_raw \
  --outdir datasets/ecoli_validation/ast_processed
```

## Caller environment, database and control

`environment/amrfinderplus.yml` pins AMRFinderPlus 4.2.7 and Python 3.12.
AMRFinderPlus 4.2.7 was installed in a temporary macOS ARM64 Conda environment
using conda-forge and bioconda, without accepting additional channel terms.
Database snapshot `2026-08-07.1` was downloaded from NCBI FTP; all 162 files in
the recorded SHA-256 manifest pass. The real binary reports software 4.2.7,
database 2026-08-07.1, and lists `Escherichia` as a supported organism option.
Package URLs and hashes, database provenance and file hashes remain in
Git-ignored `datasets/amrfinderplus-db/`; the environment was a temporary
one and has not been kept.

On a network-enabled host, create the pinned environment and download a
separate database snapshot. Keep the generated dated release directory
immutable; pass that exact directory to `--database` for controls and cohort
runs. `amrfinder --database_version` reports its software and database
versions, while the call provenance records every database file hash.

```sh
conda env create -f environment/amrfinderplus.yml
database_root="$PWD/datasets/amrfinderplus-db"
mkdir -p "$database_root"
conda run -n amrscan-amrfinder-validation amrfinder_update -d "$database_root"
find "$database_root" -mindepth 1 -maxdepth 1 -type d -print
database="<choose-the-dated-release-directory-printed-above>"
conda run -n amrscan-amrfinder-validation amrfinder --database "$database" \
  --database_version
conda run -n amrscan-amrfinder-validation amrfinder -l
```

NCBI's AMRFinderPlus installation test is a reproducible positive **software
control** from the upstream `ncbi/amr` repository at tag `amrfinder_v4.2.7`:
`test_dna.fa`, `test_prot.fa`, `test_prot.gff` and `test_both.expected`. The
expected output includes known AMR elements such as `blaTEM-156`, `blaPDC-114`
and `vanG`. NCBI identifies the software/database repository as public domain.
This is not an E. coli isolate, cohort member or AST observation. Retrieve
The pinned upstream files are [`test_dna.fa`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_dna.fa),
[`test_prot.fa`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_prot.fa),
[`test_prot.gff`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_prot.gff)
and [`test_both.expected`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_both.expected).
The [NCBI usage guide](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus)
documents this installation test, and the [NCBI AMRFinderPlus repository](https://github.com/ncbi/amr)
contains its public-domain notice. Retrieve those four upstream files into an
ignored scratch directory, preserve each source URL and SHA-256, run the
installation test with the selected immutable database, and compare output
with the expected file. Record `amrfinder --database_version` and the database
file fingerprint. This control was run with software 4.2.7 and database
2026-08-07.1. AMRFinderPlus exited successfully and 25 of 26 expected rows
match. The sole difference is the `blaEC` closest-reference annotation:
expected `WP_063610930.1` (EC-15), observed `WP_498365732.1` (EC-28), with the
same `blaEC` symbol and class. A database-content difference is plausible but
has not been confirmed; this is a partial control result, not an exact-output
pass. The tracked v4.0.23 expected-report fixture remains parser input only.
The synthetic `tests/fixtures/tiny.fna` also completed a real caller no-hit
smoke with zero data rows; it is not a biological negative control, and no-hit
must not be reported as susceptibility.

On a host with the pinned environment and database, the upstream control can
be retrieved and checked without adding sequence data to Git:

```sh
control_dir=amrfinderplus-control-4.2.7
mkdir -p "$control_dir"
for name in test_dna.fa test_prot.fa test_prot.gff test_both.expected; do
  curl -fLsS "https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/$name" \
    -o "$control_dir/$name"
done
shasum -a 256 "$control_dir"/* > "$control_dir/SHA256SUMS"
amrfinder --database "$database" --plus -n "$control_dir/test_dna.fa" \
  -p "$control_dir/test_prot.fa" -g "$control_dir/test_prot.gff" \
  -O Escherichia --print_node -o "$control_dir/test_both.actual"
diff -u "$control_dir/test_both.expected" "$control_dir/test_both.actual"

# Software no-hit smoke only; this synthetic sequence says nothing about AST.
amrfinder --database "$database" -n tests/fixtures/tiny.fna \
  -O Escherichia -o "$control_dir/no_hit.tsv"
test "$(wc -l < "$control_dir/no_hit.tsv" | tr -d ' ')" = 1
```

For every verified cohort genome, run and harmonize the caller separately;
use the exact file names below, with `pdt`, `fasta` and `database` set from the
manifest, retrieval provenance and immutable database location.

```sh
amrfinder --database_version
python3 bin/amrscan.py call --sample "$pdt" --input "$fasta" \
  --database "$database" --organism Escherichia \
  --log "datasets/ecoli_validation/calls/$pdt.log" \
  --raw "datasets/ecoli_validation/calls/$pdt.amrfinder.tsv" \
  --provenance "datasets/ecoli_validation/calls/$pdt.provenance.json"
python3 bin/amrscan.py harmonize --sample "$pdt" \
  --raw "datasets/ecoli_validation/calls/$pdt.amrfinder.tsv" \
  --provenance "datasets/ecoli_validation/calls/$pdt.provenance.json" \
  --output "datasets/ecoli_validation/calls/$pdt.harmonized.tsv"
```

## Evidence linkage and summary

`bin/amrfinder_validation.py summarize` links exact PDT → BioSample → verified
Assembly → caller result → harmonized determinant evidence → AST observations.
It verifies the downloaded assembly hashes, caller input hash, raw report hash,
sample ID, caller version and database fingerprint. AST rows are joined only
by exact versioned PDT accession, then checked against the cohort BioSample
and organism. Assembly conflicts fail closed. Genotype and AST are stored in
separate evidence tables with shared isolate keys; every AST record is kept.

The deterministic summary reports selected cohort size, valid assemblies,
successfully genotyped isolates, isolates with and without AMR determinant
calls, AST-linked isolates, exclusions and reasons, AMR determinant/gene counts,
caller element types, and AST observations grouped by the normalized submitted
phenotype, including ND/I/SDD/NS/HLAR/NA where present. No sensitivity,
specificity, F1, PPV, NPV, prediction or concordance is calculated. The
summary's provenance JSON fingerprints its manifest, genome provenance, AST
table and caller reports. Results are published to a new directory and never
overwrite an existing summary.

The 2026-10-11 cohort run used these commands (`run` is
`datasets/ecoli_validation/cohort5_20261011`, `snapshot` is
`datasets/ncbi_ast_manual/ecoli_amrfinder_cohort5_20261011/repo_export`):

```sh
python3 bin/genome_dataset.py download \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --outdir "$run/genomes"
python3 bin/ast_dataset.py build --snapshot "$snapshot" \
  --outdir "$run/ast_processed"
python3 bin/amrfinder_validation.py summarize \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes "$run/genomes" \
  --ast "$run/ast_processed/ecoli_ast.tsv" \
  --calls datasets/ecoli_validation/calls \
  --exclusions benchmark/cohorts/ecoli_amrfinder_exclusions.tsv \
  --outdir "$run/summary"
```

## Cohort results (2026-10-11)

All five assemblies passed identity and hash checks, all five were genotyped
with AMRFinderPlus 4.2.7 and database 2026-08-07.1 (fingerprint
`sha256:748425a026512fda44008b9d106a4546518cc57dac59eb6432232eb3ac2c978b`,
`--organism Escherichia`, default thresholds), and all five link to AST by
exact PDT. The `PDT000041778.1` call was made earlier on the same FASTA bytes;
its input hash, raw report hash and 162 database file hashes were rechecked and
it was reused. The other four were called in this run.

| PDT | Caller rows (AMR / STRESS) | AST rows (R / S) | Drug endpoints | Eligible | Not eligible |
| --- | --- | --- | --- | --- | --- |
| `PDT000041778.1` | 17 / 1 | 24 / 24 | 22 | 3 | 16 multiple records, 2 conflicting and multiple, 1 unknown drug |
| `PDT000045580.1` | 6 / 0 | 5 / 10 | 15 | 15 | 0 |
| `PDT000045581.1` | 6 / 0 | 4 / 11 | 15 | 15 | 0 |
| `PDT000045582.1` | 6 / 1 | 4 / 11 | 15 | 15 | 0 |
| `PDT000045583.1` | 6 / 0 | 8 / 7 | 15 | 15 | 0 |
| Total | 41 / 2 | 45 / 63 | 82 | 63 | 19 |

All 108 AST observations are kept in `ast_evidence.tsv`, including the
`PDT000041778.1` duplicates across Vitek, MicroScan and Phoenix and its two
S/R conflicts (amoxicillin-clavulanic acid, ampicillin-sulbactam). Moxifloxacin
is not in the ingestion policy and is flagged rather than mapped. No ND, I,
SDD, NS or HLAR observations occur in this cohort. The source supplies no
testing method or breakpoint version for any row, and 8 FDA rows have no
breakpoint standard. Two screened isolates are listed as exclusions. Summary
provenance: manifest
`b7a1da30958231bdd69429fd8a85312df9787e672c0cc9f42bf8a28b5631d35f`, genome
provenance `74736d95407e9f4b86d8778c6f2b6082fd48cfc94d0444f3345120f5f3a77ef8`,
AST table `149accbc8b8bfc14597c039152bcb5b5f925243821ae239d9c558cbb59f15c2b`,
summary `0ec659c8ab826ceb8ffdb85b865c4df2bc52cd5cd064661482f03995009818b8`.

Genotype and AST are reported side by side. No phenotype was inferred from
genes, no breakpoints were applied and no agreement statistic was computed.

## Validation status

| Stage | Status |
| --- | --- |
| Retrieval, linkage and summary implementation | Implemented. Identity/hash mismatches fail closed. |
| Offline tests | 73 tests pass; branch coverage of `bin/` 95% (CI threshold 90%). Mocked NCBI packages, synthetic FASTA, upstream report rows and AST fixtures; software behaviour only. |
| AST snapshot | Archived 2026-10-11: 108 AST and 10 isolate rows for five exact PDTs, plus the 21,042-row candidate list and screening log, with URLs, timestamps, headers and SHA-256. |
| NCBI assembly packages | Five packages: NCBI MD5s, accession, BioSample, organism, strain, taxid and FASTA SHA-256 verified. |
| AMRFinderPlus and database | Real binary 4.2.7 with database 2026-08-07.1; 162-file SHA-256 manifest passes. |
| Real genome calls | Five calls succeeded and were harmonized; provenance records FASTA, report and database hashes. |
| Caller controls | NCBI upstream positive software control: 25/26 rows match, with one unresolved `blaEC` reference annotation difference. Synthetic no-hit smoke: zero data rows. No biological E. coli negative control. |
| Linked cohort summary | Five isolates, five AST-linked, 63 of 82 drug endpoints eligible, two documented exclusions. |

Remaining limitations: five isolates from two submitters, four of them from
one FDA project, chosen partly for clean AST; three AST platforms disagree for
one isolate; no method or breakpoint version from the source; the positive
control is not an exact-output pass; there is no biological negative control;
the Pathogen Detection index is live and will change; and the Conda environment
was temporary and has not been kept. The cohort shows that the evidence chain works on
real data. It is too small and too selected to support any statement about
genotype–phenotype agreement.
