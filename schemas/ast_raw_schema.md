# Raw AST snapshot contract, version 1.0.0

Sources inspected on 2026-10-10:

- [NCBI AST Browser documentation](https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/)
  describes the downloadable measurement-level table.
- Its official [application configuration](https://www.ncbi.nlm.nih.gov/pathogens/static/main/dist/app.ast.js)
  identifies the export and count interfaces. These are browser backend
  interfaces, without a promised stable public API version.
- [FTP ReadMe](https://ftp.ncbi.nlm.nih.gov/pathogen/ReadMe.txt), updated
  2026-10-05, documents isolate metadata and the `AST_phenotypes` summary. That
  summary does not supply the full measurement-level assay evidence.
- [AST in BigQuery](https://www.ncbi.nlm.nih.gov/pathogens/docs/ast_gcp/) provides
  `ncbi-pathogen-detect.pdbrowser.ast`; associated metadata is in
  `ncbi-pathogen-detect.pdbrowser.isolates`. This adapter does not execute billed
  queries or need Google credentials.

## Requests and projections

Both requests use `https://www.ncbi.nlm.nih.gov/pathogens/pathogens-srv/` over
HTTPS. `download` queries the current AST table, then requests isolate metadata
for its distinct versioned PDT accessions in batches of 100. It does not scrape
HTML, download genomes, or use AMR genotypes to recover AST labels.

The count request has `action=retrieve`, `collection=ast` or `isolates`,
`filter=<SOLR query>`, `limit=0`, `fl=id`. It requires `success: true` and a
nonnegative integer `ngout.data.totalCount`.

The table request has `action=solr2txt`, `type=tsv`, `nolimit=on`, a filename,
`fields=<field>|<field>,...`, and this browser query expression:

```text
[display()].from(<collection>).usingschema(/schema/pathogen).matching(q=="<percent-encoded SOLR query>")
```

The expression is then URL-encoded as the `q` parameter. The default AST query
is `scientific_name:"Escherichia coli"`. A small integration query is
`biosample_acc:SAMN05170351`. Metadata queries look like
`target_acc:(PDT000145880.5 OR PDT...)`. Exact URLs are recorded per source file.

The explicit AST projection is:

```text
id target_acc biosample_acc scientific_name antibiotic phenotype measurement_sign
mic mic_secondary disk_diffusion standard platform vendor reagent host
collection_date geo_loc_name isolation_source epi_type bioproject_acc taxgroup_name
```

The metadata projection is:

```text
target_acc biosample_acc asm_acc Run scientific_name collection_date geo_loc_name
isolation_source host epi_type bioproject_acc species_taxid isolate_identifiers
```

All projected headers are required; individual values may be empty. The export
prefixes its first header with `#`, which the parser accepts. UTF-8 BOM and
quoted TSV cells are supported. Duplicate, missing or unexpected columns,
wrong row widths and malformed quoting fail before processed output is
published. Header-only tables are valid zero-record inputs. Additional source
columns require an adapter/schema review, not silent dropping.

Only `method` and `standard_version` are optional AST columns for explicitly
identified enriched archives. They are **not requested from the current
browser**, and the live check did not provide them. An upstream mapping for
these columns has not been validated. Browser `reagent` describes a test kit or
method version and must not be copied into breakpoint version. The platform
name likewise does not establish a testing method.

`mic` and `mic_secondary` are separate concentrations in mg/L; the latter
preserves the second component of a combination. Requesting only `mic` loses
that component. Disk diffusion is a separate field in mm. `measurement_sign`
applies to the exported measurement and is retained unchanged as well as
normalized for MIC. Source AST identifiers, phenotype, standards, optional
platform/vendor/reagent, and all projected original values remain in
`raw_record_json`.

## Snapshot directory and manifest

```text
raw/
  manifest.json
  asts.tsv
  isolates.0001.tsv
  isolates.0002.tsv      # when another accession batch is needed
```

`archive` instead uses unique `ast.0000.tsv` / `metadata.0001.tsv` names and
records the original filenames. It copies bytes from already obtained tables;
its supplied source URL/description and retrieval timestamp apply to the
whole import. For exports with different source URLs/times, retain a manifest
with a separate accurately populated entry for every file, as `download` does.

Manifest-level fields include `schema_version`, `parser_version`, `acquisition`,
`sources`, and the query for downloads. Every source entry contains:

| Field | Meaning |
| --- | --- |
| `role` | `ast` or `metadata`; exactly one AST file is required |
| `source_name` | NCBI Pathogen Detection |
| `source_url` | Exact download URL, or declared archived source description |
| `retrieved_at` | ISO timestamp with timezone; supply the actual retrieval time |
| `filename`, `original_filename` | Safe local basename and original filename |
| `sha256` | Hash of the unmodified source bytes |
| `schema_version`, `parser_version` | Raw adapter version, currently 1.0.0 |
| `expected_rows` | Count endpoint result for each downloaded table |

Snapshots and processed directories must be new paths. A failure removes only
the staging directory created by that operation. Existing snapshots are never
replaced. Build verifies every checksum, timestamp, role and schema; it rejects
path traversal, duplicate manifest filenames, symlink sources and incompatible
versions. The source files are immutable by workflow convention and checksum
validation; filesystem access controls and long-term storage remain the
researcher's responsibility.

The downloader compares each export's row count with the count request. A
mismatch aborts the snapshot. This detects truncation and many concurrent
updates, but cannot establish a transactional freeze across changing NCBI
indexes. Missing or conflicting metadata links remain visible during build.
Archive each successful snapshot and retain its manifest rather than expecting
the live service to reproduce historical bytes.
