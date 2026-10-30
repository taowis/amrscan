# E. coli AST research cohort builder

`bin/ecoli_cohort.py` builds an auditable, phenotype-only isolate and
isolate–antimicrobial census from a validated raw NCBI Pathogen Detection AST
snapshot. It reuses `bin/ast_dataset.py` for archive validation, normalization,
and source joins. It does not call resistance genes, infer phenotype from
genotype, or fit prediction models.

## Census (snapshot retrieved 2026-10-11 UTC)

### Source

The candidate population is the archived isolate metadata export for
`scientific_name:"Escherichia coli" AND AST_phenotypes:* AND asm_acc:GCA_*`
(`datasets/ncbi_ast_manual/ecoli_ast_assembly_candidate_list_20261011/isolates_candidates.tsv`,
SHA-256 `66eee1964b4777086cc32e5d0aeb914cbfa5bcf9f8f324593deae5d28ae3f9e1`,
retrieved 2026-10-11T11:28:39Z). It has 21,042 rows but 10,521 distinct
versioned PDT accessions, two pipeline rows per PDT; 21,042 is not an isolate
count.

`asm_acc` exists only in the isolates collection, so the same query returns
zero rows from the AST collection. AST and isolate metadata were therefore
exported from the official Pathogen Detection service by exact `target_acc`, in
106 batches of up to 100 PDTs in numeric PDT order, between 17:04:49Z and
17:12:08Z on 2026-10-11. Every export returned HTTP 200 on the first attempt,
matched its count request and an unchanged recount at 17:12–17:14Z, contained no
PDT outside its batch and left no PDT without AST rows. Batches do not overlap
and cover all 10,521 PDTs. Batch membership, URLs, times, row counts, headers
and SHA-256 are kept with the raw data. The AST batches are concatenated byte
for byte into the one AST table the snapshot format allows, with a per-batch
row-range trace; metadata batches are used unchanged. Snapshot manifest SHA-256:
`ab974d3e57089fcd66fe43df097d9ee1061edb58511063d370ea6aee6a50321d`.

Identity: all 162,701 AST rows match their PDT's BioSample and organism; both
metadata rows of every PDT agree on BioSample, GCA assembly, organism, taxid
(562) and strain; no BioSample or assembly maps to two PDTs; there are no
duplicate AST record ids; and no identity field changed from the candidate list.
`ast_dataset.py validate` passes and the builder marks all 10,521 isolates
`verified`. Strain is blank for 1,888 isolates. The two pipeline rows differ in
SNP cluster for 51 isolates.

### Counts

| Unit | Count |
| --- | --- |
| AST observations | 162,701 (none excluded for organism) |
| Unique PDT isolates | 10,521, all with a GenBank assembly |
| Isolate–drug endpoints | 162,629 (32 with 2 observations, 20 with 3) |
| Endpoints by status | S 98,574; R 25,785; ND 35,994; I 2,206; SDD 41; NS 21; R/S conflict 6; mixed nonbinary 2 |
| R/S endpoints with a drug outside the policy | 2,670 |
| Binary R/S (and overall-eligible) endpoints | 121,689 |
| Country / temporal / source-eligible endpoints | 119,589 / 120,087 / 60,720 |
| Isolates with at least one eligible endpoint | 8,202 overall; 8,073 country; 7,998 temporal; 4,031 source |
| Isolates with no eligible endpoint | 2,319 (2,308 ND only, 11 I only) |

The endpoint arithmetic closes: 98,574 + 25,785 − 2,670 = 121,689 eligible, and
the remaining 40,940 are nonbinary, conflicting or have an unrecognised drug.
The 72 extra observations are repeated tests that collapse into one endpoint.
The eight conflicting endpoints come from four isolates (`PDT000041778.1`,
`PDT000041780.1`, `PDT000077415.3`, `PDT001463387.1`).

Almost all ND observations carry an MIC without an interpretation; they stay ND
because no breakpoints are applied. All 716 United Kingdom isolates are ND only.
15,458 endpoints use drug names outside the ingestion policy, mostly
veterinary agents (enrofloxacin, cefpodoxime, marbofloxacin, cefovecin,
doxycycline, cephalexin, orbifloxacin, pradofloxacin); they are kept and
flagged, not mapped.

### Distributions

Antimicrobials: 92 names, 44 with at least one binary-eligible isolate. The
largest binary denominators are ampicillin (7,787 isolates, 3,621 R),
ciprofloxacin (7,707, 1,286 R), gentamicin (7,385), trimethoprim-sulfamethoxazole
(7,304), tetracycline (7,088), ceftriaxone (7,024), amoxicillin-clavulanic acid
(6,570) and meropenem (6,238, 189 R).

Countries (isolates): USA 7,401; Canada 1,334; United Kingdom 716; Netherlands
505; missing 130; then 22 countries with 91 or fewer. Collection years run
2009–2025; 2017–2020 hold 6,966 isolates; 204 have no valid year, including all
89 from Australia. Source category by isolate: human 4,096, animal 1,515, food
528, unknown 4,382. "Unknown" is a category value, so the missingness table
reports source category as complete even though 42% of isolates lack a usable
source. Host is missing for 32% of isolates and SNP cluster for 61%.

All 92 antimicrobial rows keep `HUMAN_REVIEW_REQUIRED`; no study-drug panel or
threshold is set. No genotype or caller output was used. This is a descriptive
census of NCBI-submitted phenotypes, dominated by US and Canadian surveillance
submissions; it is not a sample of E. coli in general.

## Rebuild from a complete archived AST snapshot

```bash
census=datasets/ecoli_research_v1/ast_census_20261012
python3 scripts/ecoli_research_v1/fetch_ast_batches.py \
  --candidates datasets/ncbi_ast_manual/ecoli_ast_assembly_candidate_list_20261011/isolates_candidates.tsv \
  --outdir "$census"
python3 scripts/ecoli_research_v1/assemble_snapshot.py "$census"
python3 scripts/ecoli_research_v1/audit_identity.py "$census" \
  datasets/ncbi_ast_manual/ecoli_ast_assembly_candidate_list_20261011/isolates_candidates.tsv
python3 bin/ast_dataset.py validate --snapshot "$census/snapshot"
python3 bin/ecoli_cohort.py build \
  --snapshot "$census/snapshot" \
  --config benchmark/cohorts/ecoli_research_v1/cohort_config.json \
  --outdir "$census/processed"
```

`assemble_snapshot.py` refuses to run if any batch hash, count or recount
disagrees, or if batches overlap or miss PDTs. Do not use
`ast_dataset.py download` with the candidate-list query: it returns a
header-only AST table. A fresh fetch reflects the live index on that day and
will not reproduce these counts exactly.

Output directories are immutable: choose a new path for each source snapshot or
configuration version. The output manifest records source URLs, timestamps,
hashes, parser and schema versions, configuration hash, builder version and the
code commit. `SHA256SUMS` covers the generated tables and manifests. Rebuilding
the same snapshot with the same code and configuration yields deterministic
tables; the recorded time is the latest source retrieval time, not build wall
clock time.

## Observation and eligibility rules

All normalized AST observations for E. coli are preserved in
`ast_observations.tsv`, including repeated tests, discordant results, nonbinary
categories, missing results, and unknown antimicrobials. A repeated set of the
same recognized R or S result becomes one derived isolate–drug endpoint with
`repeated_consistent` status; the original rows remain available. R/S discordance
is `CONFLICT`; mixtures of binary and nonbinary values are
`NON_BINARY_OR_CONFLICT`. Neither enters the binary endpoint denominator.
Intermediate, SDD, NS, HLAR, ND, missing, and unrecognized values are retained
and counted separately; none is silently converted to R or S.

The versioned config defines a binary endpoint as a recognized antimicrobial
with one unambiguous R or S result on a versioned PDT accession. Overall
eligibility additionally requires a valid BioSample, a versioned GenBank GCA
assembly, E. coli identity, and no identity disagreement across linked source
records. Country, collection year, and known source category are separate
analysis-specific eligibility dimensions. The output reports explicit reason
codes instead of filling unknown metadata from genotype or nearby records.
Study-drug sufficiency thresholds require human review and remain unset.

## Outputs

The builder emits normalized AST observations, candidate isolates,
isolate–drug endpoints, eligible isolates, analysis-specific exclusions,
antimicrobial/country/year/source summaries, missingness and conflict tables,
a reconciled cohort flow, and a provenance manifest. Counts distinguish
observations, isolates, and isolate–drug endpoints. Summary denominators are
documented in column names; repeated observations never increase the number of
isolates.

The separate `inventory` command can summarize an archived metadata-only
candidate list. Its result explicitly sets the AST observation count to null
and labels the endpoint census unavailable. It must not be used as a substitute
for the full `build` census.

See [the schema contract](../schemas/ecoli-research-cohort-schema.md),
[configuration](../benchmark/cohorts/ecoli_research_v1/cohort_config.json),
and [tests](../tests/README.md#e-coli-research-cohort-builder).
