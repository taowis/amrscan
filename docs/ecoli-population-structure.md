# E. coli population structure and evaluation metadata

`bin/population_structure.py` adds genome-derived population labels to the
frozen PR #5 research cohort: assembly QC, Achtman 7-gene MLST, Clermont
phylogroup and the NCBI SNP cluster, together with lineage, geography, time,
source and project summaries, grouping variables for later evaluation and a
screen of which holdout designs each antimicrobial could support. It does not
assign train/validation/test sets, call resistance determinants or fit any
model.

Population groups are defined independently of downstream predictive model
performance. Lineage, phylogroup and SNP-cluster labels are computed from the
genome and NCBI metadata only. Phenotype is joined afterwards, for description,
and no threshold in `population_config.json` was set with reference to
phenotype or model results.

## Why population structure matters here

E. coli is strongly clonal. Isolates of the same sequence type, and even more
so of the same SNP cluster, share most of their genome, including resistance
genes and plasmids carried with the lineage. A random isolate-level split puts
near-identical genomes on both sides, so a later predictor can score well by
recognising the lineage rather than the resistance mechanism. Lineage is also
tied to country, year, isolation source and the project that sequenced it, so
any of these can stand in for the others. The labels and tables here exist so
those effects can be measured and controlled before any model is evaluated.

## Relationship to PR #5

The input is the PR #5 build in
`datasets/ecoli_research_v1/ast_census_20261012/processed`. The builder checks
the SHA-256 of `candidate_isolates.tsv`, `isolate_drug_endpoints.tsv` and
`cohort_manifest.json` against the PR #5 `SHA256SUMS`, can require an expected
manifest hash (`--expected-cohort-manifest-sha256`), and refuses duplicate
isolates or an assembly shared by two isolates. Every PR #5 isolate appears
once in every per-isolate table. PR #5 endpoints are read, never changed:
phenotype columns here are counts of PR #5 eligible binary R/S endpoints.

Consumed version: cohort `ecoli_research_v1`, `cohort_manifest.json` SHA-256
`21348c7da115dc0306db07d32839bd91263328a2409c4a1fb65cce8b327d3caf`, PR #5 code
commit `a39b39e`, snapshot manifest SHA-256
`ab974d3e57089fcd66fe43df097d9ee1061edb58511063d370ea6aee6a50321d`.

## Commands

```
conda env create -f environment/population-structure.yml -p <env>

python3 bin/population_structure.py fetch-genomes \
    --cohort <PR5 processed> --cache <genome cache> --workers 3
python3 bin/population_structure.py annotate --tool mlst \
    --cohort <PR5 processed> --cache <genome cache> --env-bin <env>/bin --workers 6
python3 bin/population_structure.py annotate --tool phylogroup \
    --cohort <PR5 processed> --cache <genome cache> --env-bin <env>/bin --workers 3
python3 bin/population_structure.py build \
    --cohort <PR5 processed> --cache <genome cache> --outdir <new dir> \
    --expected-cohort-manifest-sha256 <hash> --workers 6
```

All three long steps are resumable: a rerun skips any assembly whose cached
genome or annotation still matches. The output table contract is in
[`schemas/ecoli-population-structure-schema.md`](../schemas/ecoli-population-structure-schema.md).

## Assembly retrieval

Each assembly is downloaded as an NCBI Datasets v2 genome package
(`/genome/accession/<GCA>/download`, genome FASTA and sequence report) and
checked with the PR #4 `genome_dataset.package_identity` routine: the package
must contain exactly the requested versioned GCA accession, its BioSample must
match the PR #5 isolate and the organism must be E. coli. The FASTA is written
to `<cache>/<GCA>/genome.fna` with a `provenance.json` holding the source URL,
package and genome SHA-256, BioSample, record and base counts, attempts and
retrieval time. Work happens in a staging directory that is renamed into place
only after verification, so an interrupted run leaves no partial genome.
Concurrency is capped at five workers (three were used) with retry and
backoff. A cached genome is reused only if its accession and SHA-256 still
match the provenance record.

Real run: all 10,521 assemblies were retrieved and verified between
2026-10-11T19:07:30Z and 21:14:04Z (20 from the pilot, 10,501 in the main run,
no failures, about 2 h 6 min at three workers, 51 GB on disk). Genomes, packages
and the cache are Git-ignored under `datasets/`.

## Genome QC

`genome_qc.tsv` gives total length, contig count, N50, longest contig, GC over
ACGT and the fraction of non-ACGT bases. The bounds in `population_config.json`
are broad operational limits for usable E. coli assemblies, not biological
limits: fail outside 2.5–8.0 Mb or above 5% ambiguous bases; warn outside
4.0–6.5 Mb, above 500 contigs, below 20 kb N50 or above 0.5% ambiguous bases.
A missing assembly, a checksum mismatch or malformed FASTA also fails. Fail
removes an isolate from population analysis; a warning is kept and reported.

Result: 10,500 pass, 21 warning (19 above 500 contigs, 2 also below 20 kb N50),
0 fail.

## MLST

Tool: Torsten Seemann's `mlst` 2.33.0 with blastn 2.17.0+, run with an explicit
`--scheme ecoli_achtman_4` (the Achtman 7-gene scheme: adk, fumC, gyrB, icd,
mdh, purA, recA) so that no automatic scheme choice is made. The bundled
database is dated 2025-12-29 (scheme downloaded from PubMLST 2026-03-11, scheme
last updated 2025-01-03); the profile table SHA-256 is
`8bba820e34f7f1c61568464e5b3e15558e0a85f449fc5003b9b0146e24ebfa3a`. Version 2.33.0
is pinned because 2.34 and later depend on `libxcrypt1`, which bioconda does not
build for macOS ARM64. The Pasteur scheme is not used; the two schemes give
different numbers and mixing them would make ST labels incomparable.

Every call is stored in `<cache>/<GCA>/mlst.json` with the command, raw output,
return code, input genome SHA-256 and a fingerprint of tool, database and
scheme versions; a cached call is reused only when both still match. Status is
`assigned` only when all seven alleles are exact and the profile maps to an ST.
`novel_profile` (exact alleles, no ST), `novel_or_partial_alleles` (`~` or `?`
calls), `multiple_alleles` and `missing_alleles` keep their raw profile but get
no ST.

Result: 10,398 of 10,521 genomes assigned (98.8%), 1,246 distinct STs, 691 of
them seen once. Unresolved: 95 novel or partial alleles, 14 multiple alleles,
12 novel profiles, 2 missing alleles.

MLST sequence type is an operational lineage grouping and does not capture the
complete evolutionary structure of E. coli. Distinct STs can be close relatives
(single-locus variants within a clonal complex), and a large ST such as ST131
contains sub-lineages (clades A, B, C and the fluoroquinolone-resistant C1/C2
subclades) that differ in resistance. Recombination also moves MLST loci between
backgrounds.

## Phylogroup

Tool: EzClermont 1.0.0, an in silico version of the Clermont quadruplex PCR and
its follow-up assays, run on each assembly with default settings and stored in
`<cache>/<GCA>/phylogroup.json` like MLST. A, B1, B2, C, D, E, F and G are
`assigned`; `cryptic` is a cryptic Escherichia clade; a failed trpBA control is
`control_failed`; anything else (for example `U` or `U/cryptic`) is
`unresolved`. Phylogroup is a separate field and never substitutes for ST.

Result over the cohort: B2 4,150; B1 2,108; A 1,857; D 1,315; G 300; F 270;
C 261; E 189; cryptic 28; unresolved 42; control failed 1.

Agreement with MLST is a useful check on both tools. ST131, ST73, ST95 and
ST1193 are B2 in every typed genome, ST10 is A, ST69 and ST38 are D and ST117
is G, with at most five exceptions each. The exceptions worth review are ST1914
(71 D, 20 A), ST115 (82 D, 11 A) and ST420 (33 B2, 26 unresolved).

## NCBI SNP clusters

The Pathogen Detection `erd_group` field gives the isolate's SNP cluster as a
versioned PDS accession. The version changes when NCBI rebuilds the cluster, so
the field `ncbi_snp_cluster` is the unversioned accession and
`ncbi_snp_cluster_versions_json` keeps what each of the two PR #5 pipeline rows
said. In PR #5, 51 isolates differed between rows only by version; they are
`assigned` here. An isolate with two different unversioned accessions would be
`conflicting`; there are none. Isolates outside any cluster are `missing`, which
in Pathogen Detection mostly means no other isolate within the cluster distance,
not missing data.

NCBI SNP-cluster membership is not interpreted as proof of direct
epidemiological transmission. Clusters are built by NCBI with single-linkage
thresholds over isolates submitted to the database at the time, and they change
as new genomes arrive. Here they are used only as a marker of near-identical
genomes that should not be split across evaluation sets.

Result: 4,108 isolates carry a cluster, 6,413 do not. There are 2,080 distinct
clusters, 1,139 of them with only one cohort member (the rest of the cluster is
outside the cohort). 941 clusters have two or more cohort isolates, covering
2,969 isolates. Largest: PDS000202052 (61, all ST4243, all Canada 2012–2013, one
project), PDS000117091 (50, ST131, USA and UK, 2009–2024, ten projects),
PDS000269870 (50, ST1193, mostly USA), PDS000120850 (28, ST428), PDS000003047
(25, ST131, one hospital project) and PDS000037878 (25, ST3580). Ten clusters
contain more than one ST, as expected when a cluster spans a single-locus
variant.

## Lineage definition

```
lineage_id = "ST<n>"  when mlst_status == "assigned"
lineage_id = "UNKNOWN" otherwise
```

Nothing else feeds the lineage: no phylogroup fallback, no SNP cluster, no
clonal-complex merging and no phenotype. `UNKNOWN` is a missing label
(`is_lineage = false`), not a lineage, and is excluded from lineage holdout
counts. Major lineages for reporting are those with at least 50 isolates
(`major_lineage_min_isolates`).

## Missing and ambiguous labels

123 isolates (103 of the 8,202 population-eligible) are `UNKNOWN`. They remain
in the population, in country, temporal and source groupings, and in all
counts, but are not lineage-holdout eligible. They should not be assigned to an
ST by nearest profile without review. The 42 unresolved and 1 control-failed
phylogroup calls are kept as such.

## Clonal redundancy

The cohort is not deduplicated. No isolate is dropped or chosen per ST or per
cluster. Assembly, BioSample, ST and SNP cluster are all kept per isolate in
`population_metadata.tsv` and `evaluation_groups.tsv`, so a later split can
group on any of them and the effect can be measured. Redundancy is substantial:
ST131 alone is 10.1% of the cohort (9.3% of population-eligible isolates), the
top 20 STs hold 49%, and 2,969 isolates share an NCBI SNP cluster with at least
one other cohort isolate.

## Lineage distribution

Top 20 STs (all 10,521 isolates): ST131 1,060; ST10 489; ST73 439; ST69 339;
ST95 317; ST372 314; ST127 292; ST1193 228; ST38 221; ST117 211; ST12 210;
ST58 204; ST101 134; ST155 126; ST162 118; ST115 93; ST88 91; ST1914 91;
ST4243 90; ST998 79. Among the 8,202 population-eligible isolates the order is
ST131 763, ST10 456, ST69 280, ST73 276, ST95 270. Most of the B2 extraintestinal
lineages (ST131, ST73, ST95, ST127, ST1193, ST12) come from clinical projects;
ST10, ST58, ST117 and ST101 come mostly from food and animal surveillance.
`lineage_summary.tsv` has every ST with counts, country and year range,
sources, endpoints and phylogroups.

## Geography

Population-eligible isolates are 71% USA (5,852), 16% Canada (1,292) and 6%
Netherlands (505); no other country has more than 91, and 129 have no country.
Lineage and country are entangled: ST131, ST73, ST95, ST127, ST1193, ST399,
ST420, ST428 and ST3580 are each at least 90% USA, and ST1914 and ST4243 are at
least 92% Canada. A country holdout is therefore also partly a lineage holdout,
and for most drugs only USA, Canada and the Netherlands have enough isolates of
both phenotypes. `lineage_country_summary.tsv` has the full table.

## Time

Collection years run from 2009 to 2025, with peaks in 2013 (1,123), 2017
(1,015) and 2019–2020 (1,863 and 1,570) that follow individual projects rather
than sampling effort; 204 eligible isolates have no usable year. Year follows
project closely: 1,515 of the 1,570 isolates from 2020 come from PRJNA809394,
and isolates from 2021 on come mostly from projects absent earlier. The
temporal cutoff is fixed at 2021 (see [decision 0001](decisions/0001-ecoli-evaluation-design.md));
1,197 population isolates (15% of those with a year) are recent. ST115, ST1434,
ST155, ST1914, ST355, ST3580, ST399, ST4243, ST428 and ST73 are at least 90%
before 2021. Several drugs changed panel over time: from 2021 on sulfisoxazole
has 1 endpoint against 4,068 earlier, streptomycin 16 against 4,115 and
kanamycin 21 against 1,208, and the recent cefotaxime and cefazolin endpoints
are 362 R to 2 S and 257 R to 12 S. A temporal split therefore also changes the
drug set and the case mix. `lineage_year_summary.tsv` has the full table.

## Source

Only 4,031 population-eligible isolates have a usable source category: human
3,373, food 528, animal 130; 4,171 are `unknown` (mostly NCBI records without an
isolation-source category). ST131, ST73, ST95, ST127, ST1193, ST12, ST38, ST69
and ST420 are at least 97% human where the source is known; ST1434 and ST399 are
mostly food. Source is close to a proxy for lineage and project.
`lineage_source_summary.tsv` has the full table.

## Phenotype overlay

`lineage_antimicrobial_summary.tsv` gives R and S counts for each antimicrobial
and ST from population-eligible binary endpoints. It describes the data; it was
not used to define or alter any group. Resistance is concentrated in a few
lineages for some drugs: of 1,286 ciprofloxacin-resistant isolates, ST131 has
461 and ST1193 181 (together 50%), while ST1193 is 181 R to 4 S. ST131 is also
the largest single contributor of resistance to ceftriaxone (326 of 1,459),
ampicillin (599 of 3,621), trimethoprim-sulfamethoxazole (307 of 1,382) and
gentamicin (221 of 1,237). ST4243 (90 isolates, 61 of them in one Canadian SNP
cluster) is 84 R to 6 S for gentamicin. No drug reaches the review threshold of 50% of resistant
isolates in one ST or 20% in one SNP cluster.

## Projects and submitters

There are 68 BioProjects from 53 submitters. The three largest supply 71% of
population-eligible isolates: PRJNA292663 (FDA, 2,655, all USA; 528 food, the
rest without a source category), PRJNA809394 (Walter Reed Army Institute of
Research, 2,076, all USA, all human) and PRJNA966974 (Canadian Food Inspection
Agency, 1,128, all Canada, no source category). All the large projects are
single-country; only PRJNA1322038 (seven countries) and PRJNA318589 (two) span
more than one. Some lineages come almost entirely from one project
(ST4243 99%, ST3580 97%, ST1914 92%, ST399 90%, ST127 87%, ST428 87%, ST1434 88%).
For 15 antimicrobials at least half of the resistant isolates come from one
project, for example polymyxin B (106 of 107 from PRJNA1322038),
ceftazidime-clavulanic acid (54 of 54 from PRJNA885502) and kanamycin (151 of
178 from PRJNA966974). Project-specific testing, breakpoints and sampling can
look like biology. `project_summary.tsv` has every project.

## Leakage risks

The main ways information could leak across a future split are:

- the same or near-identical genome on both sides (same SNP cluster or ST);
- lineage standing in for country, year, source or project, and the reverse;
- a resistance phenotype that is mostly one lineage or one project, so the
  model learns the label source rather than the mechanism;
- drug panels and breakpoint versions that change with time and project;
- the 6,413 isolates with no SNP cluster, which are not all unrelated: NCBI
  only clusters isolates within its distance threshold, and cohort members
  may be close to each other through genomes outside the cohort.

`confounding_summary.tsv` lists every flag raised by the thresholds in
`population_config.json` (61 in this build: 13 lineage–source, 11 lineage–country,
11 lineage–project, 10 lineages absent after the 2021 cutoff, 15
resistance–project, 1 large SNP cluster). No ST reaches the dominant-lineage threshold of 10% of
the eligible population (ST131 is 9.3%), and MLST success is above the 90%
review threshold. These flags do not exclude anything.

## Holdout eligibility and feasibility

`evaluation_groups.tsv` gives each isolate its lineage, SNP-cluster, country,
year and source group and five eligibility flags, plus the evaluation-design
columns described in the next section. It contains no split.

- population: PR #5 eligible isolate (at least one binary R/S endpoint) and
  genome QC not `fail`;
- lineage: population and `mlst_status = assigned`;
- country: population and a known country;
- temporal: population and a valid collection year;
- source: population and a source in human/animal/food/environment/other.

Counts: population 8,202; lineage 8,099; country 8,073; temporal 7,998;
source 4,031.

`holdout_feasibility.tsv` screens each antimicrobial and axis with the rules
in `feasibility_rules` and records the numbers behind the status:

- lineage: feasible with at least 50 R, 50 S, five STs with both R and S, and
  no ST above half the endpoints; marginal with 20 R, 20 S and three STs;
- country and source: feasible when at least three groups of 50 or more
  isolates each have 10 R and 10 S (country also needs the largest country
  below 80%); marginal with two such groups;
- temporal: endpoints from the recent period (2021 on) must number at least
  100 with 20 R and 20 S for feasible, or 30 with 5 R and 5 S for marginal.

Every row carries `HUMAN_REVIEW_REQUIRED`: the status says whether a holdout
could be populated, not whether it should be used. What is used is decided in
the next section.

For the 44 antimicrobials with population endpoints:

| Axis | Feasible | Marginal | Not feasible |
| --- | --- | --- | --- |
| Lineage | 30 | 3 | 11 |
| Country | 12 | 10 | 22 |
| Temporal | 24 | 8 | 12 |
| Source | 1 | 9 | 34 |

Lineage-feasible: amikacin, amoxicillin-clavulanic acid, ampicillin,
ampicillin-sulbactam, aztreonam, cefazolin, cefepime, cefotaxime,
cefotaxime-clavulanic acid, cefoxitin, ceftazidime,
ceftazidime-clavulanic acid, ceftiofur, ceftriaxone, chloramphenicol,
ciprofloxacin, ertapenem, gentamicin, imipenem, kanamycin, levofloxacin,
meropenem, nalidixic acid, piperacillin-tazobactam, streptomycin,
sulfisoxazole, tetracycline, tobramycin, trimethoprim,
trimethoprim-sulfamethoxazole. Marginal: azithromycin, nitrofurantoin,
polymyxin B.

Country-feasible: amoxicillin-clavulanic acid, ampicillin, cefoxitin,
ceftazidime, ceftiofur, ceftriaxone, chloramphenicol, ciprofloxacin,
gentamicin, nalidixic acid, tetracycline, trimethoprim-sulfamethoxazole.
Marginal: azithromycin, cefazolin, cefepime, cefotaxime, ertapenem, kanamycin,
levofloxacin, streptomycin, sulfisoxazole, trimethoprim.

Temporal-feasible: the lineage-feasible list without aztreonam, cefazolin,
cefotaxime, kanamycin, streptomycin and sulfisoxazole. Marginal: azithromycin,
aztreonam, cefazolin, ceftolozane-tazobactam, colistin, fosfomycin,
nitrofurantoin, polymyxin B.

Source: only tetracycline passes; nine more are marginal.

Drugs with too few endpoints for any holdout: amoxicillin,
ceftazidime-avibactam, doripenem, imipenem-relebactam, meropenem-vaborbactam,
minocycline, piperacillin and tigecycline.

A count screen does not see composition. Several lineage-feasible drugs come
mostly from one project (cefotaxime-clavulanic acid and
ceftazidime-clavulanic acid from PRJNA885502, kanamycin and ceftiofur from
PRJNA966974), and should be reviewed before they are used.

## Evaluation design

The screen above was reviewed and the evaluation design decided in
[decision 0001](decisions/0001-ecoli-evaluation-design.md), recorded in the
`evaluation_design` section of `population_config.json` (version 1.1.0). In
short:

- SNP clusters are never split, in any design (`split_unit`).
- Lineage holdouts use connected components of shared ST or SNP cluster
  (`lineage_component`); 1,071 components, seven of which join a large ST to a
  rare single-locus variant. UNKNOWN isolates join through their SNP cluster
  (30); the 73 without one are never test data (`lineage_cv_test_eligible`).
- Leave-one-project-out (`project_group`) is the main external-validity
  analysis. Geographic transfer is USA to all other countries pooled
  (`geography_group`). There is no source holdout.
- The temporal cutoff is 2021 (`period`), reported as a time and study shift.

`evaluation_plan.tsv` applies the drug tiers and lists, for each antimicrobial,
the numbers behind its tier and which analyses it enters:

| Tier | Drugs | Analyses |
| --- | --- | --- |
| Primary (16) | ampicillin, amoxicillin-clavulanic acid, cefoxitin, cefazolin, ceftriaxone, ceftazidime, cefepime, ciprofloxacin, levofloxacin, nalidixic acid, gentamicin, trimethoprim-sulfamethoxazole, tetracycline, chloramphenicol, streptomycin, sulfisoxazole | lineage cross-validation, within-ST by cluster, leave-one-project-out, geographic transfer; temporal for 13 (not cefazolin, streptomycin, sulfisoxazole) |
| Secondary (9) | meropenem, ertapenem, imipenem, amikacin, tobramycin, piperacillin-tazobactam, aztreonam, ceftiofur, kanamycin | lineage cross-validation only |
| Excluded (19) | all others, each with its reason in `tier_reason` | none |

## Outputs and provenance

The build writes, into a new directory, the tables listed in the schema plus
`population_summary.json`, `population_manifest.json` and `SHA256SUMS`. Output is
written to a staging directory and renamed in one step, an existing directory is
never overwritten, and a rebuild from the same inputs is byte-identical
(`generated_at` is the latest genome retrieval time, not the wall clock).
The manifest records the PR #5 manifest and file hashes, the code commit and
builder SHA-256, the config SHA-256, the genome retrieval window and a hash over
all genome hashes, every tool and database version, and annotation counts.

Real build: `datasets/ecoli_research_v1/population_20261015` (Git-ignored,
builder and config 1.1.0; the 1.0.0 build `population_20261012` is kept for
comparison; it differs only in the temporal screen, the review flags, the
evaluation columns and plan, and the summary and manifest), genome hash-list SHA-256
`402a4ca2528e545259187d78e7c96b995f19509cc5bac7a2b739bed3559d3692`. A second build
into a separate directory produced the same `SHA256SUMS`, and every file
verified. The build takes about 9 minutes with six threads (it reads and hashes
51 GB of FASTA) and peaks at about 640 MB of memory. Tool environment:
`environment/population-structure.yml` (mlst 2.33.0, EzClermont 1.0.0,
Python 3.12), with an explicit lockfile kept next to the data.

## Limitations

- MLST and phylogroup are coarse. Neither resolves sub-lineages such as ST131
  clade C, and neither replaces a core-genome phylogeny.
- NCBI SNP clusters depend on what NCBI held when the cluster was built, change
  between releases, and leave 61% of the cohort unclustered.
- Country, year and source come from submitter metadata and are often missing
  or coarse; half the population has no source category.
- The cohort is a convenience sample of public data dominated by three US and
  Canadian projects. It is not a surveillance sample of E. coli.
- Feasibility rules and review thresholds are screening values chosen before
  looking at phenotype–lineage results. The decisions built on them are in
  decision 0001 and were made before any genotype–phenotype or model result.
- EzClermont relies on primer matches in assemblies and can misplace genomes
  when a target gene is split across contigs or diverged; the ST1914 and ST115
  D/A splits may be such cases.
- `mlst` 2.33.0 with its bundled database was used rather than a fresh PubMLST
  download; STs defined after 2025-01-03 are reported as novel profiles.
