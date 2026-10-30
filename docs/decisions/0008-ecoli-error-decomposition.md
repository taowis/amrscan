# 0008. Design of the E. coli error decomposition (error_decomposition_v1)

- Status: accepted
- Date: 2026-10-26
- Applies to: `benchmark/cohorts/ecoli_research_v1/error_decomposition_v1.json`
  (config 1.0.0), `benchmark/mappings/ecoli_amr_mechanisms_v1.tsv`, read by
  `bin/error_decomposition.py` 1.0.0
- Written before: any capture status, mechanism-novelty flag or error
  category was computed. The only results read were the generalisation_v1
  outputs already reported in decision 0007 and its guide.

## Context

The generalisation benchmark (decision 0007) found that a random split is only
mildly optimistic, that lineage leakage explains little of it, and that the
larger drops under country, project, geography and temporal holdouts are
concentrated in a few drugs. The curated determinant calls tied or beat the
model under every holdout. The benchmark's features are FracMinHash 31-mers
at scaled 200, which see about one single-base change in seven, and several of
the drugs with large drops have mutation-based resistance. This analysis asks
how much of the error under domain shift comes with determinants the features
do not encode, with determinants absent from training, or with no established
determinant at all. It reuses the frozen predictions and does not change any
split, endpoint, label or model.

## Decisions

### 1. Mechanism catalogue

`ecoli_amr_mechanisms_v1.tsv` takes every rule of determinant mapping v1
(decision 0002) for the 24 study drugs, unchanged: determinant,
interpretation, `requires`, evidence source and reference. It adds:

- a mechanism class from fixed rules in the config: genes are
  `acquired_gene`, except the intrinsic blaEC family (`intrinsic_gene`, which
  mapping v1 marks `does_not_predict`); point mutations take their class from
  the gene (gyrA, gyrB, parC, parE, folP, 16S, 23S: `target_site_mutation`;
  ampC and blaTEM promoter: `promoter_mutation`; regulatory and porin genes
  listed but not in mapping v1). A mutation in an unlisted gene stops the
  build;
- for point mutations, the AMRFinderPlus 2026-08-07.1 locus (reference
  protein accession or reference DNA region, and position) and class;
- gene function (`efflux` or `other`) from the AMRFinderPlus family names.

No mechanism is inferred from model features. The catalogue's checksum is in
the config. Alternatives: a hand-written mechanism list per drug from the
literature (rejected: it would duplicate mapping v1 with weaker provenance);
CARD ontology classes (not used in mapping v1, so the two would disagree).

### 2. Known mutations

The AMRFinderPlus 4.2.7 calls behind PR #7 were made with
`--organism Escherichia` and already contain point mutations with gene
coordinates. They are read again, each checked against its provenance
(caller version, database files, genome SHA-256). `known_mutations.tsv` lists
every called point mutation with its catalogue status per drug.

### 3. Representation capture

Capture asks whether the PR #14 features hold information that separates a
genome carrying the determinant from one without it. Phenotype and model
output are never read. In each genome, after its SHA-256 is checked against
the PR #8 manifest:

- **Gene.** Every 31-mer inside the called gene is hashed with the PR #14
  sourmash parameters. The gene is `captured` if at least one hash is a
  kept PR #14 feature, `partially_represented` if some are sampled but all
  were pruned, `not_represented` if none is sampled.
- **Point mutation.** The site is located from the AMRFinderPlus gene
  coordinates and the reference protein or DNA. The changed bases are found by
  comparing the mutant codon with the nearest wild-type codon, or the
  reference DNA. The 31-mers overlapping a changed base are hashed for the
  mutant allele (from the genome) and the wild-type allele (the same sequence
  with the wild-type base). The mutation is `captured` if a mutant-allele
  k-mer is a kept feature, a direct presence signal.
  `partially_represented` means no such k-mer, but a mutant k-mer is sampled
  and pruned, or a wild-type k-mer is kept, so the mutation shows only as a
  missing wild-type feature. It is `not_represented` if neither allele has a
  sampled k-mer at the site.
- **Unknown.** The site could not be located, the call is not a listed
  substitution or insertion, or the window has ambiguous bases.

A per-task flag adds whether a captured k-mer also passes that task's training
filter (present in at least 5 and absent from at least 5 training isolates).
A captured feature that never varies in training cannot be used by the model.

### 4. Mechanism and population novelty

For each task, a determinant is seen in training when at least 5 training
endpoints carry it (any phenotype; 5 is the PR #14 feature minimum). An
endpoint's mechanism is seen when all its established determinants are, unseen
otherwise. DET- endpoints have no established mechanism. The exact
determinant combination is seen when at least one training endpoint has it.
Population novelty is the PR #14 ST, SNP cluster and country flags, plus
project. None of these flags uses predictions.

### 5. Drug groups

The groups are fixed from the PR #14 common-cohort pooled ΔAUROC against random
on the same isolates, using each drug's largest drop over its holdouts:

- `high_loss`: at least 0.05, with the interval above zero;
- `stable`: below 0.02;
- `intermediate`: everything else.

`freeze` computes the groups and writes them to the analysis manifest before
any mechanism table exists. All drugs are analysed. The case studies are every
primary drug in the high_loss or stable group, whatever the mechanism results
turn out to be.

### 6. Error categories

For a false negative (R predicted S) the category is:

- `R_KNOWN_DET_*`, `R_KNOWN_MUT_*` or `R_KNOWN_DET_AND_MUT_*` when the
  endpoint has established determinants, each with suffix `REPRESENTED` if
  any of them is captured and `NOT_REPRESENTED` otherwise;
- `R_NO_ESTABLISHED_DET` for DET-;
- `R_UNRESOLVED` for uncertain determinants, or capture unknown for all of
  them;
- `R_GENOTYPE_UNAVAILABLE` when PR #7 gave no determinant status.

A false positive is `S_KNOWN_DET_PRESENT`, `S_NO_KNOWN_DET` or
`S_UNRESOLVED`. The evidence columns (capture, novelty, QC) stay separate, so
an endpoint with both mechanisms is not forced into one.

### 7. Mutation augmentation

This is a sensitivity analysis and is labelled `mutation_augmented`; it is
never reported as the benchmark.

- **Features.** The PR #14 features plus one 0/1 column per catalogue
  point mutation with a `predicts_resistance` or `component` rule, set from
  the AMRFinderPlus calls (partial hits excluded). The columns are the same
  for every drug.
- **Which drugs.** A drug is refitted only if its catalogue has such a rule,
  at least 20 common-cohort endpoints carry the mutation, and fewer than 90%
  of those carriers have it captured.
- **How.** The common-cohort tasks are refitted with the exact frozen splits
  (checked against the PR #14 split hashes), using the same model, C grid,
  inner CV, seeds, feature filter and threshold. There is no new tuning.
- **Control.** `permuted_control` repeats the refit with the indicator rows
  permuted across isolates (seed 20261029). It keeps the prevalence but breaks
  the link to each genome, so it shows how much a change comes from adding
  columns at all.

Gap rescue is (random − holdout) ΔAUROC of the original model minus that of
the augmented model, on the same isolates. The interval comes from a paired
bootstrap of 1,000 split-unit resamples (seed 20261028).

## Corrections after freezing

Two corrections were made to how the code reads its inputs; no rule, threshold
or group changed.

- `calls` first compared a provenance field that does not exist, so the
  genome-hash check never ran. It now requires `input_sha256` to equal the
  PR #8 hash, `status` to be success and `raw_sha256` to match the report.
  All 10,521 call sets pass, as before.
- PR #14 stores a determinant status only for mapping-eligible endpoints and
  writes NA for the 2,826 mapping-uncertain ones. The first context build read
  that NA as "genotype unavailable" (805 false negatives). These endpoints have
  calls, so they now take the PR #7 status (`uncertain`) and fall under
  `R_UNRESOLVED` or `S_UNRESOLVED`.

Two later fixes in the report code changed no number: the performance helper
now returns the endpoint counts its tables list (the first report run stopped
on their absence), a rank correlation over a constant input is reported as
missing instead of NaN, and empty detail fields of the audit tables are
written as NA.

## Consequences

- PR #14 outputs are inputs. The config pins their hashes, and every stage
  stops on a difference.
- Correct model calls on R/DET- endpoints and errors without an established
  determinant are described as unexplained resistance, not as new
  mechanisms.
- Capture is a property of the features, not of the fitted model. The model
  may still not use a captured k-mer.
- The scale of the PR #14 sampling sets an expected capture rate of about 14%
  for a single-base change (31 overlapping k-mers at 1/200); observed rates
  are reported per determinant.
