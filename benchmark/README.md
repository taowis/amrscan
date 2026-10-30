# Future benchmarking and generalisability interface

The interface below was set out in PR #2, before any split was implemented.
Future evaluation takes frozen genotype/AST snapshots, an isolate mapping and a
split manifest with columns `sample_id`, `split`, `holdout_axis`, `holdout_value`,
`dataset_version`, `split_version`, `seed`. Allowed `split` values are train,
validation and test. Unknown metadata must be accounted for explicitly.

Required evaluation designs:

- **Lineage holdout:** keep all linked isolates/assemblies from a held-out lineage
  outside training and validation; specify the lineage typing method/version.
- **Country holdout:** freeze geographic definitions and handle unknown or
  conflicting provenance separately.
- **Time holdout:** set a cutoff using collection date, preserve date precision,
  and exclude missing/ambiguous dates or evaluate them as a documented cohort.

Prevent identical isolates, near duplicates and linked assemblies from crossing
splits. Fit imputation, scaling and feature selection on training data only.
Keep the test set untouched during tuning, report cohort sizes and exclusions,
per-drug class balance, calibration and uncertainty, and compare to simple
baselines. Publish split manifests and their hashes with caller/database and
AST breakpoint versions so changes in data/callers are distinguishable from
changes in model performance. Random splits alone do not establish transportability.

The phenotype-only E. coli AST census configuration and its current data
availability are documented in [the research cohort guide](../docs/ecoli-research-cohort.md).
Lineage (Achtman ST), SNP-cluster, country, year and source grouping variables,
holdout feasibility screening and confounding flags for that cohort are in
[the population structure guide](../docs/ecoli-population-structure.md); no split
is assigned there.
The established-determinant baseline that later models are compared with uses
the mapping in `mappings/` and is described in
[the genotype–phenotype baseline guide](../docs/ecoli-genotype-phenotype-baseline.md).
Mapping v1 is prespecified and is the primary baseline; v2
([decision 0003](../docs/decisions/0003-ecoli-determinant-mapping-v2.md)) was
written after seeing the v1 results and is exploratory. The primary baseline
uses PR #5 labels with R and S swapped for one project whose labels run
against its MICs
([decision 0004](../docs/decisions/0004-ecoli-phenotype-label-correction.md)).

The generalisation benchmark (`bin/generalisation_benchmark.py`, config
`cohorts/ecoli_research_v1/generalisation_v1.json`) implements these designs
for the E. coli cohort: random, cluster-grouped random, lineage, country,
project, geography and temporal splits with leakage checks, frozen by hash
before any model is fitted
([decision 0007](../docs/decisions/0007-ecoli-generalisation-benchmark.md),
[guide](../docs/ecoli-generalisation-benchmark.md)). Its summary tables are in
`cohorts/ecoli_research_v1/generalisation_v1/`.

The error decomposition (`bin/error_decomposition.py`, config
`cohorts/ecoli_research_v1/error_decomposition_v1.json`) reuses those frozen
predictions and splits. It audits whether the established determinants of
each endpoint are encoded by the benchmark's k-mer features, whether they
were seen in training, and how errors divide across these explanations; the
mechanism catalogue is `mappings/ecoli_amr_mechanisms_v1.tsv`
([decision 0008](../docs/decisions/0008-ecoli-error-decomposition.md),
[guide](../docs/ecoli-domain-shift-error-decomposition.md)). Its summary tables
are in `cohorts/ecoli_research_v1/error_decomposition_v1/`.
