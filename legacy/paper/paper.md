---
title: 'AMRScan: A hybrid R and Nextflow toolkit for rapid antimicrobial resistance gene detection from sequencing data'
tags:
  - R
  - Nextflow
  - bioinformatics
  - pathogen genomics
  - antimicrobial resistance
  - reproducibility
authors:
  - name: Kaitao Lai
    orcid: 0000-0002-9420-9352
    affiliation: 1
affiliations:
  - name: University of Sydney
    index: 1
date: 2025-07-08
bibliography: paper.bib
---

# Summary

**AMRScan** is a hybrid bioinformatics toolkit implemented in both R and [Nextflow](https://www.nextflow.io/) for the rapid and reproducible detection of antimicrobial resistance (AMR) genes from next-generation sequencing (NGS) data. The toolkit enables users to identify AMR gene hits in sequencing reads by aligning them against reference databases such as CARD using BLAST [@altschul1990blast].

The R implementation provides a concise, script-based approach suitable for single-sample analysis, teaching, and rapid prototyping. In contrast, the Nextflow implementation enables reproducible, scalable workflows for multi-sample batch processing in high-performance computing (HPC) and containerized environments. It leverages modular pipeline design with support for automated database setup, quality control, conversion, BLAST alignment, and results parsing.

AMRScan helps bridge the gap between lightweight exploratory analysis and production-ready surveillance pipelines, making it suitable for both research and public health genomics applications.

# Statement of Need

While several large-scale AMR detection platforms such as ResFinder [@zankari2012resfinder] exist, many are resource-intensive or require complex installations. AMRScan addresses the need for a minimal, transparent, and reproducible toolkit that can be used flexibly in small labs, clinical settings, or large-scale surveillance workflows.

The inclusion of a pure Nextflow implementation enables high-throughput, multi-sample analyses in cloud and HPC environments, while the standalone R script remains accessible to users in data science, microbiology, and epidemiology. Both versions use shared components (e.g., a BLAST parsing script) to ensure consistency and reproducibility of results.

# Usage Guidance

AMRScan provides two usage modes, tailored to user needs:

- **R script (`AMRScan_standalone.R`)**: Best suited for small datasets, single-sample analysis, quick local tests, educational use, and lightweight environments without workflow managers.

- **Nextflow workflow (`main.nf`)**: Designed for large-scale, automated analyses, this version excels in multi-sample settings, HPC/cloud infrastructure, and environments where reproducibility, parallelism, and containerization are priorities.

This flexible dual-mode implementation ensures that AMRScan can serve both teaching/demo scenarios and production-grade bioinformatics pipelines.

# Implementation

- The R script `scripts/AMRScan_standalone.R` encapsulates the entire pipeline in a linear script-based format.
- The Nextflow workflow `workflow/main.nf` organizes the same logic into modular processes:
  - `DownloadCARD`, `MakeBLASTdb`, `ConvertFASTQ`, `RunBLAST`, and `ParseResults`
- The shared R script `scripts/parse_blast.R` is used for post-BLAST result summarization.
- Both implementations are documented, testable, and validated using mock NGS input.

# Example Dataset and Demonstration

The example data used to validate AMRScan was obtained from a study by Munim et al. (2024) on multidrug-resistant *Klebsiella pneumoniae* isolated from poultry in Noakhali, Bangladesh [@munim2024mdr]. The assembled genome was downloaded from GenBank ([GCA_037966445.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_037966445.1)).

For antimicrobial resistance gene detection, we used the protein homolog model from the Comprehensive Antibiotic Resistance Database (CARD) [@jia2017card], version Broadstreet v4.0.1, available at:  
[https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2](https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2)

A sample output summary is shown below:

## Top AMR Hits Summary

| Query                       | Subject                                 | Identity | Length | Evalue | Bitscore | Annotation                                               |
|----------------------------|-----------------------------------------|----------|--------|--------|----------|----------------------------------------------------------|
| NZ_JBBPBW010000028.1       | gb\|BAM10414.1\|ARO:3003039\|OprA     | 40.839   | 453    | 0      | 252      | OprA \[Pseudomonas aeruginosa\]                       |
| NZ_JBBPBW010000035.1       | gb\|ABV89601.1\|ARO:3004826\|LAP-2    | 100.000  | 285    | 0      | 587      | LAP-2 \[Enterobacter cloacae\]                        |
| NZ_JBBPBW010000001.1       | gb\|CAA66729.1\|ARO:3001070\|SHV-11   | 100.000  | 286    | 0      | 581      | SHV-11 \[Klebsiella pneumoniae\]                      |
| NZ_JBBPBW010000010.1       | gb\|CCI79240.1\|ARO:3005047\|eptB     | 99.303   | 574    | 0      | 1109     | eptB \[Klebsiella pneumoniae subsp. rhinoscleromatis\]|
| NZ_JBBPBW010000104.1       | gb\|AHK10285.1\|ARO:3002859\|dfrA14   | 98.726   | 157    | 0      | 327      | dfrA14 \[Escherichia coli\]                          |
| NZ_JBBPBW010000011.1       | gb\|AAC75271.1\|ARO:3003952\|YojI     | 83.912   | 547    | 0      | 885      | YojI \[Escherichia coli str. K-12 substr. MG1655\]    |

# Software Repository

The source code for AMRScan is freely available on GitHub at:  
[https://github.com/biosciences/AMRScan](https://github.com/biosciences/AMRScan)

# Acknowledgements

The author thanks collaborators at University of Sydney for feedback on early concepts.

# References
