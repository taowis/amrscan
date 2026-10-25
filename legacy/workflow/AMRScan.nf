#!/usr/bin/env nextflow

/*
 * LEGACY CARD/BLAST workflow, retained for historical reproduction.
 * New modular entry point: ../main.nf. Highest bitscore is a display ranking,
 * not proof of biological function or phenotypic resistance.
 * AMRScan - Antimicrobial Resistance Gene Scanner
 * Nextflow implementation of the R-based AMR detection pipeline
 */

nextflow.enable.dsl = 2

// Parameters
params.input = null
params.outdir = "results"
params.card_db = "db/protein_fasta_protein_homolog_model.fasta"
params.evalue = "1e-5"
params.threads = 4
params.min_quality = 20

// Validation
if (!params.input) {
    error "Please provide an input file with --input"
}

log.info """
==============================================
AMRScan - Antimicrobial Resistance Scanner
==============================================
Input file        : ${params.input}
Output directory   : ${params.outdir}
CARD database      : ${params.card_db}
E-value threshold  : ${params.evalue}
Threads            : ${params.threads}
Min quality        : ${params.min_quality}
==============================================
"""

/*
 * Process 1: Quality control and format conversion
 */
process PREPROCESS_READS {
    tag "${sample_id}"
    publishDir "${params.outdir}/preprocessed", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file)

    output:
    tuple val(sample_id), path("*.fasta"), emit: fasta
    path("*.log"), emit: log

    script:
    def file_ext = input_file.getExtension()

    if (file_ext == 'fastq' || file_ext == 'fq') {
        """
        #!/usr/bin/env Rscript

        library(Biostrings)
        library(ShortRead)

        cat("Processing FASTQ file: ${input_file}\\n")

        # Read FASTQ
        fq <- readFastq("${input_file}")

        # Quality filtering
        fq_trimmed <- fq[quality(fq) >= ${params.min_quality}]

        # Convert to FASTA
        seqs <- sread(fq_trimmed)
        names(seqs) <- as.character(id(fq_trimmed))

        # Write FASTA
        writeXStringSet(DNAStringSet(seqs), filepath = "${sample_id}.fasta")

        # Log stats
        cat("Original reads:", length(fq), "\\n", file = "${sample_id}_preprocess.log")
        cat("Filtered reads:", length(fq_trimmed), "\\n", file = "${sample_id}_preprocess.log", append = TRUE)
        cat("Reads retained:", round(length(fq_trimmed)/length(fq)*100, 2), "%\\n", file = "${sample_id}_preprocess.log", append = TRUE)
        """
    } else if (file_ext in ['fa', 'fna', 'fasta']) {
        """
        # Input is already FASTA, just copy
        cp ${input_file} ${sample_id}.fasta
        echo "Input was already FASTA format" > ${sample_id}_preprocess.log
        """
    } else {
        error "Unsupported file format: ${file_ext}"
    }
}

/*
 * Process 2: Prepare BLAST database
 */
process MAKEBLASTDB {
    tag "CARD_database"
    storeDir "work/blast_db"

    input:
    path card_fasta

    output:
    path "card_blastdb*", emit: blast_db

    script:
    """
    makeblastdb -in ${card_fasta} -dbtype prot -out card_blastdb
    """
}

/*
 * Process 3: Run BLAST search
 */
process BLAST_SEARCH {
    tag "${sample_id}"
    publishDir "${params.outdir}/blast", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)
    path blast_db

    output:
    tuple val(sample_id), path("*_blast_results.tsv"), emit: blast_results

    script:
    """
    blastx -query ${fasta} \\
           -db card_blastdb \\
           -out ${sample_id}_blast_results.tsv \\
           -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \\
           -evalue ${params.evalue} \\
           -num_threads ${params.threads}
    """
}

/*
 * Process 4: Parse and summarize results
 */
process PARSE_RESULTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    tuple val(sample_id), path(blast_results)

    output:
    path "*_AMR_hits_summary.csv", emit: summary
    path "*_AMR_hits_detailed.csv", emit: detailed

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(magrittr)

    # Read BLAST results
    if (file.size("${blast_results}") > 0) {
        blast_data <- read.delim("${blast_results}", header = FALSE, stringsAsFactors = FALSE)
        colnames(blast_data) <- c("Query", "Subject", "Identity", "Length", "Evalue", "Bitscore", "Annotation")

        # Get top hits per query
        top_hits <- blast_data %>%
            group_by(Query) %>%
            top_n(1, Bitscore) %>%
            ungroup()

        # Write summary
        write.csv(top_hits, file = "${sample_id}_AMR_hits_summary.csv", row.names = FALSE)

        # Write detailed results
        write.csv(blast_data, file = "${sample_id}_AMR_hits_detailed.csv", row.names = FALSE)

        cat("Found", nrow(top_hits), "AMR hits in", length(unique(blast_data\$Query)), "queries\\n")

    } else {
        # No hits found
        empty_df <- data.frame(
            Query = character(0),
            Subject = character(0),
            Identity = numeric(0),
            Length = integer(0),
            Evalue = numeric(0),
            Bitscore = numeric(0),
            Annotation = character(0)
        )

        write.csv(empty_df, file = "${sample_id}_AMR_hits_summary.csv", row.names = FALSE)
        write.csv(empty_df, file = "${sample_id}_AMR_hits_detailed.csv", row.names = FALSE)

        cat("No AMR hits found\\n")
    }
    """
}

/*
 * Process 5: Generate final report
 */
process GENERATE_REPORT {
    tag "${sample_id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(summary), path(detailed), path(preprocess_log)

    output:
    path "*_AMR_report.html", emit: report

    script:
    """
    #!/usr/bin/env Rscript

    library(knitr)
    library(rmarkdown)

    # Create R Markdown report
    rmd_content <- '
---
title: "AMR Scan Report"
subtitle: "Sample: ${sample_id}"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk\$set(echo = FALSE, warning = FALSE, message = FALSE)
library(DT)
library(ggplot2)
library(dplyr)
```

## Summary

This report summarizes the antimicrobial resistance gene detection results for sample **${sample_id}**.

```{r read_data}
# Read results
summary_data <- read.csv("${summary}")
detailed_data <- read.csv("${detailed}")

# Read preprocessing log
if (file.exists("${preprocess_log}")) {
    preprocess_info <- readLines("${preprocess_log}")
} else {
    preprocess_info <- "No preprocessing information available"
}
```

## Preprocessing Information

```{r preprocessing}
cat(paste(preprocess_info, collapse = "\\n"))
```

## AMR Detection Results

### Top Hits Summary

```{r summary_table}
if (nrow(summary_data) > 0) {
    DT::datatable(summary_data, options = list(scrollX = TRUE))
} else {
    cat("No AMR genes detected in this sample.")
}
```

### Statistics

```{r stats}
if (nrow(summary_data) > 0) {
    cat("Total AMR hits:", nrow(summary_data), "\\n")
    cat("Unique queries with hits:", length(unique(summary_data\$Query)), "\\n")
    cat("Average identity:", round(mean(summary_data\$Identity), 2), "%\\n")
    cat("Average bit score:", round(mean(summary_data\$Bitscore), 2), "\\n")
} else {
    cat("No statistics available - no AMR genes detected.")
}
```

### Identity Distribution

```{r identity_plot}
if (nrow(summary_data) > 0) {
    ggplot(summary_data, aes(x = Identity)) +
        geom_histogram(binwidth = 5, fill = "steelblue", alpha = 0.7) +
        labs(title = "Distribution of Sequence Identity",
             x = "Identity (%)",
             y = "Count") +
        theme_minimal()
} else {
    cat("No plot available - no AMR genes detected.")
}
```

### Detailed Results

```{r detailed_table}
if (nrow(detailed_data) > 0) {
    DT::datatable(detailed_data, options = list(scrollX = TRUE))
} else {
    cat("No detailed results available.")
}
```

## Parameters Used

- E-value threshold: ${params.evalue}
- Minimum quality score: ${params.min_quality}
- Number of threads: ${params.threads}
- CARD database: ${params.card_db}

'

    # Write the R Markdown content to a file
    writeLines(rmd_content, "${sample_id}_report.Rmd")

    # Render the report
    rmarkdown::render("${sample_id}_report.Rmd", output_file = "${sample_id}_AMR_report.html")
    """
}

/*
 * Main workflow
 */
workflow {

    // Create input channel
    input_ch = Channel.fromPath(params.input)
        .map { file ->
            def sample_id = file.getBaseName()
            return tuple(sample_id, file)
        }

    // Preprocess reads
    PREPROCESS_READS(input_ch)

    // Prepare BLAST database
    card_db = Channel.fromPath(params.card_db)
    MAKEBLASTDB(card_db)

    // Run BLAST
    BLAST_SEARCH(
        PREPROCESS_READS.out.fasta,
        MAKEBLASTDB.out.blast_db.collect()
    )

    // Parse results
    PARSE_RESULTS(BLAST_SEARCH.out.blast_results)

    // Generate report
    report_input = PARSE_RESULTS.out.summary
        .join(PARSE_RESULTS.out.detailed)
        .join(PREPROCESS_READS.out.log)

    // GENERATE_REPORT(report_input)
}

/*
 * Completion message
 */
workflow.onComplete {
    log.info """
    ==============================================
    AMRScan workflow completed!
    ==============================================
    Results directory: ${params.outdir}
    Success: ${workflow.success}
    Duration: ${workflow.duration}
    ==============================================
    """
}
