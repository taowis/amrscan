# AMRScan - Antimicrobial Resistance Gene Scanner in R
# LEGACY demonstration: known fixed-path/output-path issues; see README.md.
# A BLAST match or top bitscore does not establish phenotypic resistance.
# Implements read QC, FASTA conversion, BLAST, and AMR hit parsing.

# Load required libraries
library(Biostrings)    # for DNA sequence handling
library(ShortRead)     # for FASTQ parsing
library(RCurl)         # for downloading reference data
library(dplyr)
library(magrittr)  # for %>% pipe

# Set paths
# Run from the repository root.
input_file <- "data/GCF_037966445.1_ASM3796644v1_genomic.fna"
sample_id <- tools::file_path_sans_ext(basename(input_file))
output_dir <- "results/"
blast_dir <- file.path(output_dir, "blast")
summary_dir <- file.path(output_dir, "final")

if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(blast_dir)) dir.create(blast_dir)
output_dir <- "results/"

# Create output directory if it doesn’t exist
if (!dir.exists(output_dir)) dir.create(output_dir)

# Step 1: Preprocess FASTQ or FASTA
# Determine input file type
file_ext <- tools::file_ext(input_file)
input_fasta <- NULL

if (file_ext == "fastq") {
  cat("Detected FASTQ input. Performing quality filtering and conversion to FASTA...\n")
  fq <- readFastq(input_file)
  fq_trimmed <- fq[quality(fq) >= 20]
  writeFastq(fq_trimmed, file = paste0(output_dir, "cleaned.fastq"))

  seqs <- sread(fq_trimmed)
  names(seqs) <- as.character(id(fq_trimmed))
  input_fasta <- paste0(output_dir, "cleaned.fasta")
  writeXStringSet(DNAStringSet(seqs), filepath = input_fasta)

} else if (file_ext %in% c("fa", "fna", "fasta")) {
  cat("Detected FASTA input. Skipping quality filtering.\n")
  input_fasta <- input_file
} else {
  stop("Unsupported input format. Please provide a .fastq or .fasta file.")
}

# Step 2: Prepare CARD BLAST database
db_fasta <- "db/protein_fasta_protein_homolog_model.fasta"  # use local copy
blast_db <- "db/card_blastdb"  # store BLAST DB files in data/

if (!file.exists(paste0(blast_db, ".pin"))) {  # check if BLAST DB is built
  cat("Preparing BLAST database from local CARD FASTA file...\n")
  system(paste("makeblastdb -in", db_fasta, "-dbtype prot -out", blast_db))
} else {
  cat("BLAST database already exists. Skipping makeblastdb.\n")
}

# Step 3: Run BLAST
cat("Running BLASTX on input FASTA:", input_fasta, "\n")
blast_output_file <- file.path(blast_dir, paste0(sample_id, "_blast_results.tsv"))
system(paste("blastx -query", input_fasta,
             "-db", blast_db,
             "-out", paste0(output_dir, "blast_results.tsv"),
             "-outfmt '6 qseqid sseqid pident length evalue bitscore stitle'",
             "-evalue 1e-5 -num_threads 4"))

# Step 4: Parse Results
blast_results <- read.delim(blast_output_file, header = FALSE)
colnames(blast_results) <- c("Query", "Subject", "Identity", "Length", "Evalue", "Bitscore", "Annotation")

top_hits <- blast_results %>%
  group_by(Query) %>%
  top_n(1, Bitscore)

summary_file <- file.path(summary_dir, paste0(sample_id, "_AMR_hits_summary.csv"))
write.csv(top_hits, file = summary_file, row.names = FALSE)

cat("AMRScan completed. Results saved to:", output_dir, "\n")
