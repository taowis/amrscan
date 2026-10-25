# Minimal test-run.R for AMRScan.R

cat("Running AMRScan.R test...\n")

# Run the standalone script
source("scripts/AMRScan.R")

# Check if output file is generated
output_file <- "results/final/GCF_037966445.1_ASM3796644v1_genomic_AMR_hits_summary.csv"
if (file.exists(output_file)) {
  cat("✅ Test passed: Output file generated at", output_file, "\n")
} else {
  stop("❌ Test failed: Output file not found.")
}
