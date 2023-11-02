#!/usr/bin/env Rscript

# Summary:
# This script performs the import of transcript quantification data, typically
# from RNA-Seq analysis tools like Salmon or Kallisto, into R using the
# tximport package. It then creates a SummarizedExperiment object containing
# count, abundance, and length data for further analysis. The script also
# provides functionality to summarize data at the gene level if a
# transcript-to-gene mapping is provided. It expects command line arguments for
# paths to data and configuration files and outputs various tables summarizing
# the imported data.

# Load required libraries
library(SummarizedExperiment)
library(tximport)

# Function to safely read a CSV file
read_csv_safe <- function(file, sep = "\t") {
  if (file.exists(file)) {
    read.csv(file, sep = sep, header = TRUE) 
  } else {
    message("File not available: ", file)
    return(NULL)
  }
}

# Function to build tables from SummarizedExperiment objects
build_table <- function(se, slot) {
  cbind(rowData(se)[, 1:2], assays(se)[[slot]])
}

# Function to get list of files based on quant_type
list_of_files <- function(path, quant_type) {
  pattern <- if (quant_type == "salmon") "quant.sf" else "abundance.h5"
  files <- list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  setNames(files, basename(dirname(files)))
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: tximport.r <coldata> <quant_out> <sample_name> <quant_type> <tx2gene>", 
       call. = FALSE)
}

# Assign arguments to variables
coldata_path <- args[1]
quant_out_path <- args[2]
sample_name <- args[3]
quant_type <- args[4]
tx2gene_path <- args[5]

# Prepare tx2gene if it exists and is not empty
tx2gene <- NULL
if (file.size(tx2gene_path) > 0) {
  tx2gene_df <- read_csv_safe(tx2gene_path)
  if (!is.null(tx2gene_df)) {
    tx2gene <- tx2gene_df[, 1:2]
    colnames(tx2gene) <- c("tx", "gene_id")
  }
}

# Process filenames based on quant_type
quant_files <- list_of_files(quant_out_path, quant_type)
coldata <- read_csv_safe(coldata_path)
if (!is.null(coldata)) {
  coldata <- cbind(files = quant_files, coldata[match(names(quant_files), coldata[,1]),])
} else {
  coldata <- data.frame(files = quant_files, names = names(quant_files))
}

# Import transcript quantifications
txi <- tximport(quant_files, type = quant_type, txOut = TRUE)

# Prepare rowData and colData for SummarizedExperiment
colData <- DataFrame(coldata)
rowData <- if (!is.null(tx2gene)) {
  tx2gene[match(rownames(txi$counts), tx2gene$tx), , drop = FALSE]
} else {
  data.frame(tx = rownames(txi$counts))
}
rownames(rowData) <- rowData$tx

# Create SummarizedExperiment object
se <- SummarizedExperiment(assays = txi, colData = colData, rowData = rowData)

# Write output files function
output_files <- function(se, prefix, suffix) {
  file_name_abundance <- sprintf("%s_%s_abundance.tsv", prefix, suffix)
  file_name_counts <- sprintf("%s_%s_counts.tsv", prefix, suffix)
  write.table(build_table(se, "abundance"), file_name_abundance, 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(build_table(se, "counts"), file_name_counts, 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Process and write gene level data if tx2gene is provided
if (!is.null(tx2gene)) {
  se_gene <- summarizeToGene(se, tx2gene = tx2gene)
  output_files(se_gene, sample_name, "gene")
}

# Write transcript level data
output_files(se, sample_name, "transcript")

# Print session information to standard out
citation("tximeta")
sessionInfo()

