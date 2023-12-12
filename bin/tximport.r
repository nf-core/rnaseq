#!/usr/bin/env Rscript

# Script for importing and processing transcript-level quantifications.
# Written by Lorena Pantano, later modified by Jonathan Manning, and released
# under the MIT license.

# Loading required libraries
library(SummarizedExperiment)
library(tximport)

# Parsing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: tximport.r <coldata_path> <path> <prefix> <quant_type> <tx2gene_path>",
        call.=FALSE)
}

# Assigning command line arguments to variables
coldata_path <- args[1]
path <- args[2]
prefix <- args[3]
quant_type <- args[4]
tx2gene_path <- args[5]

## Functions

# Build a table from a SummarizedExperiment object
build_table <- function(se.obj, slot) {
    cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]])
}

# Write a table to a file with given parameters
write_se_table <- function(params) {
    file_name <- paste0(prefix, ".", params$suffix)
    write.table(build_table(params$obj, params$slot), file_name,
                sep="\t", quote=FALSE, row.names = FALSE)
}

# Read transcript metadata from a given path
read_transcript_info <- function(tinfo_path){
    info <- file.info(tinfo_path)
    if (info$size == 0) {
        stop("tx2gene file is empty")
    }

    transcript_info <- read.csv(tinfo_path, sep="\t", header = FALSE,
                                col.names = c("tx", "gene_id", "gene_name"))

    extra <- setdiff(rownames(txi[[1]]), as.character(transcript_info[["tx"]]))
    transcript_info <- rbind(transcript_info, data.frame(tx=extra, gene_id=extra, gene_name=extra))
    transcript_info <- transcript_info[match(rownames(txi[[1]]), transcript_info[["tx"]]), ]
    rownames(transcript_info) <- transcript_info[["tx"]]

    list(transcript = transcript_info,
        gene = unique(transcript_info[,2:3]),
        tx2gene = transcript_info[,1:2])
}

# Read and process sample/column data from a given path
read_coldata <- function(coldata_path){
    if (file.exists(coldata_path)) {
        coldata <- read.csv(coldata_path, sep="\t")
        coldata <- coldata[match(names, coldata[,1]),]
        coldata <- cbind(files = fns, coldata)
    } else {
        message("ColData not available: ", coldata_path)
        coldata <- data.frame(files = fns, names = names)
    }
    rownames(coldata) <- coldata[["names"]]
}

# Create a SummarizedExperiment object with given data
create_summarized_experiment <- function(counts, abundance, length, col_data, row_data) {
    SummarizedExperiment(assays = list(counts = counts, abundance = abundance, length = length),
        colData = col_data,
        rowData = row_data)
}

# Main script starts here

# Define pattern for file names based on quantification type
pattern <- ifelse(quant_type == "kallisto", "abundance.tsv", "quant.sf")
fns <- list.files(path, pattern = pattern, recursive = T, full.names = T)
names <- basename(dirname(fns))
names(fns) <- names
dropInfReps <- quant_type == "kallisto"

# Import transcript-level quantifications
txi <- tximport(fns, type = quant_type, txOut = TRUE, dropInfReps = dropInfReps)

# Read transcript and sample data
transcript_info <- read_transcript_info(tx2gene_path)
coldata <- read_coldata(coldata_path)

# Create initial SummarizedExperiment object
se <- create_summarized_experiment(txi[["counts"]], txi[["abundance"]], txi[["length"]],
    DataFrame(coldata), transcript_info$transcript)

# Setting parameters for writing tables
params <- list(
    list(obj = se, slot = "abundance", suffix = "transcript_tpm.tsv"),
    list(obj = se, slot = "counts", suffix = "transcript_counts.tsv"),
    list(obj = se, slot = "length", suffix = "transcript_lengths.tsv")
)

# Process gene-level data if tx2gene mapping is available
if ("tx2gene" %in% names(transcript_info) && !is.null(transcript_info$tx2gene)) {
    tx2gene <- transcript_info$tx2gene
    gi <- summarizeToGene(txi, tx2gene = tx2gene)
    gi.ls <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    gi.s <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "scaledTPM")

    gene_info <- transcript_info$gene[match(rownames(gi[[1]]), transcript_info$gene[["gene_id"]]),]
    rownames(gene_info) <- gene_info[["tx"]]

    col_data_frame <- DataFrame(coldata)

    # Create gene-level SummarizedExperiment objects
    gse <- create_summarized_experiment(gi[["counts"]], gi[["abundance"]], gi[["length"]],
        col_data_frame, gene_info)
    gse.ls <- create_summarized_experiment(gi.ls[["counts"]], gi.ls[["abundance"]], gi.ls[["length"]],
        col_data_frame, gene_info)
    gse.s <- create_summarized_experiment(gi.s[["counts"]], gi.s[["abundance"]], gi.s[["length"]],
        col_data_frame, gene_info)

    params <- c(params, list(
        list(obj = gse, slot = "length", suffix = "gene_lengths.tsv"),
        list(obj = gse, slot = "abundance", suffix = "gene_tpm.tsv"),
        list(obj = gse, slot = "counts", suffix = "gene_counts.tsv"),
        list(obj = gse.ls, slot = "abundance", suffix = "gene_tpm_length_scaled.tsv"),
        list(obj = gse.ls, slot = "counts", suffix = "gene_counts_length_scaled.tsv"),
        list(obj = gse.s, slot = "abundance", suffix = "gene_tpm_scaled.tsv"),
        list(obj = gse.s, slot = "counts", suffix = "gene_counts_scaled.tsv")
    ))
}

# Writing tables for each set of parameters
done <- lapply(params, write_se_table)

# Output session information and citations
citation("tximeta")
sessionInfo()

