#!/usr/bin/env Rscript

# Written by Lorena Pantano and released under the MIT license.

library(SummarizedExperiment)
library(tximport)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("Usage: tximport.r <coldata_path> <path> <prefix> <quant_type> <tx2gene_path>", call.=FALSE)
}

coldata_path = args[1]
path = args[2]
prefix = args[3]
quant_type = args[4]
tx2gene_path = args[5]

## Functions 

# Function to build a table from a SummarizedExperiment object
build_table <- function(se.obj, slot) {
    cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]])
}

# Function to write a table to a file
write_se_table <- function(params) {
    file_name <- paste0(prefix, ".", params$suffix)
    write.table(build_table(params$obj, params$slot), file_name, sep="\t", quote=FALSE, row.names = FALSE)
}

# Read and process transcript meta
read_transcript_info <- function(tinfo_path){

    info = file.info(tinfo_path)
    if (info$size == 0) {
        stop("tx2gene file is empty")
    } 

    # Read transcript information and set column names
    transcript_info = read.csv(
        tinfo_path, 
        sep="\t", 
        header = FALSE, 
        col.names = c("tx", "gene_id", "gene_name")
    )

    # Identify and add extra transcripts
    extra = setdiff(rownames(txi[[1]]), as.character(transcript_info[["tx"]]))
    transcript_info = rbind(transcript_info, data.frame(tx=extra, gene_id=extra, gene_name=extra))

    # Reorder rows to match the order in txi and set row names
    transcript_info = transcript_info[match(rownames(txi[[1]]), transcript_info[["tx"]]), ]
    rownames(transcript_info) = transcript_info[["tx"]]

    # Return various formulations of the info for downstream use
    list(
        transcript = transcript_info,
	gene = unique(transcript_info[,2:3]),
	tx2gene = transcript_info[,1:2]
    )
}

# Read sample/column data
read_coldata <- function(coldata_path){

    if (file.exists(coldata_path)) {
        coldata = read.csv(coldata_path, sep="\t")
        coldata = coldata[match(names, coldata[,1]),]
        coldata = cbind(files = fns, coldata)
    } else {
        message("ColData not available: ", coldata_path)
        coldata = data.frame(files = fns, names = names)
    }

    rownames(coldata) = coldata[["names"]]
}

# Function to create a SummarizedExperiment object
create_summarized_experiment <- function(counts, abundance, length, col_data, row_data) {
    return(SummarizedExperiment(assays = list(counts = counts, abundance = abundance, length = length),
                                colData = col_data,
                                rowData = row_data))
}

# Read the abundance values
pattern <- ifelse(quant_type == "kallisto", "abundance.tsv", "quant.sf")
fns = list.files(path, pattern = pattern, recursive = T, full.names = T)
names = basename(dirname(fns))
names(fns) = names
dropInfReps = quant_type == "kallisto"
txi = tximport(fns, type = quant_type, txOut = TRUE, dropInfReps = dropInfReps)

# Read transcript metadata
transcript_info = read_transcript_info(tx2gene_path)

# Read sample/ column data
coldata = read_coldata(coldata_path)

# Process input data to summarized experiment objects

# Create initial SummarizedExperiment object
se <- create_summarized_experiment(txi[["counts"]], txi[["abundance"]], txi[["length"]],
                                   DataFrame(coldata), transcript_info$transcript)

# Parameters for writing tables
params <- list(
    list(obj = se, slot = "abundance", suffix = "transcript_tpm.tsv"),
    list(obj = se, slot = "counts", suffix = "transcript_counts.tsv"),
    list(obj = se, slot = "length", suffix = "transcript_lengths.tsv")
)

# Check if tx2gene is present in transcript_info
if ("tx2gene" %in% names(transcript_info) && !is.null(transcript_info$tx2gene)) {
    tx2gene <- transcript_info$tx2gene
    gi <- summarizeToGene(txi, tx2gene = tx2gene)
    gi.ls <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    gi.s <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "scaledTPM")

    # Get gene information matching the rows in gi
    gene_info <- transcript_info$gene[match(rownames(gi[[1]]), transcript_info$gene[["gene_id"]]),]
    rownames(gene_info) <- gene_info[["tx"]]

    # Create dataframe only once for reuse
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

# Apply the write_se_table function to each set of parameters
lapply(params, write_se_table)

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()

