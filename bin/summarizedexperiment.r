#!/usr/bin/env Rscript

# Written by Lorena Pantano and released under the MIT license.

library(SummarizedExperiment)

## Create SummarizedExperiment (se) object from counts

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: summarizedexperiment.r <coldata> <counts> <tpm> <tx2gene>", call. = FALSE)
}

coldata <- args[1]
counts_fn <- args[2]
tpm_fn <- args[3]
tx2gene <- args[4]

info <- file.info(tx2gene)
if (info$size == 0) {
    tx2gene <- NULL
} else {
    rowdata <- read.csv(tx2gene, sep = "\t", header = FALSE)
    colnames(rowdata) <- c("tx", "gene_id", "gene_name")
    tx2gene <- rowdata[, 1:2]
}

counts <- read.csv(counts_fn, row.names = 1, sep = "\t")
counts <- counts[, 2:ncol(counts), drop = FALSE] # remove gene_name column
tpm <- read.csv(tpm_fn, row.names = 1, sep = "\t")
tpm <- tpm[, 2:ncol(tpm), drop = FALSE] # remove gene_name column

if (length(intersect(rownames(counts), rowdata[["tx"]])) > length(intersect(rownames(counts), rowdata[["gene_id"]]))) {
    by_what <- "tx"
} else {
    by_what <- "gene_id"
    rowdata <- unique(rowdata[, 2:3])
}

if (file.exists(coldata)) {
    coldata <- read.csv(coldata, sep = "\t")
    coldata <- coldata[match(colnames(counts), coldata[, 1]), ]
    coldata <- cbind(files = fns, coldata)
} else {
    message("ColData not avaliable ", coldata)
    coldata <- data.frame(files = colnames(counts), names = colnames(counts))
}

rownames(coldata) <- coldata[["names"]]
extra <- setdiff(rownames(counts), as.character(rowdata[[by_what]]))
if (length(extra) > 0) {
    rowdata <- rbind(
        rowdata,
        data.frame(
            tx = extra,
            gene_id = extra,
            gene_name = extra
        )[, colnames(rowdata)]
    )
}

rowdata <- rowdata[match(rownames(counts), as.character(rowdata[[by_what]])), ]
rownames(rowdata) <- rowdata[[by_what]]
se <- SummarizedExperiment(
    assays = list(counts = counts, abundance = tpm),
    colData = DataFrame(coldata),
    rowData = rowdata
)

saveRDS(se, file = paste0(tools::file_path_sans_ext(counts_fn), ".rds"))
