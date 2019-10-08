#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: tximeta.r <coldata> <salmon_out>", call.=FALSE)
}

coldata = args[1]
counts_fn = args[2]
tpm_fn = args[3]

tx2gene = "tx2gene.csv"
info = file.info(tx2gene)
if (info$size == 0){
  tx2gene = NULL
}else{
  rowdata = read.csv(tx2gene, header = FALSE)
  colnames(rowdata) = c("tx", "gene_id", "gene_name")
  tx2gene = rowdata[,1:2]
}

counts = read.csv(counts_fn, row.names=1)
tpm = read.csv(tpm_fn, row.names=1)

if (length(intersect(rownames(counts), rowdata[["tx"]])) > length(intersect(rownames(counts), rowdata[["gene_id"]]))){
    by_what = "tx"
} else {
    by_what = "gene_id"
    rowdata = unique(rowdata[,2:3])
}

if (file.exists(coldata)){
    coldata = read.csv(coldata)
    coldata = coldata[match(colnames(counts), coldata[,1]),]
    coldata = cbind(files = fns, coldata)
}else{
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = colnames(counts), names = colnames(counts))
}
library(SummarizedExperiment)

rownames(coldata) = coldata[["names"]]
extra = setdiff(rownames(counts),  as.character(rowdata[[by_what]]))
if (length(extra) > 0){
    rowdata = rbind(rowdata,
                    data.frame(tx=extra,
                               gene_id=extra,
                               gene_name=extra))
}

rowdata = rowdata[match(rownames(counts), as.character(rowdata[[by_what]])),]
rownames(rowdata) = rowdata[[by_what]]
se = SummarizedExperiment(assays = list(counts = counts,
                                        abundance = tpm),
                          colData = DataFrame(coldata),
                          rowData = rowdata)

saveRDS(se, file = paste0(tools::file_path_sans_ext(counts_fn), ".rds"))
