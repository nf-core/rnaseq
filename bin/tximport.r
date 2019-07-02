#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: tximeta.r <coldata> <salmon_out>", call.=FALSE)
}

path = args[2]
coldata = args[1]

sample_name = args[3]

prefix = paste(c(sample_name, "salmon"), sep="_")

tx2gene = "tx2gene.csv"
info = file.info(tx2gene)
if (info$size == 0){
  tx2gene = NULL
}else{
  rowdata = read.csv(tx2gene, header = FALSE)
  colnames(rowdata) = c("tx", "gene_id", "gene_name")
  tx2gene = rowdata[,1:2]
}

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names = basename(dirname(fns))
names(fns) = names

if (file.exists(coldata)){
    coldata = read.csv(coldata)
    coldata = coldata[match(names, coldata[,1]),]
    coldata = cbind(files = fns, coldata)
}else{
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = fns, names = names)
}

library(SummarizedExperiment)
library(tximport)

txi = tximport(fns, type = "salmon", txOut = TRUE)
rownames(coldata) = coldata[["names"]]
rowdata = rowdata[match(rownames(txi[[1]]), rowdata[["tx"]]),]
se = SummarizedExperiment(assays = list(counts = txi[["counts"]],
                                        abundance = txi[["abundance"]],
                                        length = txi[["length"]]),
                          colData = DataFrame(coldata),
                          rowData = rowdata)
if (!is.null(tx2gene)){
  gi = summarizeToGene(txi, tx2gene = tx2gene)
  growdata = unique(rowdata[,2:3])
  growdata = growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
  gse = SummarizedExperiment(assays = list(counts = gi[["counts"]],
                                          abundance = gi[["abundance"]],
                                          length = gi[["length"]]),
                            colData = DataFrame(coldata),
                            rowData = growdata)
}

if(exists("gse")){
  saveRDS(gse, file = "gse.rds")
  write.csv(assays(gse)[["abundance"]], paste(c(prefix, "gene_tpm.csv"), collapse="_"), quote=FALSE)
  write.csv(assays(gse)[["counts"]], paste(c(prefix, "gene_counts.csv"), collapse="_"), quote=FALSE)
}

saveRDS(se, file = "se.rds")
write.csv(assays(se)[["abundance"]], paste(c(prefix, "transcript_tpm.csv"), collapse="_"), quote=FALSE)
write.csv(assays(se)[["counts"]], paste(c(prefix, "transcript_counts.csv"), collapse="_"), quote=FALSE)

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()
