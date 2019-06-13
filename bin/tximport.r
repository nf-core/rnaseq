#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: tximeta.r <coldata> <salmon_out>", call.=FALSE)
}

path = args[2]
coldata = args[1]

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
coldata = list.files(coldata, full.names = TRUE)
if (length(coldata)==0){
  coldata = "NULL"
}
if (file.exists(coldata)){
    coldata = read.csv(coldata)
    coldata = coldata[match(names, coldata[,1]),]
    coldata = cbind(files = fns, coldata)
}else{
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = fns, names = names)
}

library(SummarizedExperiment)

# if not genome version is giving
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
  write.csv(assays(se)[["abundance"]], "merged_salmon_gene_tpm.csv")
  write.csv(assays(se)[["counts"]], "merged_salmon_gene_reads.csv")
}

saveRDS(se, file = "se.rds")
write.csv(assays(se)[["abundance"]], "merged_salmon_tx_tpm.csv")
write.csv(assays(se)[["counts"]], "merged_salmon_tx_reads.csv")

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()