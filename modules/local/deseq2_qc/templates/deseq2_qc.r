#!/usr/bin/env Rscript

# Written by Harshil Patel and Gavin Kelly and released under the MIT license.

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.
## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE OPTIONS                              ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

# count_file:    count file matrix where rows are genes and columns are samples
# count_col:     first column containing sample count data
# id_col:        column containing identifiers to be used
# sample_suffix: suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'
# vst:           run vst transform instead of rlog
opt <- list(
    count_file    = '$sample.counts_gene_length_scaled',
    count_col     = 3L,
    id_col        = 1L,
    sample_suffix = '',
    outdir        = './',
    outprefix     = ifelse('$task.ext.prefix' == 'null', 'deseq2', '$task.ext.prefix'),
    vst           = FALSE,
    cores         = as.integer($task.cpus)
)
opt_types <- lapply(opt, class)

# Apply option overrides from task.ext.args
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }
    # Preserve classes from defaults
    args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    opt[[ao]] <- args_opt[[ao]]
}

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

count.table           <- read.delim(file=opt\$count_file,header=TRUE, row.names=NULL, check.names=FALSE)
rownames(count.table) <- count.table[,opt\$id_col]
count.table           <- count.table[,opt\$count_col:ncol(count.table),drop=FALSE]
colnames(count.table) <- gsub(opt\$sample_suffix,"",colnames(count.table))
colnames(count.table) <- gsub(pattern='\\\\.\$', replacement='', colnames(count.table))

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

if (file.exists(opt\$outdir) == FALSE) {
    dir.create(opt\$outdir, recursive=TRUE)
}
setwd(opt\$outdir)

samples.vec     <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
coldata         <- data.frame(sample=samples.vec, row.names=samples.vec)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct!=1 & n_distinct!=length(samples.vec)]
    if (ncol(groupings)!=0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}

DDSFile <- paste(opt\$outprefix,".dds.RData",sep="")

counts  <- count.table[,samples.vec,drop=FALSE]
# `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~1)
dds     <- estimateSizeFactors(dds)
if (min(dim(count.table))<=1)  { # No point if only one sample, or one gene
    save(dds,file=DDSFile)
    saveRDS(dds, file=sub("\\\\.dds\\\\.RData\$", ".rds", DDSFile))
    warning("Not enough samples or genes in counts file for PCA.", call.=FALSE)
    quit(save = "no", status = 0, runLast = FALSE)
}
if (!opt\$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
} else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
save(dds,file=DDSFile)
saveRDS(dds, file=sub("\\\\.dds\\\\.RData\$", ".rds", DDSFile))

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca\$sdev^2/sum(pca\$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca\$x)
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df\$PC1)-median(df\$PC1)), abs(rank(df\$PC2)-median(df\$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

PlotFile <- paste(opt\$outprefix,".plots.pdf",sep="")

pdf(file=PlotFile, onefile=TRUE, width=7, height=7)
## PCA
ntop <- c(500, Inf)
for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name, ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")\$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, label=paste0(" ", sample, " "))) +
        geom_point() +
        geom_text(check_overlap=TRUE, vjust=0.5, hjust="inward") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) +
        theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)

    if (decompose) {
        pc_names <- paste0("PC", attr(pca.data, "percentVar")\$PC)
        long_pc <- reshape(pca.data, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
        long_pc <- subset(long_pc, component<=5)
        long_pc_grp <- reshape(long_pc, varying=names(groupings), direction="long", sep="", timevar="grouper")
        long_pc_grp <- subset(long_pc_grp, grouper<=5)
        long_pc_grp\$component <- paste("PC", long_pc_grp\$component)
        long_pc_grp\$grouper <- paste0(long_pc_grp\$grouper, c("st","nd","rd","th","th")[long_pc_grp\$grouper], " prefix")
        pl <- ggplot(long_pc_grp, aes(x=Group, y=PC)) +
            geom_point() +
            stat_summary(fun=mean, geom="line", aes(group = 1)) +
            labs(x=NULL, y=NULL, subtitle = plot_subtitle, title="PCs split by sample-name prefixes") +
            facet_grid(component~grouper, scales="free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        print(pl)
    }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

## WRITE PC1 vs PC2 VALUES TO FILE
pca.vals           <- pca.data[,c("PC1","PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
write.table(pca.vals, file = paste(opt\$outprefix, ".pca.vals.txt", sep=""),
            row.names = FALSE, col.names = TRUE, sep = "\\t", quote = TRUE)

## SAMPLE CORRELATION HEATMAP
sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
)

## WRITE SAMPLE DISTANCES TO FILE
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt\$outprefix, ".sample.dists.txt", sep=""),
            row.names=FALSE, col.names=TRUE, sep="\\t", quote=FALSE)
dev.off()

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir, recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir,opt\$outprefix, ".size_factors.RData", sep="")

normFactors <- sizeFactors(dds)
save(normFactors, file=NormFactorsFile)

for (name in names(sizeFactors(dds))) {
    sizeFactorFile <- paste(SizeFactorsDir,name, ".txt", sep="")
    write(as.numeric(sizeFactors(dds)[name]), file=sizeFactorFile)
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()

################################################
################################################
################################################
################################################

################################################
################################################
## MULTIQC CUSTOM CONTENT                     ##
################################################
################################################

# The PCA values and sample distances with the headers that MultiQC reads them with, labelled by aligner
label       <- ifelse('$task.ext.args2' == 'null', '', '$task.ext.args2')
label_lower <- tolower(label)
label_upper <- toupper(label)

multiqc_content <- function(header_file, header_from, header_label, values_pattern, output_suffix) {
    header <- readLines(header_file)
    header <- gsub(header_from, paste0(label_lower, '_', header_from), header, fixed = TRUE)
    header <- gsub(header_label, paste0(label_upper, ' ', header_label), header, fixed = TRUE)
    values <- unlist(lapply(Sys.glob(values_pattern), readLines))
    writeLines(c(header, values), paste0(label_lower, output_suffix))
}

multiqc_content('$pca_header_multiqc', 'deseq2_pca', 'DESeq2 PCA', '*.pca.vals.txt', '.pca.vals_mqc.tsv')
multiqc_content('$clustering_header_multiqc', 'deseq2_clustering', 'DESeq2 sample', '*.sample.dists.txt', '.sample.dists_mqc.tsv')

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version       <- strsplit(version[['version.string']], ' ')[[1]][3]
deseq2.version  <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
'versions.yml')
