#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Usage: edgeR_heatmap_MDS.r <sample_1.bam> <sample_2.bam> <sample_3.bam> (more bam files optional)", call.=FALSE)
}

# Load / install required packages
if (!require("limma")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma", suppressUpdates=TRUE)
    library("limma")
}
if (!require("edgeR")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR", suppressUpdates=TRUE)
    library("edgeR")
}
if (!require("data.table")){
    install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("data.table")
}
if (!require("gplots")) {
    install.packages("gplots", dependencies=TRUE, repos='http://cloud.r-project.org/')
    library("gplots")
}

# Load count column from all files into a list of data frames
# Use data.tables fread as much much faster than read.table
# Row names are GeneIDs
temp <- lapply(lapply(args, fread, skip="Geneid", header=TRUE), function(x){return(as.data.frame(x)[,c(1, ncol(x))])})

# Merge into a single data frame
merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="Geneid")
}
data <- data.frame(Reduce(merge.all, temp))

# Clean sample name headers
colnames(data) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(data))

# Set GeneID as row name
rownames(data) <- data[,1]
data[,1] <- NULL

# Convert data frame to edgeR DGE object
dataDGE <- DGEList( counts=data.matrix(data) )

# Normalise counts
dataNorm <- calcNormFactors(dataDGE)

# Make MDS plot
pdf('edgeR_MDS_plot.pdf')
MDSdata <- plotMDS(dataNorm)
dev.off()

# Print distance matrix to file
write.csv(MDSdata$distance.matrix, 'edgeR_MDS_distance_matrix.csv', quote=FALSE,append=TRUE)

# Print plot x,y co-ordinates to file
MDSxy = MDSdata$cmdscale.out
colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))
write.csv(MDSxy, 'edgeR_MDS_Aplot_coordinates_mqc.csv', quote=FALSE, append=TRUE)

# Get the log counts per million values
logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

# Calculate the Pearsons correlation between samples
# Plot a heatmap of correlations
pdf('log2CPM_sample_correlation_heatmap.pdf')
hmap <- heatmap.2(as.matrix(cor(logcpm, method="pearson")),
  key.title="Pearson's Correlation", trace="none",
  dendrogram="row", margin=c(9, 9)
)
dev.off()

# Write correlation values to file
write.csv(hmap$carpet, 'log2CPM_sample_correlation_mqc.csv', quote=FALSE, append=TRUE)

# Plot the heatmap dendrogram
pdf('log2CPM_sample_distances_dendrogram.pdf')
hmap <- heatmap.2(as.matrix(dist(t(logcpm))))
plot(hmap$rowDendrogram, main="Sample Pearson's Correlation Clustering")
dev.off()

file.create("corr.done")

# Printing sessioninfo to standard out
print("Sample correlation info:")
sessionInfo()
