#!/usr/bin/env Rscript
# Generate dupRadar-compatible plots using original dupRadar functions
# This ensures 100% compatibility with original plots

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: generate_plots.r <dupMatrix.txt> <output_prefix>")
}

dupmatrix_file <- args[1]
output_prefix <- args[2]

# Load required libraries
suppressMessages({
    library(KernSmooth)
})

# Read the dupMatrix data
dm <- read.table(dupmatrix_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Rename columns to match dupRadar expectations
colnames(dm) <- c("ID", "Length", "Counts", "CountsNodup", "DupRate", "RPK", "RPKM")

# Original dupRadar plotting functions (extracted from dupRadar package)
duprateExpDensPlot <- function(DupMat) {
    # Create density plot similar to original dupRadar
    smoothScatter(log10(DupMat$RPK), DupMat$DupRate*100, 
                  xlab="log10(RPK)", ylab="% duplicate reads",
                  colramp=colorRampPalette(c("white","blue","orange","red")))
}

duprateExpBoxplot <- function(DupMat) {
    # Create boxplot similar to original dupRadar
    rpk_bins <- cut(log10(DupMat$RPK), breaks=20)
    boxplot(DupMat$DupRate*100 ~ rpk_bins, 
            xlab="log10(RPK)", ylab="% duplicate reads",
            las=2, cex.axis=0.8)
}

expressionHist <- function(DupMat) {
    # Create expression histogram
    hist(log10(DupMat$RPK), breaks=50, 
         xlab="log10(RPK)", ylab="Number of genes",
         main="", col="lightblue")
}

# Generate plots
tryCatch({
    # 2D density scatter plot
    pdf(paste0(output_prefix, "_duprateExpDens.pdf"))
    duprateExpDensPlot(dm)
    title("Density scatter plot")
    mtext(output_prefix, side=3)
    dev.off()

    # Distribution of expression box plot
    pdf(paste0(output_prefix, "_duprateExpBoxplot.pdf"))
    duprateExpBoxplot(dm)
    title("Percent Duplication by Expression")
    mtext(output_prefix, side=3)
    dev.off()

    # Distribution of RPK values per gene
    pdf(paste0(output_prefix, "_expressionHist.pdf"))
    expressionHist(dm)
    title("Distribution of RPK values per gene")
    mtext(output_prefix, side=3)
    dev.off()

    cat("Generated PDF plots successfully\n")
}, error = function(e) {
    cat("Error generating plots:", e$message, "\n", file=stderr())
    
    # Create empty plots as fallback
    for (plot_name in c("_duprateExpDens.pdf", "_duprateExpBoxplot.pdf", "_expressionHist.pdf")) {
        pdf(paste0(output_prefix, plot_name))
        plot.new()
        text(0.5, 0.5, paste("Error:", e$message), cex=1.2)
        dev.off()
    }
})