---
order: 1
---

# Differential expression analysis

This tutorial walks through differential expression (DE) analysis using output from the nf-core/rnaseq pipeline. You will learn how to run the pipeline for DE analysis, load the results into R, perform DE analysis with DESeq2, visualise the results, and run functional enrichment analysis.

## What you will learn

- How to configure nf-core/rnaseq for a DE experiment
- How to load pipeline output into R using DESeq2
- How to perform quality control (PCA, hierarchical clustering)
- How to identify differentially expressed genes
- How to visualise results (MA plots, volcano plots, heatmaps)
- How to run over-representation analysis with clusterProfiler

## Prerequisites

- A completed nf-core/rnaseq run with at least two experimental conditions and three biological replicates per condition
- R (version 4.0 or later) with RStudio or another R IDE
- The following R packages installed:
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  - [tidyverse](https://www.tidyverse.org/)
  - [pheatmap](https://cran.r-project.org/package=pheatmap)
  - [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)
  - [ggrepel](https://cran.r-project.org/package=ggrepel)
  - [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
  - [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) (or the appropriate organism annotation package)

:::note
This tutorial uses simulated data with 5 differentially expressed genes on human chromosome 21. The same workflow applies to real experimental data. It was originally developed at the nf-core hackathon Barcelona 2024 by Lorenzo Sola, Francesco Lescai, and Mariangela Santorsola.
:::

## Next steps

1. [Run the pipeline](running-the-pipeline.md) — Configure and run nf-core/rnaseq for your DE experiment
2. [DE analysis in R](de-analysis-in-r.md) — Perform differential expression analysis with DESeq2
