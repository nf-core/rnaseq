# nf-core/rnaseq: Documentation

Welcome to the nf-core/rnaseq documentation! This pipeline performs comprehensive RNA-seq analysis from raw reads through quality control, alignment, quantification, and differential expression analysis.

## 📚 Documentation Structure

The nf-core/rnaseq documentation is organized into the following main sections:

### Core Documentation

- **[Usage Guide](usage.md)** - Complete reference for running the pipeline, including all parameters, configuration options, and advanced features
- **[Output Documentation](output.md)** - Detailed descriptions of all output files, quality control metrics, and result interpretation
- **[Differential Expression Tutorial](usage/differential_expression_analysis/introduction.md)** - Hands-on workshop covering RNA-seq theory and practical differential expression analysis

---

## 🚀 Quick Start

**New to the pipeline?** Here's what you need to get started:

### 1. Prepare Your Samplesheet

Create a CSV file with your sample information ([detailed guide](usage.md#samplesheet-input)):

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_1,control_1_R1.fastq.gz,control_1_R2.fastq.gz,auto
CONTROL_2,control_2_R1.fastq.gz,control_2_R2.fastq.gz,auto
TREATMENT_1,treatment_1_R1.fastq.gz,treatment_1_R2.fastq.gz,auto
TREATMENT_2,treatment_2_R1.fastq.gz,treatment_2_R2.fastq.gz,auto
```

### 2. Run the Pipeline

```bash
nextflow run nf-core/rnaseq \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38 \
  -profile docker
```

### 3. View Your Results

Open `results/multiqc/multiqc_report.html` in your browser for a comprehensive quality control overview.

📖 **For detailed instructions**, see the complete [Usage Guide](usage.md)

---

## 🎯 Common Use Cases

### I want to...

#### Analyze Standard RNA-seq Data

- **Start here**: [Quick Start](#quick-start) above
- **Learn more**: [Usage Guide](usage.md) and [Samplesheet Preparation](usage.md#samplesheet-input)

#### Understand Differential Expression Analysis

- **Start here**: [Differential Expression Tutorial](usage/differential_expression_analysis/introduction.md)
- **Theory background**: [RNA-seq Analysis Theory](usage/differential_expression_analysis/theory.md)

#### Reprocess Existing BAM Files

- **Use case**: Skip expensive alignment steps when reanalyzing data
- **Guide**: [BAM Input for Reprocessing Workflow](usage.md#bam-input-for-reprocessing-workflow)

#### Remove Ribosomal RNA Contamination

- **Options**: SortMeRNA (default), Bowtie2, or RiboDetector
- **Guide**: [rRNA Removal Options](usage.md#rrna-removal-options)

#### Analyze Prokaryotic RNA-seq

- **Specific workflow**: Uses Bowtie2 + Salmon
- **Guide**: [Prokaryotic RNA-seq Analysis](output.md#bowtie2-and-salmon-prokaryotic)

#### Work with UMI-based Protocols

- **Features**: UMI extraction and deduplication
- **Guide**: [Unique Molecular Identifiers (UMI)](usage.md#unique-molecular-identifiers-umi)

#### Analyze 3' Digital Gene Expression Data

- **Specialized workflow**: For 3' end counting assays
- **Guide**: [3' Digital Gene Expression](usage.md#3-digital-gene-expression-assays)

#### Use My Own Reference Genome

- **Best practice**: Explicit file specification
- **Guide**: [Reference Genome Options](usage.md#reference-genome-options)

#### Screen for Contamination

- **Tools**: Kraken2/Bracken or Sylph
- **Guide**: [Contamination Screening Options](usage.md#contamination-screening-options)

#### Accelerate Analysis with GPU or Commercial Tools

- **Sentieon**: [STAR Acceleration](usage.md#sentieon-acceleration-for-star)
- **Parabricks**: [GPU Acceleration](usage.md#parabricks-gpu-acceleration-for-star)

---

## 📋 Documentation by Topic

### Input & Data Preparation

- [Samplesheet Format](usage.md#samplesheet-input) - Required columns and structure
- [Multiple Runs](usage.md#multiple-runs-of-the-same-sample) - Merging sequencing runs
- [Strandedness Detection](usage.md#strandedness-prediction) - Auto-detection and manual specification
- [BAM Input](usage.md#bam-input-for-reprocessing-workflow) - Reprocessing workflow
- [FASTQ Linting](usage.md#linting) - Quality validation

### Processing Options

- [Adapter Trimming](usage.md#adapter-trimming-options) - TrimGalore vs fastp
- [rRNA Removal](usage.md#rrna-removal-options) - SortMeRNA, Bowtie2, RiboDetector
- [FASTQ Sampling](usage.md#fastq-sampling) - Subsample reads for testing

### Alignment & Quantification

- [Alignment Tools](usage.md#alignment-options) - STAR, HISAT2, STAR+RSEM
- [Quantification Methods](usage.md#quantification-options) - Salmon, Kallisto, featureCounts
- [Pseudoalignment](output.md#pseudoalignment-and-quantification) - Lightweight quantification
- [Prokaryotic Analysis](output.md#bowtie2-and-salmon-prokaryotic) - Bowtie2 + Salmon workflow

### Quality Control & Analysis

- [MultiQC Reports](output.md#multiqc) - Comprehensive QC overview
- [RSeQC Metrics](output.md#rseqc) - Read distribution and strand specificity
- [DESeq2 Analysis](output.md#deseq2) - Differential expression and PCA
- [Contamination Screening](usage.md#contamination-screening-options) - Kraken2 or Sylph

### Configuration & Execution

- [Resource Requests](usage.md#resource-requests) - Memory and CPU allocation
- [Custom Containers](usage.md#custom-containers) - Using custom Docker/Singularity images
- [Tool Arguments](usage.md#custom-tool-arguments) - Per-tool parameter customization
- [Reproducibility](usage.md#reproducibility) - Ensuring consistent results

### Understanding Results

- [Output Files Overview](output.md#pipeline-overview) - What gets generated
- [Preprocessing Results](output.md#preprocessing) - FastQC, trimming, filtering
- [Alignment Results](output.md#alignment-and-quantification) - BAM files and counts
- [QC Metrics](output.md#quality-control) - Interpreting quality control results

---

## 🆘 Getting Help

### Documentation Resources

- **Pipeline Website**: [https://nf-co.re/rnaseq](https://nf-co.re/rnaseq)
- **Parameter Reference**: [https://nf-co.re/rnaseq/parameters](https://nf-co.re/rnaseq/parameters)
- **nf-core Documentation**: [https://nf-co.re](https://nf-co.re)
- **Nextflow Documentation**: [https://www.nextflow.io/docs/latest/](https://www.nextflow.io/docs/latest/)

### Community Support

- **Slack Channel**: [#rnaseq on nf-core Slack](https://nfcore.slack.com/channels/rnaseq)
- **GitHub Issues**: [Report bugs or request features](https://github.com/nf-core/rnaseq/issues)
- **nf-core Community**: [Join the community](https://nf-co.re/join)

### Training Resources

- **nf-core Training**: [https://nf-co.re/events/training](https://nf-co.re/events/training)
- **Nextflow Training**: [https://training.nextflow.io](https://training.nextflow.io)
- **Differential Expression Tutorial**: [Tutorial Workshop](usage/differential_expression_analysis/introduction.md)

---

## 💡 Tips for New Users

> **📖 New to Nextflow?** Start with the [Differential Expression Tutorial](usage/differential_expression_analysis/introduction.md) for a guided introduction to both RNA-seq concepts and the pipeline.

> **⚡ Quick test run?** Use the built-in test profile: `nextflow run nf-core/rnaseq -profile test,docker`

> **🔍 Parameter documentation** is auto-generated from the pipeline schema. See the [nf-core website](https://nf-co.re/rnaseq/parameters) for the complete parameter reference.

> **🎯 Not sure which aligner to use?** STAR + Salmon is the default and works well for most eukaryotic RNA-seq. See [Alignment Options](usage.md#alignment-options) for alternatives.

> **⚠️ Strandedness unknown?** Use `strandedness,auto` in your samplesheet - the pipeline will detect it automatically. See [Strandedness Prediction](usage.md#strandedness-prediction).

---

## 📊 Pipeline Overview

The pipeline performs the following key steps:

1. **Preprocessing** - Quality control (FastQC), adapter trimming (TrimGalore/fastp), rRNA removal
2. **Alignment** - Map reads to reference genome (STAR, HISAT2) or transcriptome (Salmon, Kallisto)
3. **Quantification** - Generate gene/transcript counts
4. **Quality Control** - Comprehensive metrics (RSeQC, Qualimap, dupRadar, Preseq)
5. **Differential Expression** - Optional DESeq2 analysis with PCA and clustering
6. **Contamination Screening** - Optional Kraken2/Bracken or Sylph profiling
7. **Reporting** - Integrated MultiQC report

For a detailed workflow diagram and step-by-step output descriptions, see the [Output Documentation](output.md).

---

## 🔬 Citations & Credits

If you use nf-core/rnaseq for your analysis, please cite:

- **The pipeline**: [nf-core/rnaseq publication](https://doi.org/10.1038/s41587-020-0439-x)
- **Nextflow**: [Di Tommaso et al. 2017](https://doi.org/10.1038/nbt.3820)
- **Specific tools**: Citations are provided in the [MultiQC report](output.md#multiqc) and [pipeline output](output.md#pipeline-information)

For a complete list of tools and citations, see the [Pipeline Information](output.md#pipeline-information) section.
