# nf-core/rnaseq: Usage

## :warning: Read this documentation on the nf-core website: [https://nf-co.re/rnaseq/usage](https://nf-co.re/rnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Quick start

### 1. Create a samplesheet

Describe your FASTQ files in a CSV. Each row is one library. Use `auto` for strandedness and the pipeline will detect it:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,data/control_rep1_R1.fastq.gz,data/control_rep1_R2.fastq.gz,auto
CONTROL_REP2,data/control_rep2_R1.fastq.gz,data/control_rep2_R2.fastq.gz,auto
CONTROL_REP3,data/control_rep3_R1.fastq.gz,data/control_rep3_R2.fastq.gz,auto
TREATMENT_REP1,data/treatment_rep1_R1.fastq.gz,data/treatment_rep1_R2.fastq.gz,auto
TREATMENT_REP2,data/treatment_rep2_R1.fastq.gz,data/treatment_rep2_R2.fastq.gz,auto
TREATMENT_REP3,data/treatment_rep3_R1.fastq.gz,data/treatment_rep3_R2.fastq.gz,auto
```

See the [samplesheet documentation](usage/samplesheet.md) for the full format, single-end data, and advanced options.

### 2. Run the pipeline

You need a reference genome FASTA and gene annotation GTF. Download these from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) or another provider (see [reference genomes](usage/reference-genomes.md) for guidance):

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta /path/to/genome.fa \
    --gtf /path/to/genes.gtf \
    -profile docker
```

This trims adapters, aligns reads with STAR, quantifies with Salmon, and runs comprehensive QC. Results are written to `results/`.

### 3. Check results

Open `results/multiqc/star_salmon/multiqc_report.html` for an overview of all quality metrics across your samples. Key things to look for:

- **Alignment rate**: most reads as "Uniquely mapped" in the STAR section
- **Strandedness**: check the "Strandedness Checks" table for mismatches
- **PCA plot**: samples cluster by condition, not by batch

Count matrices for downstream analysis are in `results/star_salmon/` — see the [output documentation](output.md) for details.

:::tip
If you already know the strandedness of your libraries, specify `forward`, `reverse`, or `unstranded` instead of `auto`.
:::

## Pipeline parameters

Provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters; see the [nf-core configuration documentation](https://nf-co.re/usage/configuration#custom-configuration-files).

## What would you like to do?

### Set up your input data

- [**Samplesheet format**](usage/samplesheet.md) — How to create your samplesheet, field definitions, handling multiple runs of the same sample, strandedness options, and BAM input for reprocessing.

### Configure your analysis

- [**Reference genomes**](usage/reference-genomes.md) — How to provide genome FASTA and GTF files, choosing the right annotation, GENCODE support, building and reusing indices.
- [**Alignment and quantification**](usage/alignment-and-quantification.md) — Choosing between STAR+Salmon (default), RSEM, HISAT2, or pseudoalignment-only with Salmon/Kallisto.
- [**Preprocessing**](usage/preprocessing.md) — Adapter trimming options, rRNA removal, FASTQ subsampling for test runs, and contamination screening.

### Handle specialised experiments

- [**Advanced features**](usage/advanced-features.md) — UMI handling, 3' digital gene expression assays, prokaryotic RNA-seq, GPU acceleration with Sentieon or Parabricks, and BAM reprocessing workflows.

### Analyse your results

- [**Differential expression analysis**](usage/differential_expression_analysis/) — A guided tutorial for performing differential expression analysis with DESeq2 using nf-core/rnaseq output.

### Configure Nextflow and compute resources

- [**Configuration**](usage/configuration.md) — Running the pipeline, using params files, Nextflow profiles, `-resume`, custom resource configuration, and institutional configs.
