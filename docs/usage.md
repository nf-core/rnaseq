# nf-core/rnaseq: Usage

## :warning: Read this documentation on the nf-core website: [https://nf-co.re/rnaseq/usage](https://nf-co.re/rnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Quick start

Create a samplesheet describing your samples:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,auto
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,auto
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,auto
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,auto
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,auto
```

Run the pipeline:

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    -profile docker
```

This runs the default workflow: adapter trimming with Trim Galore, STAR alignment, Salmon quantification, and comprehensive QC. Results are written to `results/`, including a MultiQC report summarising all quality metrics. See the [output documentation](output.md) for a full description of results.

:::tip
Set strandedness to `auto` and the pipeline will infer it for you. If you already know the strandedness of your libraries, specify `forward`, `reverse`, or `unstranded` instead.
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
