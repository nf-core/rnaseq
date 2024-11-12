---
order: 3
---

# The nf-core/rnaseq pipeline

In order to carry out a RNA-Seq analysis we will use the nf-core pipeline [rnaseq](https://nf-co.re/rnaseq/3.14.0).

## Overview

The pipeline is organised following the diffent blocks shown below: pre-processing, traditional alignment (or lightweight alignment) and quantification, post-processing and final QC.

<figure markdown="span">
  ![metromap](./img/nf-core-rnaseq_metro_map_grey.png){ width="1000"}
</figure>

In each process, the users can choose among a range of different options. Importantly, the users can decide to follow one of the two different routes in the alignment and quantification step:

- traditional alignment and quantification (stage 2);

- lightweight alignment and quantification (stage 3).

## Experimental Design

The number of reads and the number of biological replicates are two critical factors that researchers need to carefully consider during the design of a RNA-seq experiment. While it may seem intuitive that having a large number of reads is always desirable, an excessive number can lead to unnecessary costs and computational burdens, without providing significant improvements. Instead, it is often more beneficial to prioritise the number of biological replicates, as it allows to capture the natural biological variation of the data. Biological replicates involve collecting and sequencing RNA from distinct biological samples (e.g., different individuals, tissues, or time points), helping to detect genuine changes in gene expression.

:::warning
This concept must not be confused with technical replicates that asses the technical variability of the sequencing platform by sequencing the same RNA sample multiple time.
:::

To obtain optimal results, it is crucial to balance the number of biological replicates and the sequencing depth. While increasing the depth of sequencing enhances the ability to detect genes with low expression levels, there is a plateau beyond which no further benefits are gained. Statistical power calculations can inform experimental design by estimating the optimal number of reads and replicates required. For instance, this approach helps to establish a suitable log2 fold change threshold for the DE analysis. By incorporating multiple biological replicates into the design and optimizing sequencing depth, researchers can enhance the statistical power of the analysis, reducing the number of false positive results, and increasing the reliability of the findings.

## Library design

RNA-seq library design involves important decisions, particularly the choice between paired-end and single-end sequencing. Paired-end sequencing offers insights into structural variations and transcript isoforms, significantly improving mapping accuracy for longer transcripts and repetitive regions. In contrast, single-end sequencing—where only one end of the fragment is sequenced—can be a more cost-effective option while still delivering high-quality data for gene expression analysis. The choice between these two methods ultimately depends on the research question and experimental objectives. Paired-end sequencing is ideal for identifying novel transcripts or characterizing isoforms, whereas single-end sequencing is often sufficient for quantifying gene expression. Factors such as the type of RNA (e.g., mRNA or total RNA), read length, budget, and available computational resources also influence this decision.

## Reference genome

nf-core pipelines make use of the Illumina iGenomes collection as [reference genomes](https://nf-co.re/docs/usage/reference_genomes).

Before starting the analysis, the users might want to check whether the genome they need is part of this collection. They also might want to consider downloading the reference locally, when running on premises: this would be useful for multiple runs and to speed up the analysis. In this case the parameter `--igenomes_base` might be used to pass the root directory of the downloaded references.

One might also need to use custom files: in this case the users might either provide specific parameters at command line (`--fasta` option followed by the genome of choiche), or create a config file adding a new section to the `genome` object. See [here](https://nf-co.re/docs/usage/reference_genomes#custom-genomes) for more details.

In this tutorial we will edit the config file, since the data we will be using have been simulated on chromosome 21 of the Human GRCh38 reference, and we have prepared genome fasta and genome index containing only this chromosome locally.

The two files are `/workspace/gitpod/training/data/refs/Homo_sapiens_assembly38_chr21.fa` and `/workspace/gitpod/training/data/refs/Homo_sapiens_assembly38_chr21.fa.fai`, respectively.

## Reference annoation

The reference annotation plays a crucial role in the RNA-seq analysis. Without a high-quality reference annotation, RNA-seq analysis would result in inaccurate or incomplete results. The reference annotation provides a precise guide for aligning sequencing reads to specific genomic regions, allowing to identify genes, transcripts, and regulatory elements, as well as novel transcripts and alternative splicing events.

nf-core pipelines make use of the Illumina iGenomes collection also as [reference annotation](https://nf-co.re/docs/usage/reference_genomes).
The reference annotations are vastly out of date with respect to current annotations and miss certain features. So, the general recommendation is to download a newest annotation version compatible with the genome. A user can utilize the `--gtf` or the `--gff` options to specify the annottation files of choiche, or create a config file adding a new section to the `genome` object.

Similarly to the approach utilised for the genome, in this tutorial we will edit the config file. The annotation files include only the annotated transcripts on chromosome 21 of the Human GRCh38 reference genome and we have already prepared these files locally.

The two files are `/workspace/gitpod/training/data/refs/gencode_v29_chr21.gff` and `/workspace/gitpod/training/data/refs/gencode_v29_transcripts_chr21.fa`, respectively.

## Input files

The input data should be provided in a CSV file, according to a format that is largely common for nf-core pipelines.
The format is described in the [rnaseq usage page](https://nf-co.re/rnaseq/3.14.0/docs/usage).

The input file is `/workspace/gitpod/training/data/reads/rnaseq_samplesheet.csv`.

## Running nf-core/rnaseq

In the following sections we will:

- prepare our references;

- set our computational resources in order to be able to run the pipeline on a gitpod VM;

- edit the optional settings;

- run the pipeline.

## Reference and annotation files

Following the considerations above, we will first of all edit the `nextflow.config` file in our working directory to add a new genome.
It is sufficient to add the following code to the `parameters` directive in the config.

```groovy title="nextflow.config"
igenomes_base = '/workspace/gitpod/training/data/refs/'
genomes {
        'GRCh38chr21' {
            fasta                 = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta"
            fasta_fai             = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta.fai"
            gff                   = "${params.igenomes_base}/gencode_v29_chr21_parsed.gff"
            transcript_fasta      = "${params.igenomes_base}/gencode.v29.transcripts_chr21.fa"
            star_index            = "${params.igenomes_base}/star_index_chr21.tar.gz"
            salmon_index          = "${params.igenomes_base}/salmon_index_chr21.tar.gz"
	}
}
```

To speed up the analysis we will include the `star_index` and the `salmon_index` in the config. These files have already been created locally.

## Computing resources

Based on the choices we made when starting up the gitpod environment, we recommend to use the following additional parameters.
They can also be added to the parameters directive in the config file we just edited.

```groovy title="nextflow.config"
params {
    max_cpus      = 2
    max_memory    = '6.5GB'
    max_time      = '2.h'
}
```

## Launching the pipeline

Now we are ready to launch the pipeline, and we can use the following command line:

```bash
nextflow run nf-core/rnaseq -r 3.12.0 \
--input /workspace/gitpod/training/data/reads/rnaseq_samplesheet.csv \
--outdir ./results_star_salmon \
--genome GRCh38chr21 \
--aligner star_salmon \
--pseudo_aligner salmon \
--skip_biotype_qc \
--skip_stringtie \
--skip_bigwig \
--skip_umi_extract \
--skip_trimming \
--skip_fastqc \
--skip_markduplicates \
--skip_dupradar \
--skip_rseqc \
--skip_qualimap
```

Notice that we will run the pipeline with STAR as aligner and Salmon in alignment-based mode to quantify gene expression. We will also run the pipeline with Salmon in quasi-mapping mode to perform a lightweight alignment and quantification.

The `skip` parameters were inserted to reduce the running time.
