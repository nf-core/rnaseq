---
order: 4
---

# Preprocessing

The pipeline performs several preprocessing steps before alignment: adapter trimming, optional rRNA removal, and optional contamination screening. This page covers how to configure each.

## Trim adapters

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper tool around Cutadapt and FastQC to perform quality and adapter trimming on FastQ files. Trim Galore! will automatically detect and trim the appropriate adapter sequence. It is the default trimming tool used by this pipeline, however you can use fastp instead by specifying the `--trimmer fastp` parameter. [fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast, all-in-one preprocessing for FastQ files. It has been developed in C++ with multithreading support to achieve higher performance. You can specify additional options for Trim Galore! and fastp via the `--extra_trimgalore_args` and `--extra_fastp_args` parameters, respectively.

:::note
TrimGalore! will only run using multiple cores if you are able to use more than 5 and 6 CPUs for single- and paired-end data, respectively. The total cores available to TrimGalore! will also be capped at 4 (7 and 8 CPUs in total for single- and paired-end data, respectively) because there is no longer a run-time benefit. See [release notes](https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019) and [discussion whilst adding this logic to the nf-core/atacseq pipeline](https://github.com/nf-core/atacseq/pull/65).
:::

## rRNA removal

Ribosomal RNA (rRNA) removal can be enabled with the `--remove_ribo_rna` parameter. The pipeline supports three different tools for rRNA removal, selectable via the `--ribo_removal_tool` parameter.

:::tip
For tools that use a reference database (SortMeRNA and Bowtie2), although rRNA is the primary target, the reference database can include additional abundant contaminant sequences you wish to remove, such as tRNAs or other non-coding RNAs. Simply add the paths to your custom FASTA files in the manifest file.
:::

### SortMeRNA (default)

[SortMeRNA](https://github.com/biocore/sortmerna) uses k-mer matching against rRNA databases to identify and filter rRNA reads. This is the default option and requires an rRNA database manifest file.

```bash
nextflow run nf-core/rnaseq --remove_ribo_rna --ribo_removal_tool sortmerna ...
```

By default, [rRNA databases](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) defined in the SortMeRNA GitHub repo are used. You can see an example in the pipeline GitHub repository in `assets/rrna-db-defaults.txt` which is used by default via the `--ribo_database_manifest` parameter.

:::note
The default databases are based on SILVA 119, which requires [licensing for commercial use](https://www.arb-silva.de/silva-license-information). SILVA 138+ uses CC-BY 4.0 licensing that freely permits commercial use with attribution. If you have licensing concerns, consider using Bowtie2 with custom rRNA reference sequences via `--ribo_removal_tool bowtie2`.
:::

### Bowtie2

[Bowtie2](https://github.com/BenLangmead/bowtie2) performs alignment-based filtering against rRNA reference sequences. Reads that align to the rRNA references are filtered out, and unaligned reads are kept for downstream analysis. This option also requires an rRNA database manifest file specified via `--ribo_database_manifest`.

```bash
nextflow run nf-core/rnaseq --remove_ribo_rna --ribo_removal_tool bowtie2 ...
```

### RiboDetector

:::warning
RiboDetector has known issues with ONNX multiprocessing that can cause hangs in containerized environments (Docker, Singularity). This makes it unreliable for production use in Nextflow pipelines. We recommend using SortMeRNA or Bowtie2 for rRNA removal until these issues are resolved upstream. See [hzi-bifo/RiboDetector#61](https://github.com/hzi-bifo/RiboDetector/pull/61) for details.
:::

[RiboDetector](https://github.com/hzi-bifo/RiboDetector) uses machine learning to identify rRNA reads without requiring a reference database. This makes it particularly useful when working with organisms that lack well-characterised rRNA sequences, or when you want to avoid database licensing requirements.

```bash
nextflow run nf-core/rnaseq --remove_ribo_rna --ribo_removal_tool ribodetector ...
```

RiboDetector automatically determines read length from your data and uses its pre-trained neural network model to classify reads.

## Subsample FASTQs

To reduce the number of reads used in the analysis, for example to test pipeline operation with limited resource usage, use the FASTP trimmer option (see above). FASTP has an option to take the first `n` reads of input FASTQ file(s), reducing the reads passed to subsequent steps. For example, to pass only the first 10,000 reads for trimming, set input parameters like:

```bash
--trimmer fastp --extra_fastp_args '--reads_to_process 10000'
```

## Screen for contamination

The pipeline provides the option to scan unaligned reads for contamination from other species using either [Sylph](https://sylph-docs.github.io/) or [Kraken2](https://ccb.jhu.edu/software/kraken2/), with the possibility of applying corrections from [Bracken](https://ccb.jhu.edu/software/bracken/). Since running Bracken is not computationally expensive, we recommend always using it to refine the abundance estimates generated by Kraken2.

Sylph is a [faster and much more memory-efficient tool](https://doi.org/10.1038/s41587-024-02412-y) with about equal precision in species detection to Kraken2/Bracken. Sylph also has lower rates of false positives. However, Sylph does not assign specific reads to species; it only provides overall abundance estimates. Sylph abundance estimates also [cannot assign a certain percentage of reads as unclassified](https://github.com/bluenote-1577/sylph/issues/49).

Pre-constructed sylph databases can be found [here](https://sylph-docs.github.io/pre%E2%80%90built-databases/) and taxonomies [here](https://sylph-docs.github.io/sylph-tax/). The [documentation](https://sylph-docs.github.io/sylph-tax/) also has instructions on creating custom databases/taxonomies. As a newer tool, the effect of database choice on Sylph's performance has not been explored as thoroughly as for Kraken2 or Bracken. However, the following comments on choosing databases for Kraken2 are very likely still applicable to an extent for Sylph.

The accuracy of Kraken2 is [highly dependent on the database](https://doi.org/10.1099/mgen.0.000949) used. Specifically, it is [crucial](https://doi.org/10.1128/mbio.01607-23) to ensure that the host genome/transcriptome is included in the database. (Note that the pre-built sylph databases do _not_ appear to contain the human genome/transcriptome). If you are particularly concerned about certain contaminants, it is beneficial to use a smaller, more focused database containing primarily those contaminants instead of the full standard database. Various pre-built databases [are available for download](https://benlangmead.github.io/aws-indexes/k2), and instructions for building a custom database can be found in the [Kraken2 documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown). Additionally, genomes of contaminants detected in previous sequencing experiments are available on the [OpenContami website](https://openlooper.hgc.jp/opencontami/help/help_oct.php).

While Kraken2 is capable of detecting low-abundance contaminants in a sample, false positives can occur. Therefore, if only a very small number of reads from a contaminating species are detected, interpret these results with caution. Lastly, while Kraken2 can be used without Bracken, since running Bracken is not computationally expensive, we recommend always using it to refine the abundance estimates generated by Kraken2.
