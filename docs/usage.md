# nf-core/rnaseq: Usage

## Introduction

You will need to create a samplesheet file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple replicates

The `group` identifier is the same when you have multiple replicates from the same experimental group, just increment the `replicate` identifier appropriately. The first replicate value for any given experimental group must be 1. Below is an example for a single experimental group in triplicate:

```bash
group,replicate,fastq_1,fastq_2,strandedness
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,reverse
control,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,reverse
control,3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,reverse
```

### Multiple runs of the same library

The `group` and `replicate` identifiers are the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will concatenate the raw reads before alignment. Below is an example for two samples sequenced across multiple lanes:

```bash
group,replicate,fastq_1,fastq_2,strandedness
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,unstranded
control,1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,unstranded
treatment,1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz,unstranded
treatment,1,AEG588A4_S4_L004_R1_001.fastq.gz,AEG588A4_S4_L004_R2_001.fastq.gz,unstranded
```

### Full design

A final design file consisting of both single- and paired-end data may look something like the one below. This is for two experimental groups in triplicate, where the last replicate of the `treatment` group has been sequenced twice.

```bash
group,replicate,fastq_1,fastq_2,strandedness
control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
control,2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
control,3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward
treatment,1,AEG588A4_S4_L003_R1_001.fastq.gz,,forward
treatment,2,AEG588A5_S5_L003_R1_001.fastq.gz,,forward
treatment,3,AEG588A6_S6_L003_R1_001.fastq.gz,,forward
treatment,3,AEG588A6_S6_L004_R1_001.fastq.gz,,forward
```

| Column         | Description                                                                                                 |
|----------------|-------------------------------------------------------------------------------------------------------------|
| `group`        | Group identifier for sample. This will be identical for replicate samples from the same experimental group. |
| `replicate`    | Integer representing replicate number. Must start from `1..<number of replicates>`.                         |
| `fastq_1`      | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `fastq_2`      | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `strandedness` | Sample strand-specificity. Must be one of `unstranded`, `forward` or `reverse`.                             |

An example samplesheet has been provided with the [pipeline](https://github.com/nf-core/rnaseq/blob/master/assets/samplesheet.csv).

## featureCounts Extra Gene Names

### Default "`gene_name`" Attribute Type

By default, the pipeline uses `gene_name` as the default gene identifier group. In case you need to adjust this, specify using the option `--fc_group_features` to use a different category present in your provided GTF file. Please also take care to use a suitable attribute to categorize the `biotype` of the selected features in your GTF then, using the option `--fc_group_features_type` (default: `gene_biotype`).

### Extra Gene Names or IDs

By default, the pipeline uses `gene_names` as additional gene identifiers apart from ENSEMBL identifiers in the pipeline.
This behaviour can be modified by specifying `--fc_extra_attributes` when running the pipeline, which is passed on to featureCounts as an `--extraAttributes` parameter.
See the user guide of the [Subread package here](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf).
Note that you can also specify more than one desired value, separated by a comma:
`--fc_extra_attributes gene_id,...`

#### Default "`exon`" Type

By default, the pipeline uses `exon` as the default to assign reads. In case you need to adjust this, specify using the option `--fc_count_type` to use a different category present in your provided GTF file (3rd column). For example, for nuclear RNA-seq, one could count reads in introns in addition to exons using `--fc_count_type transcript`.

### Transcriptome mapping with Salmon

Use the `--pseudo_aligner salmon` option to perform additional quantification at the transcript- and gene-level using [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). This will be run in addition to either STAR or HiSat2 and cannot be run in isolation, mainly because it allows you to obtain QC metrics with respect to the genomic alignments. By default, the pipeline will use the genome fasta and gtf file to generate the transcript fasta file, and then to build the Salmon index. You can override these parameters using the `--transcript_fasta` and `--salmon_index`, respectively.

The default Salmon parameters and a k-mer size of 31 are used to create the index. As [discussed here](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode)), a k-mer size off 31 works well with reads that are 75bp or longer.

### Alignment tool

By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the alignment tool instead. Developed by the same group behind the popular Tophat aligner, HISAT2 has a much smaller memory footprint.

To use HISAT2, use the parameter `--aligner hisat2` or set `params.aligner = 'hisat2'` in your config file.

### Reference genome paths

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--star_index '/path/to/STAR/index' \
--hisat2_index '/path/to/HISAT2/index' \
--fasta '/path/to/reference.fasta' \
--gtf '/path/to/gene_annotation.gtf' \
--gff '/path/to/gene_annotation.gff' \
--bed12 '/path/to/gene_annotation.bed'
```

Note that only one of `--star_index` / `--hisat2_index` are needed depending on which aligner you are using (see below).

The minimum requirements are a Fasta and GTF file. Note that `--gff` and `--bed` are auto-derived from the `--gtf` where needed and are not required. If these are provided and no others, then all other reference files will be automatically generated by the pipeline. If you specify a `--gff` file, it will be converted to GTF format automatically by the pipeline. If you specify both, the GTF is preferred over the GFF by the pipeline.

### "Type" of gene

The `gene_biotype` field which is typically found in Ensembl GTF files contains a key word description regarding the type of gene e.g. `protein_coding`, `lincRNA`, `rRNA`. In GENCODE GTF files this field has been renamed to `gene_type`.

ENSEMBL version:

```bash
8       havana  transcript      70635318        70669174        .       -       .       gene_id "ENSG00000147592"; gene_version "9"; transcript_id "ENST00000522447"; transcript_version "5"; gene_name "LACTB2"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "LACTB2-203"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS6208"; tag "basic"; transcript_support_level "2";
```

GENCODE version:

```bash
chr8    HAVANA  transcript      70635318        70669174        .       -       .       gene_id "ENSG00000147592.9"; transcript_id "ENST00000522447.5"; gene_type "protein_coding"; gene_name "LACTB2"; transcript_type "protein_coding"; transcript_name "LACTB2-203"; level 2; protein_id "ENSP00000428801.1"; transcript_support_level "2"; tag "alternative_3_UTR"; tag "basic"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS6208.1"; havana_gene "OTTHUMG00000164430.2"; havana_transcript "OTTHUMT00000378747.1";
```

Therefore, for `featureCounts` to correctly count the different biotypes when using a GENCODE annotation the `fc_group_features_type` is automatically set to `gene_type` when the `--gencode` flag is specified.

### Transcript IDs in FASTA files

The transcript IDs in GENCODE fasta files are separated by vertical pipes (`|`) rather than spaces.

ENSEMBL version:

```bash
>ENST00000522447.5 cds chromosome:GRCh38:8:70635318:70669174:-1 gene:ENSG00000147592.9 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:LACTB2 description:lactamase beta 2 [Source:HGNC Symbol;Acc:HGNC:18512]
```

GENCODE version:

```bash
>ENST00000522447.5|ENSG00000147592.9|OTTHUMG00000164430.2|OTTHUMT00000378747.1|LACTB2-203|LACTB2|1034|protein_coding|
```

This [issue](https://github.com/COMBINE-lab/salmon/issues/15) can be overcome by specifying the `--gencode` flag when building the Salmon index.

### Compressed Reference File Input

By default, the pipeline assumes that the reference genome files are all uncompressed, i.e. raw fasta or gtf files. If instead you intend to use compressed or gzipped references, like directly from ENSEMBL:

```bash
nextflow run nf-core/rnaseq --input samplesheet.csv \
    --genome ftp://ftp.ensembl.org/pub/release-97/fasta/microcebus_murinus/dna_index/Microcebus_murinus.Mmur_3.0.dna.toplevel.fa.gz \
    --gtf ftp://ftp.ensembl.org/pub/release-97/gtf/microcebus_murinus/Microcebus_murinus.Mmur_3.0.97.gtf.gz
```

This assumes that ALL of the reference files are compressed, including the reference indices, e.g. for STAR, HiSat2 or Salmon. For instructions on how to create your own compressed reference files, see the instructions below. This also includes any files specified with `--additional_fasta`, which are assumed to be compressed as well when the `--fasta` file is compressed. The pipeline auto-detects `gz` input for reference files. Mixing of `gz` and non-compressed input is not possible!

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnaseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnaseq releases page](https://github.com/nf-core/rnaseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/rnaseq`](https://hub.docker.com/r/nfcore/rnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/rnaseq`](https://hub.docker.com/r/nfcore/rnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
