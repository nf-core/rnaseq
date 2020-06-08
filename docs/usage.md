# nf-core/rnaseq: Usage

## Table of contents

- [nf-core/rnaseq: Usage](#nf-corernaseq-usage)
  - [Table of contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Running the pipeline](#running-the-pipeline)
    - [Updating the pipeline](#updating-the-pipeline)
    - [Reproducibility](#reproducibility)
  - [Main arguments](#main-arguments)
    - [`-profile`](#profile)
    - [`--reads`](#reads)
    - [`--single_end`](#singleend)
    - [Library strandedness](#library-strandedness)
  - [FeatureCounts Extra Gene Names](#featurecounts-extra-gene-names)
    - [Default "`gene_name`" Attribute Type](#default-%22genename%22-attribute-type)
    - [Extra Gene Names or IDs](#extra-gene-names-or-ids)
    - [Default "`exon`" Type](#default-%22exon%22-type)
  - [Transcriptome mapping with Salmon](#transcriptome-mapping-with-salmon)
  - [Alignment tool](#alignment-tool)
  - [Reference genomes](#reference-genomes)
    - [`--genome` (using iGenomes)](#genome-using-igenomes)
    - [`--star_index`, `--hisat2_index`, `--fasta`, `--gtf`, `--bed12`](#starindex---hisat2index---fasta---gtf---bed12)
    - [`--saveReference`](#savereference)
    - [`--saveTrimmed`](#savetrimmed)
    - [`--saveUnaligned`](#saveunaligned)
    - [`--saveAlignedIntermediates`](#savealignedintermediates)
    - [`--gencode`](#gencode)
      - ["Type" of gene](#%22type%22-of-gene)
      - [Transcript IDs in FASTA files](#transcript-ids-in-fasta-files)
    - [`--skipBiotypeQC`](#skipbiotypeqc)
    - [`--additional_fasta`](#--additional_fasta)
    - [`--skipAlignment`](#skipalignment)
    - [Compressed Reference File Input](#compressed-reference-file-input)
      - [Create compressed (tar.gz) STAR indices](#create-compressed-targz-star-indices)
      - [HISAT2 indices](#hisat2-indices)
      - [Salmon index](#salmon-index)
  - [Adapter Trimming](#adapter-trimming)
    - [`--clip_r1 [int]`](#clipr1-int)
    - [`--clip_r2 [int]`](#clipr2-int)
    - [`--three_prime_clip_r1 [int]`](#threeprimeclipr1-int)
    - [`--three_prime_clip_r2 [int]`](#threeprimeclipr2-int)
    - [`--trim_nextseq [int]`](#trimnextseq-int)
    - [`--skipTrimming`](#skiptrimming)
  - [Ribosomal RNA removal](#ribosomal-rna-removal)
    - [`--removeRiboRNA`](#removeriborna)
    - [`--saveNonRiboRNAReads`](#savenonribornareads)
    - [`--rRNA_database_manifest`](#rrnadatabasemanifest)
  - [Library Prep Presets](#library-prep-presets)
    - [`--pico`](#pico)
  - [Skipping QC steps](#skipping-qc-steps)
  - [Job resources](#job-resources)
    - [Automatic resubmission](#automatic-resubmission)
    - [Custom resource requests](#custom-resource-requests)
  - [AWS Batch specific parameters](#aws-batch-specific-parameters)
    - [`--awsqueue`](#awsqueue)
    - [`--awsregion`](#awsregion)
    - [`--awscli`](#awscli)
  - [Other command line parameters](#other-command-line-parameters)
    - [`--outdir`](#outdir)
    - [`--email`](#email)
    - [`--email_on_fail`](#emailonfail)
    - [`--max_multiqc_email_size`](#maxmultiqcemailsize)
    - [`-name`](#name)
    - [`-resume`](#resume)
    - [`-c`](#c)
    - [`--custom_config_version`](#customconfigversion)
    - [`--custom_config_base`](#customconfigbase)
    - [`--max_memory`](#maxmemory)
    - [`--max_time`](#maxtime)
    - [`--max_cpus`](#maxcpus)
    - [`--hisat_build_memory`](#hisatbuildmemory)
    - [`--sampleLevel`](#samplelevel)
    - [`--percent_aln_skip`](#percentalnskip)
    - [`--plaintext_email`](#plaintextemail)
    - [`--monochrome_logs`](#monochromelogs)
    - [`--multiqc_config`](#multiqcconfig)
  - [Stand-alone scripts](#stand-alone-scripts)
    <!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

This could be added by either editing those files with `nano`, e.g. `nano ~/.bashrc`, or by appending those options to your `~/.bashrc` file with:

```bash
echo "export NXF_OPTS='-Xms1g -Xmx4g'" >> ~/.bashrc
```

Then, re-running the file with `source ~/.bashrc`. If you are on a Mac, replace `~/.bashrc` in all previous examples with `~/.bash_profile`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker
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

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
- `conda`
  - Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### Library strandedness

Three command line flags / config parameters set the library strandedness for a run:

- `--forwardStranded`
- `--reverseStranded`
- `--unStranded`

If not set, the pipeline will be run as unstranded. Specifying `--pico` makes the pipeline run in `forwardStranded` mode.

You can set a default in a cutom Nextflow configuration file such as one saved in `~/.nextflow/config` (see the [nextflow docs](https://www.nextflow.io/docs/latest/config.html) for more). For example:

```nextflow
params {
    reverseStranded = true
}
```

If you have a default strandedness set in your personal config file you can use `--unStranded` to overwrite it for a given run.

These flags affect the commands used for several steps in the pipeline - namely HISAT2, featureCounts, RSeQC (`RPKM_saturation.py`), Qualimap and StringTie:

- `--forwardStranded`
  - HISAT2: `--rna-strandness F` / `--rna-strandness FR`
  - featureCounts: `-s 1`
  - RSeQC: `-d ++,--` / `-d 1++,1--,2+-,2-+`
  - Qualimap: `-pe strand-specific-forward`
  - StringTie: `--fr`
- `--reverseStranded`
  - HISAT2: `--rna-strandness R` / `--rna-strandness RF`
  - featureCounts: `-s 2`
  - RSeQC: `-d +-,-+` / `-d 1+-,1-+,2++,2--`
  - Qualimap: `-pe strand-specific-reverse`
  - StringTie: `--rf`

## FeatureCounts Extra Gene Names

### Default "`gene_name`" Attribute Type

By default, the pipeline uses `gene_name` as the default gene identifier group. In case you need to adjust this, specify using the option `--fc_group_features` to use a different category present in your provided GTF file. Please also take care to use a suitable attribute to categorize the `biotype` of the selected features in your GTF then, using the option `--fc_group_features_type` (default: `gene_biotype`).

### Extra Gene Names or IDs

By default, the pipeline uses `gene_names` as additional gene identifiers apart from ENSEMBL identifiers in the pipeline.
This behaviour can be modified by specifying `--fc_extra_attributes` when running the pipeline, which is passed on to featureCounts as an `--extraAttributes` parameter.
See the user guide of the [Subread package here](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf).
Note that you can also specify more than one desired value, separated by a comma:
`--fc_extra_attributes gene_id,...`

### Default "`exon`" Type

By default, the pipeline uses `exon` as the default to assign reads. In case you need to adjust this, specify using the option `--fc_count_type` to use a different category present in your provided GTF file (3rd column). For example, for nuclear RNA-seq, one could count reads in introns in addition to exons using `--fc_count_type transcript`.

## Transcriptome mapping with Salmon

Use the `--pseudo_aligner salmon` option to perform additional quantification at the transcript- and gene-level using [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). This will be run in addition to either STAR or HiSat2 and cannot be run in isolation, mainly because it allows you to obtain QC metrics with respect to the genomic alignments. By default, the pipeline will use the genome fasta and gtf file to generate the transcript fasta file, and then to build the Salmon index. You can override these parameters using the `--transcript_fasta` and `--salmon_index`, respectively.

The default Salmon parameters and a k-mer size of 31 are used to create the index. As [discussed here](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode)), a k-mer size off 31 works well with reads that are 75bp or longer.

## Alignment tool

By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the alignment tool instead. Developed by the same group behind the popular Tophat aligner, HISAT2 has a much smaller memory footprint.

To use HISAT2, use the parameter `--aligner hisat2` or set `params.aligner = 'hisat2'` in your config file.

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

- Human
  - `--genome GRCh37`
- Mouse
  - `--genome GRCm38`
- _Drosophila_
  - `--genome BDGP6`
- _S. cerevisiae_
  - `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      star    = '<path to the star index folder>'
      fasta   = '<path to the genome fasta file>' // Used if no star index given
      gtf     = '<path to the genome gtf file>'
      bed12   = '<path to the genome bed file>' // Generated from GTF if not given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--star_index`, `--hisat2_index`, `--fasta`, `--gtf`, `--bed12`

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

### `--saveReference`

Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

### `--saveTrimmed`

By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveUnaligned`

By default, the pipeline doesn't export unaligned/unmapped reads to a separate file. Using this option, STAR / HISAT2 and Salmon will produce a separate BAM file or a list of reads that were not aligned in a separate output directory.

### `--saveAlignedIntermediates`

As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM files from STAR / HISAT2 and sorting steps.

### `--gencode`

If your `--gtf` file is in GENCODE format and you would like to run Salmon (`--pseudo_aligner salmon`) you will need to provide this parameter in order to build the Salmon index appropriately. The `params.fc_group_features_type=gene_type` will also be set as explained below.

[GENCODE](https://www.gencodegenes.org) gene annotations are slightly different from ENSEMBL or iGenome annotations in two ways.

#### "Type" of gene

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

#### Transcript IDs in FASTA files

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

### `--skipBiotypeQC`

This skips the BiotypeQC step in the `featureCounts` process, explicitly useful when there is no available GTF/GFF with any `biotype` or similar information that could be used before.

### `--additional_fasta`

If provided, any genes here will get concatenated to the existing genome fasta, a GTF will be automatically created using the entire sequence as the `gene`, `transcript`, and `exon` features, and the alignment index will get created off of the combined fasta and GTF. It is recommended to save the reference with `--saveReference` so you do not need to create it again.

### `--skipAlignment`

By default, the pipeline aligns the input reads to the genome using either HISAT2 or STAR and counts gene expression using featureCounts. If you prefer to skip alignment altogether and only get transcript/gene expression counts with pseudo alignment, use this flag. Note that you will also need to specify `--pseudo_aligner salmon`. If you have a custom transcriptome, supply that with `--transcript_fasta`.

### Compressed Reference File Input

By default, the pipeline assumes that the reference genome files are all uncompressed, i.e. raw fasta or gtf files. If instead you intend to use compressed or gzipped references, like directly from ENSEMBL:

```bash
nextflow run nf-core/rnaseq --reads 'data/{R1,R2}*.fastq.gz' \
    --genome ftp://ftp.ensembl.org/pub/release-97/fasta/microcebus_murinus/dna_index/Microcebus_murinus.Mmur_3.0.dna.toplevel.fa.gz \
    --gtf ftp://ftp.ensembl.org/pub/release-97/gtf/microcebus_murinus/Microcebus_murinus.Mmur_3.0.97.gtf.gz
```

This assumes that ALL of the reference files are compressed, including the reference indices, e.g. for STAR, HiSat2 or Salmon. For instructions on how to create your own compressed reference files, see the instructions below. This also includes any files specified with `--additional_fasta`, which are assumed to be compressed as well when the `--fasta` file is compressed. The pipeline auto-detects `gz` input for reference files. Mixing of `gz` and non-compressed input is not possible!

#### Create compressed (tar.gz) STAR indices

STAR indices can be created by using `--saveReference`, and then using `tar` on them:

```bash
cd results/reference_genome
tar -zcvf star.tar.gz star
```

#### HISAT2 indices

HiSAT2 indices can be created by using `--saveReference`, and then using `tar` on them:

```bash
cd results/reference_genome
tar -zcvf hisat2.tar.gz *.hisat2_*
```

#### Salmon index

Salmon indices can be created by using `--saveReference`, and then using `tar` on them:

```bash
cd results/reference_genome
tar -zcvf salmon_index.tar.gz salmon_index
```

## Adapter Trimming

If specific additional trimming is required (for example, from additional tags),
you can use any of the following command line parameters. These affect the command
used to launch TrimGalore!

### `--clip_r1 [int]`

Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

### `--clip_r2 [int]`

Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

### `--three_prime_clip_r1 [int]`

Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

### `--three_prime_clip_r2 [int]`

Instructs Trim Galore to remove bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.

### `--trim_nextseq [int]`

This enables the option --nextseq-trim=3'CUTOFF within Cutadapt in Trim Galore, which will set a quality cutoff (that is normally given with -q instead), but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases.

### `--skipTrimming`

This allows to skip the trimming process to save time when re-analyzing data that has been trimmed already.

## Ribosomal RNA removal

If rRNA removal is desired (for example, metatranscriptomics),
add the following command line parameters.
Please be adviced that by default these steps make use of the SILVA v119 database that requires [`licencing for commercial/non-academic entities`](https://www.arb-silva.de/silva-license-information).

### `--removeRiboRNA`

Instructs to use SortMeRNA to remove reads related to ribosomal RNA (or any patterns found in the sequences defined by `--rRNA_database_manifest`).

### `--saveNonRiboRNAReads`

By default, non-rRNA FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--rRNA_database_manifest`

By default, rRNA databases in github [`biocore/sortmerna/rRNA_databases`](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) are used. Here the path to a text file can be provided that contains paths to fasta files (one per line, no ' or " for file names) that will be used for database creation for SortMeRNA instead of the default ones. You can see an example in the directory `assets/rrna-default-dbs.txt`. Consequently, similar reads to these sequences will be removed.
Be aware that commercial/non-academic entities require [`licensing for SILVA`](https://www.arb-silva.de/silva-license-information) with these default databases.

## Library Prep Presets

Some command line options are available to automatically set parameters for common RNA-seq library preparation kits.

> Note that these presets override other command line arguments. So if you specify `--pico --clip_r1 0`, the `--clip_r1` bit will be ignored.

If you have a kit that you'd like a preset added for, please let us know!

### `--pico`

Sets trimming and standedness settings for the _SMARTer Stranded Total RNA-Seq Kit - Pico Input_ kit.

Equivalent to: `--forwardStranded` `--clip_r1 3` `--three_prime_clip_r2 3`

## Skipping QC steps

The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

- `--skipQC` - Skip **all QC steps**, apart from MultiQC
- `--skipFastQC` - Skip FastQC
- `--skipRseQC` - Skip RSeQC
- `--skipQualimap` - Skip Qualimap
- `--skipPreseq` - Skip Preseq
- `--skipDupRadar` - Skip dupRadar (and Picard MarkDuplicates)
- `--skipEdgeR` - Skip edgeR MDS plot and heatmap
- `--skipMultiQC` - Skip MultiQC

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `--publish_dir_mode`

Choose the publishDir mode. Default is `copy`. See [nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for more details.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--hisat_build_memory`

Required amount of memory in GB to build HISAT2 index with splice sites.
The HiSAT2 index build can proceed with or without exon / splice junction information.
To work with this, a very large amount of memory is required.
If this memory is not available, the index build will proceed without splicing information.
The `--hisat_build_memory` option changes this threshold. By default it is `200GB` - if your system
`--max_memory` is set to `128GB` but your genome is small enough to build using this, then you can
allow the exon build to proceed by supplying `--hisat_build_memory 100GB`

### `--sampleLevel`

Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `--percent_aln_skip`

The pipeline will remove any samples from further processing that receive a percentage alignment below this.
This is because downstream steps typically fail otherwise, halting the execution of the pipeline for all samples.
Default: `5` (percent reads aligned).

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.

## Stand-alone scripts

The `bin` directory contains some scripts used by the pipeline which may also be run manually:

- `gtf2bed`
  - Script used to generate the BED12 reference files used by RSeQC. Takes a `.gtf` file as input
- `dupRadar.r`
  - dupRadar script used in the _dupRadar_ pipeline process.
- `edgeR_heatmap_MDS.r`
  - edgeR script used in the _Sample Correlation_ process
