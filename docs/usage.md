# nfcore/rnaseq Usage

## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile (Swedish UPPMAX users use `-profile uppmax`). See below for more information about profiles.

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

First, go to the [nfcore/rnaseq releases page](https://github.com/nf-core/rnaseq/releases) and find the latest version number - numeric only (eg. `1.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment - follow the links below to see instructions for running on that system. Available profiles are:

* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Runs using the `local` executor and pulls software from dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
* `uppmax`, `uppmax_modules`, `uppmax_devel`
    * Designed to be used on the Swedish [UPPMAX](http://uppmax.uu.se/) clusters such as `milou`, `rackham`, `bianca` and `irma`
    * See [`docs/configuration/uppmax.md`](configuration/uppmax.md)
* `hebbe`
    * Designed to be run on the [c3se Hebbe cluster](http://www.c3se.chalmers.se/index.php/Hebbe) in Chalmers, Gothenburg.
    * See [`docs/configuration/c3se.md`](configuration/c3se.md)
* `binac`, `cfc`
    * Profiles for clusters at QBiC in Tübingen, Germany
    * See [`docs/configuration/qbic.md`](configuration/qbic.md)
* `awsbatch`
    * Profile for running on AWSBatch, specific parameters are described below
* `aws`
    * A starter configuration for running the pipeline on Amazon Web Services. Uses docker and Spark.
    * See [`docs/configuration/aws.md`](configuration/aws.md)
* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

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

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### Library strandedness
Three command line flags / config parameters set the library strandedness for a run:

* `--forward_stranded`
* `--reverse_stranded`
* `--unstranded`

If not set, the pipeline will be run as unstranded. Specifying `--pico` makes the pipeline run in `forward_stranded` mode.

You can set a default in a cutom Nextflow configuration file such as one saved in `~/.nextflow/config` (see the [nextflow docs](https://www.nextflow.io/docs/latest/config.html) for more). For example:

```groovy
params {
    reverse_stranded = true
}
```

If you have a default strandedness set in your personal config file you can use `--unstranded` to overwrite it for a given run.

These flags affect the commands used for several steps in the pipeline - namely HISAT2, featureCounts, RSeQC (`RPKM_saturation.py`) and StringTie:

* `--forward_stranded`
  * HISAT2: `--rna-strandness F` / `--rna-strandness FR`
  * featureCounts: `-s 1`
  * RSeQC: `-d ++,--` / `-d 1++,1--,2+-,2-+`
  * StringTie: `--fr`
* `--reverse_stranded`
  * HISAT2: `--rna-strandness R` / `--rna-strandness RF`
  * featureCounts: `-s 2`
  * RSeQC: `-d +-,-+` / `-d 1+-,1-+,2++,2--`
  * StringTie: `--rf`

## FeatureCounts Extra Gene Names
By default, the pipeline uses `gene_names` as additional gene identifiers apart from ENSEMBL identifiers in the pipeline.
This behaviour can be modified by specifying `--fcExtraAttributes` when running the pipeline, which is passed on to featureCounts as an `--extraAttributes` parameter.
See the user guide of the [Subread package here](http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf).
Note that you can also specify more than one desired value, separated by a comma:
``--fcExtraAttributes gene_id,...``

## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human GRCh37 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the alignment tool instead. Developed by the same group behind the popular Tophat aligner, HISAT2 has a much smaller memory footprint.

To use HISAT2, use the parameter `--aligner hisat2` or set `params.aligner = 'hisat2'` in your config file.

## Reference Genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If you are running on UPPMAX, these should work without any additional configuration. If running on AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)
There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](https://github.com/nf-core/rnaseq/blob/master/conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```groovy
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
--star_index '[path to STAR index]' \
--hisat2_index '[path to HISAT2 index]' \
--fasta '[path to Fasta reference]' \
--gtf '[path to GTF file]' \
--bed12 '[path to bed12 file]'
```

Note that only one of `--star_index` / `--hisat2_index` are needed depending on which aligner you are using (see below).

The minimum requirements are a Fasta and GTF file. If these are provided and no others, then all other reference files will be automatically generated by the pipeline.

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveAlignedIntermediates`
As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM files from STAR / HISAT2 and sorting steps.

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
Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.


## Library Prep Presets
Some command line options are available to automatically set parameters for common RNA-seq library preparation kits.

> Note that these presets override other command line arguments. So if you specify `--pico --clip_r1 0`, the `--clip_r1` bit will be ignored.

If you have a kit that you'd like a preset added for, please let us know!

### `--pico`
Sets trimming and standedness settings for the _SMARTer Stranded Total RNA-Seq Kit - Pico Input_ kit.

Equivalent to: `--forward_stranded` `--clip_r1 3` `--three_prime_clip_r2 3`


## Skipping QC steps
The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

* `--skip_qc` -                Skip **all QC steps**, apart from MultiQC
* `--skip_fastqc` -            Skip FastQC
* `--skip_rseqc` -             Skip RSeQC
* `--skip_genebody_coverage` - Skip calculating the genebody coverage
* `--skip_preseq` -            Skip Preseq
* `--skip_dupradar` -          Skip dupRadar (and Picard MarkDups)
* `--skip_edger` -             Skip edgeR MDS plot and heatmap
* `--skip_multiqc` -           Skip MultiQC

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

###
## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

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

Note - you can use this to override defaults. For example, if you don't want FastQC errors to be ignored, you can specify a config file using `-c` that contains the following:

```groovy
process.$fastqc.errorStrategy = 'terminate'
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--hisatBuildMemory`
Required amount of memory in GB to build HISAT2 index with splice sites.
The HiSAT2 index build can proceed with or without exon / splice junction information.
To work with this, a very large amount of memory is required.
If this memory is not available, the index build will proceed without splicing information.
The `--hisatBuildMemory` option changes this threshold. By default it is `200GB` - if your system
`--max_memory` is set to `128GB` but your genome is small enough to build using this, then you can
allow the exon build to proceed by supplying `--hisatBuildMemory 100GB`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.

### `--clusterOptions`
Submit arbitrary cluster scheduler options (not available for all config profiles). For instance, you could use `--clusterOptions '-p devcore'` to run on the development node (though won't work with default process time requests).

## Stand-alone scripts
The `bin` directory contains some scripts used by the pipeline which may also be run manually:

* `gtf2bed`
  * Script used to generate the BED12 reference files used by RSeQC. Takes a `.gtf` file as input
* `dupRadar.r`
  * dupRadar script used in the _dupRadar_ pipeline process.
* `edgeR_heatmap_MDS.r`
  * edgeR script used in the _Sample Correlation_ process
* `RNA-pipeline-from-BAM.sh`
  * SLURM script used to mimic pipeline QC steps, taking an aligned BAM file as input.
  * Potentially unmaintained, use at your own risk!
