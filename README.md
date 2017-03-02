# NGI-RNAseq

[![Build Status](https://travis-ci.org/SciLifeLab/NGI-RNAseq.svg?branch=master)](https://travis-ci.org/SciLifeLab/NGI-RNAseq)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg?style=flat-square)](https://gitter.im/SciLifeLab/NGI-RNAseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.22.2-brightgreen.svg
)](https://www.nextflow.io/) 



Pipeline for RNA sequencing best practice analysis at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.


## Pipeline Results
See the [pipeline documentation](docs/README.md)
for explanations of the results files.

## Installation
### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure
Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-RNAseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```
git clone https://github.com/SciLifeLab/NGI-RNAseq.git
nextflow run NGI-RNAseq/main.nf
```

## Configuration
By default, the pipeline is configured to run on the Swedish UPPMAX cluster (milou / irma).

You will need to specify your UPPMAX project ID when running a pipeline. To do this, use
the command line flag `--project <project_ID>`.

To avoid having to specify this every time you run Nextflow, you can add it to your
personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID'
```

The pipeline will exit with an error message if you try to run it pipeline with the default
UPPMAX config profile and don't set project.

### Running using Docker
First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

You can now run:
```bash
nextflow run SciLifeLab/NGI-RNAseq -profile docker --reads '<path to your reads>' --fasta <path to the genome's fasta file> --gtf <path to the genome's gtf file>
```
The fasta and GTF parameters can be specified in a configuration file (you can provide it with -c), look into the `conf/docker_test.config` for an example.

The docker image containting all the required tools will do downloaded on during the run. It can be found on [DockerHub](https://hub.docker.com/r/scilifelab/ngi-rnaseq/)

A test suite for docker has been implemented, and can be run by moving to the `tests` folder and running `./docker_test.sh`. This will download a small yeast genome and some data, and attempt to run the pipeline through docker on that small dataset.

### Running on other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up
your own config file so that the script knows where to find your reference files and how your
cluster works.

Copy the contents of [`conf/uppmax.config`](conf/uppmax.config) to your own config file somewhere
and then reference it with `-c` when running the pipeline.

If you think that there are other people using the pipeline who would benefit from your configuration
(eg. other common cluster setups), please let us know. It should be easy to create a new config file
in `conf` and reference this as a named profile in [`nextflow.config`](nextflow.config). Then these
configuration options can be used by specifying `-profile <name>` when running the pipeline.

### Running on Amazon EC2
There are multiple ways of running this pipeline over Amazon's EC2 service.

The simplest way consists of creating an EC2 instance and running the docker flavour of this pipeline over the machine.

A slightly more complex way is to use our prebuilt AMI to create a ready-to-go virtual machine. The AMI is called `scilifelab/ngi-rnaseq`, id : `ami-f23c6081`. It is available in the Ireland region (eu-west-1). This AMI comes with all the tools installed, including docker and nextflow.  The pipeline can then be run by usint the following : 
```
nextflow run SciLifeLab/NGI-RNAseq -profile base --reads '<path to your reads>' --fasta <path to the genome's fasta file> --gtf <path to the genome's gtf file>
```
As for Docker, the two last parameters can be entered in a configuration file that would be supplied with -c. An example would be `conf/amazon_test.config`

A third approach will be to use nextflow to generate a cluster of machines, and run the pipeline there. This is currently not implemented.


## Running the pipeline
The typical command for running the pipeline is as follows:
```
nextflow run SciLifeLab/NGI-RNAseq --reads '*_R{1,2}.fastq.gz'
```

Note that the pipeline will create files in your working directory:
```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### `--reads`
Location of the input FastQ files:
```
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

**NB: Must be enclosed in quotes!**

Note that the `{1,2}` parentheses are required to specify paired end data. Running `--reads '*.fastq'` will treat
all files as single end. The file path should be in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the working directory.

## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads
to the reference genome. STAR is fast and common, but requires a great deal of RAM to run, typically
around 38GB for the Human GRCh37 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the
alignment tool instead. Thought of as the successor to Tophat by many, HISAT2 has a much smaller
memory footprint.

To use HISAT2, use the parameter `--aligner hisat2` or set `params.aligner = 'hisat2'` in your config file.

## Reference Genomes

### `--genome`
The reference genome to use of the analysis, needs to be one of the genome specified in the config file.
The human `GRCh37` genome is used by default.

See [`conf/uppmax.config`](conf/uppmax.config) for a list of the supported reference genomes
and their keys. Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* Drosophila
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

If you're not running on UPPMAX (the default profile), you can create your own config
file with paths to your reference genomes.
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)
for instructions on where to add this.

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

### `--star_index`, `--fasta`, `--gtf`, `--bed12`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:
```
--star_index [path to STAR index] \
--fasta [path to Fasta reference] \
--gtf [path to GTF file] \
--bed12 [path to bed12 file]
```

### `--downloadFasta`, `--downloadGTF`
If no STAR / Fasta reference is supplied, a URL can be supplied to download a Fasta file
at the start of the pipeline. The same with a GTF reference file. A required STAR index
and BED12 files will then be generated from these downloaded files.

### `--saveReference`
Supply this parameter to save any generated reference genome files to your results folder.
These can then be used for future pipeline runs, reducing processing times.

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs,
memory and time. For most of the steps in the pipeline, if the job exits
on UPPMAX with an error code of `143` (exceeded requested resources) it will
automatically resubmit with higher requests (2*original, then 3*original).
If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default
value can be overwritten with config variables. These can be set on the command
line or in a config file. The names are set as `[process name]_[resource type]`,
for example `star_memory`.

So, to override the defaults for STAR, you can do run the pipeline as follows:
```bash
nextflow run SciLifeLab/NGI-RNAseq --star_cpus 1 --star_memory '10 GB' --star_time '24h'
```

Alternative, these can be set in a config file:
```groovy
params {
  star_cpus = 1
  star_memory = '10 GB'
  star_time = '24h'
}
```

## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `--strandRule`
Some RSeQC jobs need to know the stranded nature of the library. By default, the pipeline will use
`++,--` for single end libraries and `1+-,1-+,2++,2--` for paired end libraries. These codes are for
strand specific libraries (antisense). `1+-,1-+,2++,2--` decodes as:

*  Reads 1 mapped to `+` => parental gene on `+`
*  Reads 1 mapped to `-` => parental gene on `-`
*  Reads 2 mapped to `+` => parental gene on `-`
*  Reads 2 mapped to `-` => parental gene on `+`

Use this parameter to override these defaults. For example, if your data is paired end and strand specific,
but same-sense to the reference, you could run:
```
nextflow run NGI-RNAseq/main.nf --strandRule '1++,1--,2+-,2-+'
```
Use `--strandRule 'none'` if your data is not strand specific.

### `--rlocation`
Some steps in the pipeline run R with required modules. By default, the pipeline will install
these modules to `~/R/nxtflow_libs/` if not present. You can specify what path to use with this
command line flag.

### `--clusterOptions`
Submit arbitrary SLURM options (UPPMAX profile only). For instance, you could use `--clusterOptions '-p devcore'`
to run on the development node (though won't work with default process time requests).

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes. **NB:** one hyphen only (core Nextflow parameter).

## Stand-alone scripts
There is a folder with some unmaintained standalone scripts that you can use:
[stand-alone scripts](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/stand-alone-scripts).
Currently it only contains one file, a `sbatch` script that starts after the alignment step with BAM files.

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn).

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a>
<a href="stand_alone/https://www.scilifelab.se/platforms/genomics/" target= _blank><img src="docs/images/NGI-final-small.png" title="NGI"></a>
</p>

