# NGI-RNAseq
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
nextflow run SciLifeLab/NGI-RNAseq -profile docker -c <your_config>  --reads '<path to your reads>'
```

Note that you will need to configure the reference genomes to use. See below
for more instructions.

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

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
