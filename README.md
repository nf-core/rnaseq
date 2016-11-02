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
We are in the process of setting up a Docker image with the pipeline requirements, for easier
use and better reproducibility. Check back soon for more information on this!

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

### `--index`, `--gtf`, `--bed12`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:
```
--index [path to STAR index]   --gtf [path to GTF file]   --bed12 [path to bed12 file]
```

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

## Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
Written by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn).

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
