# RNA-BP
Pipeline for RNA sequencing best practice analysis at the NGI at Scilifelab Stockholm, Sweden

Written by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn)

> **See the [pipeline documentation](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/Docs/README.md) for explanations of the results files.**

## Installation
### NextFlow installation
To use this pipeline, you need to have a working version of NextFlow installed. You can find more
information about this pipeline tool at [nextflow.io](http://www.nextflow.io/). The typical installation
of NextFlow looks like this:

```
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```
Note that if you're running on the Swedish UPPMAX cluster (Milou) you can load NextFlow as an
environment module:
```
module load nextflow
```

### NextFlow configuration
Next, you need to set up a config file so that NextFlow knows how to run and where to find reference
indexes. You can find an example configuration file for UPPMAX (milou) with this repository:
[`example_uppmax_config`](https://github.com/SciLifeLab/NGI-RNAseq/blob/master/example_uppmax_config).

Copy this file to `~/.nextflow/config` and edit the line `'-A b2013064'` to contain your own UPPMAX project
identifier instead.

It is entirely possible to run this pipeline on other clusters - just note that you may need to customise
the `process` environment (eg. if you're using a cluster system other than SLURM) and the paths to reference
files.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-RNAseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:
```
git clone https://github.com/SciLifeLab/NGI-RNAseq.git
nextflow run NGI-RNAseq/main.nf
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```
nextflow run SciLifeLab/NGI-RNAseq --reads '*_R{1,2}.fastq.gz' --genome 'GRCm38'
```
or using a more manual approach ( require you to clone the git repository)

```
nextflow path_to_NGI-RNAseq/main.nf -c path_to_your_nextflow_config --reads '*_R{1,2}.fastq.gz' --genome 'GRCm38'
```

Note that the pipeline will create files in your working directory:
```bash
work            # Directory containing the nextflow working files
results         # Finished results for each sample, one directory per pipeline step
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
The human `GRCh37` genome is set as default.
```
--genome 'GRCm38'
```
The `example_uppmax_config` file currently has the location of references for `GRCh37` (Human), `GRCm38` (Mouse)
and `sacCer2` (Yeast).

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap, which require at least three samples to work. I.e use this when
running on one or two sample only.

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


### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.

## Credits
These scripts were written for use at the
[National Genomics Infrastructure](https://portal.scilifelab.se/genomics/)
at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
For more information, please get in touch with
[Phil Ewels](https://github.com/ewels).

<p align="center"><a href="stand_alone/http://www.scilifelab.se/" target="_blank"><img src="Docs/images/SciLifeLab_logo.png" title="SciLifeLab"></a></p>
