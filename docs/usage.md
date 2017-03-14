# NGI-RNAseq Usage

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
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

```bash
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
The reference genome to use for the analysis, needs to be one of the genome specified in the config file. This is `False` by default and needs to be specified (unless index files are supplied, see below).

See [`conf/uppmax.config`](conf/uppmax.config) for a list of the supported reference genomes
and their keys. Common genomes that are supported are:

* Human
  * `--genome GRCh37`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
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

```bash
--star_index '[path to STAR index]' \
--fasta '[path to Fasta reference]' \
--gtf '[path to GTF file]' \
--bed12 '[path to bed12 file]'
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
automatically resubmit with higher requests (2 x original, then 3 x original).
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

```bash
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


-----------------------------------------------------------------------------------------

<p align="center"><a href="http://www.scilifelab.se/" target="_blank"><img src="images/SciLifeLab_logo.png" title="SciLifeLab"></a>
<a href="https://ngisweden.scilifelab.se/" target= _blank><img src="images/NGI-final-small.png" title="NGI" style="height:100px;"></a>
</p>
