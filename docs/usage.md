# NGI-RNAseq Usage

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run SciLifeLab/NGI-RNAseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37
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

Note that the `{1,2}` parentheses are required to specify paired end data. The file path should be
in quotation marks to prevent shell glob expansion.

If left unspecified, the pipeline will assume that the data is in a directory called `data` in the
working directory (`data/*{1,2}.fastq.gz`).

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, specify `--singleEnd`
on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks,
can then be used for `--reads`. For example: `--singleEnd --reads '*.fastq'`

It is not possible to run a mixture of single-end and paired-end files in one run.

### Library strandedness
Three command line flags / config parameters set the library strandedness for a run:

* `--forward_stranded`
* `--reverse_stranded`
* `--unstranded`

If not set, the pipeline will be run as unstranded. The UPPMAX configuration file sets `reverse_stranded` to true by default. Use `--unstranded` or `--forward_stranded` to overwrite this. Specifying `--pico` makes the pipeline run in `forward_stranded` mode.

These flags affect the commands used for several steps in the pipeline - namely HISAT2, featureCounts, RSeQC (`RPKM_saturation.py`)
and StringTie:
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

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveAlignedIntermediates`
As above, by default intermediate BAM files will not be saved. The final BAM files created
after the Picard MarkDuplicates step are always saved. Set to true to also copy out BAM
files from STAR / HISAT2 and sorting steps.

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


## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits on UPPMAX with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). Useful if using different UPPMAX
projects or different sets of reference genomes.

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, we run on UPPMAX but don't want to use the MultiQC
environment module as is the default. So we specify a config file using `-c` that contains the following:

```groovy
process.$multiqc.module = []
```

### `--rlocation`
Some steps in the pipeline run R with required modules. By default, the pipeline will install
these modules to `~/R/nxtflow_libs/` if not present. You can specify what path to use with this
command line flag.

### `--multiqc_config`
If you would like to supply a custom config file to MultiQC, you can specify a path with `--multiqc_config`. This is used instead of the config file specific to the pipeline.

### `--clusterOptions`
Submit arbitrary SLURM options (UPPMAX profile only). For instance, you could use `--clusterOptions '-p devcore'`
to run on the development node (though won't work with default process time requests).

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


---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
