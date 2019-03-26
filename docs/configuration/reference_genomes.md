# nf-core/rnaseq: Reference Genomes Configuration

The nf-core/rnaseq pipeline needs a reference genome for alignment and annotation.

These paths can be supplied on the command line at run time (see the [usage docs](../usage.md)),
but for convenience it's often better to save these paths in a nextflow config file.
See below for instructions on how to do this.
Read [Adding your own system](adding_your_own.md) to find out how to set up custom config files.

The minimal requirements are a FASTA file and a GTF file. If STAR and BED12 references are also specified, the pipeline won't have to generate them and will run faster. Use the command line option `--saveReference` to keep the generated references so that they can be added to your config and used again in the future.

Reference genome paths can be specified on the command line each time you run with `--star_index`, `--hisat_index`, `--fasta`, `--gtf` and `--bed12`. Fasta is only required if building a STAR or HISAT2 index.

## Adding paths to a config file
Specifying long paths every time you run the pipeline is a pain.
To make this easier, the pipeline comes configured to understand reference genome keywords which correspond to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline.

Note that this genome key can also be specified in a config file if you always use the same genome.

To use this system, add paths to your config file using the following template:

```nextflow
params {
  genomes {
    'YOUR-ID' {
      bed12  = '<PATH TO BED FILE>/genes.bed'
      fasta  = '<PATH TO FASTA FILE>/genome.fa'
      gtf    = '<PATH TO GTF FILE>/genes.gtf'
      star    = '<PATH TO STAR INDEX>/STARIndex/'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```

You can add as many genomes as you like as long as they have unique IDs.

## illumina iGenomes
To make the use of reference genomes easier, illumina has developed a centralised resource called [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).
Multiple reference index types are held together with consistent structure for multiple genomes.

We have put a copy of iGenomes up onto AWS S3 hosting and this pipeline is configured to use this by default.
The hosting fees for AWS iGenomes are currently kindly funded by a grant from Amazon.
The pipeline will automatically download the required reference files when you run the pipeline.

For more information about AWS iGenomes, see [https://ewels.github.io/AWS-iGenomes/](https://ewels.github.io/AWS-iGenomes/).

Downloading the files takes time and bandwidth, so we recommend making a local copy of the iGenomes resource.
Once downloaded, you can customise the variable `params.igenomes_base` in your custom configuration file to point to the reference location.
For example:

```nextflow
params.igenomes_base = '/path/to/data/igenomes/'
```
