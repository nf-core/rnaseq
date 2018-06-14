# nfcore/rnaseq: Reference Genomes Configuration

The nfcore/rnaseq pipeline needs a reference genome for alignment and annotation. If not already available, start by downloading the relevant reference, for example from [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

The minimal requirements are a FASTA file and a GTF file. If STAR and BED12 references are also specified, the pipeline won't have to generate them and will run faster. Use the command line option `--saveReference` to keep the generated references so that they can be added to your config and used again in the future.

Reference genome paths can be specified on the command line each time you run with `--star_index`, `--hisat_index`, `--fasta`, `--gtf` and `--bed12`. Fasta is only required if building a STAR or HISAT2 index.

## Adding paths to a config file
Specifying long paths every time you run the pipeline is a pain. To make this easier, the pipeline comes configured to understand reference genome keywords which correspond to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline. This method is used on known systems such as UPPMAX.

Note that this genome key can also be specified in a config file if you always use the same genome.

To use this system, add paths to your config file using the following template:

```groovy
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
To make the use of reference genomes easier, illumina has developed a centralised resource called [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Multiple reference index types are held together with consistent structure for multiple genomes.

If possible, we recommend making this resource available on your cluster. We have put a copy of iGenomes up onto AWS S3 hosting and this pipeline is configured to use this for some profiles (`uppmax`, `docker`, `aws`). These profiles will automatically pull the required reference files when you run the pipeline.

To add iGenomes to your config file, add the following line to the end of your config file:

```groovy
includeConfig '/path/to/nfcore-rnaseq/conf/igenomes.config'
```

This works best when you have a `profile` set up in the pipeline - see [`nextflow.config`](../../nextflow.config).

The hosting fees for AWS iGenomes are currently funded by a grant from Amazon. We hope that this work will be extended past the end of the grant expiry date (mid 2018), but we can't be sure at this point.

For more information about the AWS iGenomes, see https://ewels.github.io/AWS-iGenomes/
