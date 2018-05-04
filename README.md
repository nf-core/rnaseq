# ![nfcore/rnaseq](docs/images/nfcore-rnaseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/rnaseq.svg?branch=master)](https://travis-ci.org/nf-core/rnaseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/rnaseq.svg)](https://hub.docker.com/r/nfcore/rnaseq/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)


### Introduction

**nfcore/rnaseq** is a bioinformatics analysis pipeline used for RNA sequencing data.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)), aligns the reads ([STAR](https://github.com/alexdobin/STAR) or [HiSAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)), generates gene counts ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/), [StringTie](https://ccb.jhu.edu/software/stringtie/)) and performs extensive quality-control on the results ([RSeQC](http://rseqc.sourceforge.net/), [dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html), [Preseq](http://smithlabresearch.org/software/preseq/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [MultiQC](http://multiqc.info/)). See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Documentation
The nfcore/rnaseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Amazon Web Services (aws)](docs/configuration/aws.md)
    * [Swedish UPPMAX clusters](docs/configuration/uppmax.md)
    * [Swedish cs3e Hebbe cluster](docs/configuration/c3se.md)
    * [Tübingen QBiC](docs/configuration/qbic.md)
    * [CCGA Kiel](docs/configuration/ccga.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
These scripts were written at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
The pipeline was initially developed by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

Many thanks to other who have helped out along the way too, including (but not limited to):
[@Galithil](https://github.com/Galithil),
[@pditommaso](https://github.com/pditommaso),
[@orzechoj](https://github.com/orzechoj),
[@apeltzer](https://github.com/apeltzer),
[@colindaven](https://github.com/colindaven).

### Participating Institutes
nfcore/rnaseq is now used by a number of core sequencing and bioinformatics facilities. Some of these are listed below. If you use this pipeline too, please let us know in an issue and we will add you to the list.

<table>
  <tr>
    <td><img src="docs/images/NGI_logo.png" width="200"></td>
    <td>National Genomics Infrastructure (NGI), Sweden</td>
    <td>https://ngisweden.scilifelab.se/</td>
  </tr>
  <tr>
    <td><img src="docs/images/QBiC_logo.png" width="200"></td>
    <td>Quantitative Biology Center (QBiC), Germany</td>
    <td>https://portal.qbic.uni-tuebingen.de/portal/</td>
  </tr>
</table>
