# ![nf-core/rnaseq](docs/images/nf-core-rnaseq_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/rnaseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rnaseq/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnaseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rnaseq/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnaseq/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1400710-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1400710)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnaseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnaseq)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnaseq** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation.

On release, automated continuous integration tests run the pipeline on a [full-sized dataset](https://github.com/nf-core/test-datasets/tree/rnaseq#full-test-dataset-origin) obtained from the ENCODE Project Consortium on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from running the full-sized tests individually for each `--aligner` option can be viewed on the [nf-core website](https://nf-co.re/rnaseq/results) e.g. the results for running the pipeline with `--aligner star_salmon` will be in a folder called `aligner_star_salmon` and so on.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Download FastQ files via SRA, ENA or GEO ids and auto-create input samplesheet ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html); *if required*)
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
5. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
6. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
7. Choice of multiple alignment and quantification routes:
    1. [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
    2. [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
    3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
8. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
9. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
10. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
11. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
12. Create bigWig coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
13. Extensive quality control:
    1. [`RSeQC`](http://rseqc.sourceforge.net/)
    2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
    3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
    4. [`Preseq`](http://smithlabresearch.org/software/preseq/)
    5. [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
14. Pseudo-alignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/); *optional*)
15. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

> **NB:** Quantification isn't performed if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments. However, you can use this route if you have a preference for the alignment, QC and other types of downstream analysis compatible with the output of HISAT2.
> **NB:** The `--aligner star_rsem` option will require STAR indices built from version 2.7.6a or later. However, in order to support legacy usage of genomes hosted on AWS iGenomes the `--aligner star_salmon` option requires indices built with STAR 2.6.1d or earlier. Please refer to this [issue](https://github.com/nf-core/rnaseq/issues/498) for further details.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`).

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_. Note: This pipeline does not currently support running with Conda on macOS if the `--remove_ribo_rna` parameter is used because the latest version of the SortMeRNA package is not available for this platform.

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/rnaseq -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    * Typical command for RNA-seq analysis:

        ```bash
        nextflow run nf-core/rnaseq \
            --input samplesheet.csv \
            --genome GRCh37 \
            -profile <docker/singularity/podman/conda/institute>
        ```

    * Typical command for downloading public data:

        ```bash
        nextflow run nf-core/rnaseq \
            --public_data_ids ids.txt \
            -profile <docker/singularity/podman/conda/institute>
        ```

    > **NB:** The commands to obtain public data and to run the main arm of the pipeline are completely independent. This is intentional because it allows you to download all of the raw data in an initial pipeline run (`results/public_data/`) and then to curate the auto-created samplesheet based on the available sample metadata before you run the pipeline again properly.

See [usage](https://nf-co.re/rnaseq/usage) and [parameter](https://nf-co.re/rnaseq/parameters) docs for all of the available options when running the pipeline.

## Documentation

The nf-core/rnaseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/rnaseq/usage) and [output](https://nf-co.re/rnaseq/output).

## Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

The pipeline was re-written in Nextflow DSL2 by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

Many thanks to other who have helped out along the way too, including (but not limited to):
[@Galithil](https://github.com/Galithil),
[@pditommaso](https://github.com/pditommaso),
[@orzechoj](https://github.com/orzechoj),
[@apeltzer](https://github.com/apeltzer),
[@colindaven](https://github.com/colindaven),
[@lpantano](https://github.com/lpantano),
[@olgabot](https://github.com/olgabot),
[@jburos](https://github.com/jburos).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rnaseq` channel](https://nfcore.slack.com/channels/rnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/rnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1400710](https://doi.org/10.5281/zenodo.1400710)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
