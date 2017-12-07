# NGI-RNAseq Installation

To start using the NGI-RNAseq pipeline, there are three steps described below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)
3. Configure the pipeline
    * [Swedish UPPMAX System](#31-configuration-uppmax)
    * [Swedish Hebbe (C3SE) System](#32-configuration-hebbe-c3se)
    * [Other Clusters](#33-configuration-other-clusters)
    * [Docker](#34-configuration-docker)
    * [Amazon AWS](#35-configuration-amazon-ec2)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

**You need NextFlow version >= 0.24 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `SciLifeLab/NGI-RNAseq` is specified as the pipeline name.

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/SciLifeLab/NGI-RNAseq/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/NGI-RNAseq-master
```

## 3.1) Configuration: UPPMAX
The pipeline comes bundled with configurations to use the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (tested on `milou`, `rackham`, `bianca` and `irma`). As such, you shouldn't need to add any custom configuration - everything _should_ work out of the box.

To use the pipeline on UPPMAX, you **must** specificy `-profile uppmax` when running the pipeline (as of Nov 2017).

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```

## 3.2) Configuration: Hebbe (C3SE)
This pipeline has been successfully used on the [Hebbe cluster](http://www.c3se.chalmers.se/index.php/Hebbe) in Gothenburg, though it requires significantly more setup work than at UPPMAX. This is mainly due to the fact that none of the required software is pre-installed.

To use, follow the steps described below ([3.3) Configuration: Other clusters](#33-configuration-other-clusters)) to install and configure the required software. Create a config file in your home directory (`~/.nextflow/config`) with paths to your reference genome indices (see below for [instructions](#reference-genomes)). Finally, run the pipeline with `-profile hebbe --project [project-id]`. This will launch the [hebbe config](../conf/hebbe.config) which has been pre-configured with a setup suitable for the Hebbe cluster. Note that to date it has only been tested on Yeast data - if jobs are failing due to insufficient resources, please [let us know](https://github.com/SciLifeLab/NGI-RNAseq/issues) and we will update it accordingly.

## 3.3) Configuration: Other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the script knows where to find your reference files and how your cluster works.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config`.

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile with [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

### Cluster Environment
By default, pipeline uses the `local` Nextflow executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```groovy
process {
  executor = 'YOUR_SYSTEM_TYPE'
}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```groovy
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```

### Reference Genomes
The NGI-RNAseq pipeline needs a reference genome for alignment and annotation. If not already available, start by downloading the relevant reference, for example from [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

The minimal requirements are a FASTA file and a GTF file. If STAR and BED12 references are also specified, the pipeline won't have to generate them and will run faster. Use the command line option `--saveReference` to keep the generated references so that they can be added to your config and used again in the future.

Reference genome paths can be specified on the command line each time you run with `--star_index`, `--fasta`, `--gtf` and `--bed12`. Alternatively, add the paths to the config under a relevant id and just specify this id with `--genome ID` when you run the pipeline _(this can also be set as a default in your config)_:

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


### Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system.

#### Docker Image
Docker is a software tool that allows you to run your analysis inside a "container" - basically a miniature self-contained software environment. We've already made a docker image for you, so if you can run docker and nextflow then you don't need to worry about any other software dependencies.

The pipeline comes with a script to build a docker image. This runs automatically and creates a hosted docker image that you can find here: https://hub.docker.com/r/scilifelab/ngi-rnaseq/

If you run the pipeline with `-profile docker` or `-with-docker 'scilifelab/ngi-rnaseq'` then nextflow will download this docker image automatically and run using this.

Note that the docker images are tagged with version as well as the code, so this is a great way to ensure reproducibility. You can specify pipeline version when running with `-r`, for example `-r v1.3`. This uses pipeline code and docker image from this tagged version.

#### Singularity image
Many HPC environments are not able to run Docker due to problems with needing administrator privileges. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub. The UPPMAX configuration uses Singularity by default, meaning no problems with software dependencies and great reproducibility.

To use the singularity image with a different config, use `-with-singularity 'docker://scilifelab/ngi-rnaseq'` when running the pipeline.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name ngi-rnaseq.img docker://scilifelab/ngi-rnaseq
```

Then transfer this file and run the pipeline with this path:

```
nextflow run /path/to/NGI-RNAseq -with-singularity /path/to/ngi-rnaseq.img
```

#### Environment Modules
If your cluster uses environment modules, you can use the pipeline with these. There is a bundled config file to use these on UPPMAX (as was done in earlier versions of this pipeline). To use this, run the pipeline with `-profile uppmax_modules`.

If running on another system, add lines to your custom config file as follows _(customise module names and versions as appropriate)_:

```groovy
process {
  $makeSTARindex.module = ['star']
  $makeHisatSplicesites.module = ['HISAT2']
  $makeHISATindex.module = ['HISAT2']
  $fastqc.module = ['FastQC']
  $trim_galore.module = ['FastQC', 'TrimGalore']
  $star.module = ['star']
  $hisat2Align.module = ['samtools/1.3', 'HISAT2']
  $hisat2_sortOutput.module = ['samtools/1.3']
  $rseqc.module = ['rseqc', 'samtools/1.3']
  $preseq.module = ['preseq']
  $markDuplicates.module = ['picard/2.0.1']
  $dupradar.module = ['R/3.2.3']
  $featureCounts.module = ['subread']
  $stringtieFPKM.module = ['StringTie/1.2.0']
  $sample_correlation.module = ['R/3.2.3']
  $multiqc.module = ['MultiQC']
}
```

##### R Package Location
If you are using a central installation of R, you may not have write permissions for installing custom modules. If this is the case, add the following to your Nextflow configuration file to specify where these files should be saved:

```groovy
params {
  rlocation = "$HOME/R/nxtflow_libs/" // or any path
}
```

#### Manual Installation
If the software is not already available, you will need to install it.

If you are able to use [Docker](https://www.docker.com/), you can use the [sclifelab/ngi-rnaseq](https://hub.docker.com/r/scilifelab/ngi-rnaseq/) image which comes with all requirements. This is pulled by Nextflow automatically if you use `-profile docker` (see below for [further instructions](#34-configuration-docker)).

We recommend using [Bioconda](https://bioconda.github.io/) to install the required software as the process is quite easy in our experience. This can be done as follows:

##### 1) Install miniconda in your home directory
``` bash
cd
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

##### 2) Add the bioconda conda channel (and others)
```bash
conda config --add channels anaconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels salilab
```

##### 2) Create a conda environment, with all necessary packages:
```bash
conda create --name rna_seq_py2.7 python=2.7
source activate rna_seq_py2.7
conda install --yes \
    bioconductor-dupradar=1.2.2 \
    r-essentials \
    samtools \
    star=2.5.3a \
    bedtools=2.26.0 \
    trim-galore=0.4.1 \
    fastqc=0.11.5 \
    rseqc=2.6.4 \
    multiqc \
    gsl=1.16 \
    preseq=2.0.2 \
    picard=2.9.0 \
    subread=1.5.0.post3 \
    nextflow=0.23.4 \
    pysam=0.10.0 \
    stringtie=1.3.3 \
    graphviz=2.38.0 \
    hisat2=2.0.5
```
_(Feel free to adjust versions as required.)_

##### 3) Set up Picard
Picard requires the `PICARD_HOME` environment variable to be set. In some cases, a temporary directory must also be specified (if the default does not have enough available space). To automatically set and unset these when you activate and deactivate your conda environment, do the following:

```bash
cd ~/miniconda3/envs/rna_seq_py2.7 # Or path to your conda environment
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```

Put in `./etc/conda/activate.d/env_vars.sh`:
```bash
#!/bin/sh
export PICARD_HOME="$HOME/miniconda3/envs/rna_seq_py2.7/share/picard-2.9.0-0/"
export _JAVA_OPTIONS=-Djava.io.tmpdir='/path/to/tmp'
```

Then add to `./etc/conda/deactivate.d/env_vars.sh`:
```bash
#!/bin/sh
unset PICARD_HOME
unset _JAVA_OPTIONS
```

##### 4) Usage
Once created, the conda environment can be activated before running the pipeline and deactivated afterwards:

```bash
source activate rna_seq_py2.7
# run pipeline
source deactivate
```

## 3.4) Configuration: Docker
Docker is a great way to run NGI-RNAseq, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required.

First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run SciLifeLab/NGI-RNAseq -profile docker --reads '<path to your reads>' --fasta '<path to fasta ref>' --gtf '<path to gtf>'
```

Nextflow will recognise `SciLifeLab/NGI-RNAseq` and download the pipeline from GitHub. The `-profile docker` configuration lists the [sclifelab/ngi-rnaseq](https://hub.docker.com/r/scilifelab/ngi-rnaseq/) image that we have created and is hosted at dockerhub, and this is downloaded.

A reference genome is still required by the pipeline. Specifying paths to FASTA and GTF files is the minimum requirement, STAR and BED12 references will automatically be generated. See the above [Reference Genomes](#reference-genomes) documentation for instructions on how to configure Nextflow with preset paths to make this easier.

A test suite for docker comes with the pipeline, and can be run by moving to the [`tests` directory](https://github.com/ewels/NGI-RNAseq/tree/master/tests) and running `./docker_test.sh`. This will download a small yeast genome and some data, and attempt to run the pipeline through docker on that small dataset. This is automatically run using [Travis](https://travis-ci.org/SciLifeLab/NGI-RNAseq/) whenever changes are made to the pipeline.

## 3.5) Configuration: Amazon EC2
There are multiple ways of running this pipeline over Amazon's EC2 service. Please see the [NGI-RNAseq with AWS](amazon_web_services.md) documentation for more information.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
