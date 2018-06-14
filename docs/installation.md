# nfcore/rnaseq Installation

To start using the nfcore/rnaseq pipeline, there are three steps described below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)
3. Configure the pipeline
    * [Local configuration](configuration/local.md)
    * [Amazon Web Services (aws)](configuration/aws.md)
    * [Swedish UPPMAX clusters](configuration/uppmax.md)
    * [Swedish cs3e Hebbe cluster](configuration/c3se.md)
    * [Tübingen QBiC](configuration/qbic.md)
    * [Adding your own system](configuration/adding_your_own.md)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

**You need NextFlow version >= 0.24 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `nf-core/rnaseq` is specified as the pipeline name.

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/nf-core/rnaseq/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/rnaseq-master
```
