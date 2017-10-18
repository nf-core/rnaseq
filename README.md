# ![NGI-RNAseq](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/NGI-RNAseq_logo.png)

[![Build Status](https://travis-ci.org/SciLifeLab/NGI-RNAseq.svg?branch=master)](https://travis-ci.org/SciLifeLab/NGI-RNAseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg?style=flat-square)](https://gitter.im/SciLifeLab/NGI-RNAseq)

### Introduction

NGI-RNAseq is a bioinformatics best-practice analysis pipeline used for RNA sequencing data at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline is primarily used with a SLURM cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se). However, the pipeline should be able to run on any system that Nextflow supports. We have done some limited testing using Docker and AWS, and the pipeline comes with some configuration for these systems. See the [installation docs](docs/installation.md) for more information.

### Documentation
The NGI-RNAseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)
4. [Troubleshooting](Troubleshooting.md)

If you're interested in running the pipeline in the cloud, please read the docs about using the pipeline with [Amazon Web Services](docs/amazon_web_services.md).

### Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
The pipeline was developed by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)). Docker and AWS integration was led by Denis Moreno ([@Galithil](https://github.com/Galithil)) and Phil Ewels ([@ewels](https://github.com/ewels)).

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
