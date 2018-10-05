# nfcore/rnaseq: Configuration for other clusters

It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the pipeline knows how to work with your cluster.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config` (see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more).

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile - see [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

## Reference genomes
Remember that you will need to define a reference genome to use. See [reference_genomes.md](reference_genomes.md) for instructions.

## Cluster Environment
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


## Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system. If possible, we _highly_ recommend using either Docker or Singularity.

### Docker
Docker is a great way to run nfcore/rnaseq, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required.

First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run nf-core/rnaseq -profile docker --reads '<path to your reads>' --fasta '<path to fasta ref>' --gtf '<path to gtf>'
```

Nextflow will recognise `nf-core/rnaseq` and download the pipeline from GitHub. The `-profile docker` configuration lists the [nfcore/rnaseq](https://hub.docker.com/r/nfcore/rnaseq/) image that we have created and is hosted at dockerhub, and this is downloaded.

The public docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r 1.0`. This uses pipeline code and docker image from this tagged version.

To add docker support to your own config file (instead of using the `docker` profile, which runs locally), add the following:

```groovy
docker {
  enabled = true
}
process {
  container = wf_container
}
```

The variable `wf_container` is defined dynamically and automatically specifies the image tag if Nextflow is running with `-r`.

A test suite for docker comes with the pipeline, and can be run by moving to the [`tests` directory](https://github.com/nf-core/rnaseq/tree/master/tests) and running `./docker_test.sh`. This will download a small yeast genome and some data, and attempt to run the pipeline through docker on that small dataset. This is automatically run using [Travis](https://travis-ci.org/nf-core/rnaseq/) whenever changes are made to the pipeline.

### Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

To specify singularity usage in your pipeline config file, add the following:

```groovy
singularity {
  enabled = true
}
process {
  container = "docker://$wf_container"
}
```

The variable `wf_container` is defined dynamically and automatically specifies the image tag if Nextflow is running with `-r`.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

> NB: The "tag" at the end of this command corresponds to the pipeline version.
> Here, we're pulling the docker image for version 1.0 of the nfcore/rnaseq pipeline
> Make sure that this tag corresponds to the version of the pipeline that you're using

```bash
singularity pull --name nfcore-rnaseq-1.0.img docker://nfcore/rnaseq:1.0
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/nfcore-rnaseq -with-singularity /path/to/nfcore-rnaseq-1.0.img
```

### Bioconda
The workflow comes with a conda environment definition - a file called
[`environment.yml`](../environment.yml) which lists conda channels and package names / versions.

To use conda for this pipeline, first make sure that you have conda installed. We recommend miniconda:
https://conda.io/miniconda.html

Next, create a new environment using the [`environment.yml`](../environment.yml) file:

```bash
# Download the environment.yml file
curl https://raw.githubusercontent.com/nf-core/rnaseq/master/environment.yml -o environment.yml

# Create a new conda environment using it
conda env create -f environment.yml

# Activate the new conda environment
source activate nfcore-rnaseq
```

> NB: The above link grabs the latest version of `environment.yml`. It's best to be running
> a tagged release of the pipeline - if so, make sure that you use the corresponding conda env file.


### Environment Modules
If you can't use docker or singularity, but your cluster uses environment modules, you can use the pipeline with these. There is a bundled config file to use these on UPPMAX (as was done in earlier versions of this pipeline) that can be used with `-profile uppmax_modules`.

To use environment modules in your own config, add lines to your custom config file as follows _(customise module names and versions as appropriate)_:

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
