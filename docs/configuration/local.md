# nfcore/rnaseq: Local Configuration

If running the pipeline in a local environment, we highly recommend using either Docker or Singularity.

## Docker
Docker is a great way to run nfcore/rnaseq, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required. The nfcore/rnaseq profile comes with a configuration profile for docker, making it very easy to use. This also comes with the required presets to use the AWS iGenomes resource, meaning that if using common reference genomes you just specify the reference ID and it will be autaomtically downloaded from AWS S3.

First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run nf-core/rnaseq -profile docker --genome '<genome ID>' --reads '<path to your reads>'
```

Nextflow will recognise `nf-core/rnaseq` and download the pipeline from GitHub. The `-profile docker` configuration lists the [nfcore/rnaseq](https://hub.docker.com/r/nfcore/rnaseq/) image that we have created and is hosted at dockerhub, and this is downloaded.

For more information about how to work with reference genomes, see [`docs/configuration/reference_genomes.md`](docs/configuration/reference_genomes.md).

### Pipeline versions
The public docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r 1.0`. This uses pipeline code and docker image from this tagged version.


## Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

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
