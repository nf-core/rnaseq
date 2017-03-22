# NGI-RNAseq Installation

### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure
Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when run if
`SciLifeLab/NGI-RNAseq` is specified as the pipeline name.

If you prefer, you can download the files yourself from GitHub and run them directly:

```bash
git clone https://github.com/SciLifeLab/NGI-RNAseq.git
nextflow run NGI-RNAseq/main.nf
```

## Configuration
By default, the pipeline is configured to run on the Swedish UPPMAX cluster (milou / irma).

You will need to specify your UPPMAX project ID when running a pipeline. To do this, use
the command line flag `--project <project_ID>`.

To avoid having to specify this every time you run Nextflow, you can add it to your
personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID'
```

The pipeline will exit with an error message if you try to run it pipeline with the default
UPPMAX config profile and don't set project.

### Running using Docker
First, install docker on your system : [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

You can now run:
```bash
nextflow run SciLifeLab/NGI-RNAseq -profile docker --reads '<path to your reads>' --fasta <path to the genome's fasta file> --gtf <path to the genome's gtf file>
```
The fasta and GTF parameters can be specified in a configuration file (you can provide it with -c), look into the `conf/docker_test.config` for an example.

The docker image containting all the required tools will do downloaded on during the run. It can be found on [DockerHub](https://hub.docker.com/r/scilifelab/ngi-rnaseq/)

A test suite for docker has been implemented, and can be run by moving to the `tests` folder and running `./docker_test.sh`. This will download a small yeast genome and some data, and attempt to run the pipeline through docker on that small dataset.

### Running on other clusters
It is entirely possible to run this pipeline on other clusters, though you will need to set up
your own config file so that the script knows where to find your reference files and how your
cluster works.

Copy the contents of [`conf/uppmax.config`](conf/uppmax.config) to your own config file somewhere
and then reference it with `-c` when running the pipeline.

If you think that there are other people using the pipeline who would benefit from your configuration
(eg. other common cluster setups), please let us know. It should be easy to create a new config file
in `conf` and reference this as a named profile in [`nextflow.config`](nextflow.config). Then these
configuration options can be used by specifying `-profile <name>` when running the pipeline.

### Running on Amazon EC2
There are multiple ways of running this pipeline over Amazon's EC2 service.

The simplest way consists of creating an EC2 instance and running the docker flavour of this pipeline over the machine.

A slightly more complex way is to use our prebuilt AMI to create a ready-to-go virtual machine. The AMI is called `scilifelab/ngi-rnaseq`, id: `ami-f23c6081`. It is available in the Ireland region (`eu-west-1`). This AMI comes with all the tools installed, including docker and nextflow.  The pipeline can then be run by creating an instance using this AMI, logging in and using the following command:

```bash
nextflow run SciLifeLab/NGI-RNAseq -profile base --reads 'path/to/data/sample_*_{1,2}.fastq' --fasta 'path/to/fasta.fz' --gtf 'path/to/genes.gtf'
```

As for Docker, the two last parameters can be entered in a configuration file that would be supplied with `-c`. An example would be `conf/amazon_test.config`

A third approach will be to use nextflow to generate a cluster of machines, and run the pipeline there. This is currently not implemented.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---