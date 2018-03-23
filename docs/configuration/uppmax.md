# NGI-RNAseq: UPPMAX Configuration

The pipeline comes bundled with configurations to use the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (tested on `milou`, `rackham`, `bianca` and `irma`). As such, you shouldn't need to add any custom configuration - everything _should_ work out of the box.

To use the pipeline on UPPMAX, you **must** specificy `-profile uppmax` when running the pipeline (as of Nov 2017).

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```

## Running offline
If you are running the pipeline on Bianca or Irma, you will not have an active internet connection and some automated features will not be able to function. Specifically, you'll need to transfer the pipeline files and the singularity image manually.

First, to generate the singularity image, run the following command. Note that you need singularity installed - this is available on the other UPPMAX clusters (Milou and Rackham):

First, pull the image file where you have an internet connection:

> NB: The "tag" at the end of this command corresponds to the pipeline version.
> Here, we're pulling the docker image for version 1.4 of the NGI-RNAseq pipeline
> Make sure that this tag corresponds to the version of the pipeline that you're using

```bash
singularity pull --name ngi-rnaseq-1.4.img docker://scilifelab/ngi-rnaseq:1.4
pwd # Prints path to your singularity container
```

The NGI-RNAseq pipeline files can be downloaded from https://github.com/SciLifeLab/NGI-RNAseq/releases

Download the pipeline files and transfer the compressed archive (the `.zip`
or `.tar.gz` file). Once transferred, extract the pipeline files.
For example, with a `.zip` file:

```bash
unzip 1.3.1.zip
mv NGI-RNAseq-1.3.1 NGI-RNAseq # rename the folder
cd NGI-RNAseq
pwd # Prints full path to your pipeline
```

Finally, move to the directory where you want to run the pipeline
and execute Nextflow with the path to the pipeline, as so:

```bash
cd /path/to/my/data/analysis
nextflow run /path/to/NGI-RNAseq -with-singularity /path/to/singularity/ngi-rnaseq-1.4.img
```

(Note that you'll need the other common flags such as `--reads` and `--genome` in addition to this command).

## Environment modules and development
If you would prefer to use environment modules instead of singularity, you can use the old version of the configuration by specifying `-profile uppmax_modules` (we don't recommend this).

For pipeline development work on `milou`, use `-profile uppmax_devel` - this uses the milou [devel partition](http://www.uppmax.uu.se/support/user-guides/slurm-user-guide/#tocjump_030509106905141747_8) for testing the pipeline quickly. Please note that this is _not_ suitable for proper analysis runs - only tiny test datasets.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
