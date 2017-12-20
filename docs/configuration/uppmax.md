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

```bash
singularity pull --name ngi-rnaseq.img docker://scilifelab/ngi-rnaseq
```

The NGI-RNAseq pipeline files can be downloaded from https://github.com/SciLifeLab/NGI-RNAseq/releases

Once transferred, run the pipeline with the specific paths, as so:

```bash
nextflow run /path/to/NGI-RNAseq -with-singularity /path/to/ngi-rnaseq.img
```

(Note that you'll need the other common flags such as `--reads` and `--genome` in addition to this).

## Environment modules and development
If you would prefer to use environment modules instead of singularity, you can use the old version of the configuration by specifying `-profile uppmax_modules` (we don't recommend this).

For pipeline development work on `milou`, use `-profile uppmax_devel` - this uses the milou [devel partition](http://www.uppmax.uu.se/support/user-guides/slurm-user-guide/#tocjump_030509106905141747_8) for testing the pipeline quickly. Please note that this is _not_ suitable for proper analysis runs - only tiny test datasets.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
