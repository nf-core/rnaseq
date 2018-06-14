# nfcore/rnaseq: CCGA Configuration

The IKMB/CCGA Kiel currently maintains one HPC Cluster at Kiel University.

Please follow these steps to run the pipeline:

First-time only: Create the conda environment:

```
module load miniconda2
conda env create -f /path/to/git/environment.yml
```

Each time you wish to run the pipeline, activate the conda environment, specifying the desired version (VERSION, e.g. 1.5):

```
module load miniconda2
source activate nfcore-rnaseq-VERSION
```

Finally, load the Nextflow module:

```
module load IKMB Java/1.8 Nextflow
```

To run the pipeline with the pre-configured environment, use `-profile ccga`. A full command line call will then look as follows:

```
nextflow run nf-core/rnaseq -profile ccga --reads '*_R{1,2}.fastq.gz' --genome GRCh37
```
