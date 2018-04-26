# nfcore/rnaseq: QBiC Configuration


QBiC has two clusters available to run and analyze data directly. For both systems, an appropriate configuration profile is existing, enabling a direct usability of the pipeline on the respective cluster environment. 

# BINAC

You may use the pipeline with the `-profile binac` switch when starting the pipeline. A typical call could work like this for example
```
nextflow run -profile binac nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker
``` 
# Core Facility Cluster (CFC)

You may use the pipeline with the `-profile cfc` switch when starting the pipeline. A typical call could work like this for example
```
nextflow run -profile cfc nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile docker
``` 

---

[![QBiC](images/QBiC_logo.png)](https://portal.qbic.uni-tuebingen.de/portal/)

---
