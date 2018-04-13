# nfcore/RNAseq: CCGA Configuration

The IKMB/CCGA Kiel currently maintains one HPC Cluster at Kiel University. 

# CCGA

To run the pipeline with the pre-configured environment, use `-profile ccga`. A full command line call will then look as follows:
```
nextflow run -profile ccga nf-core/RNAseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 
``` 
---

[![IKMB](images/IKMB_logo.png)](https://http://www.ikmb.uni-kiel.de/)

---
