# Troubleshooting 

## Generic problems

### Problems with input specification
This is already covered in the [usage README](usage.md#--reads) but it's worth mentioning again.

If only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration.

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample. 
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.

Note that if your sample name is "messy" then you have to be very perticular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will. 

### Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your fastq files are scattered in different paths then we recomend that you generate a direcotry with symlinked files. If running in paried end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point. 

## UPPMAX specific options

### Bioinfo-tools not loaded
For some unkown reason the bioinfo-tools module is sometimes not loaded properly which means that no other modules are loaded and the pipeline chrashes. Fortunately this is easily fixed by simply making shure to already have it loaded when starting the pipeline. 

I.e. before you start the pipeline run:

```
module load bioinfo-tools
```
Example error message:

```
Command wrapper:
  
  Lmod Warning: Did not find: FastQC/0.11.5
  
  Try: "module spider FastQC/0.11.5"

```


## Extra resources and getting help
If you still have an issue with runing the pipeline then feel free to contact us. 
Have look at the [issue tracker for our repo](https://github.com/SciLifeLab/NGI-RNAseq/issues). Maybe someone has already had the same problem?

Gitter is a chatt client connected to Github, feel free to come in and chat with us;
[NGI-RNAseq Gitter]((https://gitter.im/SciLifeLab/NGI-RNAseq)) 

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow). 



