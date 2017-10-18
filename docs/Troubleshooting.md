# Troubleshooting 

## Generic problems

### Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample. 
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4.  If you are running Single end data make sure to specify `--singleEnd`

If the pipeline can't find your files then you will get the following error

```
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will. 

The above information is also covered in the [usage README](usage.md#--reads).



### Some of my samples are missing!
If the pipeline ran without crashing and some samples are missing then there are two possible explanations:

1. You made a mistake when declaring the input files - in that case see above. 
2. Your sample didn't reach 5 % read alignment.

Samples that receive less than 5% alignment are skipped for further analysis. You will thus see STAR/Hisat2 output for such samples but nothing more. This limit was set for two reasons, one it doesn't really make any sense to waste CPU hours running on data that bad. Secondly, and most importantly we found that some of the downstream processes would occasionally crash if we allowed it. 

### Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your fastq files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point. 

## UPPMAX specific options

### Bioinfo-tools not loaded
For some unknown reason the bioinfo-tools module is sometimes not loaded properly which means that no other modules are loaded and the pipeline crashes. Fortunately this is easily fixed by simply making sure to already have it loaded when starting the pipeline. 
Example error message:

```
Command wrapper:
  
  Lmod Warning: Did not find: FastQC/0.11.5
  
  Try: "module spider FastQC/0.11.5"

```

To fix this before you start the pipeline run:

```
module load bioinfo-tools
```


## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us. 
Have look at the [issue tracker for our repo](https://github.com/SciLifeLab/NGI-RNAseq/issues). Maybe someone has already had the same problem?

Gitter is a chatt client connected to Github, feel free to come in and chat with us;
[NGI-RNAseq Gitter]((https://gitter.im/SciLifeLab/NGI-RNAseq)) 

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow). 



