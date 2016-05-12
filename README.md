# RNA-BP
Pipeline for RNA sequencing best practice abalysis at the NGI at Scilifelab Stockholm, Sweden

Authors:
Phil Ewels (@ewels)

Rickard Hammarén (@Hammarn)



## Pipeline parameters

#### `--reads` 
Points to the read files 

```
 --reads '/home/dataset/sample_*_{1,2}.fastq'
```
By default RNA-BP assumes the reads to be in the folder data in the working directory i.e.
`data/*.fastq.gz`

### `--genome`
The reference genome to use of the anaylys, needs to be one of the genome specified in the config file. 
The human `GRCh37` genome is set as default. 
```
--genome 'GRCm38'
```

### '-c'
Path to config file - nextflow command


An example run:
```
$ nextflow main.nf --reads '*R{1,2}.fastq.gz' --genome 'GRCm38'
```
