#Manual RNA BP 2.0
This document describes the different programs and analyses that are part of the NGI-RNAseq best practice 2.0 analysis pipeline. 

The programs used in this pipeline all have their own outputs, here we will mostly describe what the output looks like in the MultiQC report that is generated at the end of the pipeline.

##Pipeline overview:

* FastQC - read quility control
* Cutadapt - trimming
* STAR - alignment
* RSeQC
   - bam stat
   - infer experiment
   - junction saturation
   - RPKM saturation
   - read duplication
   - inner distance
   - gene body coverage
   - read distribution
   - junction annotation
* dupRadar - investigate duplication problems
* preseq - 
* subread featureCounts - gene counts, biotype counts, rRNA estimation.
* String Tie - FPKMs for genes and transcripts
* edgeR - create MDS plot and sample pairwise distance heatmap / dendrogram
* MultiQC - aggregate report

## FastQC
FastQC gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get inormarmation wheter there is adapter contamination in your reads, or if there are other overrepresented sequences present.

For further reading and documentation see the [FastQC website](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

Note that the FastQC plots displayed in the MultiQC report is for the untrimmed reads. They will thus contain adapeter sequence and potenitally regions with low quality. FastQC is ran another time after trimming, that output is included in the results dir under FastQC. 

## Cutadapt
The NGI RNA BP 2.0 pipeline uses TrimGalore for removal of adapter contamination and trimming of low quality regions. TrimGalore is a wrapper for FastQC and Cutadapt. Where FastQC is used for indetifying  adapter contamination and regions with too low quality. Cutadapt is the used to remove these regions of the reads. 
[TrimGalore Website](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[Cutadapt Documentation](http://cutadapt.readthedocs.io/en/stable/guide.html)
## STAR
STAR is a read aligner designed for RNA sequencing.  STAR stands for Spliced Transcripts Alignment to a Reference, and it is designed with speed in mind. 

The STAR section of the MultiQC report is rather self explanatory the more blue (Uniquely mapped or Mapped to multiple loci) the better. Conversely, the more pink or burgundy (Unmapped), the worse the alignment.

<img src= "images/star_alignment_plot.png">

##RSeQC
RSeQC is a package of scripts designed to evaluate the quality of RNA seq data. [RSeQC manual](http://rseqc.sourceforge.net/)

###bam_stat.py

Determines the amount of uniquely mapped reads based on the mapping quality. An example output:

```
#Output (all numbers are read count)
#==================================================
Total records:                                 41465027
QC failed:                                     0
Optical/PCR duplicate:                         0
Non Primary Hits                               8720455
Unmapped reads:                                0

mapq < mapq_cut (non-unique):                  3127757
mapq >= mapq_cut (unique):                     29616815
Read-1:                                        14841738
Read-2:                                        14775077
Reads map to '+':                              14805391
Reads map to '-':                              14811424
Non-splice reads:                              25455360
Splice reads:                                  4161455
Reads mapped in proper pairs:                  21856264
Proper-paired reads map to different chrom:    7648
```
MultiQC plots each of these statistics in a dot plot (each sample is a dot - hover to see the sample highlighted across all fields).
###infer_experiment.py

This script guesses how the reads were stranded for strand specific RNA-seq data. 
Example output from an unstranded(~50% sense/antisense) library of paired end data:

**From MultiQC report:**
<img src="images/rseqc_infer_experiment_plot.png">

**From the infer_experiment.txt file:**

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0409
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4839
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4752
```



###junction_saturation.py
This script aims to determine if the used sequencing depth is enough to perform  alternative splicing analysis on the data. 

<img src="images/rseqc_junction_saturation_plot.png">

A line that has plateaued indicates that all junctions have been detected and that further secuencing will not yield more results. Insteads the same fragment would be sequenced again and again. None of the lines in the example above have plataued and thus these samples could reveal more alternative splicing information if they were sequenced deeper. 

###RPKM_saturation.py
Resamples a subset of the total RNA reads and calculates the RPKM value for each subset. We use the default subsets of every 5% of the total reads. I.e 5%,10%...95%,100%. 
A percent relative error is then calcuated based on the subsamples, this is the y-axis in the graph. 
As of this text being written there isn't a MultiQC submodule for displaying the output of RPKM_saturation.py. Below is instead the pdfs generated by the script itself:

<img src= "images/saturation.png" width="500" height="500" >


###read_duplication.py
This plot shows the number of reads (y-axis) with a given number of exact duplicates (x-axis). Most reads in an RNA-seq library should have a low number of exact duplicates. Samples which have many reads with many duplicates (a large area under the curve) may be suffering excessive technical duplication.
<img src= "images/rseqc_read_dups_plot.png" >

###inner_distance.py
inner_distance.py  tries to calculate the inner distance (or insert size) between two paired RNA reads. Calculated as _start of reads 2_ - _end of read 2_. 
Note that values can be negative if the reads overlap. One of the p
<img src= "images/rseqc_inner_distance_plot.png">

###gene\_body_coverage.py

gene\_body_coverage.py calculates the reads coverage over the gene body. This makes it easy to identify 3 or 5' skewness in the libraries. 
<img src= "images/rseqc_gene_body_coverage_plot.png">


###read_distribution.py
read\_distribution.py calculates how mapped reads are distributed over genome features. A good result for a standard RNA seq experiments is genreally to have as many exonic reads as possible - CDS_Exons should be close to 100%. A large amount of intronic reads could be indicitative of DNA contamination in your sample or some other problem.

<img src="images/rseqc_read_distribution_plot.png">


###junction_annotation.py
Junction annotation compares detected splice junctions to a reference gene model. An RNA read can be spliced 2 or more times, each time is called a splicing event.
<img src= "images/rseqc_junction_annotation_junctions_plot.png">

##dupRadar
dupRadar is a Bioconductor library for R. 
It plots the duplication rate per gene. A good sample (i.e. without PCR duplication?) will have a sigmodial shape as seen in the example output below:  

##preseq
Preseq estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing the total read count. A shallow curve indicates that the libary has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads. 
<img src="images/preseq_plot.png">

##subread featureCounts 
The featureCounts program from the subread package summarises the read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations. 

In addition we also plots the Biotype information, which includes miRNA, rRNA, diffferent kinds of pseudogenes to mention a few. 
In both plot you generally want the blue bar to be as big ass possible. That means assigned feature in the general featureCounts plot and protein coding in the Biotype plot. The biotype plot is a useful tool for detecting rRNA contamination in your samples. 

<img src="images/featureCounts_assignment_plot.png">
<img src="images/featureCounts_biotype_plot.png">
 
##String Tie
StringTie assembles RNA-Seq alignments into potential transcripts. It assembles and quantitates full-length transcripts representing multiple splice variants for each gene locus.

StringTie outputs FPKM metrics for genes and transcripts. The analysis generates the following files:
_\<sample\>_ligned.sortedByCoord.out.bam_transcripts.gtf
This `.gtf` file contains all of the assembled transcipts from StringTie 

_\<sample\>_Aligned.sortedByCoord.out.bam.gene_abund.txt 
Gene aboundances, FPKM values 

_\<sample\>_Aligned.sortedByCoord.out.bam.cov_refs.gtf
This `.gtf` file contains the transcripts that are fully covered by reads.

## edgeR
edgeR is a bioconductor package for R used for RNAseq data analysis. The script included in the pipeline uses edgeR to normalise read counts and create a heatmap / dendrogram showing pairwise euclidean distance (sample similarity). It also creates a 2D MDS scatter plot showing sample grouping.
[edgeR Bioconductor page](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

**Heatmap:**
<img src= "images/mqc_hcplot_hocmzpdjsq.png" >
**MDS plot:**
<img src= "images/mqc_hcplot_ltqchiyxfz.png" >
The heat map and MDS plot together try and give a feeling for which of your sample libraries that are most simmilar to eah other. 

##MultiQC
[MultiQC](http://multiqc.info/) is a visualation tool and is what we at NGI use to generate our hmtl reports. It collects all the reports for your samples and genearate clear and simple reports for all your samples for an easy overview of your project. Almost all of the output results of the RNA seq BP pipeline 2.0 is available in the `multiqc_report.html` file.  

