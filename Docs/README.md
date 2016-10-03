#Manual RNA BP 2.0
This is an an attempt at descripting the different programs and analyses that consitutes the NGI RNA BP 2.0 analysis pipeline and to inform the use of what to look for when analysing their samples. 

The programs used in this pipeline all have their own outputs, here we will mostly describe what the output looks like in the MultiQC report that is genrated at the end of the pipeline.

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

For further reading and documentation see the [FastQC website](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or the 
[FastQC Manual](https://biof-edu.colorado.edu/videos/dowell-short-read-class/day-4/fastqc-manual)

## Cutadapt
The NGI RNA BP 2.0 pipeline uses TrimGalore for removal of adapter contamination and trimming of low quality regions. TrimGalore is a wrapper for FastQC and Cutadapt. Where FastQC is used for indetifying  adapter contamination and regions with too low quality. Cutadapt is the used to remove these regions of the reads. 
[TrimGalore Website](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[Cutadapt Documentation](http://cutadapt.readthedocs.io/en/stable/guide.html)
## STAR
STAR is a read aligner designed for RNA sequencing.  STAR stands for Spliced Transcripts Alignment to a Reference, and it is designed with speed in mind. 

The STAR section of the MultiQC report is rather self explanatory the more blue (Uniquely mapped or Mapped to multiple loci) the better conversely, the more pink or burgundy (Unmapped), the worse the alignment.

<img src= "images/star_alignment_plot.png">

##RSeQC
RSeQC is a package of scripts designed to evaluate the quality of RNA seq data. [RSeQC manual](http://rseqc.sourceforge.net/)

###bam_stat.py

Determines the amount of uniquely mapped reads based on the mapping quality. An example output:

```
bam_stat.py  -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam

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
The MultiQC plot of this information can be a bit overwhelming at first but it displayes the above information simultaneously for all the samples in your report. With each sample represented as a dot. 
###infer_experiment.py

This scripts guesses how the reads were stranded for strand specific RNA-seq data. 
Example output from an unstranded(~50% sense/antisense) library of paired end data:

**From MultiQC report:**
<img src= "images/infer_experiment.png">

**From the infer_experiment.txt file:**

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0409
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4839
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4752
```



###junction_saturation.py
This scripts aims to determine if the used sequencing depth is enough to perform  alternative splicing analysis on the data. 

![Junction Saturation](images/Junction_saturation.png)


A line that has plateaued indicates that all junctions have been detected and that further secuencing will not yield more results. Insteads the same fragment would be sequenced again and again. None of the lines in the example above have plataued and thus these samples could reveal more alternative splicing information if they were sequenced deeper. 

###RPKM_saturation.py
Resamples a subset of the total RNA reads and calculates the RPKM value for each subset. We use the default subsets of every 5% of the total reads. I.e 5%,10%...95%,100%. 
A percent relative error is then calcuated based on the subsamples, this is the y-axis in the graph. 
As of this text being written there isn't a MultiQC submodule for displaying the output of RPKM_saturation.py. Below is instead the pdfs generated by the script itself:

<img src= "images/saturation.png" width="500" height="500" >


###read_duplication.py
This script investigates the read duplication rate. A good library will generally have a slightly declining graph. The red library in the example graph really stands out and needs investigating. 

<img src= "images/rseqc_read_dups_plot.png" >

###inner_distance.py
inner_distance.py  tries to calculate the inner distance (or insert size) between two paired RNA reads. Calculated as _start of reads 2_ - _end of read 2_. 
Note that values can be negative if the reads overlap, this is indicative of a poor library, and fragmented RNA. As in the example picture below. 

<img src= "images/rseqc_inner_distance_plot.png">

###gene\_body_coverage.py


###read_distribution.py
read_distribution.py calculates how mapped reads are distributed over genome features. A good result for a standard RNA seq experiments is genreally to have as much exonic reads as possible, i.e CDS_Exons near 100%. A large amount of intronic reads could be indicitative of DNA contamination in your sample or some other problem.

<img src="images/rseqc_read_distribution_plot.png">


###junction_annotation.py
MultiQC module under way. Rscript included. Source it to generate a pie chart over the known and novel splicing junctions. 
##dupRadar
dupRadar is A Bioconductor library for R. 
It plots the duplication rate per gene. A good sample (i.e. without PCR duplication?) will have a sigmodial shape as seen in the example output below:  

##preseq
Preseq estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing the total read count. A shallow curve indicates that the libary has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads. 
<img src="images/preseq_complexity_curve.png">

##subread featureCounts 
The featureCounts program from the subread package summarises the read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations. 

In addition we also plots the Biotype information, which includes miRNA, rRNA, diffferent kinds of pseudogenes to mention a few. 
In both plot you generally want the blue bar to be as big ass possible. That means assigned feature in the general featureCounts plot and protein coding in the Biotype plot. The biotype plot is a useful tool for detecting rRNA contamination in your samples. 

<img src="images/featureCounts_assignment_plot.png">
<img src="images/featureCounts_biotype_plot.png">
 
##String Tie
StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus.
The output from stringtie is FPKMs for genes and transcripts. The analysis generates the following files: 
_\<sample\>_ligned.sortedByCoord.out.bam_transcripts.gtf
_\<sample\>_Aligned.sortedByCoord.out.bam_transcripts.gtf
_\<sample\>_Aligned.sortedByCoord.out.bam.cov_refs.gtf

## edgeR
edgR is a bioconductor package for R and is used to do differential extression analysis on RNAseq data.
Rcreate MDS plot and sample pairwise distance heatmap / dendrogram
[edgeR Bioconductor page](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
The output generated from edgeR process is a heatmap and a MDS-plot.
The heatmap may look something like this;

**Heatmap:**
<img src= "images/mqc_hcplot_hocmzpdjsq.png" >
**MDS plot:**
<img src= "images/mqc_hcplot_ltqchiyxfz.png" >
The hatmap and MDS plot together try and give you a feeling for which of your sample libraries that are most simmilar to eah other. 

##MultiQC
[MultiQC](http://multiqc.info/) is a visualation tool and is what we at NGI use to generate our hmtl reports. It collects all the reports for your samples and genearate clear and simple reports for all your samples for an easy overview of your project. Almost all of the output results of the RNA seq BP pipeline 2.0 is available in the `multiqc_report.html` file.  

