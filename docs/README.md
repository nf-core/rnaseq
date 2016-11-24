# NGI-RNAseq

NGI-RNAseq is the new RNA-seq Best Practice pipeline used by the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab](https://www.scilifelab.se/platforms/ngi/) in Stockholm, Sweden.

## Pipeline overview:

* [FastQC](#fastqc) - read quality control
* [TrimGalore](#trimgalore) - adapter trimming
* [STAR](#star) - alignment
* [RSeQC](#rseqc) - RNA quality control metrics
   - [BAM stat](#bam-stat)
   - [Infer experiment](#infer-experiment)
   - [Junction saturation](#junction-saturation)
   - [RPKM saturation](#rpkm-saturation)
   - [Read duplication](#read-duplication)
   - [Inner distance](#inner-distance)
   - [Gene body coverage](#gene-body-coverage)
   - [Read distribution](#read-distribution)
   - [Junction annotation](#junction-annotation)
* [dupRadar](#dupradar) - technical / biological read duplication
* [Preseq](#preseq) - library complexity
* [featureCounts](#featurecounts) - gene counts, biotype counts, rRNA estimation.
* [StringTie](#stringtie) - FPKMs for genes and transcripts
* [Sample_correlation](#Sample_correlation) - create MDS plot and sample pairwise distance heatmap / dendrogram
* [MultiQC](#multiqc) - aggregate report

## Output files provided
Each process of the pipeline generates it's own directory with it's output files contained within. 

* fastqc

Contains two files for each of your samples:
`sample_fastqc.zip` and `sample_fastqc.html` . The html file is the fastqc report which contains quality metrics for your untrimmed raw fastq files. The zip directory also contains the fastqc report alongside it's images. 

* trim_galore

Contains two (four for paired end) files for each of your samples:
`sample_trimmed.fq.gz` wich is a fastq file containing quality trimmed reads and `sample_trimmed.fq.gz_triming_report.txt` which contains **trim\_galore**´s triming report wherein you can find out which parameters that were used. 
	    
* STAR

In this directory you will find 5 files for each of your samples. 
`Sample_Aligned.sortedByCoord.out.bam` this is the aligned BAM file and probably of most interest. 
`Sample_Log.final.out` - contains the finished STAR alignment report, mapping information 
`Sample_Log.out` main log file with a lot of detailed information about the run. This file is most useful
for troubleshooting and debugging.
`Sample_Log.progress.out` is written coninously by STAR and is mostly useful for debugging purposes, if the run finsihed succesfully then all info in here is found in the Log.final.out file displayed in a more compact manor. 
`Sample_SJ.out.tab` - contains filtered *splice junctions* detected in the mapping. 

* RSeQC

This directory contains the following files for each sample:

 `Sample_infer_experiment.txt`             
 `Sample_junction_annotation_log.txt`        
 `Sample_read_distribution.txt`       
 `Sample_read_duplication.DupRate_plot.pdf`  
 `Sample_read_duplication.DupRate_plot.r`    
 `Sample_read_duplication.pos.DupRate.xls`   
 `Sample_read_duplication.seq.DupRate.xls`   
 `Sample_RPKM_saturation.eRPKM.xls`   
 `Sample_RPKM_saturation.rawCount.xls`  
 `Sample_RPKM_saturation.saturation.pdf`     
 `Sample_RPKM_saturation.saturation.r`       
 `Sample_rseqc.geneBodyCoverage.curves.pdf`
 `Sample_rseqc.geneBodyCoverage.r`
 `Sample_rseqc.geneBodyCoverage.txt`
 `Sample_rseqc.inner_distance.txt`
 `Sample_rseqc.inner_distance_freq.txt`
 `Sample_rseqc.inner_distance_plot.r`
 `Sample_rseqc.junction.xls`
 `Sample_rseqc.junction_plot.r`
 `Sample_rseqc.junctionSaturation_plot.pdf`
 `Sample_rseqc.junctionSaturation_plot.r`
 `Sample_rseqc.splice_events.pdf`
 `Sample_rseqc.splice_junction.pdf`
 `Sample_bam_stat.txt`
 
 These are all quality metrics files and contains the basis for the plots in the MultiQC report. 
 In general; the `.r` files are R scripts for generating the figures, the `.txt` are summary files, the `.xls` are data tables and the `.pdf` files are the actual figures. 

* dupRadar

Sample_markDups.bam_dupMatrix.txt            Sample_markDups.bam_duprateExpDens.pdf
Sample_markDups.bam_duprateExpBoxplot.pdf    Sample_markDups.bam_expressionHist.pdf
Sample_markDups.bam_duprateExpDensCurve.txt  Sample_markDups.bam_intercept_slope.txt

* preseq

`sample_ccurve.txt` is the only output file from preseq. This file contains plot values for the complexity curve (the bundled MultiQC report uses this for it's plot). 

* featureCounts 
 In here you will find three files: `Sample.bam_biotype_counts.txt` This file contains the read counts for the different bioitypes that featureCounts distinguishes.   `Sample.featureCounts.txt` contains a list of the counts for each gene id in the provided `gtf` file `Sample.featureCounts.txt.summary` is the summary file wich contains statistics about the mapping. 

* StringTie

These output files are covered in the explanation of the program below: 
[StringTie](#stringtie)

* Sample correlation

In here you will find six files; `edgeR_MDS_distance_matrix.txt` and `edgeR_MDS_plot_coordinates.txt` conatins the data for the MDS plot, which you willl find in the `edgeR_MDS_plot.pdf` file. `log2CPM_sample_distances_dendrogram.pdf` contains a dendrogram of the euclidian distance of your samples.    `log2CPM_sample_distances.txt` contains the data used to generate the heatmap which can be found in `log2CPM_sample_distances_heatmap.pdf`.

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report is for the untrimmed reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

## TrimGalore
The NGI-RNAseq BP 2.0 pipeline uses [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for removal of adapter contamination and trimming of low quality regions. TrimGalore uses [Cutadapt](https://github.com/marcelm/cutadapt) for adapter trimming and runs FastQC after it finishes.

MultiQC reports the percentage of bases removed by TrimGalore in the _General Statistics_ table and produces a line plot showing where reads were trimmed.

## STAR
STAR is a read aligner designed for RNA sequencing.  STAR stands for Spliced Transcripts Alignment to a Reference, it produces results comparable to TopHat (the aligned previously used by NGI for RNA alignments) but is much faster.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

![STAR](images/star_alignment_plot.png)

## RSeQC
RSeQC is a package of scripts designed to evaluate the quality of RNA seq data. You can find out more about the package at the [RSeQC website](http://rseqc.sourceforge.net/).

This pipeline runs several, but not all RSeQC scripts. Results can be found in the `rseqc` directory.

### BAM stat

This script gives numerous statistics about the aligned BAM files produced by STAR. A typical output looks as follows:

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

MultiQC plots each of these statistics in a dot plot. Each sample in the project is a dot - hover to see the sample highlighted across all fields.

RSeQC documentation - [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py)

### Infer experiment

This script predicts the mode of library preparation (sense-stranded or antisense-stranded) according to how aligned reads overlay gene features in the reference genome.
Example output from an unstranded (~50% sense/antisense) library of paired end data:

**From MultiQC report:**
![infer_experiment](images/rseqc_infer_experiment_plot.png)

**From the `infer_experiment.txt` file:**

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0409
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4839
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4752
```
RSeQC documentation - [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)


### Junction saturation

This script shows the number of splice sites detected at the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield more observations. A good sample should approach such a plateau of _Known junctions_, very deep sequencing is typically requires to saturate all _Novel Junctions_ in a sample.

None of the lines in this example have plateaued and thus these samples could reveal more alternative splicing information if they were sequenced deeper.

![Junction saturation](images/rseqc_junction_saturation_plot.png)

RSeQC documentation - [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)

### RPKM saturation

This tool resamples a subset of the total RNA reads and calculates the RPKM value for each subset. We use the default subsets of every 5% of the total reads. I.e 5%,10%...95%,100%.
A percent relative error is then calculated based on the subsamples, this is the y-axis in the graph. A typical PDF figure looks as follows:

![RPKM saturation](images/saturation.png)

A complex library will have low resampling error in well expressed genes. This data is not yet reported in the MultiQC report.

RSeQC documentation - [RPKM_saturation.py](http://rseqc.sourceforge.net/#rpkm-saturation-py)


### Read duplication
This plot shows the number of reads (y-axis) with a given number of exact duplicates (x-axis). Most reads in an RNA-seq library should have a low number of exact duplicates. Samples which have many reads with many duplicates (a large area under the curve) may be suffering excessive technical duplication.

![Read duplication](images/rseqc_read_dups_plot.png)

RSeQC documentation - [read_duplication.py](http://rseqc.sourceforge.net/#read-duplication-py)

### Inner distance
inner_distance.py  tries to calculate the inner distance (or insert size) between two paired RNA reads. This is calculated as the start of read 1 to the end of read 2:
![inner distance concept](images/inner_distance_concept.png)
> _Credit: RSeQC documentation._

Note that values can be negative if the reads overlap. A typical set of samples may look like this:
![Inner distance](images/rseqc_inner_distance_plot.png)

Note that this plot will not be generated for single-end data. Very short insert sizes are often seen in old or degraded samples (_eg._ FFPE).

RSeQC documentation - [inner_distance.py](http://rseqc.sourceforge.net/#inner-distance-py)

### Gene body coverage

This script calculates the reads coverage across gene bodies. This makes it easy to identify 3' or 5' skew in libraries, which can happen in degraded samples prepared with poly-A selection.

A typical set of libraries with little or no bias will look as follows:

![Gene body coverage](images/rseqc_gene_body_coverage_plot.png)

RSeQC documentation - [gene\_body_coverage.py](http://rseqc.sourceforge.net/#genebody-coverage-py)


### Read distribution
This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA seq experiments is generally to have as many exonic reads as possible (_CDS_Exons_). A large amount of intronic reads could be indicative of DNA contamination in your sample or some other problem.

![Read distribution](images/rseqc_read_distribution_plot.png)

RSeQC documentation - [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)


### Junction annotation
Junction annotation compares detected splice junctions to a reference gene model. An RNA read can be spliced 2 or more times, each time is called a splicing event.

![Junction annotation](images/rseqc_junction_annotation_junctions_plot.png)

RSeQC documentation - [junction_annotation.py](http://rseqc.sourceforge.net/#junction-annotation-py)

## dupRadar
[dupRadar](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html) is a Bioconductor library for R. It plots the duplication rate against expression (RPKM) for every gene. A good sample with little technical duplication will only show high numbers of duplicates for highly expressed genes. Samples with technical duplication will have high duplication for all genes, irrespective of transcription level.

![dupRadar](images/dupRadar_plot.png)
> _Credit: [dupRadar documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)_

## Preseq
[Preseq](http://smithlabresearch.org/software/preseq/) estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing the total read count. A shallow curve indicates that the library has reached complexity saturation and further sequencing would likely not add further unique reads. The dashed line shows a perfectly complex library where total reads = unique reads.

Note that these are predictive numbers only, not absolute. The MultiQC plot can sometimes give extreme sequencing depth on the X axis - click and drag from the left side of the plot to zoom in on more realistic numbers.

![preseq](images/preseq_plot.png)

## featureCounts
[featureCounts](http://bioinf.wehi.edu.au/featureCounts/) from the subread package summarises the read distribution over genomic features such as genes, exons, promotors, gene bodies, genomic bins and chromosomal locations.
RNA reads should mostly overlap genes, so be assigned.

![featureCounts](images/featureCounts_assignment_plot.png)

We also use featureCounts to count overlaps with different classes of features. This gives a good idea of where aligned reads are ending up and can show potential problems such as rRNA contamination.
![biotypes](images/featureCounts_biotype_plot.png)
 
## StringTie
[StringTie](https://ccb.jhu.edu/software/stringtie/) assembles RNA-Seq alignments into potential transcripts. It assembles and quantitates full-length transcripts representing multiple splice variants for each gene locus.

StringTie outputs FPKM metrics for genes and transcripts as well as the transcript features that it generates. The analysis generates the following files:

* `<sample>_Aligned.sortedByCoord.out.bam.gene_abund.txt`
  * Gene aboundances, FPKM values
* `<sample>_Aligned.sortedByCoord.out.bam_transcripts.gtf`
  * This `.gtf` file contains all of the assembled transcipts from StringTie
* `<sample>_Aligned.sortedByCoord.out.bam.cov_refs.gtf`
  * This `.gtf` file contains the transcripts that are fully covered by reads.

## Sample_correlation
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) is a Bioconductor package for R used for RNA-seq data analysis. The script included in the pipeline uses edgeR to normalise read counts and create a heatmap / dendrogram showing pairwise euclidean distance (sample similarity). It also creates a 2D MDS scatter plot showing sample grouping. These help to show sample similarity and can reveal batch effects and sample groupings.

**Heatmap:**
![heatmap](images/mqc_hcplot_hocmzpdjsq.png)
**MDS plot:**
![mds_plot](images/mqc_hcplot_ltqchiyxfz.png)


## MultiQC
[MultiQC](http://multiqc.info/) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.
