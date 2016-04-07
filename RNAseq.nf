#!/usr/bin/env nextflow

/*
========================================================================================
                    R N A - S E Q    T W O    P O I N T    Z E R O
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 @Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow rnaseq.nf
 
 Pipeline variables can be configured with the following command line options:
 --genome [GRCh37 | GRCm38]
 --index [path to STAR index]
 --gtf [path to GTF file]
 --files [path to input files]
 --mode [single | paired]
 
 For example:
 $ nextflow rnaseq.nf --files path/to/data/*fq.gz
----------------------------------------------------------------------------------------
 Pipeline overview:
 - cutadapt - trimming
 - STAR - align
 - RSeQC
   - bam_stat
   - infer_experiment
   - splice junction saturation
   - RPKM saturation
   - read duplication
   - inner distance
   - gene body coverage
   - read distribution
   - junction annotation
 - dupRadar
 - preseq
 - subread featureCounts - gene counts. rRNA estimation.
 - String Tie - FPKMs for genes and transcripts
 - MultiQC
----------------------------------------------------------------------------------------
 GA project GA_14_20 RNA-Seq Pipeline. See planning document:
 https://docs.google.com/document/d/1_I4r-yYLl_nA5SzMKtABjDKxQxHSb5N9FMWyomVSWVU/edit#heading=h.uc2543wvne80
----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Reference genome index
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].star
params.gtf   = params.genomes[ params.genome ].gtf
params.bed12 = params.genomes[ params.genome ].bed12

// Input files
params.read1 = "data/*_1.fastq.gz"
params.read2 = "data/*_2.fastq.gz"
read1 = file(params.read1)
read2 = file(params.read2)
// Output path
params.out = "$PWD"

log.info "===================================="
log.info " RNAbp : RNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Read 1       : ${params.read1}"
log.info "Read 2       : ${params.read2}"
log.info "Genome       : ${params.genome}"
log.info "Index        : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.out}"
log.info "===================================="

//Setup R libary directory
params.rlocation = "/home/richam/R/nxtflow_libs/"
nxtflow_libs=file('/home/richam/R/nxtflow_libs')
nxtflow_libs.mkdirs()

// Set up nextflow objects
index = file(params.index)
gtf   = file(params.gtf)
bed12 = file(params.bed12)

// Validate inputs
if( !index.exists() ) exit 1, "Missing STAR index: ${index}"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: ${gtf}"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: ${bed12}"



/*
 * STEP 1 - FastQC
 */

process fastqc {
    
    module 'bioinfo-tools'
    module 'FastQC'
    
    memory '2 GB'
    time '1h'
    
    input:
    file read1 from read1
    file read2 from read2

    output:
    file '*_fastqc.html' into fastqc_html
    file '*_fastqc.zip' into fastqc_zip

    """
    fastqc -q ${read1} ${read2}
    """
}



/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
    
    module 'bioinfo-tools'
    module 'FastQC'
    module 'cutadapt'
    module 'TrimGalore'
    
    cpus 3
    memory '3 GB'
    time '8h'

    input:
    file read1 from read1
    file read2 from read2
    
    output:
    file '*_val_1.fq.gz' into trimmed_read1
    file '*_val_2.fq.gz' into trimmed_read2

    """
    trim_galore --paired --gzip --fastqc_args "-q" $read1 $read2
    """
}



/*
 * STEP 3 - align with STAR
 * Inspired by https://github.com/AveraSD/nextflow-rnastar
 */

process star {
    
    module 'bioinfo-tools'
    module 'star'
    
    cpus 8
    memory '64 GB'
    time '5h'

    input:
    file index
    file gtf
    file trimmed_read1 from trimmed_read1
    file trimmed_read2 from trimmed_read2

    output:
    file '*.Aligned.sortedByCoord.out.bam' into bam4, bam5, bam6, bam7, bam8, bam9
    file '*.Log.final.out' into results
    file '*.Log.out' into results
    file '*.Log.progress.out' into results
    file '*.SJ.out.tab' into results

    """
    prefix=\$(echo $trimmed_read1 | sed 's/_.*/./')
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $trimmed_read1 $trimmed_read2 \\
         --runThreadN ${task.cpus} \\
         --twopassMode Basic \\
         --outWigType bedGraph \\
         --outSAMtype BAM SortedByCoordinate\\
         --readFilesCommand zcat\\
         --outFileNamePrefix \$prefix
    """
}



/*
 * STEP 4 - RNASeQC analysis
 */

process rnaseqc {
    
    module 'bioinfo-tools'
    module 'rseqc'
    
    memory '64 GB'
    time '2h'
   
    errorStrategy 'ignore' 
    
    input:
    file bam4
    file bed12 from bed12
   
     
    output:
    file '*.bam_stat.txt' into results                          // bam_stat
    file '*.splice_events.{txt,pdf}' into results               // junction_annotation
    file '*.splice_junction.{txt,pdf}' into results             // junction_annotation
    file '*.junctionSaturation_plot.{txt,pdf}' into results     // junction_saturation
    file '*.inner_distance.{txt,pdf}' into results              // inner_distance
    file '*.curves.{txt,pdf}' into results                      // geneBody_coverage
//    file '*.heatMap.{txt,pdf}' into results                     // geneBody_coverage
    file '*.infer_experiment.txt' into results                  // infer_experiment
    file '*.read_distribution.txt' into results                 // read_distribution
    file '*DupRate.xls' into results                            // read_duplication
    file '*DupRate_plot.pdf' into results                      // read_duplication
    file '*.saturation.{txt,pdf}' into results             // RPKM_saturation
    
    """
    bam_stat.py -i $bam4 2> ${bam4}.bam_stat.txt
    junction_annotation.py -i $bam4 -o ${bam4}.rseqc -r $bed12
    junction_saturation.py -i $bam4 -o ${bam4}.rseqc -r $bed12
    inner_distance.py -i $bam4 -o ${bam4}.rseqc -r $bed12
    geneBody_coverage.py -i $bam4 -o ${bam4}.rseqc -r $bed12
    infer_experiment.py -i $bam4 -r $bed12 > ${bam4}.infer_experiment.txt
    read_distribution.py -i $bam4 -r $bed12 > ${bam4}.read_distribution.txt
    read_duplication.py -i $bam4 -o ${bam4}.read_duplication
    RPKM_saturation.py -i $bam4 -r $bed12 -d '1+-,1-+,2++,2--' -o ${bam4}.RPKM_saturation
    """
}






/*
 * STEP 5 - preseq analysis
 */

process preseq {
    
    module 'bioinfo-tools'
    module 'preseq'
    
    memory '4 GB'
    time '2h'
    
    input:
    file bam5
    
    output:
    file '*.ccurve.txt' into results
    
    """
    preseq lc_extrap -v -B $bam5 -o ${bam5}.ccurve.txt
    """
}


/*
* STEP 6 Mark duplicates
*/

process markDuplicates {
    module 'bioinfo-tools'
    module 'picard/2.0.1'
    
    memory '16GB'
    time '2h'
    
    input:
    file bam6
    
    output: 
    file '*.markDups.bam' into bam_md
    file '*markDups_metrics.txt' into results

    """
    echo \$PICARD_HOME
    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates INPUT=${bam6} OUTPUT=${bam6}.markDups.bam METRICS_FILE=${bam6}.markDups_metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true PROGRAM_RECORD_ID='null' VALIDATION_STRINGENCY=LENIENT 
    """
    }



/*
STEP 7 - dupRadar
 */

process dupradar {
    
    module 'bioinfo-tools'
    module 'R/3.2.3'
    
    memory '16 GB'
    time '2h'
    
    input:
    file bam_md 
    file gtf from gtf
    
    output:
    file '*_duprateExpDens.pdf' into results
    file '*_intercept_slope.txt' into results
    file '*_expressionHist.pdf' into results
    
    shell
    """
    #!/usr/bin/env Rscript
    if (!("dupRadar" %in% installed.packages()[,"Package"])){
        .libPaths( c( "/home/richam/R/nxtflow_libs", .libPaths() ) )
        source("https://bioconductor.org/biocLite.R")
        biocLite("dupRadar")
    }
    library("dupRadar")
          
    # Duplicate stats
    stranded <- 2
    paired <- TRUE
    threads <- 8
    dm <- analyzeDuprates("${bam_md}", "${gtf}", stranded, paired, threads)
    write.table(dm, file=paste("${bam_md}", "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\t")
    
    # 2D density scatter plot
    pdf(paste0("${bam_md}", "_duprateExpDens.pdf"))
    duprateExpDensPlot(DupMat=dm)
    title("Density scatter plot")
    dev.off()
    fit <- duprateExpFit(DupMat=dm)
    cat("duprate at low read counts: ", fit\$intercept, "progression of the duplication rate: ", fit\$slope, "\n", fill=TRUE, labels="${bam_md}", file=paste0("${bam_md}", "_intercept_slope.txt"), append=FALSE)
   
    # Distribution of RPK values per gene
    pdf(paste0("${bam_md}", "_expressionHist.pdf"))
                         expressionHist(DupMat=dm)
    title("Distribution of RPK values per gene")
    dev.off()
    """
    }

 /*
 * STEP 8 Feature counts
 */


process featureCounts {
    
    module 'bioinfo-tools'
    module 'subread'
    
    memory '4 GB'
    time '2h'
    
    input:
    file bam8
    file gtf from gtf
    
    output:
    file '*_gene.featureCounts.txt' into results
    file '*_biotype.featureCounts.txt' into results
    file '*_rRNA_counts.txt' into results
    
    """
    featureCounts -a $gtf -g gene_id -o ${bam8}_gene.featureCounts.txt -p -s 2 $bam8
    featureCounts -a $gtf -g gene_biotype -o ${bam8}_biotype.featureCounts.txt -p -s 2 $bam8
    cut -f 1,7 ${bam8}_biotype.featureCounts.txt | sed '1,2d' | grep 'rRNA' > ${bam8}_rRNA_counts.txt
    """
}



/*
 * STEP 9 - stringtie FPKM
 */

process stringtieFPKM {
    
    module 'bioinfo-tools'
    module 'StringTie'
    
    memory '4 GB'
    time '2h'
    
    input:
    file bam9
    file gtf from gtf
    
    output:
    file '*_transcripts.gtf' into results
    file '*.gene_abund.txt' into results
    file '*.cov_refs.gtf' into results
    
    """
    stringtie $bam9 -o ${bam9}_transcripts.gtf -v -G $gtf -A ${bam9}.gene_abund.txt -C ${bam9}.cov_refs.gtf -e -b ${bam9}_ballgown
    """
}



