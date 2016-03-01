#!/usr/bin/env nextflow

/* ========================================================================================
 *                      R N A - S E Q    T W O    P O I N T    Z E R O
 * ========================================================================================
 * New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 * @Authors
 * Phil Ewels <phil.ewels@scilifelab.se>
 * ----------------------------------------------------------------------------------------
 * Basic command:
 * $ nextflow rnaseq.nf
 *
 * Pipeline variables can be configured with the following command line options:
 * --genome [GRCh37 | GRCm38]
 * --index [path to STAR index]
 * --gtf [path to GTF file]
 * --files [path to input files]
 * --mode [single | paired]
 *
 * For example:
 * $ nextflow rnaseq.nf --files path/to/data/*fq.gz
 * ----------------------------------------------------------------------------------------
 * Pipeline overview:
 * - cutadapt - trimming
 * - STAR - align
 * - RSeQC
 *   - bam_stat
 *   - infer_experiment
 *   - splice junction saturation
 *   - RPKM saturation
 *   - read duplication
 *   - inner distance
 *   - gene body coverage
 *   - read distribution
 *   - junction annotation
 * - dupRadar
 * - preseq
 * - subread featureCounts - gene counts. rRNA estimation.
 * - String Tie - FPKMs for genes and transcripts
 * ----------------------------------------------------------------------------------------
 * GA project GA_14_20 RNA-Seq Pipeline. See planning document:
 * https://docs.google.com/document/d/1_I4r-yYLl_nA5SzMKtABjDKxQxHSb5N9FMWyomVSWVU/edit#heading=h.uc2543wvne80
 * ----------------------------------------------------------------------------------------
 */



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Reference genome index
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].star
params.gtf = params.genomes[ params.genome ].gtf

// Input files
params.read1 = file("data/*_1.fastq.gz")
params.read2 = file("data/*_2.fastq.gz")

// Output path
params.out = "$PWD"

log.info "===================================="
log.info " RNAbp : RNA-Seq Best Practice v${version}"
log.info "===================================="
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

// Set up nextflow objects
gtf = file(params.gtf)
index = file(params.index)


/*
 * STEP 1 - FastQC
 */

process fastqc {
    
    module 'bioinfo-tools'
    module 'FastQC'
    
    memory '2 GB'
    time '1h'
    
    input:
    file read1 from params.read1
    file read2 from params.read2

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
    memory '3GB'
    time '8h'

    input:
    file read1 from params.read1
    file read2 from params.read2
    
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
    
    cpus 6
    memory '32GB'
    time '5h'

    input:
    file index
    file gtf
    file trimmed_read1
    file trimmed_read2

    output:
    file '*.Aligned.sortedByCoord.out.bam' into bam
    file '*.Log.final.out' into results
    file '*.Log.out' into results
    file '*.Log.progress.out' into results
    file '*.SJ.out.tab' into results

    """
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $trimmed_read1 $trimmed_read2 \\
         --runThreadN ${task.cpus} \\
		 --twopassMode Basic \\
         --outWigType bedGraph \\
		 --outSAMtype BAM SortedByCoordinate
    """
}

