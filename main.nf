#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                R N A - S E Q    T W O    P O I N T    Z E R O
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-RNAseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.2

// Configurable variables
params.genome = 'GRCh37'
if( params.genomes ) {
    params.star_index = params.genomes[ params.genome ].star
    params.gtf   = params.genomes[ params.genome ].gtf
    params.bed12 = params.genomes[ params.genome ].bed12
}
params.reads = "data/*{_1,_2}*.fastq.gz"
params.outdir = './results'

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs = file(params.rlocation)
nxtflow_libs.mkdirs()

def single
params.sampleLevel = false
params.strandRule = false

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Validate inputs
if( !params.star_index && !params.fasta && !params.download_fasta ){
    exit 1, "No reference genome specified!"
}
if( params.star_index ){
    star_index = file(params.star_index)
    if( !star_index.exists() ) exit 1, "STAR index not found: $star_index"
} else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: $star_index"
}
if( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "Missing GTF annotation: $gtf"
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// Header log info
log.info "========================================="
log.info " NGI-RNAseq : RNA-Seq Best Practice v${version}"
log.info "========================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
if(params.star_index)          log.info "STAR Index     : ${params.star_index}"
else if(params.fasta)          log.info "Fasta Ref      : ${params.fasta}"
else if(params.download_fasta) log.info "Fasta URL      : ${params.download_fasta}"
if(params.gtf)                 log.info "GTF Ref        : ${params.gtf}"
else if(params.download_gtf)   log.info "GTF URL        : ${params.download_gtf}"
if(params.bed12)               log.info "BED12 Ref      : ${params.bed12}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "R libraries    : ${params.rlocation}"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
if( params.clip_r1 > 0) log.info "Trim R1        : ${params.clip_r1}"
if( params.clip_r2 > 0) log.info "Trim R2        : ${params.clip_r2}"
if( params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if( params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "Config Profile : ${workflow.profile}"
log.info "========================================="

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }

/*
 * PREPROCESSING - Download Fasta
 */
if( !params.star_index && !params.fasta && params.downloadFasta ) {
    process downloadFASTA {
        publishDir { params.saveReference ? "${params.outdir}/reference_genome" : null }, mode: 'copy'

        input:
        set url from params.downloadFasta

        output:
        file "*.{fa,fasta}" into fasta
        
        script:
        """
        curl -O -L $url
        """
    }
}
/*
 * PREPROCESSING - Download GTF
 */
if( !params.gtf && params.downloadGTF ) {
    process downloadGTF {
        publishDir { params.saveReference ? "${params.outdir}/reference_genome" : null }, mode: 'copy'

        input:
        set url from params.downloadGTF

        output:
        file "*.gtf" into gtf
        
        script:
        """
        curl -O -L $url
        """
    }
}
/*
 * PREPROCESSING - Build STAR index
 */
if( fasta && !params.star_index) {
    process makeSTARindex {
        publishDir { params.saveReference ? "${params.outdir}/reference_genome" : null }, mode: 'copy'

        cpus { $params.makeSTARindex.cpus ?: 12 }
        memory { $params.makeSTARindex.memory ?: 30.GB }
        time { $params.makeSTARindex.time ?: 5.h }
        errorStrategy = 'terminate'

        input:
        file fasta

        output:
        file 'star' into star_index
        
        script:
        """
        STAR \
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}
/*
 * PREPROCESSING - Build BED12 file
 */
if ( !params.bed12 ){
    bed12 = file(params.bed12)
} else {
    process makeBED12 {
        publishDir { params.saveReference ? "${params.outdir}/reference_genome" : null }, mode: 'copy'

        time { $params.makeBED12.time ?: 2.h }
        errorStrategy = 'terminate'

        input:
        file gtf

        output:
        file '*.bed12' into bed12
        
        script:
        """
        convert2bed \\
            --input=gtf \\
            --output=bed \\
            --max-mem=${task.mem} \\
            > ${gtf}.bed12
        """
    }
}
if( !bed12.exists() ) exit 1, "Missing BED12 annotation: $bed12"


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    memory { $params.fastqc.memory ?: 2.GB * task.attempt }
    time { $params.fastqc.time ?: 4.h * task.attempt }
    errorStrategy = task.exitStatus == 143 ? 'retry' : 'ignore'

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    cpus { $params.trim_galore.cpus ?: 2 * task.attempt }
    memory { $params.trim_galore.memory ?: 4.GB * task.attempt }
    time { $params.trim_galore.time ?: 8.h * task.attempt }
    errorStrategy = task.exitStatus == 143 ? 'retry' : 'terminate'

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results

    script:
    single = reads instanceof Path
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (single) {
        """
        trim_galore --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}


/*
 * STEP 3 - align with STAR
 * Originally inspired by https://github.com/AveraSD/nextflow-rnastar
 */
process star {
    tag "$reads"
    publishDir "${params.outdir}/STAR", mode: 'copy'

    cpus { $params.star.cpus ?: 10 * task.attempt }
    memory { $params.star.memory ?: 80.GB * task.attempt }
    time { $params.star.time ?: 5.h * task.attempt }
    errorStrategy = task.exitStatus == 143 ? 'retry' : 'terminate'

    input:
    file star_index
    file gtf
    file reads from trimmed_reads

    output:
    set file('*Log.final.out'), file ('*.bam') into aligned
    file '*.out' into star_logs
    file '*SJ.out.tab'

    script:
    """
    #Getting STAR prefix
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
    STAR --genomeDir $star_index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outFileNamePrefix \$f
    """
}


// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    if(percent_aligned.toFloat() <='10'.toFloat() ){
        println "#################### VERY POOR ALIGNMENT RATE ONLY ${percent_aligned}%! FOR ${logs}"
        false
    } else {
        println "Passed aligment with ${percent_aligned}%! FOR ${logs}"
        true
    }
}
// Filter removes all 'aligned' channels that fail the check
aligned
    .filter { logs, bams -> check_log(logs) }
    .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }


/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    tag "$bam_rseqc"
    publishDir "${params.outdir}/rseqc" , mode: 'copy'

    memory { $params.rseqc.memory ?: 32.GB * task.attempt }
    time { $params.rseqc.time ?: 7.h * task.attempt }

    input:
    file bam_rseqc
    file bed12 from bed12

    output:
    file '*.{txt,pdf,r,xls}' into rseqc_results
    /*  The following files are being generated by this process:
        .bam_stat.txt                         // bam_stat
        .splice_events.{txt,pdf}              // junction_annotation
        .splice_junction.{txt,pdf}            // junction_annotation
        .junctionSaturation_plot.{txt,pdf,r}  // junction_saturation
        .inner_distance.{txt,pdf}             // inner_distance
        .curves.{txt,pdf}                     // geneBody_coverage
        .geneBodyCoverage.txt                 // geneBody_coverage
        .heatMap.{txt,pdf}                    // geneBody_coverage
        .infer_experiment.txt                 // infer_experiment
        .read_distribution.txt                // read_distribution
        DupRate.xls                           // read_duplication
        DupRate_plot.pdf                      // read_duplication
        .saturation.{txt,pdf}                 // RPKM_saturation
    */

    script:
    def strandRule = params.strandRule ?: (single ? '++,--' : '1+-,1-+,2++,2--')

    """
    samtools index $bam_rseqc
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc}.infer_experiment.txt
    RPKM_saturation.py -i $bam_rseqc -r $bed12 -d $strandRule -o ${bam_rseqc}.RPKM_saturation
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12 2> ${bam_rseqc}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12
    geneBody_coverage.py -i ${bam_rseqc} -o ${bam_rseqc}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc}.read_duplication
    echo "Filename $bam_rseqc RseQC version: "\$(read_duplication.py --version)
    """
}



/*
 * STEP 5 - preseq analysis
 */
process preseq {
    tag "$bam_preseq"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    memory { $params.preseq.memory ?: 4.GB * task.attempt }
    time { $params.preseq.time ?: 2.h * task.attempt }

    input:
    file bam_preseq

    output:
    file '*.ccurve.txt' into preseq_results

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq}.ccurve.txt
    echo "File name: $bam_preseq  preseq version: "\$(preseq)
    """
}


/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {
    tag "$bam_markduplicates"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy'

    memory { $params.markDuplicates.memory ?: 16.GB * task.attempt }
    time { $params.markDuplicates.time ?: 2.h * task.attempt }

    input:
    file bam_markduplicates

    output:
    file '*.markDups.bam' into bam_md
    file '*markDups_metrics.txt' into picard_results

    script:
    """
    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam_markduplicates \\
        OUTPUT=${bam_markduplicates}.markDups.bam \\
        METRICS_FILE=${bam_markduplicates}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT

    #Printing out version number to standard out
    echo "File name: $bam_markduplicates Picard version "\$(java -Xmx2g -jar \$PICARD_HOME/picard.jar  MarkDuplicates --version 2>&1)
    """
}


/*
 * STEP 7 - dupRadar
 */
process dupradar {
    tag "$bam_md"
    publishDir "${params.outdir}/dupradar", pattern: '*.{pdf,txt}', mode: 'copy'

    memory { $params.dupradar.memory ?: 16.GB * task.attempt }
    time { $params.dupradar.time ?: 2.h * task.attempt }

    input:
    file bam_md
    file gtf from gtf

    output:
    file '*.{pdf,txt}' into dupradar_results

    script:
    def paired = single ? 'FALSE' :  'TRUE'

    """
    #!/usr/bin/env Rscript

    # Load / install dupRadar package
    .libPaths( c( "${params.rlocation}", .libPaths() ) )
    if (!require("dupRadar")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("dupRadar", suppressUpdates=TRUE, lib="${params.rlocation}")
        library("dupRadar")
    }

    # Duplicate stats
    stranded <- 2
    threads <- 8
    dm <- analyzeDuprates("$bam_md", "$gtf", stranded, $paired, threads)
    write.table(dm, file=paste("$bam_md", "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\\t")

    # 2D density scatter plot
    pdf(paste0("$bam_md", "_duprateExpDens.pdf"))
    duprateExpDensPlot(DupMat=dm)
    title("Density scatter plot")
    mtext("$bam_md", side=3)
    dev.off()
    fit <- duprateExpFit(DupMat=dm)
    cat(
      paste("- dupRadar Int (duprate at low read counts):", fit\$intercept),
      paste("- dupRadar Sl (progression of the duplication rate):", fit\$slope),
      fill=TRUE, labels="$bam_md",
      file=paste0("$bam_md", "_intercept_slope.txt"), append=FALSE
    )
    
    # Get numbers from dupRadar GLM
    curve_x <- sort(log10(dm\$RPK))
    curve_y = 100*predict(fit\$glm,data.frame(x=curve_x),type="response")
    # Remove all of the infinite values
    infs = which(curve_x %in% c(-Inf,Inf))
    curve_x = curve_x[-infs]
    curve_y = curve_y[-infs]
    # Reduce number of data points
    curve_x <- curve_x[seq(1, length(curve_x), 10)]
    curve_y <- curve_y[seq(1, length(curve_y), 10)]
    # Convert x values back to real counts
    curve_x = 10^curve_x
    # Write to file
    write.table(
      cbind(curve_x, curve_y),
      file=paste0("$bam_md", "_duprateExpDensCurve.txt"),
      quote=FALSE, row.names=FALSE
    )
    
    # Distribution of expression box plot
    pdf(paste0("$bam_md", "_duprateExpBoxplot.pdf"))
    duprateExpBoxplot(DupMat=dm)
    title("Percent Duplication by Expression")
    mtext("$bam_md", side=3)
    dev.off()
    
    # Distribution of RPK values per gene
    pdf(paste0("$bam_md", "_expressionHist.pdf"))
    expressionHist(DupMat=dm)
    title("Distribution of RPK values per gene")
    mtext("$bam_md", side=3)
    dev.off()

    # Printing sessioninfo to standard out
    print("$bam_md")
    citation("dupRadar")
    sessionInfo()
    """
}


/*
 * STEP 8 Feature counts
 */
process featureCounts {
    tag "$bam_featurecounts"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    memory { $params.dupradar.memory ?: 4.GB * task.attempt }
    time { $params.dupradar.time ?: 2.h * task.attempt }

    input:
    file bam_featurecounts
    file gtf from gtf

    output:
    file '*_gene.featureCounts.txt' into geneCounts
    file '*_gene.featureCounts.txt.summary' into featureCounts_logs
    file '*_biotype_counts.txt' into featureCounts_biotype

    script:
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts}_gene.featureCounts.txt -p -s 2 $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts}_biotype.featureCounts.txt -p -s 2 $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts}_biotype.featureCounts.txt > ${bam_featurecounts}_biotype_counts.txt
    """
}


/*
 * STEP 9 - stringtie FPKM
 */
process stringtieFPKM {
    tag "$bam_stringtieFPKM"
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy'

    memory { $params.dupradar.memory ?: 4.GB * task.attempt }
    time { $params.dupradar.time ?: 2.h * task.attempt }

    input:
    file bam_stringtieFPKM
    file gtf from gtf

    output:
    file '*_transcripts.gtf'
    file '*.gene_abund.txt'
    file '*.cov_refs.gtf'
    stdout into stringtie_log

    script:
    """
    stringtie $bam_stringtieFPKM \\
        -o ${bam_stringtieFPKM}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM}_ballgown

    echo "File name: $bam_stringtieFPKM Stringtie version "\$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }



/*
 * STEP 10 - edgeR MDS and heatmap
 */
process sample_correlation {
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    memory { $params.dupradar.memory ?: 16.GB * task.attempt }
    time { $params.dupradar.time ?: 2.h * task.attempt }

    input:
    file input_files from geneCounts.toList()
    bam_count

    output:
    file '*.{txt,pdf}' into sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script:
    """
    #!/usr/bin/env Rscript

    # Load / install required packages
    .libPaths( c( "${params.rlocation}", .libPaths() ) )
    if (!require("limma")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("limma", suppressUpdates=TRUE, lib="${params.rlocation}")
        library("limma")
    }

    if (!require("edgeR")){
        source("http://bioconductor.org/biocLite.R")
        biocLite("edgeR", suppressUpdates=TRUE, lib="${params.rlocation}")
        library("edgeR")
    }

    if (!require("data.table")){
        install.packages("data.table", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("data.table")
    }

    if (!require("gplots")) {
        install.packages("gplots", dependencies=TRUE, repos='http://cloud.r-project.org/', lib="${params.rlocation}")
        library("gplots")
    }

    # Load input counts data
    datafiles = c( "${(input_files as List).join('", "')}" )

    # Load count column from all files into a list of data frames
    # Use data.tables fread as much much faster than read.table
    # Row names are GeneIDs
    temp <- lapply(datafiles, fread, skip="Geneid", header=TRUE, colClasses=c(NA, rep("NULL", 5), NA))

    # Merge into a single data frame
    merge.all <- function(x, y) {
        merge(x, y, all=TRUE, by="Geneid")
    }
    data <- data.frame(Reduce(merge.all, temp))

    # Clean sample name headers
    colnames(data) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(data))

    # Set GeneID as row name
    rownames(data) <- data[,1]
    data[,1] <- NULL

    # Convert data frame to edgeR DGE object
    dataDGE <- DGEList( counts=data.matrix(data) )

    # Normalise counts
    dataNorm <- calcNormFactors(dataDGE)

    # Make MDS plot
    pdf('edgeR_MDS_plot.pdf')
    MDSdata <- plotMDS(dataNorm)
    dev.off()

    # Print distance matrix to file
    write.table(MDSdata\$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE, sep="\\t")

    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))
    write.table(MDSxy, 'edgeR_MDS_plot_coordinates.txt', quote=FALSE, sep="\\t")

    # Get the log counts per million values
    logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

    # Calculate the euclidean distances between samples
    dists = dist(t(logcpm))

    # Plot a heatmap of correlations
    pdf('log2CPM_sample_distances_heatmap.pdf')
    hmap <- heatmap.2(as.matrix(dists),
      main="Sample Correlations", key.title="Distance", trace="none",
      dendrogram="row", margin=c(9, 9)
    )
    dev.off()

    # Plot the heatmap dendrogram
    pdf('log2CPM_sample_distances_dendrogram.pdf')
    plot(hmap\$rowDendrogram, main="Sample Dendrogram")
    dev.off()

    # Write clustered distance values to file
    write.table(hmap\$carpet, 'log2CPM_sample_distances.txt', quote=FALSE, sep="\\t")

    file.create("corr.done")

    # Printing sessioninfo to standard out
    print("Sample correlation info:")
    sessionInfo()
    """
}


/*
 * STEP 11 MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    memory { $params.multiqc.memory ?: 4.GB * task.attempt }
    time { $params.multiqc.time ?: 4.h * task.attempt }
    errorStrategy = 'ignore'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('star/*') from star_logs.flatten().toList()
    file ('rseqc/*') from rseqc_results.flatten().toList()
    file ('preseq/*') from preseq_results.flatten().toList()
    file ('dupradar/*') from dupradar_results.flatten().toList()
    file ('featureCounts/*') from featureCounts_logs.flatten().toList()
    file ('featureCounts_biotype/*') from featureCounts_biotype.flatten().toList()
    file ('stringtie/*') from stringtie_log.flatten().toList()
    file ('sample_correlation_results/*') from sample_correlation_results.flatten().toList()

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
    multiqc -f .
    """
}


