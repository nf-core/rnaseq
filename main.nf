#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                R N A - S E Q    T W O    P O I N T    Z E R O
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 @Authors
 Phil Ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow main.nf
 
 Pipeline variables can be configured with the following command line options:
 --genome [ID]
 --index [path to STAR index]
 --gtf [path to GTF file]
 --reads [path to input files]
 --sampleLevel [set to true to run on sample and not project level, i.e skipping MDS plot]
 --strandRule [overwrite default strandRule used by RSeQC]
 
 For example:
 $ nextflow main.nf --reads 'path/to/data/sample_*_{1,2}.fq.gz'
---------------------------------------------------------------------------------------
The pipeline can determine whether the input data is single or paired end. This relies on
specifying the input files correctly. For paired en data us the example above, i.e.
'sample_*_{1,2}.fastq.gz'. Without the glob {1,2} (or similiar) the data will be treated
as single end.
----------------------------------------------------------------------------------------
 Pipeline overview:
 - FastQC - read quility control
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
 - subread featureCounts - gene counts, biotype counts, rRNA estimation.
 - String Tie - FPKMs for genes and transcripts
 - edgeR - create MDS plot and sample pairwise distance heatmap / dendrogram
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

// Configurable variables
params.genome = 'GRCh37'
params.index = params.genomes[ params.genome ].star
params.gtf   = params.genomes[ params.genome ].gtf
params.bed12 = params.genomes[ params.genome ].bed12
params.name = "RNA-Seq Best practice"
params.reads = "data/*{_1,_2}*.fastq.gz"
params.outdir = './results'

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs = file(params.rlocation)
nxtflow_libs.mkdirs()

def single
params.sampleLevel = false
params.strandRule = false

log.info "===================================="
log.info " RNAbp : RNA-Seq Best Practice v${version}"
log.info "===================================="
log.info "Reads       : ${params.reads}"
log.info "Genome      : ${params.genome}"
log.info "Index       : ${params.index}"
log.info "Annotation   : ${params.gtf}"
log.info "Current home : $HOME"
log.info "Current user : $USER"
log.info "Current path : $PWD"
log.info "R libraries  : ${params.rlocation}"
log.info "Script dir   : $baseDir"
log.info "Working dir  : $workDir"
log.info "Output dir   : ${params.outdir}"
log.info "===================================="

// Validate inputs
index = file(params.index)
gtf   = file(params.gtf)
bed12 = file(params.bed12)
if( !index.exists() ) exit 1, "Missing STAR index: $index"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: $gtf"
if( !bed12.exists() ) exit 2, "Missing BED12 annotation: $bed12"

/*
 * Create a channel for input read files
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path ->
        def prefix = readPrefix(path, params.reads)
        tuple(prefix, path)
    }
    .groupTuple(sort: true)
    .set { read_files }
 
read_files.into { read_files_fastqc; read_files_trimming; name_for_star }


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$prefix"
    
    module 'bioinfo-tools'
    module 'FastQC'
    
    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    set val(prefix), file(reads:'*') from read_files_fastqc
    
    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$prefix"
    
    module 'bioinfo-tools'
    module 'TrimGalore'
    
    cpus 3
    memory { 3.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/trim_galore", mode: 'copy'
    
    input:
    set val(prefix), file(reads:'*') from read_files_trimming
    
    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results
    
    script:
    single = reads instanceof Path
    if(single) {
        """
        trim_galore --gzip $reads
        """
    } else {
        """
        trim_galore --paired --gzip $reads
        """
    }
}



/*
 * STEP 3 - align with STAR
 * Inspired by https://github.com/AveraSD/nextflow-rnastar
 */
process star {
    tag "$prefix"
    
    module 'bioinfo-tools'
    module 'star'
    
    cpus 10
    memory '80GB'
    time  { 5.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/STAR", mode: 'copy'
    
    input:
    file index
    file gtf
    file (reads:'*') from trimmed_reads
    
    output:
    set file('*Log.final.out'), file ('*.bam') into aligned
    file '*.out' into star_logs
    file '*SJ.out.tab'
    
    script:
    """
    #Getting STAR prefix
    f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1}
    prefix=\$f
    STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outFileNamePrefix \$prefix
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
    .set { SPLIT_BAMS }
SPLIT_BAMS.into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }


/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    tag "$bam_rseqc"
    
    module 'bioinfo-tools'
    module 'rseqc'
    module 'samtools'
    memory { 32.GB * task.attempt }
    time  {7.h * task.attempt }
    
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/rseqc" , mode: 'copy'
    
    input:
    file bam_rseqc
    file bed12 from bed12
    
    output:
    file '*.{txt,pdf,r,xls}' into rseqc_results
    /*  The following files are being generated by this process:
        .bam_stat.txt                     // bam_stat
        .splice_events.{txt,pdf}            // junction_annotation
        .splice_junction.{txt,pdf}           // junction_annotation
        .junctionSaturation_plot.{txt,pdf,r}   // junction_saturation
        .inner_distance.{txt,pdf}            // inner_distance
        .curves.{txt,pdf}                  // geneBody_coverage
        .geneBodyCoverage.txt               // geneBody_coverage
        .heatMap.{txt,pdf}                 // geneBody_coverage
        .infer_experiment.txt               // infer_experiment
        .read_distribution.txt              // read_distribution
        DupRate.xls                       // read_duplication
        DupRate_plot.pdf                   // read_duplication
        .saturation.{txt,pdf}               // RPKM_saturation
    */
    
    script:
    if (!params.strandRule){
        if (single){
            strandRule ='++,--'
        } else {
            strandRule = '1+-,1-+,2++,2--'
        }
    } else {
        strandRule = params.strandRule
    }
    
    """
    samtools index $bam_rseqc
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc}.infer_experiment.txt
    RPKM_saturation.py -i $bam_rseqc -r $bed12 -d $strandRule -o ${bam_rseqc}.RPKM_saturation
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc}.rseqc -r $bed12
    geneBody_coverage.py -i ${bam_rseqc} -o ${bam_rseqc}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc}.read_duplication
    echo "Filename $bam_rseqc RseQC version: "$(read_duplication.py --version)
    """
}



/*
 * STEP 5 - preseq analysis
 */
process preseq {
    tag "$bam_preseq"
    
    module 'bioinfo-tools'
    module 'preseq'
    
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/preseq", mode: 'copy'
    
    input:
    file bam_preseq
    
    output:
    file '*.txt' into preseq_results
    
    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq}.ccurve.txt
    echo "File name: $bam_preseq ----- preseq version: "$(preseq)
    """
}


/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {
    tag "$bam_markduplicates"
    
    module 'bioinfo-tools'
    module 'picard/2.0.1'
    
    memory { 16.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/markDuplicates", mode: 'copy'
    
    input:
    file bam_markduplicates
    
    output:
    file '*.markDups.bam' into bam_md
    file '*markDups_metrics.txt' into picard_results
    file 'versioninfo.txt'
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
    echo "File name: $bam_markduplicates MarkDuplicates version "$(java -Xmx2g -jar $PICARD_HOME/picard.jar  MarkDuplicates --version 2>&1)
    """
}


/*
 * STEP 7 - dupRadar
 */
process dupradar {
    tag "$bam_md"
    
    module 'bioinfo-tools'
    module 'R/3.2.3'
    
    memory { 16.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/dupradar", pattern: '*.{pdf,txt}', mode: 'copy'
    
    input:
    file bam_md
    file gtf from gtf
    
    output:
    file '*.{pdf,txt}' into dupradar_results
    
    script:
    def paired
    if(single) {
        paired= 'FALSE'
    } else {
        paired = 'TRUE'
    }
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
    write.table(dm, file=paste("$bam_md", "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\t")
    
    # 2D density scatter plot
    pdf(paste0("$bam_md", "_duprateExpDens.pdf"))
    duprateExpDensPlot(DupMat=dm)
    title("Density scatter plot")
    dev.off()
    fit <- duprateExpFit(DupMat=dm)
    cat("duprate at low read counts: ", fit\$intercept, "progression of the duplication rate: ", fit\$slope, "\n",
        fill=TRUE, labels="${bam_md}", file=paste0("$bam_md", "_intercept_slope.txt"), append=FALSE )
    
    # Distribution of RPK values per gene
    pdf(paste0("$bam_md", "_expressionHist.pdf"))
    expressionHist(DupMat=dm)
    title("Distribution of RPK values per gene")
    dev.off()
    
    #Printing sessioninfo to standard out
    print("$bam_md")
    sessionInfo()
    """

}



/*
 * STEP 8 Feature counts
 */
process featureCounts {
    tag "$bam_featurecounts"
    
    module 'bioinfo-tools'
    module 'subread'
    
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/featureCounts", mode: 'copy'
    
    input:
    file bam_featurecounts
    file gtf from gtf
    
    output:
    file '*_gene.featureCounts.txt' into geneCounts
    file '*_biotype.featureCounts.txt'
    file '*_rRNA_counts.txt'
    file '*.summary' into featureCounts_logs
    
    script:
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts}_gene.featureCounts.txt -p -s 2 $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts}_biotype.featureCounts.txt -p -s 2 $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts}_biotype.featureCounts.txt | sed '1,2d' | grep 'rRNA' > ${bam_featurecounts}_rRNA_counts.txt
    """
}


/*
 * STEP 9 - stringtie FPKM
 */
process stringtieFPKM {
    tag "$bam_stringtieFPKM"
    
    module 'bioinfo-tools'
    module 'StringTie'
    
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy'
    
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

    echo "Stringtie version "$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }



/*
 * STEP 10 - edgeR MDS and heatmap
 */
process sample_correlation {
    module 'bioinfo-tools'
    module 'R/3.2.3'
    
    memory { 16.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors '-1'
    
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'
    
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
    write.table(MDSdata\$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE, sep="\t")
    
    # Print plot x,y co-ordinates to file
    MDSxy = MDSdata\$cmdscale.out
    colnames(MDSxy) = c(paste(MDSdata\$axislabel, '1'), paste(MDSdata\$axislabel, '2'))
    write.table(MDSxy, 'edgeR_MDS_plot_coordinates.txt', quote=FALSE, sep="\t")
    
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
    write.table(hmap\$carpet, 'log2CPM_sample_distances.txt', quote=FALSE, sep="\t")
    
    file.create("corr.done")

    #Printing sessioninfo to standard out
    sessionInfo()
    """
}


/*
 * STEP 11 MultiQC
 */
process multiqc {
    module 'bioinfo-tools'
    
    memory '4GB'
    time '4h'
    
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    errorStrategy 'ignore'
    
    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('star/*') from star_logs.toList()
    file ('rseqc/*') from rseqc_results.toList()
    file ('preseq/*') from preseq_results.toList()
    file ('dupradar/*') from dupradar_results.toList()
    file ('featureCounts/*') from featureCounts_logs.toList()
    file ('stringtie/*') from stringtie_log.toList()
    file ('sample_correlation_results/*') from sample_correlation_results.toList()
    
    output:
    file '*multiqc_report.html'
    file '*multiqc_data'
    
    script:
    """
    # Load MultiQC with environment module if not already in PATH
    type multiqc >/dev/null 2>&1 || { module load multiqc; };
    multiqc -f -t ngi .
    """
}


/*
 * Helper function, given a file Path
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 *
 * For example:
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 *
 * Returns:
 *   'file_alpha'
 */
def readPrefix( Path actual, template ) {
    
    final fileName = actual.getFileName().toString()
    
    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') )
        filePattern = '*' + filePattern
    
    def regex = filePattern
        .replace('.','\\.')
        .replace('*','(.*)')
        .replace('?','(.?)')
        .replace('{','(?:')
        .replace('}',')')
        .replace(',','|')
    
    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {
        def end = matcher.end(matcher.groupCount() )
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
            prefix=prefix[0..-2]
        return prefix
    }
    return fileName
}
