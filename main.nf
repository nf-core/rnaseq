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
params.star_index = false
params.fasta = false
params.gtf = false
params.bed12 = false
if( params.genomes && !params.star_index && !params.fasta && !params.gtf && !params.bed12 ) {
    params.star_index = params.genomes[ params.genome ].star
    params.fasta = params.genomes[ params.genome ].fasta
    params.gtf   = params.genomes[ params.genome ].gtf
    params.bed12 = params.genomes[ params.genome ].bed12
}
params.hisat_index = false
params.hisatBuildMemory = 200
params.reads = "data/*{1,2}.fastq.gz"
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

// Choose aligner
params.aligner = 'star'
if (params.aligner != 'star' && params.aligner != 'hisat2'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}

// Validate inputs
def star_index, hisat_index, fasta, gtf, bed12, indexing_splicesites, alignment_splicesites

if( params.star_index && params.aligner == 'star' ){
    star_index = file(params.star_index)
    if( !star_index.exists() ) exit 1, "STAR index not found: $star_index"
}
else if ( params.hisat_index && params.aligner == 'hisat2' ){
    hisat_index = file(params.hisat_index)
    if( !star_index.exists() ) exit 1, "HISAT2 index not found: $hisat_index"
}
else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: $fasta"
}
else if ( ( params.aligner == 'hisat2' && !params.download_hisat2index ) && !params.download_fasta ){
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF annotation file not found: $gtf"
}
if( params.bed12 ){
    bed12 = file(params.bed12)
    if( !bed12.exists() ) exit 1, "BED12 annotation file not found: $bed12"
}
if( params.aligner == 'hisat2' && params.splicesites ){
    indexing_splicesites = file(params.splicesites)
    alignment_splicesites = file(params.splicesites)
    if( !alignment_splicesites.exists() ) exit 1, "HISAT2 splice sites file not found: $alignment_splicesites"
}
if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project"

// Header log info
log.info "========================================="
log.info " NGI-RNAseq : RNA-Seq Best Practice v${version}"
log.info "========================================="
log.info "Reads          : ${params.reads}"
log.info "Genome         : ${params.genome}"
if(params.aligner == 'star'){
    log.info "Aligner        : STAR"
    if(params.star_index)          log.info "STAR Index     : ${params.star_index}"
    else if(params.fasta)          log.info "Fasta Ref      : ${params.fasta}"
    else if(params.download_fasta) log.info "Fasta URL      : ${params.download_fasta}"
} else if(params.aligner == 'hisat2') {
    log.info "Aligner        : HISAT2"
    if(params.hisat2_index)        log.info "HISAT2 Index   : ${params.hisat2_index}"
    else if(params.download_hisat2index) log.info "HISAT2 Index   : ${params.download_hisat2index}"
    else if(params.fasta)          log.info "Fasta Ref      : ${params.fasta}"
    else if(params.download_fasta) log.info "Fasta URL      : ${params.download_fasta}"
    if(params.splicesites)         log.info "Splice Sites   : ${params.splicesites}"
}
if(params.gtf)                 log.info "GTF Annotation : ${params.gtf}"
else if(params.download_gtf)   log.info "GTF URL        : ${params.download_gtf}"
if(params.bed12)               log.info "BED Annotation : ${params.bed12}"
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
log.info "Config Profile : " + (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
log.info "========================================="

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }

/*
 * PREPROCESSING - Download Fasta
 */
if(!params.star_index && !params.fasta && params.download_fasta){
    process downloadFASTA {
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        val dl_url from params.download_fasta

        output:
        file "*.{fa,fasta}" into fasta

        script:
        """
        curl -O -L $dl_url
        """
    }
}
/*
 * PREPROCESSING - Download GTF
 */
if(!params.gtf && params.downloadGTF){
    process downloadGTF {
        tag params.downloadGTF
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        val url from params.downloadGTF

        output:
        file "*.gtf" into gtf

        script:
        """
        curl -O -L $url
        """
    }
}
/*
 * PREPROCESSING - Download HISAT2 Index
 */
 if( params.aligner == 'hisat2' && params.download_hisat2index && !params.hisat_index){
    process downloadGTF {
        tag params.downloadGTF
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        val url from params.download_hisat2index

        output:
        file "*/*.1.ht2" into hisat2_index // Use single file as a placeholder for the base
        file "*/*.ht2" into hs2_indices

        script:
        """
        curl -O -L $url
        tar xzf *.tar.gz
        """
    }
}
/*
 * PREPROCESSING - Build STAR index
 */
params.makeSTARindex_cpus = 12
params.makeSTARindex_memory = 30.GB
params.makeSTARindex_time = 5.h
if(params.aligner == 'star' && !params.star_index && fasta){
    process makeSTARindex {
        tag fasta
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        cpus params.makeSTARindex_cpus
        memory params.makeSTARindex_memory
        time params.makeSTARindex_time
        errorStrategy 'terminate'

        input:
        file fasta from fasta
        file gtf from gtf

        output:
        file "star" into star_index
        
        script:
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang 149 \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}
/*
 * PREPROCESSING - Build HISAT2 splice sites file
 */
params.makeHisatSplicesites_time = 2.h
if(params.aligner == 'hisat2' && !params.splicesites){
    process makeHisatSplicesites {
        tag gtf
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        time params.makeHisatSplicesites_time
        errorStrategy 'terminate'

        input:
        file gtf from gtf

        output:
        file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

        script:
        """
        hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
        """
    }
}
/*
 * PREPROCESSING - Build HISAT2 index
 */
params.makeHISATindex_cpus = 10
params.makeHISATindex_memory = 16.GB
params.makeHISATindex_time = 5.h
if(params.aligner == 'hisat2' && !params.hisat_index && !params.download_hisat2index && fasta){
    process makeHISATindex {
        tag fasta
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        cpus params.makeHISATindex_cpus
        memory params.makeHISATindex_memory
        time params.makeHISATindex_time
        errorStrategy 'terminate'

        input:
        file fasta from fasta
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf

        output:
        file "${fasta.baseName}.hisat2_index.1.ht2" into hisat2_index // Use a fake file as a placeholder for the base file
        file "${fasta.baseName}.*.ht2" into hs2_indices

        script:
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        if( task.memory.toGiga() > params.hisatBuildMemory ){
            log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
            extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
            ss = "--ss $indexing_splicesites"
            exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
        } else {
            log.info "[HISAT2 index build] Less than ${params.hisatBuildMemory} GB available, so NOT using splice sites and exons in HISAT2 index."
            log.info "[HISAT2 index build] Use --hisatBuildMemory [small number] and/or --makeHISATindex_memory [big number] to override."
            extract_exons = ''
            ss = ''
            exon = ''
        }
        """
        $extract_exons
        hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
        touch ${fasta.baseName}.hisat2_index
        """
    }
}
/*
 * PREPROCESSING - Build BED12 file
 */
params.makeBED12_time = 2.h
if(!params.bed12){
    process makeBED12 {
        tag gtf
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        time params.makeBED12_time
        errorStrategy 'terminate'

        input:
        file gtf from gtf

        output:
        file "${gtf.baseName}.bed" into bed12

        script:
        """
        convert2bed \\
            --input=gtf \\
            --output=bed \\
            --max-mem=${task.memory.toGiga()}G \\
            < $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * STEP 1 - FastQC
 */
params.fastqc_memory = 2.GB
params.fastqc_time = 4.h
process fastqc {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    memory { params.fastqc_memory * task.attempt }
    time { params.fastqc_time * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
params.trim_galore_cpus = 2
params.trim_galore_memory = 4.GB
params.trim_galore_time = 8.h
process trim_galore {
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    cpus params.trim_galore_cpus
    memory { params.trim_galore_memory * task.attempt }
    time { params.trim_galore_time * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    file "*fq.gz" into trimmed_reads
    file "*trimming_report.txt" into trimgalore_results

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
params.star_cpus = 10
params.star_memory = 80.GB
params.star_time = 5.h
if(params.aligner == 'star'){
    process star {
        tag "$reads"
        publishDir "${params.outdir}/STAR", mode: 'copy'

        cpus params.star_cpus
        memory { params.star_memory * task.attempt }
        time { params.star_time * task.attempt }
        errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }

        input:
        file star_index from star_index
        file gtf from gtf
        file reads from trimmed_reads

        output:
        set file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"

        script:
        """
        #Getting STAR prefix
        f=($reads);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
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
    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }
}


/*
 * STEP 3 - align with HISAT2
 */
params.star_cpus = 10
params.star_memory = 80.GB
params.star_time = 5.h
if(params.aligner == 'hisat2'){
    process hisat2 {
        tag "$reads"
        publishDir "${params.outdir}/HISAT2", mode: 'copy'

        cpus params.star_cpus
        memory { params.star_memory * task.attempt }
        time { params.star_time * task.attempt }
        errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }

        input:
        file index from hisat2_index // placeholder filename stub
        file hs2_indices
        file alignment_splicesites from alignment_splicesites
        file reads from trimmed_reads

        output:
        file "*.bam" into bam_count, bam_rseqc, bam_preseq, bam_markduplicates, bam_featurecounts, bam_stringtieFPKM
        file "*.hisat2_log.txt" into alignment_logs

        script:
        index_base = index - ~/.1.ht2$/
        if (single) {
            """
            f='$reads';f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
            hisat2 -x $index_base \\
                   -U $reads \\
                   --known-splicesite-infile $alignment_splicesites \\
                   -p ${task.cpus} \\
                   --met-stderr \\
                   | samtools view -bS -F 4 -F 256 - > \${f}.bam
                   2> \${f}.hisat2_log.txt
            """
        } else {
            """
            f=($reads);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_1};f=\${f%_R1}
            hisat2 -x $index_base \\
                   -1 $reads[0] \\
                   -2 $reads[0] \\
                   --known-splicesite-infile $alignment_splicesites \\
                   --no-mixed \\
                   --no-discordant \\
                   -p ${task.cpus} \\
                   --met-stderr \\
                   | samtools view -bS -F 4 -F 8 -F 256 - > \${f}.bam
                   2> \${f}.hisat2_log.txt
            """
        }
    }
}




/*
 * STEP 4 - RSeQC analysis
 */
params.rseqc_memory = 32.GB
params.rseqc_time = 7.h
process rseqc {
    tag "$bam_rseqc"
    publishDir "${params.outdir}/rseqc" , mode: 'copy'

    memory { params.rseqc_memory * task.attempt }
    time { params.rseqc_time * task.attempt }

    input:
    file bam_rseqc
    file bed12 from bed12

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results
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
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    RPKM_saturation.py -i $bam_rseqc -r $bed12 -d $strandRule -o ${bam_rseqc.baseName}.RPKM_saturation
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    geneBody_coverage.py -i ${bam_rseqc.baseName} -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    echo "Filename $bam_rseqc RseQC version: "\$(read_duplication.py --version)
    """
}



/*
 * STEP 5 - preseq analysis
 */
params.preseq_memory = 4.GB
params.preseq_time = 2.h
process preseq {
    tag "$bam_preseq"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    memory { params.preseq_memory * task.attempt }
    time { params.preseq_time * task.attempt }

    input:
    file bam_preseq

    output:
    file "${bam_preseq.baseName}.ccurve.txt" into preseq_results

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
    echo "File name: $bam_preseq  preseq version: "\$(preseq)
    """
}


/*
 * STEP 6 Mark duplicates
 */
params.markDuplicates_memory = 16.GB
params.markDuplicates_time = 2.h
process markDuplicates {
    tag "$bam_markduplicates"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy'

    memory { params.markDuplicates_memory * task.attempt }
    time { params.markDuplicates_time * task.attempt }

    input:
    file bam_markduplicates

    output:
    file "${bam_markduplicates.baseName}.markDups.bam" into bam_md
    file "${bam_markduplicates.baseName}.markDups_metrics.txt" into picard_results

    script:
    """
    java -Xmx2g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
        INPUT=$bam_markduplicates \\
        OUTPUT=${bam_markduplicates.baseName}.markDups.bam \\
        METRICS_FILE=${bam_markduplicates.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT

    # Print version number to standard out
    echo "File name: $bam_markduplicates Picard version "\$(java -Xmx2g -jar \$PICARD_HOME/picard.jar  MarkDuplicates --version 2>&1)
    """
}


/*
 * STEP 7 - dupRadar
 */
params.dupradar_memory = 16.GB
params.dupradar_time = 2.h
process dupradar {
    tag "${bam_md.baseName}"
    publishDir "${params.outdir}/dupradar", pattern: '*.{pdf,txt}', mode: 'copy'

    memory { params.dupradar_memory * task.attempt }
    time { params.dupradar_time * task.attempt }

    input:
    file bam_md
    file gtf from gtf

    output:
    file "*.{pdf,txt}" into dupradar_results

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
    write.table(dm, file=paste("${bam_md.baseName}", "_dupMatrix.txt", sep=""), quote=F, row.name=F, sep="\\t")

    # 2D density scatter plot
    pdf(paste0("${bam_md.baseName}", "_duprateExpDens.pdf"))
    duprateExpDensPlot(DupMat=dm)
    title("Density scatter plot")
    mtext("${bam_md.baseName}", side=3)
    dev.off()
    fit <- duprateExpFit(DupMat=dm)
    cat(
      paste("- dupRadar Int (duprate at low read counts):", fit\$intercept),
      paste("- dupRadar Sl (progression of the duplication rate):", fit\$slope),
      fill=TRUE, labels="${bam_md.baseName}",
      file=paste0("${bam_md.baseName}", "_intercept_slope.txt"), append=FALSE
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
      file=paste0("${bam_md.baseName}", "_duprateExpDensCurve.txt"),
      quote=FALSE, row.names=FALSE
    )
    
    # Distribution of expression box plot
    pdf(paste0("${bam_md.baseName}", "_duprateExpBoxplot.pdf"))
    duprateExpBoxplot(DupMat=dm)
    title("Percent Duplication by Expression")
    mtext("${bam_md.baseName}", side=3)
    dev.off()
    
    # Distribution of RPK values per gene
    pdf(paste0("${bam_md.baseName}", "_expressionHist.pdf"))
    expressionHist(DupMat=dm)
    title("Distribution of RPK values per gene")
    mtext("${bam_md.baseName}", side=3)
    dev.off()

    # Printing sessioninfo to standard out
    print("${bam_md.baseName}")
    citation("dupRadar")
    sessionInfo()
    """
}


/*
 * STEP 8 Feature counts
 */
params.dupradar_memory = 4.GB
params.dupradar_time = 2.h
process featureCounts {
    tag "$bam_featurecounts"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    memory { params.dupradar_memory * task.attempt }
    time { params.dupradar_time * task.attempt }

    input:
    file bam_featurecounts
    file gtf from gtf

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${bam_featurecounts.baseName}_biotype_counts.txt" into featureCounts_biotype

    script:
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s 2 $bam_featurecounts
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s 2 $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt > ${bam_featurecounts.baseName}_biotype_counts.txt
    """
}


/*
 * STEP 9 - stringtie FPKM
 */
params.dupradar_memory = 4.GB
params.dupradar_time = 2.h
process stringtieFPKM {
    tag "$bam_stringtieFPKM"
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy'

    memory { params.dupradar_memory * task.attempt }
    time { params.dupradar_time * task.attempt }

    input:
    file bam_stringtieFPKM
    file gtf from gtf

    output:
    file "${bam_stringtieFPKM.baseName}_transcripts.gtf"
    file "${bam_stringtieFPKM.baseName}.gene_abund.txt"
    file "${bam_stringtieFPKM}.cov_refs.gtf"
    stdout into stringtie_log

    script:
    """
    stringtie $bam_stringtieFPKM \\
        -o ${bam_stringtieFPKM.baseName}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM.baseName}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM.baseName}_ballgown

    echo "File name: $bam_stringtieFPKM Stringtie version "\$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }



/*
 * STEP 10 - edgeR MDS and heatmap
 */
params.dupradar_memory = 16.GB
params.dupradar_time = 2.h
process sample_correlation {
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    memory { params.dupradar_memory * task.attempt }
    time { params.dupradar_time * task.attempt }

    input:
    file input_files from geneCounts.toList()
    bam_count

    output:
    file "*.{txt,pdf}" into sample_correlation_results

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
params.multiqc_memory = 4.GB
params.multiqc_time = 4.h
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    memory { params.multiqc_memory * task.attempt }
    time { params.multiqc_time * task.attempt }
    errorStrategy 'ignore'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('alignment/*') from alignment_logs.flatten().toList()
    file ('rseqc/*') from rseqc_results.flatten().toList()
    file ('preseq/*') from preseq_results.flatten().toList()
    file ('dupradar/*') from dupradar_results.flatten().toList()
    file ('featureCounts/*') from featureCounts_logs.flatten().toList()
    file ('featureCounts_biotype/*') from featureCounts_biotype.flatten().toList()
    file ('stringtie/*') from stringtie_log.flatten().toList()
    file ('sample_correlation_results/*') from sample_correlation_results.flatten().toList()

    output:
    file "*multiqc_report.html"
    file "*multiqc_data"

    script:
    """
    multiqc -f .
    """
}


