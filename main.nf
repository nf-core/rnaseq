#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               N G I - R N A S E Q    B E S T    P R A C T I C E
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started March 2016.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-RNAseq
 #### Authors
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén @Hammarn  <rickard.hammaren@scilifelab.se>
 Docker and AWS integration by
 Denis Moreno @Galithil <denis.moreno@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.0

// Configurable variables
params.project = false
params.genome = false
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.splicesites = false
params.download_hisat2index = false
params.download_fasta = false
params.download_gtf = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.reads = "data/*{1,2}.fastq.gz"
params.outdir = './results'

// R library locations
params.rlocation = false
if (params.rlocation){
    nxtflow_libs = file(params.rlocation)
    nxtflow_libs.mkdirs()
}


def single
params.sampleLevel = false

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Preset trimming options
params.pico = false
if (params.pico){
  params.clip_r1 = 3
  params.clip_r2 = 0
  params.three_prime_clip_r1 = 0
  params.three_prime_clip_r2 = 3
}

// Choose aligner
params.aligner = 'star'
if (params.aligner != 'star' && params.aligner != 'hisat2'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}

// Validate inputs
if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}
else if ( params.hisat2_index && params.aligner == 'hisat2' ){
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}
else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
else if ( ( params.aligner == 'hisat2' && !params.download_hisat2index ) && !params.download_fasta ){
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeBED12;
              gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
}
else if ( !params.download_gtf ){
    exit 1, "No GTF annotation specified!"
}
if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
}
if( params.aligner == 'hisat2' && params.splicesites ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .into { indexing_splicesites; alignment_splicesites }
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
log.info "Strandedness   : ${params.unstranded ? 'None' : params.forward_stranded ? 'Forward' : params.reverse_stranded ? 'Reverse' : 'None'}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "R libraries    : ${params.rlocation}"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Output dir     : ${params.outdir}"
log.info "Save Reference : ${params.saveReference}"
log.info "Save Trimmed   : ${params.saveTrimmed}"
log.info "Save Intermeds : ${params.saveAlignedIntermediates}"
if( params.pico       ) log.info "Trim Profile   : SMARTer Stranded Total RNA-Seq Kit - Pico Input"
if( params.clip_r1 > 0) log.info "Trim R1        : ${params.clip_r1}"
if( params.clip_r2 > 0) log.info "Trim R2        : ${params.clip_r2}"
if( params.three_prime_clip_r1 > 0) log.info "Trim 3' R1     : ${params.three_prime_clip_r1}"
if( params.three_prime_clip_r2 > 0) log.info "Trim 3' R2     : ${params.three_prime_clip_r2}"
log.info "Config Profile : " + (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
if(params.project) log.info "UPPMAX Project : ${params.project}"
log.info "========================================="

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_trimming }

/*
 * PREPROCESSING - Download Fasta
 */
if(((params.aligner == 'star' && !params.star_index) || (params.aligner == 'hisat2' && !params.hisat2_index)) && !params.fasta && params.download_fasta){
    process downloadFASTA {
        tag "${params.download_fasta}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        output:
        file "*.{fa,fasta}" into fasta

        script:
        """
        curl -O -L ${params.download_fasta}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}
/*
 * PREPROCESSING - Download GTF
 */
if(!params.gtf && params.download_gtf){
    process downloadGTF {
        tag "${params.download_gtf}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        output:
        file "*.gtf" into gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts, gtf_stringtieFPKM

        script:
        """
        curl -O -L ${params.download_gtf}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}
/*
 * PREPROCESSING - Download HISAT2 Index
 */
 if( params.aligner == 'hisat2' && params.download_hisat2index && !params.hisat2_index){
    process downloadHS2Index {
        tag "${params.download_hisat2index}"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        output:
        file "*/*.ht2" into hs2_indices

        script:
        """
        curl -O -L ${params.download_hisat2index}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}
/*
 * PREPROCESSING - Build STAR index
 */
if(params.aligner == 'star' && !params.star_index && fasta){
    process makeSTARindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta
        file gtf from gtf_makeSTARindex

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
if(params.aligner == 'hisat2' && !params.splicesites){
    process makeHisatSplicesites {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeHisatSplicesites

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
if(params.aligner == 'hisat2' && !params.hisat2_index && !params.download_hisat2index && fasta){
    process makeHISATindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf_makeHISATindex

        output:
        file "${fasta.baseName}.*.ht2" into hs2_indices

        script:
        if( task.memory == null ){
            log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
            avail_mem = 0
        } else {
            log.info "[HISAT2 index build] Available memory: ${task.memory}"
            avail_mem = task.memory.toGiga()
        }
        if( avail_mem > params.hisatBuildMemory ){
            log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
            extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
            ss = "--ss $indexing_splicesites"
            exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
        } else {
            log.info "[HISAT2 index build] Less than ${params.hisatBuildMemory} GB available, so NOT using splice sites and exons in HISAT2 index."
            log.info "[HISAT2 index build] Use --hisatBuildMemory [small number] to skip this check."
            extract_exons = ''
            ss = ''
            exon = ''
        }
        """
        $extract_exons
        hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
        """
    }
}
/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed12

        script: // This script is bundled with the pipeline, in NGI-RNAseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

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
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else params.saveTrimmed ? filename : null
        }

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    file "*fq.gz" into trimmed_reads
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    single = reads instanceof Path
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (single) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}


/*
 * STEP 3 - align with STAR
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
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if(params.aligner == 'star'){
    process star {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".out") > 0) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }

        input:
        file reads from trimmed_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()

        output:
        set file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix
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
if(params.aligner == 'hisat2'){
    process hisat2Align {
        tag "$prefix"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_log.txt") > 0) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }

        input:
        file reads from trimmed_reads
        file hs2_indices from hs2_indices.collect()
        file alignment_splicesites from alignment_splicesites.collect()

        output:
        file "${prefix}.bam" into hisat2_bam
        file "${prefix}.hisat2_log.txt" into alignment_logs

        script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        def rnastrandness = ''
        if (params.forward_stranded && !params.unstranded){
            rnastrandness = single ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (params.reverse_stranded && !params.unstranded){
            rnastrandness = single ? '--rna-strandness R' : '--rna-strandness RF'
        }
        if (single) {
            """
            set -o pipefail   # Capture exit codes from HISAT2, not samtools
            hisat2 -x $index_base \\
                   -U $reads \\
                   $rnastrandness \\
                   --known-splicesite-infile $alignment_splicesites \\
                   -p ${task.cpus} \\
                   --met-stderr \\
                   | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
                   2> ${prefix}.hisat2_log.txt
            """
        } else {
            """
            set -o pipefail   # Capture exit codes from HISAT2, not samtools
            hisat2 -x $index_base \\
                   -1 ${reads[0]} \\
                   -2 ${reads[1]} \\
                   $rnastrandness \\
                   --known-splicesite-infile $alignment_splicesites \\
                   --no-mixed \\
                   --no-discordant \\
                   -p ${task.cpus} \\
                   --met-stderr \\
                   | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
                   2> ${prefix}.hisat2_log.txt
            """
        }
    }

    process hisat2_sortOutput {
        tag "${hisat2_bam.baseName}"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename -> params.saveAlignedIntermediates ? "aligned_sorted/$filename" : null }

        input:
        file hisat2_bam

        output:
        file "${hisat2_bam.baseName}.sorted.bam" into bam_count, bam_rseqc, bam_preseq, bam_markduplicates, bam_featurecounts, bam_stringtieFPKM

        script:
        def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"
        """
        samtools sort \\
            $hisat2_bam \\
            -@ ${task.cpus} $avail_mem \\
            -o ${hisat2_bam.baseName}.sorted.bam
        """
    }
}


/*
 * STEP 4 - RSeQC analysis
 */
process rseqc {
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)                "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)              "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else "$filename"
        }

    input:
    file bam_rseqc
    file bed12 from bed12.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    def strandRule = ''
    if (params.forward_stranded && !params.unstranded){
        strandRule = single ? '-d ++,--' : '-d 1++,1--,2+-,2-+'
    } else if (params.reverse_stranded && !params.unstranded){
        strandRule = single ? '-d +-,-+' : '-d 1+-,1-+,2++,2--'
    }
    """
    samtools index $bam_rseqc
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    RPKM_saturation.py -i $bam_rseqc -r $bed12  $strandRule -o ${bam_rseqc.baseName}.RPKM_saturation
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
process preseq {
    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

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
process markDuplicates {
    tag "${bam_markduplicates.baseName - '.sorted'}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    input:
    file bam_markduplicates

    output:
    file "${bam_markduplicates.baseName}.markDups.bam" into bam_md
    file "${bam_markduplicates.baseName}.markDups_metrics.txt" into picard_results

    script:
    if( task.memory == null ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    java -Xmx${avail_mem}g -jar \$PICARD_HOME/picard.jar MarkDuplicates \\
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
process dupradar {
    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    input:
    file bam_md
    file gtf from gtf_dupradar.collect()

    output:
    file "*.{pdf,txt}" into dupradar_results

    script: // This script is bundled with the pipeline, in NGI-RNAseq/bin/
    def paired = single ? 'FALSE' :  'TRUE'
    def rlocation = params.rlocation ?: ''
    """
    dupRadar.r $bam_md $gtf $paired $rlocation
    """
}


/*
 * STEP 8 Feature counts
 */
process featureCounts {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_biotype_counts.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${bam_featurecounts.baseName}_biotype_counts.txt" into featureCounts_biotype

    script:
    def featureCounts_direction = 0
    if (params.forward_stranded && !params.unstranded) {
        featureCounts_direction = 1
    } else if (params.reverse_stranded && !params.unstranded){
        featureCounts_direction = 2
    }
    """
    featureCounts -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
    featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt > ${bam_featurecounts.baseName}_biotype_counts.txt
    """
}


/*
 * STEP 9 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()

    output:
    file 'merged_gene_counts.txt'

    script:
    """
    merge_featurecounts.py -o merged_gene_counts.txt -i $input_files
    """
}


/*
 * STEP 10 - stringtie FPKM
 */
process stringtieFPKM {
    tag "${bam_stringtieFPKM.baseName - '.sorted'}"
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
            else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
            else "$filename"
        }

    input:
    file bam_stringtieFPKM
    file gtf from gtf_stringtieFPKM.collect()

    output:
    file "${bam_stringtieFPKM.baseName}_transcripts.gtf"
    file "${bam_stringtieFPKM.baseName}.gene_abund.txt"
    file "${bam_stringtieFPKM}.cov_refs.gtf"
    stdout into stringtie_log

    script:
    def StringTie_direction = ''
    if (params.forward_stranded && !params.unstranded){
        StringTie_direction = "--fr"
    } else if (!params.reverse_stranded && !params.unstranded){
        StringTie_direction = "--rf"
    }
    """
    stringtie $bam_stringtieFPKM \\
        $StringTie_direction \\
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
 * STEP 11 - edgeR MDS and heatmap
 */
process sample_correlation {
    tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    input:
    file input_files from geneCounts.collect()
    bam_count

    output:
    file "*.{txt,pdf}" into sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in NGI-RNAseq/bin/
    def rlocation = params.rlocation ?: ''
    """
    edgeR_heatmap_MDS.r "rlocation=$rlocation" $input_files
    """
}


/*
 * STEP 12 MultiQC
 */
process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    echo true

    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect()
    file ('preseq/*') from preseq_results.collect()
    file ('dupradar/*') from dupradar_results.collect()
    file ('featureCounts/*') from featureCounts_logs.collect()
    file ('featureCounts_biotype/*') from featureCounts_biotype.collect()
    file ('stringtie/*') from stringtie_log.collect()
    file ('sample_correlation_results/*') from sample_correlation_results.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*multiqc_data"
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    """
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -f . 2>&1
    """
}


/*
 * STEP 13 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    val prefix from multiqc_prefix

    output:
    file "results_description.html"

    script:
    def rlocation = params.rlocation ?: ''
    """
    markdown_to_html.r $baseDir/docs/output.md results_description.html $rlocation
    """
}
