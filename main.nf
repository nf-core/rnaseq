#!/usr/bin/env nextflow
/*
===============================================================
 nf-core/rnaseq
===============================================================
 RNA-Seq Analysis Pipeline. Started March 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
 #### Authors
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén @Hammarn  <rickard.hammaren@scilifelab.se>
---------------------------------------------------------------
*/

nextflow.enable.modules = true
include 'modules/help' as _ 

// Show help emssage
if (params.help){
    print_usage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  }

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false


ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt")
ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt")
ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt")
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
Channel.fromPath("$baseDir/assets/where_are_my_files.txt").set{ ch_where }

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
if (params.pico){
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}

// Validate inputs
if (params.aligner != 'star' && params.aligner != 'hisat2'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2'"
}
if( params.star_index && params.aligner == 'star' ){
    Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
        .set { ch_star_index }
}
else if ( params.hisat2_index && params.aligner == 'hisat2' ){
    Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
        .set { ch_hs2_indices }
}
else if ( params.fasta ){
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
           .set { ch_fasta }
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set { ch_gtf_make }

} else if( params.gff ){
  Channel
        .fromPath(params.gff)
        .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
        .set { ch_gtf_file }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}

if( params.bed12 ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .set { ch_bed12 }
}
if( params.aligner == 'hisat2' && params.splicesites ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .set { ch_splicesites }
}
if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set { ch_raw_reads }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set { ch_raw_reads }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set { ch_raw_reads }
}

include 'modules/rnaseq' as _ params(params)


log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
 nf-core/rnaseq : RNA-Seq Best Practice v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
if( params.pico ) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness'] = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
summary['Trim R1'] = clip_r1
summary['Trim R2'] = clip_r2
summary["Trim 3' R1"] = three_prime_clip_r1
summary["Trim 3' R2"] = three_prime_clip_r2
if(params.aligner == 'star'){
    summary['Aligner'] = "STAR"
    if(params.star_index)          summary['STAR Index']   = params.star_index
    else if(params.fasta)          summary['Fasta Ref']    = params.fasta
} else if(params.aligner == 'hisat2') {
    summary['Aligner'] = "HISAT2"
    if(params.hisat2_index)        summary['HISAT2 Index'] = params.hisat2_index
    else if(params.fasta)          summary['Fasta Ref']    = params.fasta
    if(params.splicesites)         summary['Splice Sites'] = params.splicesites
}
if(params.gtf)                 summary['GTF Annotation']  = params.gtf
if(params.gff)                 summary['GFF3 Annotation']  = params.gff
if(params.bed12)               summary['BED Annotation']  = params.bed12
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) {
    summary['E-mail Address'] = params.email
    summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard'){
    if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
        log.error "====================================================\n" +
                  "  WARNING! You are running with the default 'standard'\n" +
                  "  pipeline config profile, which runs on the head node\n" +
                  "  and assumes all software is on the PATH.\n" +
                  "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                  "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                  "============================================================"
    }
}

/*
 * PREPROCESSING - Build STAR index
 */
if(params.aligner == 'star' && !params.star_index && params.fasta){
    makeSTARindex(ch_fasta, ch_gtf_make) .set{ ch_star_index }
}

/*
 * PREPROCESSING - Build HISAT2 splice sites file
 */
if(params.aligner == 'hisat2' && !params.splicesites){

    makeHisatSplicesites(ch_gtf_make).set { ch_splicesites } 

}
/*
 * PREPROCESSING - Build HISAT2 index
 */
if(params.aligner == 'hisat2' && !params.hisat2_index && params.fasta){

    makeHISATindex(
        ch_fasta, 
        ch_splicesites, 
        ch_gtf_make) 
        .set { ch_hs2_indices} 

}

/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if(params.gff){
  convertGFFtoGTF(ch_gtf_file).set { ch_gtf_make }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    makeBED12(ch_gtf_make).set { ch_bed12 }
}


/*
 * STEP 1 - FastQC
 */


 fastqc(ch_raw_reads)


/*
 * STEP 2 - Trim Galore!
 */
 trim_galore(ch_raw_reads, ch_where.collect())

 trim_galore.output.first.set { trimmed_reads }
 trim_galore.output.second.set { trimgalore_results }
 trim_galore.output.third.set { trimgalore_fastqc_reports }

/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
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
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if(params.aligner == 'star'){
    Channel.from(false).set { hisat_stdout }

    star( trimmed_reads, ch_star_index.collect(), ch_gtf_make.collect(), ch_where.collect() )

    // Filter removes all 'aligned' channels that fail the check
    star.output
        .first 
        .filter { logs, bams -> rnaseq.check_log(logs) }
        .flatMap {  logs, bams -> bams }
        .set { ch_bam }
    
    star.output.second.set { alignment_logs }
    star.output.third.set { star_log }
}

/*
 * STEP 3 - align with HISAT2
 */
if(params.aligner == 'hisat2'){
    Channel.from(false).set { star_log }

    hisat2Align( 
            trimmed_reads, 
            ch_hs2_indices.collect(), 
            ch_splicesites.collect(), 
            ch_where.collect() ) 

    hisat2Align.output.first.set { hisat2_bam } 
    hisat2Align.output.second.set { alignment_logs } 

    hisat2_sort( hisat2_bam, ch_where.collect() )  
    hisat2_sort.output.first.set { ch_bam }
    hisat2_sort.output.second.set { ch_bam_index }
}

/*
 * STEP 4 - RSeQC analysis
 */

 rseqc( bam_ch, ch_bam_index, ch_bed12.collect() )

/*
 * Step 4.1 Rseqc create BigWig coverage
 */

createBigWig(ch_bam, ch_bam_index) 

/*
 * Step 4.2 Rseqc genebody_coverage
 */

genebody_coverage(bigwig.output, ch_bed12.collect()) 


/*
 * STEP 5 - preseq analysis
 */
 preseq( ch_bam )


/*
 * STEP 6 Mark duplicates
 */

 markDuplicates( ch_bam )
 markDuplicates.output.first.set { bam_md }
 markDuplicates.output.second.set { picard_results }


/*
 * STEP 7 - dupRadar
 */

 dupradar( markDuplicates.output.first, ch_gtf_make.collect() ) 


/*
 * STEP 8 Feature counts
 */
 featureCounts( ch_bam, ch_gtf_make.collect(), ch_biotypes_header.collect() )
 featureCounts.output.first.set {ch_feature_counts }
 featureCounts.output.second { featureCounts_logs }
 featureCounts.output.third { featureCounts_biotype }

/*
 * STEP 9 - Merge featurecounts
 */
 merge_featureCounts( ch_feature_counts.collect() )


/*
 * STEP 10 - stringtie FPKM
 */
 stringtieFPKM(ch_bam, ch_gtf_make.collect())
 stringtieFPKM.output.fourth.set { stringtie_log }


/*
 * STEP 11 - edgeR MDS and heatmap
 */
 sample_correlation(
     ch_feature_counts.collect(),
     ch_bam.count(), 
     ch_mdsplot_header, 
     ch_heatmap_header
 ) 

/*
 * Parse software version numbers
 */
get_software_versions() 


/*
 * Pipeline parameters to go into MultiQC report
 */
 workflow_summary_mqc() 


/*
 * STEP 12 MultiQC
 */
 multiqc(
    ch_multiqc_config,
    fastqc.output.collect().ifEmpty([]),
    trimgalore_results.collect(),
    alignment_logs.collect(),
    rseqc.output.collect().ifEmpty([]),
    genebody_coverage.output.collect().ifEmpty([]),
    preseq.output.collect().ifEmpty([]),
    dupradar.output.collect().ifEmpty([]),
    featureCounts_logs.collect(),
    featureCounts_biotype.collect(),
    stringtie_log.collect(),
    sample_correlation.output.collect().ifEmpty([]),
    get_software_versions.output.collect(),
    workflow_summary_mqc.output.collect()
 ) 


/*
 * STEP 13 - Output Description HTML
 */
 output_documentation( ch_output_docs )


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nfcore/rnaseq] Successful: $workflow.runName"
    if(skipped_poor_alignment.size() > 0){
        subject = "[nfcore/rnaseq] Partially Successful (${skipped_poor_alignment.size()} skipped): $workflow.runName"
    }
    if(!workflow.success){
      subject = "[nfcore/rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = skipped_poor_alignment

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skip_multiqc) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nfcore/rnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/rnaseq] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    ngirnaseqlogo = new File("$baseDir/assets/nfcore-rnaseq_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    if(skipped_poor_alignment.size() > 0){
        log.info "[nfcore/rnaseq] WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor alignment scores!"
    }

    log.info "[nfcore/rnaseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }

}
