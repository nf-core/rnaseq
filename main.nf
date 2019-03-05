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

summary = create_summary()

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
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
        .set { ch_bam }
    
    star.output.second.set { alignment_logs }
    star.output.fourth.set { star_log }
    star.output.sixth.set { ch_bam_index }
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

 rseqc( ch_bam, ch_bam_index, ch_bed12.collect() )

/*
 * Step 4.1 Rseqc create BigWig coverage
 */

createBigWig(ch_bam, ch_bam_index) 

/*
 * Step 4.2 Rseqc genebody_coverage
 */

genebody_coverage(createBigWig.output, ch_bed12.collect()) 

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
 featureCounts.output.second.set { featureCounts_logs }
 featureCounts.output.third.set { featureCounts_biotype }

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

    send_email(summary, skipped_poor_alignment)

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
