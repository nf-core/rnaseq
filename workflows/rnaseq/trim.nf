//
// Read QC, UMI extraction and trimming
//

include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'
include { TRIMGALORE            } from '../../../modules/nf-core/trimgalore/main'

workflow FASTQ_FASTQC_UMITOOLS_TRIM {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    trimmer           // fastp or trimgalore
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    skip_trimming     // boolean: true/false
    adapter_fasta     // file: adapter.fasta
    save_trimmed_fail // boolean: true/false
    save_merged       // boolean: true/false
    min_trimmed_reads // integer: > 0

    emit:
    reads.map { meta, fastq ->
        def fastqc_raw = !skip_fastqc
            ? FASTQC_RAW (meta, fastq)
            : null

        def umi = with_umi && !skip_umi_extract
            ? UMITOOLS_EXTRACT ( meta, fastq )
            : null
        
        // Discard R1 / R2 if required
        if (umi && umi_discard_read in [1,2] && !meta.single_end) {
            meta.single_end = true
            umi.reads = [ umi.reads[umi_discard_read % 2] ]
        }

        if (skip_trimming) {
            return [ meta, fastqc_raw, umi, null ]
        }

        if (trimmer == Trimmer.FASTP) {
            def trim = FASTP ( meta, fastq, null, save_trimmed_fail, save_merged )
            trim.num_reads = getFastpReadsAfterFiltering(trim.json)
            def adapter_seq = getFastpAdapterSequence(trim.json)
            def fastqc_trim = !skip_fastqc
                ? FASTQC_TRIM ( meta, fastq )
                : null

            return [ meta, fastqc_raw, umi, [ trim, adapter_seq, fastqc_trim ] ]
        }

        if (trimmer == Trimmer.TRIMGALORE) {
            def trim = TRIMGALORE ( meta, fastq )
            trim.num_reads = trim.log
                ? getTrimGaloreReadsAfterFiltering(meta.single_end ? trim.log : trim.log[-1])
                : null

            return [ meta, fastqc_raw, umi, [ trim ] ]
        }

        // throw new IllegalStateException()
    }
}

//
// Function that parses fastp json output file to get total number of reads after trimming
//
def getFastpReadsAfterFiltering(json_file) {
    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
}

def getFastpAdapterSequence(json_file){
    def json = new groovy.json.JsonSlurper().parseText(json_file.text)
    try{
        adapter = json['adapter_cutting']['read1_adapter_sequence']
    } catch(Exception ex){
        adapter = ""
    }
    return adapter
}

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
}

enum Trimmer {
    TRIMGALORE,
    FASTP
}
