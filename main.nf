#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnaseq
========================================================================================
 nf-core/rnaseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

if (params.help) {
    def command = "nextflow run nf-core/rnaseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

////////////////////////////////////////////////////
/* --          PARAMETER SUMMARY               -- */
////////////////////////////////////////////////////

summary  = Schema.params_summary(workflow, params)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(26)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Workflow to run the main pipeline
include { RNASEQ                  } from './rnaseq' addParams( summary: summary )

// Workflow to auto-create the input samplesheet for the pipeline from public database ids
// include { SRA_TO_SAMPLESHEET } from './sra_to_samplesheet'

include { SRA_IDS_TO_RUNINFO } from './modules/local/process/sra_ids_to_runinfo' addParams( options: [:] )
include { SRA_RUNINFO_TO_FTP } from './modules/local/process/sra_runinfo_to_ftp' addParams( options: [:] )
include { SRA_FASTQ_FTP      } from './modules/local/process/sra_fastq_ftp'      addParams( options: [:] )

if (params.public_data_ids) { 
    Channel
        .from(file(params.public_data_ids).readLines())
        .set { ch_public_data_ids }
} else { 
    exit 1, 'Input file with public database ids not specified!' 
}

workflow {
    if (params.public_data_ids) {
        SRA_IDS_TO_RUNINFO (
            ch_public_data_ids
        )

        SRA_RUNINFO_TO_FTP (
            SRA_IDS_TO_RUNINFO.out.tsv.collect()
        )

        SRA_RUNINFO_TO_FTP
            .out
            .tsv
            .splitCsv(header:true, sep:',')
            .map { 
                row -> [
                    [ id:row.id, single_end:row.single_end.toBoolean(), is_ftp:row.is_ftp.toBoolean(), md5_1:row.md5_1, md5_2:row.md5_2 ],
                    [ row.fastq_1, row.fastq_2 ],
                ]
            }
            .set { ch_sra_reads }

        SRA_FASTQ_FTP (
            ch_sra_reads
        )
        

    } else {
        RNASEQ ()
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
