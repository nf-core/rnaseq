////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

if (params.public_data_ids) { 
    Channel
        .from(file(params.public_data_ids, checkIfExists: true))
        .splitCsv(header:false, sep:'', strip:true)
        .map { it[0] }
        .unique()
        .set { ch_public_data_ids }
} else { 
    exit 1, 'Input file with public database ids not specified!' 
}

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { SRA_IDS_TO_RUNINFO    } from './modules/local/process/sra_ids_to_runinfo'    addParams( options: modules['sra_ids_to_runinfo']    )
include { SRA_RUNINFO_TO_FTP    } from './modules/local/process/sra_runinfo_to_ftp'    addParams( options: modules['sra_runinfo_to_ftp']    )
include { SRA_FASTQ_FTP         } from './modules/local/process/sra_fastq_ftp'         addParams( options: modules['sra_fastq_ftp']         )
include { SRA_TO_SAMPLESHEET    } from './modules/local/process/sra_to_samplesheet'    addParams( options: modules['sra_to_samplesheet'], results_dir: modules['sra_fastq_ftp'].publish_dir )
include { SRA_MERGE_SAMPLESHEET } from './modules/local/process/sra_merge_samplesheet' addParams( options: modules['sra_merge_samplesheet'] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow SRA_DOWNLOAD {

    /*
     * MODULE: Get SRA run information for public database ids
     */
    SRA_IDS_TO_RUNINFO (
        ch_public_data_ids
    )

    /*
     * MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
     */
    SRA_RUNINFO_TO_FTP (
        SRA_IDS_TO_RUNINFO.out.tsv
    )

    SRA_RUNINFO_TO_FTP
        .out
        .tsv
        .splitCsv(header:true, sep:'\t')
        .map { 
            meta -> 
                meta.single_end = meta.single_end.toBoolean()
                [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
        }
        .set { ch_sra_reads }

    if (!params.skip_sra_fastq_download) {
        /*
         * MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
         */
        SRA_FASTQ_FTP (
            ch_sra_reads.map { meta, reads -> if (meta.fastq_1)  [ meta, reads ] }
        )

        /*
         * MODULE: Stage FastQ files downloaded by SRA together and auto-create a samplesheet for the pipeline
         */
        SRA_TO_SAMPLESHEET (
            SRA_FASTQ_FTP.out.fastq
        )

        /*
         * MODULE: Create a merged samplesheet across all samples for the pipeline
         */
        SRA_MERGE_SAMPLESHEET (
            SRA_TO_SAMPLESHEET.out.csv.collect{it[1]}
        )

        /*
         * If ids don't have a direct FTP download link write them to file for download outside of the pipeline
         */
        def no_ids_file = ["${params.outdir}", "${modules['sra_fastq_ftp'].publish_dir}", "IDS_NOT_DOWNLOADED.txt" ].join(File.separator)
        ch_sra_reads
            .map { meta, reads -> if (!meta.fastq_1) "${meta.id.split('_')[0..-2].join('_')}" }
            .unique()
            .collectFile(name: no_ids_file, sort: true, newLine: true)

    }
}

////////////////////////////////////////////////////
/* --            COMPLETION EMAIL              -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Checks.sra_download(log)
    Completion.email(workflow, params, params.summary_params, baseDir, log)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

