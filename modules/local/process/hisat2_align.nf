// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/hisat2:2.2.0--py37hfa133b6_4"
    //container "https://depot.galaxyproject.org/singularity/hisat2:2.2.0--py37hfa133b6_4"

    conda (params.conda ? "bioconda::hisat2=2.2.0" : null)

    input:
    tuple val(meta), path(reads)
    path index
    path splicesites
    val options

    output:
    path "hisat2", emit: index
    path "*.version.txt", emit: version

    script:
    def avail_mem = 0
    if (!task.memory) {
        log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
    } else {
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        avail_mem = task.memory.toGiga()
    }

    def extract_exons = ''
    def ss = ''
    def exon = ''
    if (avail_mem > params.hisat_build_memory) {
        log.info "[HISAT2 index build] Over ${params.hisat_build_memory} GB available, so using splice sites and exons in HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.exons.txt"
        ss = "--ss $splicesites"
        exon = "--exon ${gtf.baseName}.exons.txt"
    } else {
        log.info "[HISAT2 index build] Less than ${params.hisat_build_memory} GB available, so NOT using splice sites and exons in HISAT2 index."
        log.info "[HISAT2 index build] Use --hisat_build_memory [small number] to skip this check."
    }

    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    mkdir hisat2
    $extract_exons
    hisat2-build \\
        -p $task.cpus \\
        $ss \\
        $exon \\
        $ioptions.args \\
        $fasta \\
        hisat2/${fasta.baseName}
    echo \$(hisat2 --version 2>&1) | sed 's/^.*version //; s/64.*\$//' > ${software}.version.txt
    """
}

//         process HISAT2_ALIGN {
//             tag "$name"
//             label 'high_memory'
//             publishDir "${params.outdir}/hisat2", mode: params.publish_dir_mode,
//                 saveAs: { filename ->
//                     if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
//                     else if (!params.save_align_intermeds && filename == "where_are_my_files.txt") filename
//                     else if (params.save_align_intermeds && filename != "where_are_my_files.txt") filename
//                     else null
//                 }
//
//             input:
//             tuple val(name), path(reads) from trimmed_reads_alignment
//             path index from ch_hisat2_index
//             path splicesites from ch_splicesites
//             path wherearemyfiles from ch_where_are_my_files
//
//             output:
//             path "${prefix}.bam" into hisat2_bam
//             path "${prefix}.hisat2_summary.txt" into alignment_logs
//             path "where_are_my_files.txt"
//             path "unmapped.hisat2*" optional true
//
//             script:
//             index_base = index[0].toString() - ~/.\d.ht2l?/
//             prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
//             seq_center = params.seq_center ? "--rg-id ${prefix} --rg CN:${params.seq_center.replaceAll('\\s','_')} SM:$prefix" : "--rg-id ${prefix} --rg SM:$prefix"
//             def rnastrandness = ''
//             if (forward_stranded && !unstranded) {
//                 rnastrandness = params.single_end ? '--rna-strandness F' : '--rna-strandness FR'
//             } else if (reverse_stranded && !unstranded) {
//                 rnastrandness = params.single_end ? '--rna-strandness R' : '--rna-strandness RF'
//             }
//             if (params.single_end) {
//                 unaligned = params.save_unaligned ? "--un-gz unmapped.hisat2.gz" : ''
//                 """
//                 hisat2 \\
//                     -x $index_base \\
//                     -U $reads \\
//                     $rnastrandness \\
//                     --known-splicesite-infile $splicesites \\
//                     -p $task.cpus $unaligned \\
//                     --met-stderr \\
//                     --new-summary \\
//                     --dta \\
//                     $params.hisat2_align_options \\
//                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
//                     | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
//                 """
//             } else {
//                 unaligned = params.save_unaligned ? "--un-conc-gz unmapped.hisat2.gz" : ''
//                 """
//                 hisat2 \\
//                     -x $index_base \\
//                     -1 ${reads[0]} \\
//                     -2 ${reads[1]} \\
//                     $rnastrandness \\
//                     --known-splicesite-infile $splicesites \\
//                     --no-mixed \\
//                     --no-discordant \\
//                     -p $task.cpus $unaligned \\
//                     --met-stderr \\
//                     --new-summary \\
//                     $params.hisat2_align_options \\
//                     --summary-file ${prefix}.hisat2_summary.txt $seq_center \\
//                     | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
//                 """
//             }
//         }
