// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process STRINGTIE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/stringtie:2.1.2--h7e0af3c_1"
    //container "https://depot.galaxyproject.org/singularity/stringtie:2.1.2--h7e0af3c_1"

    conda (params.conda ? "bioconda::stringtie=2.1.2" : null)

    input:
    tuple val(meta), path(bam)
    path  gtf
    val   options

    // output:
    // path "${bam.baseName}" into qualimap_results
    // tuple val(meta), path("*.bam"), emit: bam
    // path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    // Figure out strandedness from pipeline parameters
    def unstranded       = params.unstranded
    def forward_stranded = params.forward_stranded
    def reverse_stranded = params.reverse_stranded
    if (params.pico) {
        unstranded       = false
        forward_stranded = true
        reverse_stranded = false
    }
    def strandedness = 'non-strand-specific'
    if (forward_stranded) {
        strandedness = 'strand-specific-forward'
    } else if (reverse_stranded) {
        strandedness = 'strand-specific-reverse'
    }
    def paired_end = meta.single_end ? '' : '-pe'
    def memory     = task.memory.toGiga() + "G"
    """
    unset DISPLAY
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        $ioptions.args \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        $paired_end \\
        -outdir $prefix
    """
}

//     process STRINGTIE {
//         tag "${bam.baseName - '.sorted'}"
//         publishDir "${params.outdir}/stringtie", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
//                 else if (filename.indexOf("cov_refs.gtf") > 0) "cov_refs/$filename"
//                 else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
//                 else "$filename"
//             }
//
//         input:
//         path bam from bam_stringtieFPKM
//         path gtf from ch_gtf
//
//         output:
//         path "${bam.baseName}_transcripts.gtf"
//         path "${bam.baseName}.gene_abund.txt"
//         path "${bam}.cov_refs.gtf"
//         path "${bam.baseName}_ballgown"
//
//         script:
//         def st_direction = ''
//         if (forward_stranded && !unstranded) {
//             st_direction = "--fr"
//         } else if (reverse_stranded && !unstranded) {
//             st_direction = "--rf"
//         }
//         def ignore_gtf = params.stringtie_ignore_gtf ? "" : "-e"
//         """
//         stringtie $bam \\
//             $st_direction \\
//             -o ${bam.baseName}_transcripts.gtf \\
//             -v \\
//             -G $gtf \\
//             -A ${bam.baseName}.gene_abund.txt \\
//             -C ${bam}.cov_refs.gtf \\
//             -b ${bam.baseName}_ballgown \\
//             $ignore_gtf
//         """
//     }
