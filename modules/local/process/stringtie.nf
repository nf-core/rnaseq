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
    // path bam from bam_stringtieFPKM
    // path gtf from ch_gtf

    // output:
    // path "${bam.baseName}" into qualimap_results
    // tuple val(meta), path("*.bam"), emit: bam
    // path  "*.version.txt"         , emit: version

    // output:
    // path "${bam.baseName}_transcripts.gtf"
    // path "${bam.baseName}.gene_abund.txt"
    // path "${bam}.cov_refs.gtf"
    // path "${bam.baseName}_ballgown"

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandedness \\
        -G $gtf \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene_abundance.txt \\
        -C ${prefix}.coverage.gtf \\
        -b ${prefix}.ballgown \\
        $ioptions.args

    stringtie --version > ${software}.version.txt
    """
}
// saveAs: { filename ->
//     if (filename.indexOf("transcripts.gtf") > 0) "transcripts/$filename"
//     else if (filename.indexOf("coverage.gtf") > 0) "cov_refs/$filename"
//     else if (filename.indexOf("ballgown") > 0) "ballgown/$filename"
//     else "$filename"
// }
