// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DUPRADAR {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/bioconductor-dupradar:1.18.0--r40_0"
    //container "https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.18.0--r40_0"

    conda (params.conda ? "bioconda::bioconductor-dupradar=1.18.0" : null)

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
//     process DUPRADAR {
//         tag "${bam.baseName - '.sorted.markDups'}"
//         label 'high_time'
//         publishDir "${params.outdir}/dupradar", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
//                 else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
//                 else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
//                 else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
//                 else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
//                 else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
//                 else "$filename"
//             }
//
//         input:
//         path bam from bam_md
//         path gtf from ch_gtf
//
//         output:
//         path "*.{pdf,txt}" into dupradar_results
//
//         script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
//         def dupradar_direction = 0
//         if (forward_stranded && !unstranded) {
//             dupradar_direction = 1
//         } else if (reverse_stranded && !unstranded) {
//             dupradar_direction = 2
//         }
//         def paired = params.single_end ? 'single' :  'paired'
//         """
//         dupRadar.r $bam $gtf $dupradar_direction $paired $task.cpus
//         """
//     }
