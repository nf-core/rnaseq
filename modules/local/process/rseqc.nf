// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process RSEQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    //container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"

    conda (params.conda ? "bioconda::rseqc=3.0.1" : null)

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed
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

//     process RSEQC {
//         tag "${bam.baseName - '.sorted'}"
//         label 'mid_memory'
//         publishDir "${params.outdir}/rseqc" , mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                 if (filename.indexOf("bam_stat.txt") > 0)                           "bam_stat/$filename"
//                 else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
//                 else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
//                 else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
//                 else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
//                 else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
//                 else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
//                 else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
//                 else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
//                 else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
//                 else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
//                 else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
//                 else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
//                 else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
//                 else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
//                 else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
//                 else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
//                 else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
//                 else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
//                 else if (filename.indexOf("junction_annotation_log.txt") > 0)       "junction_annotation/$filename"
//                 else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
//                 else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
//                 else filename
//             }
//
//         when:
//         !params.skip_qc && !params.skip_rseqc
//
//         input:
//         path bam from bam_rseqc
//         path bai from bam_index_rseqc
//         path bed12 from ch_bed12
//
//         output:
//         path "*.{txt,pdf,r,xls}" into rseqc_results
//
//         script:
//         """
//         infer_experiment.py -i $bam -r $bed12 > ${bam.baseName}.infer_experiment.txt
//         junction_annotation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12 2> ${bam.baseName}.junction_annotation_log.txt
//         bam_stat.py -i $bam 2> ${bam.baseName}.bam_stat.txt
//         junction_saturation.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
//         inner_distance.py -i $bam -o ${bam.baseName}.rseqc -r $bed12
//         read_distribution.py -i $bam -r $bed12 > ${bam.baseName}.read_distribution.txt
//         read_duplication.py -i $bam -o ${bam.baseName}.read_duplication
//         """
//     }
