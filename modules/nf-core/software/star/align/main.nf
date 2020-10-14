// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // Don't upgrade me - 2.7X indices incompatible with iGenomes.
    container "quay.io/biocontainers/star:2.6.1d--0"
    //container "https://depot.galaxyproject.org/singularity/star:2.6.1d--0"

    conda (params.conda ? "bioconda::star=2.6.1d" : null)

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    
    output:
    tuple val(meta), path("*Aligned.out.bam") , emit: bam
    tuple val(meta), path("*Log.final.out")   , emit: log_final
    tuple val(meta), path("*Log.out")         , emit: log_out
    tuple val(meta), path("*Log.progress.out"), emit: log_progress
    path  "*.version.txt"                     , emit: version

    tuple val(meta), path("*sortedByCoord.out.bam")  , optional:true, emit: bam_sorted
    tuple val(meta), path("*toTranscriptome.out.bam"), optional:true, emit: bam_transcript
    tuple val(meta), path("*fastq.gz")               , optional:true, emit: fastq
    tuple val(meta), path("*.tab")                   , optional:true, emit: tab

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $ignore_gtf \\
        $seq_center \\
        $options.args

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    STAR --version | sed -e "s/STAR_//g" > ${software}.version.txt
    """
}
