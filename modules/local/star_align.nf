// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.6.1d' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/star:2.6.1d--0'
    } else {
        container 'quay.io/biocontainers/star:2.6.1d--0'
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path "versions.yml"                       , emit: version

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    script:
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def out_sam_type = (options.args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (options.args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $seq_center \\
        $options.args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
