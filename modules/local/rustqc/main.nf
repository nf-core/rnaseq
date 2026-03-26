process RUSTQC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // TODO: pin to a release tag before merge
    container "ghcr.io/seqeralabs/rustqc:dev"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    tuple val(meta), path("rustqc")                                                , emit: results
    tuple val(meta), path("rustqc/rseqc/infer_experiment/*.infer_experiment.txt")  , emit: inferexperiment_txt, optional: true
    tuple val("${task.process}"), val('rustqc'), eval("rustqc --version | sed 's/rustqc //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strandedness = meta.strandedness ?: 'unstranded'
    def paired       = meta.single_end ? '' : '--paired'
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        --stranded ${strandedness} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir rustqc \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rustqc/{dupradar,featurecounts,preseq,samtools} \
            rustqc/rseqc/{bam_stat,infer_experiment,read_duplication,read_distribution,junction_annotation,junction_saturation,inner_distance,tin} \
            rustqc/qualimap/{raw_data_qualimapReport,images_qualimapReport}
    touch rustqc/rseqc/infer_experiment/${prefix}.infer_experiment.txt
    """
}
