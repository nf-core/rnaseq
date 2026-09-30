nextflow.enable.types = true

include { StringtieInput } from '../../types'

process STRINGTIE_STRINGTIE {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3facd74a0f728c9bb9e9a731b58c343895d2dbdfeb812ce5747f701103fc61cf/data' :
        'community.wave.seqera.io/library/stringtie:3.0.3--e8043d00caecd051' }"

    input:
    sample: StringtieInput
    mode: List<String>
    annotation_gtf: Path?

    output:
    record(
        id:             sample.id,
        meta:           sample.meta,
        transcript_gtf: file("${prefix}.transcripts.gtf"),
        abundance:      file("${prefix}.gene.abundance.txt"),
        coverage_gtf:   file("${prefix}.coverage.gtf", optional: true),
        ballgown:       files("${prefix}.ballgown/*.ctab", optional: true)
    )

    topic:
    tuple(task.process, 'stringtie', eval('stringtie --version')) >> 'versions'

    script:
    def args      = task.ext.args ?: ''
    def args2     = task.ext.args2 ?: ''
    prefix        = task.ext.prefix ?: "${sample.meta.id}"
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    def ballgown  = annotation_gtf ? "-b ${prefix}.ballgown" : ""
    def coverage  = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""

    // atleast one bam must be provided
    if (!sample.bam && !sample.lrbam) {
        error "At least one of bam or lrbam must be provided for ${sample.meta.id}"
    }

    // check for mode validity and required inputs for each mode
    def run_mode = ''
    if (mode) {
        def valid_modes = ['expression-estimation', 'long-reads-assembly', 'mix-reads-assembly', 'nascent-aware-assembly']
        def modes = mode
        modes.each { m ->
            if (!(m in valid_modes)) {
                error "Invalid mode: ${m}. Valid options are: ${valid_modes.join(', ')}"
            }
        }

        // check for required inputs based on modes
        if (modes.contains('mix-reads-assembly') && !(sample.bam && sample.lrbam)) {
            error "mode 'mix-reads-assembly' requires both bam and lrbam to be provided for ${sample.meta.id}"
        }
        if (modes.contains('long-reads-assembly') && !sample.lrbam) {
            error "mode 'long-reads-assembly' requires lrbam to be provided for ${sample.meta.id}"
        }
        if (modes.contains('expression-estimation') && !annotation_gtf) {
            error "mode 'expression-estimation' (-e) requires annotation_gtf to be provided for ${sample.meta.id}"
        }

        // add mode flags based on the provided modes
        def mode_flags = []
        if (modes.contains('expression-estimation')) {
            mode_flags += (sample.lrbam && !sample.bam) ? ['-L', '-e'] : ['-e']
        }
        if (modes.contains('long-reads-assembly') && !modes.contains('expression-estimation')) {
            mode_flags += ['-L']
        }
        if (modes.contains('mix-reads-assembly')) {
            mode_flags += ['--mix']
        }
        if (modes.contains('nascent-aware-assembly')) {
            mode_flags += ['-N', '--nasc']
        }

        run_mode = mode_flags.join(' ')
    }

    // --mix requires the short-read alignments first, long-read alignments second
    def bam_inputs = (sample.bam && sample.lrbam) ? "${sample.bam} ${sample.lrbam}" : (sample.bam ? "${sample.bam}" : "${sample.lrbam}")
    """
    stringtie \\
        -o ${prefix}.transcripts.gtf \\
        ${run_mode} \\
        ${args2} \\
        ${reference} \\
        -A ${prefix}.gene.abundance.txt \\
        ${coverage} \\
        ${ballgown} \\
        -p ${task.cpus} \\
        ${bam_inputs} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def has_annotation = annotation_gtf ? true : false

    """
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.gene.abundance.txt
    ${has_annotation ? "touch ${prefix}.coverage.gtf" : ''}
    ${has_annotation ? "mkdir -p ${prefix}.ballgown" : ''}
    ${has_annotation ? "touch ${prefix}.ballgown/e_data.ctab ${prefix}.ballgown/i_data.ctab ${prefix}.ballgown/t_data.ctab ${prefix}.ballgown/e2t.ctab ${prefix}.ballgown/i2t.ctab" : ''}
    """
}
