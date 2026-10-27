process STRINGTIE_STRINGTIE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3facd74a0f728c9bb9e9a731b58c343895d2dbdfeb812ce5747f701103fc61cf/data' :
        'community.wave.seqera.io/library/stringtie:3.0.3--e8043d00caecd051' }"

    input:
    tuple val(meta), path(srbam), path(lrbam)
    val(mode)
    path(annotation_gtf)

    output:
    tuple val(meta), path("${prefix}.transcripts.gtf")   , emit: transcript_gtf
    tuple val(meta), path("${prefix}.gene.abundance.txt"), emit: abundance
    tuple val(meta), path("${prefix}.coverage.gtf")      , optional: true, emit: coverage_gtf
    tuple val(meta), path("${prefix}.ballgown/*.ctab")   , optional: true, emit: ballgown
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args2     = task.ext.args2 ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    def ballgown  = annotation_gtf ? "-b ${prefix}.ballgown" : ""
    def coverage  = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""

    // atleast one bam must be provided
    if (!srbam && !lrbam) {
        error "At least one of srbam or lrbam must be provided for ${meta.id}"
    }

    // check for mode validity and required inputs for each mode
    def run_mode = ''
    if (mode) {
        def valid_modes = ['expression-estimation', 'long-reads-assembly', 'mix-reads-assembly', 'nascent-aware-assembly']
        def modes = (mode instanceof List) ? mode :
                    (mode instanceof String) ? mode.toString().split(',').collect { x -> x.trim() } : []
        modes.each { m ->
            if (!(m in valid_modes)) {
                error "Invalid mode: ${m}. Valid options are: ${valid_modes.join(', ')}"
            }
        }

        // check for required inputs based on modes
        if (modes.contains('mix-reads-assembly') && !(srbam && lrbam)) {
            error "mode 'mix-reads-assembly' requires both srbam and lrbam to be provided for ${meta.id}"
        }
        if (modes.contains('long-reads-assembly') && !lrbam) {
            error "mode 'long-reads-assembly' requires lrbam to be provided for ${meta.id}"
        }
        if (modes.contains('expression-estimation') && !annotation_gtf) {
            error "mode 'expression-estimation' (-e) requires annotation_gtf to be provided for ${meta.id}"
        }

        // add mode flags based on the provided modes
        def mode_flags = []
        if (modes.contains('expression-estimation')) {
            mode_flags += (lrbam && !srbam) ? ['-L', '-e'] : ['-e']
        }
        if (modes.contains('long-reads-assembly') && !modes.contains('expression-estimation')) {
            mode_flags += '-L'
        }
        if (modes.contains('mix-reads-assembly')) {
            mode_flags += '--mix'
        }
        if (modes.contains('nascent-aware-assembly')) {
            mode_flags += ['-N', '--nasc']
        }

        run_mode = mode_flags.join(' ')
    }

    // --mix requires the short-read alignments first, long-read alignments second
    def bam_inputs = (srbam && lrbam) ? "$srbam $lrbam" : (srbam ? "$srbam" : "$lrbam")
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
    prefix = task.ext.prefix ?: "${meta.id}"
    def has_annotation = annotation_gtf ? true : false

    """
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.gene.abundance.txt
    ${has_annotation ? "touch ${prefix}.coverage.gtf" : ''}
    ${has_annotation ? "mkdir -p ${prefix}.ballgown" : ''}
    ${has_annotation ? "touch ${prefix}.ballgown/e_data.ctab ${prefix}.ballgown/i_data.ctab ${prefix}.ballgown/t_data.ctab ${prefix}.ballgown/e2t.ctab ${prefix}.ballgown/i2t.ctab" : ''}
    """
}
