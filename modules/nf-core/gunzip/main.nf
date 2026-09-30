nextflow.enable.types = true

include { ArchiveInput } from '../types'

record GunzipResult {
    id:   String
    meta: Map
    file: Path
}

process GUNZIP {
    tag "${sample.archive}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    sample: ArchiveInput

    output:
    record(id: sample.id, meta: sample.meta, file: file("${gunzip}")) as GunzipResult

    topic:
    tuple(task.process, 'gunzip', eval('gunzip --version 2>&1 | head -1 | sed "s/^.*(gzip) //; s/ Copyright.*//"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def nameWithoutGz = sample.archive.extension == 'gz' ? sample.archive.baseName : sample.archive.name
	def extension = file(nameWithoutGz).extension
	def name = file(nameWithoutGz).baseName
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"
    """
    # Not calling gunzip itself because it creates files
    # with the original group ownership rather than the
    # default one for that user / the work directory
    gzip \\
        -cd \\
        ${args} \\
        ${sample.archive} \\
        > ${gunzip}
    """

    stub:
    def nameWithoutGz = sample.archive.extension == 'gz' ? sample.archive.baseName : sample.archive.name
	def extension = file(nameWithoutGz).extension
	def name = file(nameWithoutGz).baseName
    def prefix = task.ext.prefix ?: name
    gunzip = prefix + ".${extension}"
    """
    touch ${gunzip}
    """
}
