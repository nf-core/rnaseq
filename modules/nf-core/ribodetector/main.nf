process RIBODETECTOR {
	tag "$meta.id"
	label 'process_medium'

	conda "${moduleDir}/environment.yml"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4de8fe74d21198e6fc8218cb3209d929b3d7dab750678501b096b0ccc324307b/data' :
        'community.wave.seqera.io/library/ribodetector:0.3.2--cbe1c77fa14eeb53' }"

	input:
	tuple val(meta), path(fastq)
	val length

	output:
	tuple val(meta), path("*.nonrna*.fastq.gz"), emit: fastq
	tuple val(meta), path("*.log")             , emit: log
	tuple val("${task.process}"), val('ribodetector'), eval('ribodetector --version | sed "s/ribodetector //"'), emit: versions_ribodetector, topic: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"
	ribodetector_bin = task.accelerator ? "ribodetector" : "ribodetector_cpu"
	ribodetector_mem = task.accelerator ? "-m $task.memory.toGiga()" : ""
	output = meta.single_end ? "${prefix}.nonrna.fastq.gz" : "${prefix}.nonrna.1.fastq.gz ${prefix}.nonrna.2.fastq.gz"

	"""
	${ribodetector_bin} \\
		-i ${fastq} \\
		-o ${output} \\
		-l ${length} \\
		-t ${task.cpus} \\
		--log ${prefix}.log \\
		${ribodetector_mem} \\
		${args}
	"""

	stub:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"

	"""
	echo $args

	echo | gzip > ${prefix}.nonrna.1.fastq.gz
	echo | gzip > ${prefix}.nonrna.2.fastq.gz
	touch ${prefix}.log
	"""
}
