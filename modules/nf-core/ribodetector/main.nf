nextflow.enable.types = true

record RibodetectorInput {
    id:     String
    meta:   Map
    reads:  List<Path>
    length: Integer
}

process RIBODETECTOR {
	tag "$sample.meta.id"
	label 'process_medium'

	conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
	container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/51100097fc2e31d7b78074bc954774f77726099220e820382e939278817a66da/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/463b8ad941e7f1f2decef20844d666c1c8ac233e166d2bc766164c4a93905a3c/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/ribodetector_pytorch-gpu_cuda-version:fa9183da731515ea' : 'community.wave.seqera.io/library/ribodetector:0.3.3--ad3d7071e408b502') }"

	input:
	sample: RibodetectorInput

	output:
	record(
		id:    sample.id,
		meta:  sample.meta,
		reads: files('*.nonrna*.fastq.gz').toSorted { f -> f.name },
		log:   file('*.log')
	)

	topic:
	tuple(task.process, 'ribodetector', eval('ribodetector --version | sed "s/ribodetector //"')) >> 'versions'
	tuple(task.process, 'cuda', eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"')) >> 'versions'

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${sample.meta.id}"
	ribodetector_bin = task.accelerator ? "ribodetector" : "ribodetector_cpu"
	ribodetector_mem = task.accelerator ? "-m ${task.memory.toGiga()}" : ""
	output = sample.meta.single_end ? "${prefix}.nonrna.fastq.gz" : "${prefix}.nonrna.1.fastq.gz ${prefix}.nonrna.2.fastq.gz"

	"""
	${ribodetector_bin} \\
		-i ${sample.reads} \\
		-o ${output} \\
		-l ${sample.length} \\
		-t ${task.cpus} \\
		--log ${prefix}.log \\
		${ribodetector_mem} \\
		${args}
	"""

	stub:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${sample.meta.id}"

	"""
	echo $args

	echo "" | gzip > ${prefix}.nonrna.1.fastq.gz
	echo "" | gzip > ${prefix}.nonrna.2.fastq.gz
	touch ${prefix}.log
	"""
}
