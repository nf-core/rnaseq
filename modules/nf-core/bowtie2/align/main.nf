process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val   save_unaligned
    val   sort_bam

    output:
    tuple val(meta), path("*.sam")      , emit: sam     , optional:true
    tuple val(meta), path("*.bam")      , emit: bam     , optional:true
    tuple val(meta), path("*.cram")     , emit: cram    , optional:true
    tuple val(meta), path("*.csi")      , emit: csi     , optional:true
    tuple val(meta), path("*.crai")     , emit: crai    , optional:true
    tuple val(meta), path("*.log")      , emit: log
    tuple val(meta), path("*fastq.gz")  , emit: fastq   , optional:true
    tuple val("${task.process}"), val('bowtie2'), eval("bowtie2 --version 2>&1 | sed -n '1s/.*bowtie2-align-s version //p'"), emit: versions_bowtie2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //'"), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        -x \$INDEX \\
        $reads_args \\
        --threads $task.cpus \\
        $unaligned \\
        $rg \\
        $args \\
        2>| >(tee ${prefix}.bowtie2.log >&2) \\
        | samtools $samtools_command $args2 --threads $task.cpus ${reference} -o ${prefix}.${extension} -

    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi

    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi
    """

    stub:
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"
    def create_unmapped = ""
    if (meta.single_end) {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped.fastq.gz" : ""
    } else {
        create_unmapped = save_unaligned ? "touch ${prefix}.unmapped_1.fastq.gz && touch ${prefix}.unmapped_2.fastq.gz" : ""
    }
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }

    """
    touch ${prefix}.${extension}
    ${create_index}
    touch ${prefix}.bowtie2.log
    ${create_unmapped}
    """

}
