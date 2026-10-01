nextflow.enable.types = true

include { Bowtie2Logs } from '../../types'

record Bowtie2AlignInput {
    id:    String
    meta:  Map
    reads: List<Path>
    args:  String?
}

record Bowtie2AlignResult {
    id:       String
    meta:     Map
    sam:      Path?
    raw_bams: List<Path>
    cram:     Path?
    csi:      Path?
    crai:     Path?
    unmapped: List<Path>
    bowtie2:  Bowtie2Logs
}

process BOWTIE2_ALIGN {
    tag "${sample.meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    sample: Bowtie2AlignInput
    index: Path
    fasta: Path?
    save_unaligned: Boolean
    sort_bam: Boolean

    output:
    record(
        id:       sample.id,
        meta:     sample.meta,
        sam:      file('*.sam',  optional: true),
        raw_bams: files('*.bam', optional: true).toSorted { f -> f.name },
        cram:     file('*.cram', optional: true),
        csi:      file('*.csi',  optional: true),
        crai:     file('*.crai', optional: true),
        unmapped: files('*fastq.gz', optional: true).toSorted { f -> f.name },
        bowtie2:  record(log: file('*.log'))
    ) as Bowtie2AlignResult

    topic:
    tuple(task.process, 'bowtie2', eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'")) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'
    tuple(task.process, 'pigz', eval("pigz --version 2>&1 | sed 's/pigz //'")) >> 'versions'

    script:
    def args = task.ext.args ?: sample.args ?: ''
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"

    def unaligned = ""
    def reads_args = ""
    if (sample.meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${sample.reads.join(' ')}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${sample.reads[0]} -2 ${sample.reads[1]}"
    }

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? (extension_matcher[0][2] as String).toLowerCase() : "bam"
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
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? ((args2 =~ extension_pattern)[0][2] as String).toLowerCase() : "bam"
    def create_unmapped = ""
    if (sample.meta.single_end) {
        create_unmapped = save_unaligned ? "echo | gzip > ${prefix}.unmapped.fastq.gz" : ""
    } else {
        create_unmapped = save_unaligned ? "echo | gzip > ${prefix}.unmapped_1.fastq.gz && echo | gzip > ${prefix}.unmapped_2.fastq.gz" : ""
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
