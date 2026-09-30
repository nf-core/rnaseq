nextflow.enable.types = true

include { FastpReads } from '../types'

record FastpResult {
    id:           String
    meta:         Map
    reads:        List<Path>
    json:         Path
    html:         Path
    log:          Path
    reads_fail:   List<Path>
    reads_merged: Path?
}

process FASTP {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013aad5427d824afe472e6607ea47685ff0181f1fb09e52a179e0ec39e43e88/data'
:         'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde' }"

    input:
    sample: FastpReads
    discard_trimmed_pass: Boolean
    save_trimmed_fail: Boolean
    save_merged: Boolean

    output:
    record(
        id:           sample.id,
        meta:         sample.meta,
        reads:        files('*.fastp.fastq.gz', optional: true).toSorted { f -> f.name },
        json:         file('*.json'),
        html:         file('*.html'),
        log:          file('*.log'),
        reads_fail:   files('*.fail.fastq.gz', optional: true).toSorted { f -> f.name },
        reads_merged: file('*.merged.fastq.gz', optional: true)
    ) as FastpResult

    topic:
    tuple(task.process, 'fastp', eval('fastp --version 2>&1 | sed -e "s/fastp //g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def adapter_list = sample.adapter_fasta ? "--adapter_fasta ${sample.adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail && sample.meta.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !sample.meta.single_end ? "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_R1.fail.fastq.gz --unpaired2 ${prefix}_R2.fail.fastq.gz" : ''
    def out_fq1 = discard_trimmed_pass ? 'true' : ( sample.meta.single_end ? "--out1 ${prefix}.fastp.fastq.gz" : "--out1 ${prefix}_R1.fastp.fastq.gz" )
    def out_fq2 = discard_trimmed_pass ? 'true' : "--out2 ${prefix}_R2.fastp.fastq.gz"
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf ${sample.reads[0]} ${prefix}.fastq.gz

        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2) \\
        | gzip -c > ${prefix}.fastp.fastq.gz
        """
    } else if (sample.meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf ${sample.reads[0]} ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            $out_fq1 \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_R1.fastq.gz ] && ln -sf ${sample.reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f  ${prefix}_R2.fastq.gz ] && ln -sf ${sample.reads[1]} ${prefix}_R2.fastq.gz
        fastp \\
            --in1 ${prefix}_R1.fastq.gz \\
            --in2 ${prefix}_R2.fastq.gz \\
            $out_fq1 \\
            $out_fq2 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    }

    stub:
    def prefix              = task.ext.prefix ?: "${sample.meta.id}"
    def is_single_output    = task.ext.args?.contains('--interleaved_in') || sample.meta.single_end
    def touch_reads         = (discard_trimmed_pass) ? "" : (is_single_output) ? "echo '' | gzip > ${prefix}.fastp.fastq.gz" : "echo '' | gzip > ${prefix}_R1.fastp.fastq.gz ; echo '' | gzip > ${prefix}_R2.fastp.fastq.gz"
    def touch_merged        = (!is_single_output && save_merged) ? "echo '' | gzip >  ${prefix}.merged.fastq.gz" : ""
    def touch_fail_fastq    = (!save_trimmed_fail) ? "" : sample.meta.single_end ? "echo '' | gzip > ${prefix}.fail.fastq.gz" : "echo '' | gzip > ${prefix}.paired.fail.fastq.gz ; echo '' | gzip > ${prefix}_R1.fail.fastq.gz ; echo '' | gzip > ${prefix}_R2.fail.fastq.gz"
    """
    $touch_reads
    $touch_fail_fastq
    $touch_merged
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"
    """
}
