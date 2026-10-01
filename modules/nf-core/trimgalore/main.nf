nextflow.enable.types = true


record TrimgaloreInput {
    id:    String
    meta:  Map
    reads: List<Path>
    args:  String?
}

record TrimgaloreResult {
    id:       String
    meta:     Map
    reads:    List<Path>
    log:      List<Path>
    json:     List<Path>
    unpaired: List<Path>
    html:     List<Path>
    zip:      List<Path>
}

process TRIMGALORE {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e00369598bd6b7b34a7c83d5496c381104bf8b885c31a4b65b92e6ea2059fbb3/data' :
        'community.wave.seqera.io/library/trim-galore:2.3.0--6a38a479b4972363'}"

    input:
    sample: TrimgaloreInput

    output:
    record(
        id:       sample.id,
        meta:     sample.meta,
        reads:    files("*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz").toSorted { f -> f.name },
        log:      files("*report.txt", optional: true).toSorted { f -> f.name },
        json:     files("*report.json", optional: true).toSorted { f -> f.name },
        unpaired: files("*unpaired{,_1,_2}.fq.gz", optional: true).toSorted { f -> f.name },
        html:     files("*.html", optional: true).toSorted { f -> f.name },
        zip:      files("*.zip", optional: true).toSorted { f -> f.name }
    ) as TrimgaloreResult

    topic:
    tuple(task.process, "trimgalore", eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"')) >> 'versions'

    script:
    // sample.args replaces the options of task.ext.args that it repeats
    def base_args  = task.ext.args ? task.ext.args.replaceAll('\\s(?=--)', '\u0001').tokenize('\u0001') : []
    def extra_args = sample.args   ? sample.args.replaceAll('\\s(?=--)', '\u0001').tokenize('\u0001') : []
    def extra_options = extra_args.collect { a -> a.trim().tokenize(' ')[0] }
    def kept_args = base_args.findAll { a -> !extra_options.contains(a.trim().tokenize(' ')[0]) }.toList()
    def args = (kept_args + extra_args).collect { a -> a.trim() }.join(' ')
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (sample.meta.single_end) {
            cores = (task.cpus as int) - 3
        }
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    if (sample.meta.single_end) {
        def args_list = args.replaceAll('\\s(?=--)', '\u0001').tokenize('\u0001').findAll { arg -> !arg.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s ${sample.reads[0]} ${prefix}.fastq.gz
        trim_galore \\
            ${args_list.join(' ')} \\
            --cores ${cores} \\
            --gzip \\
            ${prefix}.fastq.gz
        """
    }
    else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${sample.reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${sample.reads[1]} ${prefix}_2.fastq.gz
        trim_galore \\
            ${args} \\
            --cores ${cores} \\
            --paired \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    if (sample.meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}_trimmed.fq.gz ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.txt ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.json"
    }
    else {
        output_command = "echo '' | gzip > ${prefix}_1_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.json ;"
        output_command += "echo '' | gzip > ${prefix}_2_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.json"
    }
    """
    ${output_command}
    """
}
