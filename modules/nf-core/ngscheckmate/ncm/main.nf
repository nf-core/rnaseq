process NGSCHECKMATE_NCM {
    label 'process_low'

    conda "bioconda::ngscheckmate=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.0--py27r41hdfd78af_3':
        'quay.io/biocontainers/ngscheckmate:1.0.0--py27r41hdfd78af_3' }"

    input:
    path files
    path snp_bed
    path fasta

    output:
    path "*.pdf"            , emit: pdf, optional: true
    path "*_corr_matrix.txt", emit: corr_matrix
    path "*_matched.txt"    , emit: matched
    path "*_all.txt"        , emit: all
    path "*.vcf"            , emit: vcf, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "output"
    def unzip = files.any { it.toString().endsWith(".vcf.gz") }
    """
    if $unzip
    then
        for VCFGZ in *.vcf.gz; do
            gunzip -cdf \$VCFGZ > \$( basename \$VCFGZ .gz );
        done
    fi

    NCM_REF="./"${fasta} ncm_edited.py -d . -bed ${snp_bed} -O . -N ${prefix} $args

    if $unzip
    then
        rm -f *.vcf  # clean up decompressed vcfs
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}
