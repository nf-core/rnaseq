process RSEM_MERGE_COUNTS {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path ('genes/*')
    path ('isoforms/*')

    output:
    path "rsem.merged.gene_counts.tsv"      , emit: counts_gene
    path "rsem.merged.gene_tpm.tsv"         , emit: tpm_gene
    path "rsem.merged.transcript_counts.tsv", emit: counts_transcript
    path "rsem.merged.transcript_tpm.tsv"   , emit: tpm_transcript
    path "rsem.merged.genes_long.tsv"       , emit: genes_long
    path "rsem.merged.isoforms_long.tsv"    , emit: isoforms_long
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tmp/genes
    cut -f 1,2 `ls ./genes/* | head -n 1` > gene_ids.txt
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        echo \$samplename > tmp/genes/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.counts.txt
        echo \$samplename > tmp/genes/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.tpm.txt
    done

    mkdir -p tmp/isoforms
    cut -f 1,2 `ls ./isoforms/* | head -n 1` > transcript_ids.txt
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        echo \$samplename > tmp/isoforms/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.counts.txt
        echo \$samplename > tmp/isoforms/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.tpm.txt
    done

    paste gene_ids.txt tmp/genes/*.counts.txt > rsem.merged.gene_counts.tsv
    paste gene_ids.txt tmp/genes/*.tpm.txt > rsem.merged.gene_tpm.tsv
    paste transcript_ids.txt tmp/isoforms/*.counts.txt > rsem.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms/*.tpm.txt > rsem.merged.transcript_tpm.tsv

    # Create long format for genes (idx=1-4, concat columns 5-7)
    echo -e "sample_name\tgene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM" > rsem.merged.genes_long.tsv
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        tail -n+2 \$fileid | awk -v sample=\$samplename 'BEGIN{OFS="\t"}{print sample,\$1,\$2,\$3,\$4,\$5,\$6,\$7}' >> rsem.merged.genes_long.tsv
    done

    # Create long format for isoforms (idx=1-4, concat columns 5-8)
    echo -e "sample_name\ttranscript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct" > rsem.merged.isoforms_long.tsv
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        tail -n+2 \$fileid | awk -v sample=\$samplename 'BEGIN{OFS="\t"}{print sample,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' >> rsem.merged.isoforms_long.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch rsem.merged.gene_counts.tsv
    touch rsem.merged.gene_tpm.tsv
    touch rsem.merged.transcript_counts.tsv
    touch rsem.merged.transcript_tpm.tsv
    touch rsem.merged.genes_long.tsv
    touch rsem.merged.isoforms_long.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
