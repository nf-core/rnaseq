nextflow.enable.types = true

record CustomRsemmergecountsInput {
    id:       String
    meta:     Map
    genes:    List<Path>
    isoforms: List<Path>
}

record RsemMerge {
    counts_gene:       Path
    tpm_gene:          Path
    counts_transcript: Path
    tpm_transcript:    Path
    genes_long:        Path
    isoforms_long:     Path
}

// The single 'all_samples' row of merged RSEM tables; CUSTOM_RSEMMERGECOUNTS output without meta
record RsemMergeSample {
    id:         String
    rsem_merge: RsemMerge
}

process CUSTOM_RSEMMERGECOUNTS {
    tag "${sample.meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    sample: CustomRsemmergecountsInput

    stage:
    stageAs sample.genes, 'genes/*'
    stageAs sample.isoforms, 'isoforms/*'

    output:
    record(
        id: sample.id,
        meta: sample.meta,
        rsem_merge: record(
            counts_gene:       file("${prefix}.gene_counts.tsv"),
            tpm_gene:          file("${prefix}.gene_tpm.tsv"),
            counts_transcript: file("${prefix}.transcript_counts.tsv"),
            tpm_transcript:    file("${prefix}.transcript_tpm.tsv"),
            genes_long:        file("${prefix}.genes_long.tsv"),
            isoforms_long:     file("${prefix}.isoforms_long.tsv")
        )
    )

    topic:
    tuple(task.process, 'sed', eval("sed --version 2>&1 | sed '1!d;s/^.*) //'")) >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
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

    paste gene_ids.txt tmp/genes/*.counts.txt > ${prefix}.gene_counts.tsv
    paste gene_ids.txt tmp/genes/*.tpm.txt > ${prefix}.gene_tpm.tsv
    paste transcript_ids.txt tmp/isoforms/*.counts.txt > ${prefix}.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms/*.tpm.txt > ${prefix}.transcript_tpm.tsv

    # Create long format for genes (idx=1-4, concat columns 5-7)
    echo -e "sample_name\\tgene_id\\ttranscript_id(s)\\tlength\\teffective_length\\texpected_count\\tTPM\\tFPKM" > ${prefix}.genes_long.tsv
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        tail -n+2 \$fileid | awk -v sample=\$samplename 'BEGIN{OFS="\\t"}{print sample,\$1,\$2,\$3,\$4,\$5,\$6,\$7}' >> ${prefix}.genes_long.tsv
    done

    # Create long format for isoforms (idx=1-4, concat columns 5-8)
    echo -e "sample_name\\ttranscript_id\\tgene_id\\tlength\\teffective_length\\texpected_count\\tTPM\\tFPKM\\tIsoPct" > ${prefix}.isoforms_long.tsv
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        tail -n+2 \$fileid | awk -v sample=\$samplename 'BEGIN{OFS="\\t"}{print sample,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' >> ${prefix}.isoforms_long.tsv
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.gene_counts.tsv
    touch ${prefix}.gene_tpm.tsv
    touch ${prefix}.transcript_counts.tsv
    touch ${prefix}.transcript_tpm.tsv
    touch ${prefix}.genes_long.tsv
    touch ${prefix}.isoforms_long.tsv
    """
}
