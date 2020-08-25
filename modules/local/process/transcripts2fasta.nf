// Import generic module functions
include { saveFiles } from './functions'

process TRANSCRIPTS2FASTA {
    tag "$fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'genome', publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path fasta
    path gtf
    val options

    output:
    path "*.fa", emit: fasta
    
    script: // filter_gtf_for_genes_in_genome.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf_for_genes_in_genome.py --gtf $gtf --fasta $fasta -o ${gtf.baseName}_in_${fasta.baseName}.gtf
    gffread -F -w transcripts.fa -g $fasta ${gtf.baseName}_in_${fasta.baseName}.gtf
    """
}
