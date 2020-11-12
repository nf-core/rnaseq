// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process GTF_GENE_FILTER {
    tag "$fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path fasta
    path gtf
    
    output:
    path "*.gtf"
    
    script: // filter_gtf_for_genes_in_genome.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf_for_genes_in_genome.py --gtf $gtf --fasta $fasta -o ${fasta.baseName}_genes.gtf
    """
}
