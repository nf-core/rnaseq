// Import generic module functions
include { saveFiles } from './functions'

process GTF_GENE_FILTER {
    tag "$fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'genome', publish_id:'') }

    container "quay.io/biocontainers/python:3.8.3"
    //container  https://depot.galaxyproject.org/singularity/python:3.8.3

    conda (params.conda ? "conda-forge::python=3.8.3" : null)

    input:
    path fasta
    path gtf
    val  options

    output:
    path "*.gtf"
    
    script: // filter_gtf_for_genes_in_genome.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf_for_genes_in_genome.py --gtf $gtf --fasta $fasta -o ${fasta.baseName}_genes.gtf
    """
}
