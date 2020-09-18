// Import generic module functions
include { saveFiles } from './functions'

/*
 * Concatenate additional fasta file e.g. ERCC spike-ins, GTF etc to primary fasta and gtf genome annotation
 */
process CAT_ADDITIONAL_FASTA {
    tag "$add_fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'genome', publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path fasta
    path gtf
    path add_fasta
    val  options

    output:
    path "${name}.fasta", emit: fasta
    path "${name}.gtf"  , emit: gtf

    script:
    def genome_name = params.genome ? params.genome : fasta.getBaseName()
    def add_name    = add_fasta.getBaseName()
    name            = "${genome_name}_${add_name}"
    """
    fasta2gtf.py -o ${add_fasta.baseName}.gtf $add_fasta
    cat $fasta $add_fasta > ${name}.fasta
    cat $gtf ${add_fasta.baseName}.gtf > ${name}.gtf
    """
}
