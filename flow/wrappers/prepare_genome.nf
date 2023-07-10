// Test this workflow using:
// nextflow run ./flow/wrappers/prepare_genome.nf -profile docker,test -c ./nextflow.config -c ./flow/conf/test_wrapper.config --outdir ./results

include { PREPARE_GENOME } from '../../subworkflows/local/prepare_genome.nf'

// If additional fasta file is provided biotype value to use when appending entries to GTF file
def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type

workflow {
    // Run wrapper
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        null, // GFF
        params.additional_fasta,
        params.transcript_fasta,
        null, // gene_bed
        null, // splicesites
        null, // bbsplit_fasta_list
        null, // star_index
        null, // rsem_index
        null, // salmon_index
        null, // hisat2_index
        null, // bbsplit_index
        params.gencode,
        false, // is_aws_igenome
        biotype,
        ["star_salmon", "star_rsem", "hisat2", "salmon"]
    )
}
