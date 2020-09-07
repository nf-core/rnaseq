// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process EDGER_CORRELATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path counts
    path mdsplot_header
    path heatmap_header
    val  options
    // path input_files from geneCounts.collect()
    // val num_bams from bam_count.count()

    //output:
    //tuple val(meta), path("*.pdf")    , emit: pdf
    //tuple val(meta), path("*.txt")    , emit: txt
    //tuple val(meta), path("*_mqc.txt"), emit: multiqc
    //path  "*.version.txt"             , emit: version
    //         path "*.{txt,pdf,csv}" into sample_correlation_results

    //         when:
    //         num_bams > 2 && (!params.sample_level)

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def software = getSoftwareName(task.process)
    """
    edgeR_heatmap_MDS.r $counts
    cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
    mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
    cat $heatmap_header log2CPM_sample_correlation_mqc.csv >> tmp_file
    mv tmp_file log2CPM_sample_correlation_mqc.csv

    Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='${software}.version.txt')"
    """
}
