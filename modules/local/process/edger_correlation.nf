// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process EDGER_CORRELATION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path counts
    path mdsplot_header
    path heatmap_header
    
    output:
    path "*.pdf"         , optional:true, emit: pdf
    path "*_matrix.csv"  , optional:true, emit: matrix
    path "*_mqc.csv"     , optional:true, emit: multiqc
    path  "*.version.txt", emit: version

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def software = getSoftwareName(task.process)
    if (counts.size() > 2) {
        """
        edger_heatmap_mds.r $counts
        cat $mdsplot_header edger_mds_aplot_coordinates_mqc.csv >> tmp_file
        mv tmp_file edger_mds_aplot_coordinates_mqc.csv
        cat $heatmap_header log2cpm_sample_correlation_mqc.csv >> tmp_file
        mv tmp_file log2cpm_sample_correlation_mqc.csv

        Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='${software}.version.txt')"
        """
    } else {
        """
        Rscript -e "library(edgeR); write(x=as.character(packageVersion('edgeR')), file='${software}.version.txt')"
        """
    }
}
