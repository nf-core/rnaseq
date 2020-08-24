// Import generic module functions
include { saveFiles } from './functions'

/*
 * Output Markdown documentation to HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:"pipeline_info", publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path output_docs
    path images
    val options

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
