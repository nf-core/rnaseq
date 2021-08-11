// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]

process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pyaml==15.8.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pyaml%3A15.8.2--py36_0"
    } else {
        container "quay.io/biocontainers/pyaml:15.8.2--py36_0"
    }

    cache false

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path 'software_versions_mqc.yml', emit: mqc_yaml

    script:
    """
    #!/usr/bin/env python

    import yaml
    from textwrap import dedent

    def _make_versions_html(versions):
        html = ["""
        <style>
        tbody:nth-child(even) {
            background-color: #f2f2f2;;
        }
        </style>
        <table style="width:100%">
            <thead>
                <tr>
                    <th> Process Name </th>
                    <th> Software </th>
                    <th> Version  </th>
                </tr>
            </thead>
        """]
        for process, tmp_versions in versions.items():
            html.append("<tbody>")
            print_module = 0
            for tool, version in tmp_versions.items():
                print_module += 1
                html.append(
                    dedent(
                        f'''\
                        <tr>
                            <td><samp>{process if (print_module == 1) else ''}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "\\n".join(html)

    with open("$versions") as f:
        versions = yaml.safe_load(f)

    versions["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version"
    }

    versions_mqc = {
        'id': 'software_versions',
        'section_name': '${workflow.manifest.name} Software Versions',
        'section_href': 'https://github.com/${workflow.manifest.name}',
        'plot_type': 'html',
        'description': 'are collected at run time from the software output.',
        'data': _make_versions_html(versions)
    }

    with open("software_versions.yml", 'w') as f:
        yaml.dump(versions, f, default_flow_style=False)
    with open("software_versions_mqc.yml", 'w') as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)
    """
}
