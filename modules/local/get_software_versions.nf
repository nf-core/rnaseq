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
        html = ["<dl>"]
        for process, tmp_versions in versions.items():
            for tool, version in tmp_versions.items():
                html.append(
                    dedent(
                        f'''\\
                        <dl>
                            <dt>{process}</dt>
                            <dt>{tool}</dt>
                            <dd>{version}</dd>
                        </dl>
                        '''
                    )
                )
        html.append("</dl>")
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
