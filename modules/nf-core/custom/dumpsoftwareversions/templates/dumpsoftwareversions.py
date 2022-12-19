#!/usr/bin/env python


"""Provide functions to merge multiple versions.yml files."""


import platform
from textwrap import dedent

import yaml


def _make_versions_html(versions):
    """Generate a tabular HTML output of all versions for MultiQC."""
    html = [
        dedent(
            """\\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(
                    f"""\\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\\n".join(html)


def main():
    """Load all version files and generate merged output."""
    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }

    with open("$versions") as f:
        versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) | versions_this_module

    # aggregate versions by the module name (derived from fully-qualified process name)
    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        try:
            if versions_by_module[module] != process_versions:
                raise AssertionError(
                    "We assume that software versions are the same between all modules. "
                    "If you see this error-message it means you discovered an edge-case "
                    "and should open an issue in nf-core/tools. "
                )
        except KeyError:
            versions_by_module[module] = process_versions

    versions_by_module["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version",
    }

    versions_mqc = {
        "id": "software_versions",
        "section_name": "${workflow.manifest.name} Software Versions",
        "section_href": "https://github.com/${workflow.manifest.name}",
        "plot_type": "html",
        "description": "are collected at run time from the software output.",
        "data": _make_versions_html(versions_by_module),
    }

    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    with open("software_versions_mqc.yml", "w") as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)


if __name__ == "__main__":
    main()
