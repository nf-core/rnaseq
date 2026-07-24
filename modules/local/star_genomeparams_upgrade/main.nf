process STAR_GENOMEPARAMS_UPGRADE {
    tag "${meta.id ?: index.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(index, stageAs: 'input_index')

    output:
    tuple val(meta), path('star'), emit: index
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p star
    for f in input_index/*; do
        name=\$(basename "\$f")
        if [ "\$name" = "genomeParameters.txt" ]; then
            continue
        fi
        ln -s "\$(readlink -f "\$f")" "star/\$name"
    done

    awk -F'\\t' -v OFS='\\t' '
        \$1 == "versionGenome" && \$2 == "20201" {
            print "versionGenome", "2.7.4a"
            seen_upgraded = 1
            next
        }
        \$1 == "genomeType"          { seen_genomeType = 1 }
        \$1 == "genomeTransformType" { seen_transformType = 1 }
        \$1 == "genomeTransformVCF"  { seen_transformVCF = 1 }
        { print }
        END {
            if (seen_upgraded) {
                if (!seen_genomeType)    print "genomeType", "Full"
                if (!seen_transformType) print "genomeTransformType", "None"
                if (!seen_transformVCF)  print "genomeTransformVCF", "-"
            }
        }
    ' "input_index/genomeParameters.txt" > star/genomeParameters.txt
    """

    stub:
    """
    mkdir -p star
    for f in input_index/*; do
        ln -s "\$(readlink -f "\$f")" "star/\$(basename "\$f")"
    done
    """
}
