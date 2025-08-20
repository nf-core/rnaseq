process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${
        def use_fast_mode = task.ext.use_fast_dupradar ?: false
        use_fast_mode ? "${moduleDir}/environment_fast.yml" : "${moduleDir}/environment.yml"
    }"
    container "${
        def use_fast_mode = task.ext.use_fast_dupradar ?: false
        // For fast mode, use pre-built Wave container with featureCounts
        use_fast_mode ? 'community.wave.seqera.io/library/subread_python_matplotlib_numpy_pruned:70e1e046570ad3a3' : (
            workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/bioconductor-dupradar:1.32.0--r43hdfd78af_0' :
            'biocontainers/bioconductor-dupradar:1.32.0--r43hdfd78af_0'
        )
    }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*_duprateExpDens.pdf")   , emit: scatter2d
    tuple val(meta), path("*_duprateExpBoxplot.pdf"), emit: boxplot
    tuple val(meta), path("*_expressionHist.pdf")   , emit: hist
    tuple val(meta), path("*_dupMatrix.txt")        , emit: dupmatrix
    tuple val(meta), path("*_intercept_slope.txt")  , emit: intercept_slope
    tuple val(meta), path("*_mqc.txt")              , emit: multiqc
    tuple val(meta), path("*.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def use_fast_mode = task.ext.use_fast_dupradar ?: false
    if (!use_fast_mode) {
        template 'dupradar.r'
    } else {
        """
        # Parse inputs from Nextflow template variables
        input_bam='$bam'
        output_prefix='$meta.id'
        if [[ '$task.ext.prefix' != 'null' ]]; then
            output_prefix='$task.ext.prefix'
        fi
        annotation_gtf='$gtf'
        threads=$task.cpus

        # Parse strandedness
        stranded=0
        if [[ '${meta.strandedness}' == 'forward' ]]; then
            stranded=1
        elif [[ '${meta.strandedness}' == 'reverse' ]]; then
            stranded=2
        fi

        # Parse paired-end
        if [[ '${meta.single_end}' != 'true' ]]; then
            paired_flag="-p"
        else
            paired_flag=""
        fi

        # Feature type (default: exon)
        feature_type="exon"
        # Parse task.ext.args for feature_type if provided
        if [[ '$task.ext.args' == *"feature_type"* ]]; then
            feature_type=\$(echo '$task.ext.args' | grep -o 'feature_type[[:space:]]*[^[:space:]]*' | cut -d' ' -f2 || echo "exon")
        fi

        echo "=== Fast dupRadar Analysis ==="
        echo "Input BAM: \$input_bam"
        echo "Input GTF: \$annotation_gtf"
        echo "Strandness: \$(echo 'unstranded forward reverse' | cut -d' ' -f\$((stranded+1)))"
        echo "Library type: \$(if [[ '${meta.single_end}' == 'true' ]]; then echo 'single-end'; else echo 'paired-end'; fi)"
        echo "Feature type: \$feature_type"
        echo "Threads: \$threads"
        echo "Output prefix: \$output_prefix"

        # Step 1: Run featureCounts with duplicates included (no --ignoreDup flag)
        echo "Running featureCounts with duplicates..."
        if [[ -n "\$paired_flag" ]]; then
            featureCounts \\
                -T \$threads \\
                -s \$stranded \\
                \$paired_flag \\
                -t \$feature_type \\
                -g gene_id \\
                -a \$annotation_gtf \\
                -o "\${output_prefix}_with_dups.txt" \\
                \$input_bam
        else
            featureCounts \\
                -T \$threads \\
                -s \$stranded \\
                -t \$feature_type \\
                -g gene_id \\
                -a \$annotation_gtf \\
                -o "\${output_prefix}_with_dups.txt" \\
                \$input_bam
        fi

        # Step 2: Run featureCounts without duplicates (with --ignoreDup flag)
        echo "Running featureCounts without duplicates..."
        if [[ -n "\$paired_flag" ]]; then
            featureCounts \\
                -T \$threads \\
                -s \$stranded \\
                \$paired_flag \\
                -t \$feature_type \\
                -g gene_id \\
                -a \$annotation_gtf \\
                -o "\${output_prefix}_no_dups.txt" \\
                --ignoreDup \\
                \$input_bam
        else
            featureCounts \\
                -T \$threads \\
                -s \$stranded \\
                -t \$feature_type \\
                -g gene_id \\
                -a \$annotation_gtf \\
                -o "\${output_prefix}_no_dups.txt" \\
                --ignoreDup \\
                \$input_bam
        fi

        # Step 3: Process results and calculate duplication rates using Python
        echo "Running dupRadar analysis..."
        dupradar_fast_analysis.py \\
            --with-dups "\${output_prefix}_with_dups.txt" \\
            --no-dups "\${output_prefix}_no_dups.txt" \\
            --prefix "\$output_prefix"

        # Step 4: Create session info log (simplified)
        cat > "\${output_prefix}.R_sessionInfo.log" << EOF
Fast dupRadar replacement (Python implementation)
Date: \$(date)
featureCounts version: \$(featureCounts -v 2>&1 | head -1)
Python version: \$(python3 --version)
System: \$(uname -a)

This analysis was performed using native featureCounts calls instead of R/dupRadar
to improve performance with Fusion S3 filesystem while maintaining MultiQC compatibility.
EOF

        # Step 5: Create versions.yml file
        cat > versions.yml << EOF
"${task.process}":
    subread: \$(featureCounts -v 2>&1 | head -1 | sed 's/featureCounts //' | sed 's/ .*//')
    python: \$(python3 --version | sed 's/Python //')
EOF

        echo "=== Fast dupRadar Analysis Complete ==="
        """
    }

    stub:
    def use_fast_mode = task.ext.use_fast_dupradar ?: false
    if (use_fast_mode) {
        """
        touch ${meta.id}_duprateExpDens.pdf
        touch ${meta.id}_duprateExpBoxplot.pdf
        touch ${meta.id}_expressionHist.pdf
        touch ${meta.id}_dupMatrix.txt
        touch ${meta.id}_intercept_slope.txt
        touch ${meta.id}_dup_intercept_mqc.txt
        touch ${meta.id}_duprateExpDensCurve_mqc.txt
        touch ${meta.id}.R_sessionInfo.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            subread: \$(featureCounts -v 2>&1 | head -1 | sed 's/featureCounts //')
            python: \$(python3 --version | sed 's/Python //')
        END_VERSIONS
        """
    } else {
        """
        touch ${meta.id}_duprateExpDens.pdf
        touch ${meta.id}_duprateExpBoxplot.pdf
        touch ${meta.id}_expressionHist.pdf
        touch ${meta.id}_dupMatrix.txt
        touch ${meta.id}_intercept_slope.txt
        touch ${meta.id}_dup_intercept_mqc.txt
        touch ${meta.id}_duprateExpDensCurve_mqc.txt
        touch ${meta.id}.R_sessionInfo.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bioconductor-dupradar: \$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
        END_VERSIONS
        """
    }
}
