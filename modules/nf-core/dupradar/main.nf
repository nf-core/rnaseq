process DUPRADAR {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/subread_python_matplotlib_numpy_pruned:70e1e046570ad3a3'

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
    def prefix = task.ext.prefix ?: meta.id
    def strandedness = meta.strandedness == 'forward' ? 1 : meta.strandedness == 'reverse' ? 2 : 0
    def paired_flag = meta.single_end ? '' : '-p'
    """
    # Feature type (default: exon)
    feature_type="exon"
    # Parse task.ext.args for feature_type if provided
    if [[ '$task.ext.args' == *"feature_type"* ]]; then
        feature_type=\$(echo '$task.ext.args' | grep -o 'feature_type[[:space:]]*[^[:space:]]*' | cut -d' ' -f2 || echo "exon")
    fi

    echo "=== dupRadar Analysis ==="
    echo "Input BAM: $bam"
    echo "Input GTF: $gtf"
    echo "Strandedness: $strandedness"
    echo "Library type: ${meta.single_end ? 'single-end' : 'paired-end'}"
    echo "Feature type: \$feature_type"
    echo "Threads: $task.cpus"
    echo "Output prefix: $prefix"

    # Step 1: Run featureCounts with duplicates included
    echo "Running featureCounts with duplicates..."
    featureCounts \\
        -T $task.cpus \\
        -s $strandedness \\
        $paired_flag \\
        -t \$feature_type \\
        -g gene_id \\
        -a $gtf \\
        -o "${prefix}_with_dups.txt" \\
        $bam

    # Step 2: Run featureCounts without duplicates
    echo "Running featureCounts without duplicates..."
    featureCounts \\
        -T $task.cpus \\
        -s $strandedness \\
        $paired_flag \\
        -t \$feature_type \\
        -g gene_id \\
        -a $gtf \\
        -o "${prefix}_no_dups.txt" \\
        --ignoreDup \\
        $bam

    # Step 3: Process results and calculate duplication rates
    echo "Running dupRadar analysis..."
    dupradar_analysis.py \\
        --with-dups "${prefix}_with_dups.txt" \\
        --no-dups "${prefix}_no_dups.txt" \\
        --prefix "$prefix"

    # Step 4: Create session info log
    cat > "${prefix}.R_sessionInfo.log" << EOF
dupRadar Analysis (featureCounts + Python implementation)
Date: \$(date)
featureCounts version: \$(featureCounts -v 2>&1 | head -1)
Python version: \$(python3 --version)
System: \$(uname -a)

This analysis uses native featureCounts for improved performance with cloud storage.
EOF

    # Step 5: Create versions.yml
    cat > versions.yml << EOF
"${task.process}":
    subread: \$(featureCounts -v 2>&1 | head -1 | sed 's/featureCounts //' | sed 's/ .*//')
    python: \$(python3 --version | sed 's/Python //')
EOF

    echo "=== dupRadar Analysis Complete ==="
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_duprateExpDens.pdf
    touch ${prefix}_duprateExpBoxplot.pdf
    touch ${prefix}_expressionHist.pdf
    touch ${prefix}_dupMatrix.txt
    touch ${prefix}_intercept_slope.txt
    touch ${prefix}_dup_intercept_mqc.txt
    touch ${prefix}_duprateExpDensCurve_mqc.txt
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | head -1 | sed 's/featureCounts //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
