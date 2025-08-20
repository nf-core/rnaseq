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
        #!/usr/bin/env bash
        # Fast dupRadar replacement using native featureCounts
        # Written for nf-core/rnaseq to fix Fusion S3 performance issues
        # Maintains MultiQC compatibility while bypassing R overhead

        set -euo pipefail

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
        echo "Calculating duplication rates..."
        python3 << PYEOF
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Read featureCounts outputs
with_dups = pd.read_csv("\${output_prefix}_with_dups.txt", sep='\t', comment='#', skiprows=1)
no_dups = pd.read_csv("\${output_prefix}_no_dups.txt", sep='\t', comment='#', skiprows=1)

# Extract count columns (last column is the BAM file)
with_dups_counts = with_dups.iloc[:, -1]
no_dups_counts = no_dups.iloc[:, -1]

# Calculate RPK (Reads per Kilobase)
gene_lengths = with_dups['Length'] / 1000  # Convert to kb
with_dups_rpk = with_dups_counts / gene_lengths
no_dups_rpk = no_dups_counts / gene_lengths

# Calculate duplication rate
duplication_rate = np.where(with_dups_rpk > 0, 1 - (no_dups_rpk / with_dups_rpk), 0)

# Create dupMatrix
dup_matrix = pd.DataFrame({
    'gene_id': with_dups['Geneid'],
    'rpkm': with_dups_rpk,  # Using RPK instead of RPKM for simplicity
    'duplication_rate': duplication_rate
})

# Filter out genes with zero expression
dup_matrix_filtered = dup_matrix[dup_matrix['rpkm'] > 0]

# Save dupMatrix
dup_matrix_filtered.to_csv("\${output_prefix}_dupMatrix.txt", sep='\t', index=False)

# Calculate intercept and slope
if len(dup_matrix_filtered) > 1:
    log_rpkm = np.log10(dup_matrix_filtered['rpkm'])
    slope, intercept, r_value, p_value, std_err = linregress(log_rpkm, dup_matrix_filtered['duplication_rate'])
    
    # Save intercept and slope
    with open("\${output_prefix}_intercept_slope.txt", 'w') as f:
        f.write(f"intercept\\t{intercept}\\n")
        f.write(f"slope\\t{slope}\\n")
        f.write(f"r_squared\\t{r_value**2}\\n")
    
    # Generate MultiQC files
    with open("\${output_prefix}_dup_intercept_mqc.txt", 'w') as f:
        f.write("# plot_type: 'generalstats'\\n")
        f.write("# pconfig:\\n")
        f.write("#     dupradar_intercept:\\n")
        f.write("#         title: 'dupRadar Intercept'\\n")
        f.write("#         description: 'Intercept of the dupRadar fitted curve'\\n")
        f.write("#         format: '{:.3f}'\\n")
        f.write("Sample\\tdupradar_intercept\\n")
        f.write(f"\${output_prefix}\\t{intercept:.3f}\\n")
    
    with open("\${output_prefix}_duprateExpDensCurve_mqc.txt", 'w') as f:
        f.write("# plot_type: 'linegraph'\\n")
        f.write("# section_name: 'dupRadar'\\n")
        f.write("# description: 'Expression vs Duplication Rate'\\n")
        f.write("# pconfig:\\n")
        f.write("#     title: 'dupRadar: Expression vs Duplication Rate'\\n")
        f.write("#     xlab: 'Expression (log10 RPK)'\\n")
        f.write("#     ylab: 'Duplication Rate'\\n")
        
        # Sample data points for the curve
        x_points = np.linspace(log_rpkm.min(), log_rpkm.max(), 100)
        y_points = slope * x_points + intercept
        
        f.write("log10_rpk\\tduplication_rate\\n")
        for x, y in zip(x_points, y_points):
            f.write(f"{x:.3f}\\t{y:.3f}\\n")

    # Generate basic plots using matplotlib
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Scatter plot
    ax1.scatter(log_rpkm, dup_matrix_filtered['duplication_rate'], alpha=0.6, s=1)
    ax1.plot(x_points, y_points, 'r-', linewidth=2)
    ax1.set_xlabel('Expression (log10 RPK)')
    ax1.set_ylabel('Duplication Rate')
    ax1.set_title('dupRadar: Expression vs Duplication Rate')
    
    # Box plot (simplified)
    ax2.boxplot([dup_matrix_filtered['duplication_rate']])
    ax2.set_ylabel('Duplication Rate')
    ax2.set_title('Duplication Rate Distribution')
    
    # Expression histogram
    ax3.hist(log_rpkm, bins=50, alpha=0.7)
    ax3.set_xlabel('Expression (log10 RPK)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Expression Distribution')
    
    # Duplication rate histogram
    ax4.hist(dup_matrix_filtered['duplication_rate'], bins=50, alpha=0.7)
    ax4.set_xlabel('Duplication Rate')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Duplication Rate Distribution')
    
    plt.tight_layout()
    plt.savefig("\${output_prefix}_duprateExpDens.pdf")
    plt.close()
    
    # Individual plots for compatibility
    plt.figure(figsize=(8, 6))
    plt.boxplot([dup_matrix_filtered['duplication_rate']])
    plt.ylabel('Duplication Rate')
    plt.title('Duplication Rate Distribution')
    plt.savefig("\${output_prefix}_duprateExpBoxplot.pdf")
    plt.close()
    
    plt.figure(figsize=(8, 6))
    plt.hist(log_rpkm, bins=50, alpha=0.7)
    plt.xlabel('Expression (log10 RPK)')
    plt.ylabel('Frequency')
    plt.title('Expression Distribution')
    plt.savefig("\${output_prefix}_expressionHist.pdf")
    plt.close()

else:
    # Handle case with no/insufficient data
    with open("\${output_prefix}_intercept_slope.txt", 'w') as f:
        f.write("intercept\\tNA\\n")
        f.write("slope\\tNA\\n")
        f.write("r_squared\\tNA\\n")
    
    # Generate empty MultiQC files
    with open("\${output_prefix}_dup_intercept_mqc.txt", 'w') as f:
        f.write("# plot_type: 'generalstats'\\n")
        f.write("# pconfig:\\n")
        f.write("#     dupradar_intercept:\\n")
        f.write("#         title: 'dupRadar Intercept'\\n")
        f.write("#         description: 'Intercept of the dupRadar fitted curve'\\n")
        f.write("#         format: '{:.3f}'\\n")
        f.write("Sample\\tdupradar_intercept\\n")
        f.write(f"\${output_prefix}\\tNA\\n")
    
    with open("\${output_prefix}_duprateExpDensCurve_mqc.txt", 'w') as f:
        f.write("# plot_type: 'linegraph'\\n")
        f.write("# section_name: 'dupRadar'\\n")
        f.write("# description: 'Expression vs Duplication Rate'\\n")
        f.write("# pconfig:\\n")
        f.write("#     title: 'dupRadar: Expression vs Duplication Rate'\\n")
        f.write("#     xlab: 'Expression (log10 RPK)'\\n")
        f.write("#     ylab: 'Duplication Rate'\\n")
        f.write("log10_rpk\\tduplication_rate\\n")
        f.write("0\\t0\\n")
        f.write("1\\t0\\n")
    
    # Create empty plots
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, 'Insufficient data for analysis', 
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    plt.savefig("\${output_prefix}_duprateExpDens.pdf")
    plt.savefig("\${output_prefix}_duprateExpBoxplot.pdf")
    plt.savefig("\${output_prefix}_expressionHist.pdf")
    plt.close()

print("Analysis complete")
PYEOF

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

        # Step 6: Clean up intermediate files
        rm -f "\${output_prefix}_with_dups.txt" "\${output_prefix}_no_dups.txt"
        rm -f "\${output_prefix}_with_dups.txt.summary" "\${output_prefix}_no_dups.txt.summary"

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
