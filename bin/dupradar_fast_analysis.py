#!/usr/bin/env python3
"""
Fast dupRadar analysis script
Replaces R/dupRadar with Python for high-performance processing
Maintains full MultiQC compatibility
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import argparse

def main():
    parser = argparse.ArgumentParser(description='Fast dupRadar analysis')
    parser.add_argument('--with-dups', required=True, help='featureCounts output with duplicates')
    parser.add_argument('--no-dups', required=True, help='featureCounts output without duplicates')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    args = parser.parse_args()

    # Read featureCounts outputs
    with_dups = pd.read_csv(args.with_dups, sep='\t', comment='#', skiprows=1)
    no_dups = pd.read_csv(args.no_dups, sep='\t', comment='#', skiprows=1)

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
    dup_matrix_filtered.to_csv(f"{args.prefix}_dupMatrix.txt", sep='\t', index=False)

    # Calculate intercept and slope
    if len(dup_matrix_filtered) > 1:
        log_rpkm = np.log10(dup_matrix_filtered['rpkm'])
        slope, intercept, r_value, p_value, std_err = linregress(log_rpkm, dup_matrix_filtered['duplication_rate'])
        
        # Save intercept and slope
        with open(f"{args.prefix}_intercept_slope.txt", 'w') as f:
            f.write(f"intercept\t{intercept}\n")
            f.write(f"slope\t{slope}\n")
            f.write(f"r_squared\t{r_value**2}\n")
        
        # Generate MultiQC files
        with open(f"{args.prefix}_dup_intercept_mqc.txt", 'w') as f:
            f.write("# plot_type: 'generalstats'\n")
            f.write("# pconfig:\n")
            f.write("#     dupradar_intercept:\n")
            f.write("#         title: 'dupRadar Intercept'\n")
            f.write("#         description: 'Intercept of the dupRadar fitted curve'\n")
            f.write("#         format: '{:.3f}'\n")
            f.write("Sample\tdupradar_intercept\n")
            f.write(f"{args.prefix}\t{intercept:.3f}\n")
        
        with open(f"{args.prefix}_duprateExpDensCurve_mqc.txt", 'w') as f:
            f.write("# plot_type: 'linegraph'\n")
            f.write("# section_name: 'dupRadar'\n")
            f.write("# description: 'Expression vs Duplication Rate'\n")
            f.write("# pconfig:\n")
            f.write("#     title: 'dupRadar: Expression vs Duplication Rate'\n")
            f.write("#     xlab: 'Expression (log10 RPK)'\n")
            f.write("#     ylab: 'Duplication Rate'\n")
            
            # Sample data points for the curve
            x_points = np.linspace(log_rpkm.min(), log_rpkm.max(), 100)
            y_points = slope * x_points + intercept
            
            f.write("log10_rpk\tduplication_rate\n")
            for x, y in zip(x_points, y_points):
                f.write(f"{x:.3f}\t{y:.3f}\n")

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
        plt.savefig(f"{args.prefix}_duprateExpDens.pdf")
        plt.close()
        
        # Individual plots for compatibility
        plt.figure(figsize=(8, 6))
        plt.boxplot([dup_matrix_filtered['duplication_rate']])
        plt.ylabel('Duplication Rate')
        plt.title('Duplication Rate Distribution')
        plt.savefig(f"{args.prefix}_duprateExpBoxplot.pdf")
        plt.close()
        
        plt.figure(figsize=(8, 6))
        plt.hist(log_rpkm, bins=50, alpha=0.7)
        plt.xlabel('Expression (log10 RPK)')
        plt.ylabel('Frequency')
        plt.title('Expression Distribution')
        plt.savefig(f"{args.prefix}_expressionHist.pdf")
        plt.close()

    else:
        # Handle case with no/insufficient data
        with open(f"{args.prefix}_intercept_slope.txt", 'w') as f:
            f.write("intercept\tNA\n")
            f.write("slope\tNA\n")
            f.write("r_squared\tNA\n")
        
        # Generate empty MultiQC files
        with open(f"{args.prefix}_dup_intercept_mqc.txt", 'w') as f:
            f.write("# plot_type: 'generalstats'\n")
            f.write("# pconfig:\n")
            f.write("#     dupradar_intercept:\n")
            f.write("#         title: 'dupRadar Intercept'\n")
            f.write("#         description: 'Intercept of the dupRadar fitted curve'\n")
            f.write("#         format: '{:.3f}'\n")
            f.write("Sample\tdupradar_intercept\n")
            f.write(f"{args.prefix}\tNA\n")
        
        with open(f"{args.prefix}_duprateExpDensCurve_mqc.txt", 'w') as f:
            f.write("# plot_type: 'linegraph'\n")
            f.write("# section_name: 'dupRadar'\n")
            f.write("# description: 'Expression vs Duplication Rate'\n")
            f.write("# pconfig:\n")
            f.write("#     title: 'dupRadar: Expression vs Duplication Rate'\n")
            f.write("#     xlab: 'Expression (log10 RPK)'\n")
            f.write("#     ylab: 'Duplication Rate'\n")
            f.write("log10_rpk\tduplication_rate\n")
            f.write("0\t0\n")
            f.write("1\t0\n")
        
        # Create empty plots
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'Insufficient data for analysis', 
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plt.savefig(f"{args.prefix}_duprateExpDens.pdf")
        plt.savefig(f"{args.prefix}_duprateExpBoxplot.pdf")
        plt.savefig(f"{args.prefix}_expressionHist.pdf")
        plt.close()

    print("Analysis complete")

if __name__ == "__main__":
    main()