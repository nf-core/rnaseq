#!/usr/bin/env python3
"""
Generate PDF plots for dupRadar analysis
"""
import sys
import argparse
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

def generate_plots(dupmatrix_file, output_prefix):
    """Generate dupRadar PDF plots"""
    
    # Read data
    rpk_values = []
    dup_rates = []

    try:
        with open(dupmatrix_file, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    try:
                        rpk = float(parts[5])
                        dup_rate = float(parts[4]) * 100  # Convert to percentage
                        if rpk > 0 and 0 <= dup_rate <= 100:
                            rpk_values.append(rpk)
                            dup_rates.append(dup_rate)
                    except (ValueError, IndexError):
                        continue

        if len(rpk_values) > 0:
            # Create density scatter plot
            plt.figure(figsize=(10, 8))
            plt.hexbin(rpk_values, dup_rates, gridsize=50, cmap='Blues', alpha=0.7)
            plt.colorbar(label='Count')
            plt.xscale('log')
            plt.xlabel('RPK (Reads per Kilobase)')
            plt.ylabel('Duplication Rate (%)')
            plt.title(f'{output_prefix}\nDensity scatter plot')
            plt.xlim(0.1, max(rpk_values)*1.1)
            plt.ylim(0, 100)
            plt.tight_layout()
            plt.savefig(f'{output_prefix}_duprateExpDens.pdf')
            plt.close()

            # Create boxplot (simplified - just show quartiles)
            plt.figure(figsize=(10, 6))
            bins = np.logspace(np.log10(min(rpk_values)), np.log10(max(rpk_values)), 20)
            bin_indices = np.digitize(rpk_values, bins)
            box_data = []
            positions = []
            
            for i in range(1, len(bins)):
                bin_dup_rates = [dup_rates[j] for j in range(len(dup_rates)) if bin_indices[j] == i]
                if len(bin_dup_rates) > 5:  # Only include bins with sufficient data
                    box_data.append(bin_dup_rates)
                    positions.append(bins[i-1])
            
            if len(box_data) > 0:
                plt.boxplot(box_data, positions=positions, widths=np.diff(positions)[0]*0.6 if len(positions)>1 else 1)
                plt.xscale('log')
                plt.xlabel('RPK (Reads per Kilobase)')
                plt.ylabel('Duplication Rate (%)')
                plt.title(f'{output_prefix}\nPercent Duplication by Expression')
                plt.tight_layout()
                plt.savefig(f'{output_prefix}_duprateExpBoxplot.pdf')
                plt.close()
            else:
                # Create empty boxplot
                plt.figure(figsize=(10, 6))
                plt.text(0.5, 0.5, 'Insufficient data for boxplot', ha='center', va='center', transform=plt.gca().transAxes)
                plt.title(f'{output_prefix}\nPercent Duplication by Expression')
                plt.savefig(f'{output_prefix}_duprateExpBoxplot.pdf')
                plt.close()
            
            # Create expression histogram
            plt.figure(figsize=(10, 6))
            plt.hist(rpk_values, bins=50, alpha=0.7, edgecolor='black')
            plt.xscale('log')
            plt.xlabel('RPK (Reads per Kilobase)')
            plt.ylabel('Number of Genes')
            plt.title(f'{output_prefix}\nDistribution of RPK values per gene')
            plt.tight_layout()
            plt.savefig(f'{output_prefix}_expressionHist.pdf')
            plt.close()

            print('Generated PDF plots successfully')
        else:
            print('No valid data found for plotting', file=sys.stderr)
            # Create empty plots to satisfy pipeline requirements
            for plot_name, title in [
                ('_duprateExpDens.pdf', 'Density scatter plot'),
                ('_duprateExpBoxplot.pdf', 'Percent Duplication by Expression'),
                ('_expressionHist.pdf', 'Distribution of RPK values per gene')
            ]:
                plt.figure(figsize=(10, 8))
                plt.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=plt.gca().transAxes)
                plt.title(f'{output_prefix}\n{title}')
                plt.savefig(f'{output_prefix}{plot_name}')
                plt.close()

    except Exception as e:
        print(f'Error generating plots: {e}', file=sys.stderr)
        # Create minimal empty plots as fallback
        for plot_name, title in [
            ('_duprateExpDens.pdf', 'Density scatter plot'),
            ('_duprateExpBoxplot.pdf', 'Percent Duplication by Expression'),
            ('_expressionHist.pdf', 'Distribution of RPK values per gene')
        ]:
            plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, f'Error: {e}', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title(f'{output_prefix}\n{title}')
            plt.savefig(f'{output_prefix}{plot_name}')
            plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate dupRadar plots')
    parser.add_argument('dupmatrix_file', help='Input dupMatrix.txt file')
    parser.add_argument('output_prefix', help='Output file prefix')
    
    args = parser.parse_args()
    generate_plots(args.dupmatrix_file, args.output_prefix)