#!/usr/bin/env python3
"""
Generate MultiQC compatible files for dupRadar
"""
import sys
import argparse

def generate_multiqc_files(dupmatrix_file, intercept_slope_file, output_prefix):
    """Generate MultiQC compatible files"""
    
    # Extract intercept from intercept_slope file
    try:
        with open(intercept_slope_file, 'r') as f:
            for line in f:
                if 'dupRadar Int' in line:
                    intercept = float(line.split(':')[1].strip())
                    break
        else:
            intercept = 0.1  # Default
    except (FileNotFoundError, ValueError):
        intercept = 0.1

    # Create MultiQC intercept file
    sample_name = output_prefix.replace('Aligned.sortedByCoord.out.markDups', '')
    intercept_file = f"{output_prefix}_dup_intercept_mqc.txt"
    
    with open(intercept_file, 'w') as f:
        f.write("#id: DupInt\n")
        f.write("#plot_type: 'generalstats'\n")
        f.write("#pconfig:\n")
        f.write("#    dupRadar_intercept:\n")
        f.write("#        title: 'dupInt'\n")
        f.write("#        namespace: 'DupRadar'\n")
        f.write("#        description: 'Intercept value from DupRadar'\n")
        f.write("#        max: 100\n")
        f.write("#        min: 0\n")
        f.write("#        scale: 'RdYlGn-rev'\n")
        f.write("Sample dupRadar_intercept\n")
        f.write(f"{sample_name} {intercept}\n")

    # Create MultiQC curve data file
    curve_file = f"{output_prefix}_duprateExpDensCurve_mqc.txt"
    
    with open(curve_file, 'w') as f:
        f.write("#id: dupradar\n")
        f.write("#plot_type: 'linegraph'\n")
        f.write("#section_name: 'DupRadar'\n")
        f.write("#section_href: 'bioconductor.org/packages/release/bioc/html/dupRadar.html'\n")
        f.write('#description: "provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.\n')
        f.write('#    This plot shows the general linear models - a summary of the gene duplication distributions. "\n')
        f.write("#pconfig:\n")
        f.write("#    title: 'DupRadar General Linear Model'\n")
        f.write("#    xlog: True\n")
        f.write("#    xlab: 'expression (reads/kbp)'\n")
        f.write("#    ylab: '% duplicate reads'\n")
        f.write("#    ymax: 100\n")
        f.write("#    ymin: 0\n")
        f.write("#    tt_label: '<b>{point.x:.1f} reads/kbp</b>: {point.y:,.2f}% duplicates'\n")
        f.write("#    x_lines:\n")
        f.write("#        - color: 'green'\n")
        f.write("#          dash: 'LongDash'\n")
        f.write("#          label:\n")
        f.write("#                text: '0.5 RPKM'\n")
        f.write("#          value: 0.5\n")
        f.write("#          width: 1\n")
        f.write("#        - color: 'red'\n")
        f.write("#          dash: 'LongDash'\n")
        f.write("#          label:\n")
        f.write("#                text: '1 read/bp'\n")
        f.write("#          value: 1000\n")
        f.write("#          width: 1\n")

        # Generate curve data points (sample every 10th point for efficiency)
        try:
            with open(dupmatrix_file, 'r') as infile:
                header = infile.readline()  # Skip header
                data_points = []
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) >= 7:
                        try:
                            rpk = float(parts[5])
                            dup_rate = float(parts[4]) * 100  # Convert to percentage
                            if rpk > 0 and 0 <= dup_rate <= 100:
                                data_points.append((rpk, dup_rate))
                        except (ValueError, IndexError):
                            continue
                
                # Sort by RPK and sample every 10th point
                data_points.sort(key=lambda x: x[0])
                for i in range(0, len(data_points), 10):
                    rpk, dup_rate = data_points[i]
                    f.write(f"{rpk} {dup_rate}\n")
                    
        except FileNotFoundError:
            print(f"Warning: Could not read {dupmatrix_file} for curve data", file=sys.stderr)

    print(f"Generated MultiQC files: {intercept_file}, {curve_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate MultiQC files for dupRadar')
    parser.add_argument('dupmatrix_file', help='Input dupMatrix.txt file')
    parser.add_argument('intercept_slope_file', help='Input intercept_slope.txt file')
    parser.add_argument('output_prefix', help='Output file prefix')
    
    args = parser.parse_args()
    generate_multiqc_files(args.dupmatrix_file, args.intercept_slope_file, args.output_prefix)