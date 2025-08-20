#!/usr/bin/env python3
"""
Calculate intercept and slope for dupRadar GLM from duplication matrix
"""
import sys
import math
import argparse

def calculate_intercept_slope(dupmatrix_file, output_prefix):
    """Calculate intercept and slope using simple linear regression"""
    data = []
    
    try:
        with open(dupmatrix_file, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    try:
                        rpk = float(parts[5])
                        dup_rate = float(parts[4])
                        if rpk > 0 and dup_rate >= 0 and dup_rate <= 1:
                            data.append((math.log10(rpk), dup_rate * 100))  # Convert to percentage
                    except (ValueError, OverflowError):
                        continue
    except FileNotFoundError:
        print(f"Error: Could not find file {dupmatrix_file}", file=sys.stderr)
        sys.exit(1)

    if len(data) < 10:
        print('Warning: Insufficient data for regression', file=sys.stderr)
        intercept, slope = 0.1, 0.1  # Default values
    else:
        # Simple linear regression: y = mx + b
        n = len(data)
        sum_x = sum(x for x, y in data)
        sum_y = sum(y for x, y in data)
        sum_xy = sum(x*y for x, y in data)
        sum_xx = sum(x*x for x, y in data)
        
        if n * sum_xx - sum_x * sum_x != 0:
            slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x)
            intercept = (sum_y - slope * sum_x) / n
        else:
            intercept, slope = 0.1, 0.1

    # Convert intercept to the expected dupRadar format (fraction)
    intercept_frac = intercept / 100.0 if intercept > 1 else intercept

    # Write to output file
    output_file = f"{output_prefix}_intercept_slope.txt"
    with open(output_file, 'w') as f:
        f.write(f'- dupRadar Int (duprate at low read counts): {intercept_frac}\n')
        f.write(f'- dupRadar Sl (progression of the duplication rate): {slope}\n')

    print(f'Intercept: {intercept_frac}, Slope: {slope}')
    return intercept_frac, slope

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate dupRadar intercept and slope')
    parser.add_argument('dupmatrix_file', help='Input dupMatrix.txt file')
    parser.add_argument('output_prefix', help='Output file prefix')
    
    args = parser.parse_args()
    calculate_intercept_slope(args.dupmatrix_file, args.output_prefix)