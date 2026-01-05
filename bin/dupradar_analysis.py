#!/usr/bin/env python3
"""
dupRadar analysis script
Calculates duplication rates from featureCounts output
Maintains full MultiQC compatibility
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import argparse


def calculate_duplication_rates(with_dups_counts, no_dups_counts, gene_lengths):
    """Calculate RPK and duplication rates from count data."""
    gene_lengths_kb = gene_lengths / 1000  # Convert to kb
    with_dups_rpk = with_dups_counts / gene_lengths_kb
    no_dups_rpk = no_dups_counts / gene_lengths_kb
    duplication_rate = np.where(with_dups_rpk > 0, 1 - (no_dups_rpk / with_dups_rpk), 0)
    return with_dups_rpk, duplication_rate


def create_dup_matrix(gene_ids, rpk, duplication_rate):
    """Create dupMatrix dataframe."""
    return pd.DataFrame(
        {
            "gene_id": gene_ids,
            "rpkm": rpk,  # Using RPK instead of RPKM for simplicity
            "duplication_rate": duplication_rate,
        }
    )


def calculate_regression(log_rpkm, duplication_rate):
    """Calculate linear regression for duplication rate vs expression."""
    slope, intercept, r_value, p_value, std_err = linregress(log_rpkm, duplication_rate)
    return slope, intercept, r_value**2


def write_intercept_slope(prefix, intercept, slope, r_squared):
    """Write intercept and slope file."""
    with open(f"{prefix}_intercept_slope.txt", "w") as f:
        f.write(f"intercept\t{intercept}\n")
        f.write(f"slope\t{slope}\n")
        f.write(f"r_squared\t{r_squared}\n")


def write_multiqc_files(prefix, intercept, slope, log_rpkm):
    """Generate MultiQC-compatible output files."""
    # General stats file
    with open(f"{prefix}_dup_intercept_mqc.txt", "w") as f:
        f.write("# plot_type: 'generalstats'\n")
        f.write("# pconfig:\n")
        f.write("#     dupradar_intercept:\n")
        f.write("#         title: 'dupRadar Intercept'\n")
        f.write("#         description: 'Intercept of the dupRadar fitted curve'\n")
        f.write("#         format: '{:.3f}'\n")
        f.write("Sample\tdupradar_intercept\n")
        f.write(f"{prefix}\t{intercept:.3f}\n")

    # Line graph data
    with open(f"{prefix}_duprateExpDensCurve_mqc.txt", "w") as f:
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


def generate_plots(prefix, log_rpkm, duplication_rate, slope, intercept):
    """Generate PDF plots for dupRadar analysis."""
    x_points = np.linspace(log_rpkm.min(), log_rpkm.max(), 100)
    y_points = slope * x_points + intercept

    # Main scatter plot (density plot equivalent)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    ax1.scatter(log_rpkm, duplication_rate, alpha=0.6, s=1)
    ax1.plot(x_points, y_points, "r-", linewidth=2)
    ax1.set_xlabel("Expression (log10 RPK)")
    ax1.set_ylabel("Duplication Rate")
    ax1.set_title("dupRadar: Expression vs Duplication Rate")

    ax2.boxplot([duplication_rate])
    ax2.set_ylabel("Duplication Rate")
    ax2.set_title("Duplication Rate Distribution")

    ax3.hist(log_rpkm, bins=50, alpha=0.7)
    ax3.set_xlabel("Expression (log10 RPK)")
    ax3.set_ylabel("Frequency")
    ax3.set_title("Expression Distribution")

    ax4.hist(duplication_rate, bins=50, alpha=0.7)
    ax4.set_xlabel("Duplication Rate")
    ax4.set_ylabel("Frequency")
    ax4.set_title("Duplication Rate Distribution")

    plt.tight_layout()
    plt.savefig(f"{prefix}_duprateExpDens.pdf")
    plt.close()

    # Individual plots for compatibility
    plt.figure(figsize=(8, 6))
    plt.boxplot([duplication_rate])
    plt.ylabel("Duplication Rate")
    plt.title("Duplication Rate Distribution")
    plt.savefig(f"{prefix}_duprateExpBoxplot.pdf")
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.hist(log_rpkm, bins=50, alpha=0.7)
    plt.xlabel("Expression (log10 RPK)")
    plt.ylabel("Frequency")
    plt.title("Expression Distribution")
    plt.savefig(f"{prefix}_expressionHist.pdf")
    plt.close()


def generate_empty_outputs(prefix):
    """Generate placeholder outputs when insufficient data."""
    with open(f"{prefix}_intercept_slope.txt", "w") as f:
        f.write("intercept\tNA\n")
        f.write("slope\tNA\n")
        f.write("r_squared\tNA\n")

    with open(f"{prefix}_dup_intercept_mqc.txt", "w") as f:
        f.write("# plot_type: 'generalstats'\n")
        f.write("# pconfig:\n")
        f.write("#     dupradar_intercept:\n")
        f.write("#         title: 'dupRadar Intercept'\n")
        f.write("#         description: 'Intercept of the dupRadar fitted curve'\n")
        f.write("#         format: '{:.3f}'\n")
        f.write("Sample\tdupradar_intercept\n")
        f.write(f"{prefix}\tNA\n")

    with open(f"{prefix}_duprateExpDensCurve_mqc.txt", "w") as f:
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

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(
        0.5,
        0.5,
        "Insufficient data for analysis",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    plt.savefig(f"{prefix}_duprateExpDens.pdf")
    plt.savefig(f"{prefix}_duprateExpBoxplot.pdf")
    plt.savefig(f"{prefix}_expressionHist.pdf")
    plt.close()


def run_tests():
    """Run baked-in tests to verify functionality."""
    import tempfile
    import os

    print("Running baked-in tests...")

    # Test 1: calculate_duplication_rates
    with_dups = np.array([100, 200, 50, 0])
    no_dups = np.array([80, 150, 50, 0])
    lengths = np.array([1000, 2000, 1000, 1000])
    rpk, dup_rate = calculate_duplication_rates(with_dups, no_dups, lengths)

    assert np.isclose(rpk[0], 100.0), (
        f"RPK calculation failed: expected 100, got {rpk[0]}"
    )
    assert np.isclose(dup_rate[0], 0.2), (
        f"Dup rate failed: expected 0.2, got {dup_rate[0]}"
    )
    assert np.isclose(dup_rate[2], 0.0), (
        f"Zero dup rate failed: expected 0, got {dup_rate[2]}"
    )
    assert np.isclose(dup_rate[3], 0.0), f"Zero expression dup rate failed"
    print("  ✓ calculate_duplication_rates")

    # Test 2: create_dup_matrix
    gene_ids = ["gene1", "gene2"]
    rpk = np.array([100.0, 50.0])
    dup_rate = np.array([0.2, 0.1])
    df = create_dup_matrix(gene_ids, rpk, dup_rate)

    assert len(df) == 2, "dup_matrix row count wrong"
    assert list(df.columns) == ["gene_id", "rpkm", "duplication_rate"], (
        "dup_matrix columns wrong"
    )
    assert df["gene_id"].iloc[0] == "gene1", "gene_id column wrong"
    print("  ✓ create_dup_matrix")

    # Test 3: calculate_regression
    log_rpkm = np.array([1.0, 2.0, 3.0, 4.0])
    dup_rate = np.array([0.1, 0.2, 0.3, 0.4])
    slope, intercept, r_sq = calculate_regression(log_rpkm, dup_rate)

    assert np.isclose(slope, 0.1, atol=0.01), f"Slope wrong: {slope}"
    assert np.isclose(intercept, 0.0, atol=0.01), f"Intercept wrong: {intercept}"
    assert np.isclose(r_sq, 1.0, atol=0.01), f"R-squared wrong: {r_sq}"
    print("  ✓ calculate_regression")

    # Test 4: File output (write_intercept_slope)
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = os.path.join(tmpdir, "test")
        write_intercept_slope(prefix, 0.5, 0.1, 0.95)

        with open(f"{prefix}_intercept_slope.txt") as f:
            content = f.read()
        assert "intercept\t0.5" in content, "intercept not written correctly"
        assert "slope\t0.1" in content, "slope not written correctly"
        assert "r_squared\t0.95" in content, "r_squared not written correctly"
        print("  ✓ write_intercept_slope")

    print("All tests passed!")
    return True


def main():
    parser = argparse.ArgumentParser(description="dupRadar analysis")
    parser.add_argument(
        "--with-dups", required=False, help="featureCounts output with duplicates"
    )
    parser.add_argument(
        "--no-dups", required=False, help="featureCounts output without duplicates"
    )
    parser.add_argument("--prefix", required=False, help="Output prefix")
    parser.add_argument("--test", action="store_true", help="Run baked-in tests")
    args = parser.parse_args()

    if args.test:
        success = run_tests()
        sys.exit(0 if success else 1)

    if not all([args.with_dups, args.no_dups, args.prefix]):
        parser.error(
            "--with-dups, --no-dups, and --prefix are required unless --test is specified"
        )

    # Read featureCounts outputs
    with_dups = pd.read_csv(args.with_dups, sep="\t", comment="#", skiprows=1)
    no_dups = pd.read_csv(args.no_dups, sep="\t", comment="#", skiprows=1)

    # Extract count columns (last column is the BAM file)
    with_dups_counts = with_dups.iloc[:, -1].values
    no_dups_counts = no_dups.iloc[:, -1].values
    gene_lengths = with_dups["Length"].values

    # Calculate duplication rates
    rpk, duplication_rate = calculate_duplication_rates(
        with_dups_counts, no_dups_counts, gene_lengths
    )

    # Create and filter dupMatrix
    dup_matrix = create_dup_matrix(with_dups["Geneid"], rpk, duplication_rate)
    dup_matrix_filtered = dup_matrix[dup_matrix["rpkm"] > 0]

    # Save dupMatrix
    dup_matrix_filtered.to_csv(f"{args.prefix}_dupMatrix.txt", sep="\t", index=False)

    # Calculate and output results
    if len(dup_matrix_filtered) > 1:
        log_rpkm = np.log10(dup_matrix_filtered["rpkm"].values)
        dup_rates = dup_matrix_filtered["duplication_rate"].values

        slope, intercept, r_squared = calculate_regression(log_rpkm, dup_rates)

        write_intercept_slope(args.prefix, intercept, slope, r_squared)
        write_multiqc_files(args.prefix, intercept, slope, log_rpkm)
        generate_plots(args.prefix, log_rpkm, dup_rates, slope, intercept)
    else:
        generate_empty_outputs(args.prefix)

    print("Analysis complete")


if __name__ == "__main__":
    main()
