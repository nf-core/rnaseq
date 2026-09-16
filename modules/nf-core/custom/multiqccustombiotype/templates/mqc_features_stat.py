#!/usr/bin/env python

# Combine featureCounts biotype counts with a header for MultiQC,
# then calculate feature percentages for the general stats table.
#
# Written by Senthilkumar Panneerselvam and released under the MIT license.
# Adapted for nf-core/modules by Jonathan Manning.

import argparse
import logging
import os
import platform
import shlex
import sys

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

# Template variables from Nextflow
count_file = "${count}"
header_file = "${header}"
prefix = "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${meta.id}"
sample_name = "${meta.id}"


def parse_ext_args(args_string):
    """Parse arguments supplied via Nextflow's ext.args."""
    if args_string in (None, "", "null"):
        args_string = ""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max_biotypes",
        type=int,
        default=100,
        help="Maximum number of biotype rows to write to the MultiQC TSV. "
        "Above this, the TSV is dropped so MultiQC does not try to plot a "
        "category per gene.",
    )
    return parser.parse_args(shlex.split(args_string))


ext_args = parse_ext_args("${task.ext.args}")
max_biotypes = ext_args.max_biotypes


def prepare_biotype_counts(count_file, header_file, prefix):
    """Cut gene ID and count columns from featureCounts output, prepend header."""
    outfile = f"{prefix}.biotype_counts_mqc.tsv"

    with open(header_file) as hf:
        header = hf.read()

    with open(count_file) as cf:
        lines = cf.readlines()[2:]  # skip featureCounts header lines

    n_rows = 0
    with open(outfile, "w") as of:
        of.write(header)
        for line in lines:
            fields = line.strip().split("\\t")
            if len(fields) >= 7:
                of.write(f"{fields[0]}\\t{fields[6]}\\n")
                n_rows += 1

    return outfile, n_rows


mqc_main = """#id: 'biotype-gs'
#plot_type: 'generalstats'
#pconfig:"""

mqc_pconf = """#    percent_{ft}:
#        title: '% {ft}'
#        namespace: 'Biotype Counts'
#        description: '% reads overlapping {ft} features'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'"""


def mqc_feature_stat(bfile, features, outfile, sname=None):
    """Calculate features percentage for biotype counts."""
    if not sname:
        sname = os.path.splitext(os.path.basename(bfile))[0]

    fcounts = {}
    try:
        with open(bfile) as bfl:
            for ln in bfl:
                if ln.startswith("#"):
                    continue
                parts = ln.strip().split("\\t")
                if len(parts) == 2:
                    ft, cn = parts
                    fcounts[ft] = float(cn)
    except Exception:
        logger.error(f"Trouble reading the biocount file {bfile}")
        return

    total_count = sum(fcounts.values())
    if total_count == 0:
        logger.error("No biocounts found, exiting")
        return

    fpercent = {f: (fcounts[f] / total_count) * 100 if f in fcounts else 0 for f in features}
    if len(fpercent) == 0:
        logger.error(f"None of given features '{', '.join(features)}' found in the biocount file {bfile}")
        return

    out_head, out_value, out_mqc = ("Sample", f"'{sname}'", mqc_main)
    for ft, pt in fpercent.items():
        out_head = f"{out_head}\\tpercent_{ft}"
        out_value = f"{out_value}\\t{pt}"
        out_mqc = f"{out_mqc}\\n{mqc_pconf.format(ft=ft)}"

    with open(outfile, "w") as ofl:
        out_final = "\\n".join([out_mqc, out_head, out_value]).strip()
        ofl.write(out_final + "\\n")


if __name__ == "__main__":
    # Step 1: Prepare biotype counts TSV from featureCounts output
    biotype_file, n_biotypes = prepare_biotype_counts(count_file, header_file, prefix)

    # Step 2: Fail loudly if the first column of the featureCounts output has
    # too many unique values. MultiQC cannot render a bar plot with thousands
    # of categories, so a high-cardinality group attribute (e.g. a per-gene
    # identifier) passed to `featureCounts -g` is almost certainly a config
    # mistake upstream.
    if n_biotypes > max_biotypes:
        logger.error(
            f"Too many categories in '{count_file}': {n_biotypes} > "
            f"--max_biotypes={max_biotypes}. Column 1 of this file reflects "
            "the attribute passed to `featureCounts -g`; if a high-cardinality "
            "attribute was used (e.g. a per-gene identifier), MultiQC cannot "
            "plot it. Re-run featureCounts with a lower-cardinality attribute "
            "(e.g. gene_biotype, gene_type), or raise the limit via "
            "`ext.args = '--max_biotypes N'` if this is expected."
        )
        sys.exit(1)

    # Step 3: Calculate rRNA percentage for MultiQC general stats
    mqc_feature_stat(
        biotype_file,
        ["rRNA"],
        f"{prefix}.biotype_counts_rrna_mqc.tsv",
        sname=sample_name,
    )

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
