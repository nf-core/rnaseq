#!/usr/bin/env python

# Written by Senthilkumar Panneerselvam and released under the MIT license.

import argparse
import logging
import os

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

mqc_main = """#id: 'biotype-gs'
#plot_type: 'generalstats'
#pconfig:"""

mqc_pconf = """#    percent_{ft}:
#        title: '% {ft}'
#        namespace: 'Biotype Counts'
#        description: '% reads overlapping {ft} features'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'
#        format: '{{:.2f}}%'"""


def mqc_feature_stat(bfile, features, outfile, sname=None):
    # If sample name not given use file name
    if not sname:
        sname = os.path.splitext(os.path.basename(bfile))[0]

    # Try to parse and read biocount file
    fcounts = {}
    try:
        with open(bfile, "r") as bfl:
            for ln in bfl:
                if ln.startswith("#"):
                    continue
                ft, cn = ln.strip().split("\t")
                fcounts[ft] = float(cn)
    except:
        logger.error("Trouble reading the biocount file {}".format(bfile))
        return

    total_count = sum(fcounts.values())
    if total_count == 0:
        logger.error("No biocounts found, exiting")
        return

    # Calculate percentage for each requested feature
    fpercent = {f: (fcounts[f] / total_count) * 100 if f in fcounts else 0 for f in features}
    if len(fpercent) == 0:
        logger.error("Any of given features '{}' not found in the biocount file".format(", ".join(features), bfile))
        return

    # Prepare the output strings
    out_head, out_value, out_mqc = ("Sample", "'{}'".format(sname), mqc_main)
    for ft, pt in fpercent.items():
        out_head = "{}\tpercent_{}".format(out_head, ft)
        out_value = "{}\t{}".format(out_value, pt)
        out_mqc = "{}\n{}".format(out_mqc, mqc_pconf.format(ft=ft))

    # Write the output to a file
    with open(outfile, "w") as ofl:
        out_final = "\n".join([out_mqc, out_head, out_value]).strip()
        ofl.write(out_final + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Calculate features percentage for biotype counts""")
    parser.add_argument("biocount", type=str, help="File with all biocounts")
    parser.add_argument(
        "-f",
        "--features",
        dest="features",
        required=True,
        nargs="+",
        help="Features to count",
    )
    parser.add_argument("-s", "--sample", dest="sample", type=str, help="Sample Name")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="biocount_percent.tsv",
        type=str,
        help="Sample Name",
    )
    args = parser.parse_args()
    mqc_feature_stat(args.biocount, args.features, args.output, args.sample)
