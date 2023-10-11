#!/usr/bin/env python
from __future__ import print_function
import logging
from itertools import groupby
import argparse

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def is_header(line):
    return line[0] == ">"


def extract_fasta_seq_names(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    from https://www.biostars.org/p/710/
    """
    # first open the file outside
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, is_header))

    for i, header in enumerate(faiter):
        line = next(header)
        if is_header(line):
            # drop the ">"
            headerStr = line[1:].strip().split()[0]
        yield headerStr


def extract_genes_in_genome(fasta, gtf_in, gtf_out):
    seq_names_in_genome = set(extract_fasta_seq_names(fasta))
    logger.info("Extracted chromosome sequence names from : %s" % fasta)
    logger.info("All chromosome names: " + ", ".join(sorted(x for x in seq_names_in_genome)))
    seq_names_in_gtf = set([])

    n_total_lines = 0
    n_lines_in_genome = 0
    with open(gtf_out, "w") as f:
        with open(gtf_in) as g:
            for line in g.readlines():
                n_total_lines += 1
                seq_name_gtf = line.split("\t")[0]
                seq_names_in_gtf.add(seq_name_gtf)
                if seq_name_gtf in seq_names_in_genome:
                    n_lines_in_genome += 1
                    f.write(line)
    logger.info(
        "Extracted %d / %d lines from %s matching sequences in %s" % (n_lines_in_genome, n_total_lines, gtf_in, fasta)
    )
    logger.info("All sequence IDs from GTF: " + ", ".join(sorted(x for x in seq_name_gtf)))

    logger.info("Wrote matching lines to %s" % gtf_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Filter GTF only for features in the genome""")
    parser.add_argument("--gtf", type=str, help="GTF file")
    parser.add_argument("--fasta", type=str, help="Genome fasta file")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="genes_in_genome.gtf",
        type=str,
        help="GTF features on fasta genome sequences",
    )

    args = parser.parse_args()
    extract_genes_in_genome(args.fasta, args.gtf, args.output)
