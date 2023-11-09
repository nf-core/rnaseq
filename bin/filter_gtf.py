#!/usr/bin/env python
from __future__ import print_function
import logging
from itertools import groupby
import argparse

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def is_header(line: str) -> bool:
    """Returns True if the given line is a header line in a FASTA file."""
    return line[0] == ">"


def extract_fasta_seq_names(fasta_name: str) -> set:
    """Extracts the sequence names from a FASTA file.

    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    from https://www.biostars.org/p/710/

    Args:
      fasta_name: The path to the FASTA file.

    Returns:
      A set of the sequence names in the FASTA file.
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


def extract_genes_in_genome(fasta: str, gtf_in: str, prefix: str) -> None:
    """Extracts the genes in the genome from a GTF file.

    Args:
      fasta: The path to the FASTA file.
      gtf_in: The path to the input GTF file.
      prefix: Prefix for output GTF
    """
    gtf_out = prefix + "_in_genome.gtf"
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


def remove_features_without_transcript_id(gtf_in, prefix):
    """
    Removes gene rows with absent or empty transcript_id attributes from a GTF file.

    Args:
      gtf_in: Path to the input GTF file.
      prefix: Path to the output GTF file.
    """
    gtf_out = prefix + "_with_transcript_ids.gtf"

    with open(gtf_in, "r") as f_in, open(gtf_out, "w") as f_out:
        for line in f_in:
            transcript_id = line.split("\t")[8].split(" transcript_id ")[1].split(";")[0].replace('"', "")
            if transcript_id and transcript_id.isalnum():
                f_out.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Filter GTF for various reasons""")
    parser.add_argument("--gtf", type=str, help="GTF file")
    parser.add_argument("--fasta", type=str, help="Genome fasta file")
    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        default="genes",
        type=str,
        help="Prefix for output GTF files",
    )

    args = parser.parse_args()
    extract_genes_in_genome(args.fasta, args.gtf, args.prefix)
    remove_features_without_transcript_id(args.gtf, args.prefix)
