#!/usr/bin/env python
from __future__ import print_function
import logging
from itertools import groupby
import argparse
import re

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
    seqnames = set()
    with open(fasta_name) as fasta:
        for line in fasta:
            if line[0] == ">":
                seqnames.add(line[1:].split(None, 1)[0])
    return seqnames


def filter_gtf(fasta: str, gtf_in: str, gtf_in_genome_out: str, gtf_transcript_out: str) -> None:
    """Extracts the genes in the genome from a GTF file.

    Args:
      fasta: The path to the FASTA file.
      gtf_in: The path to the input GTF file.
      gtf_in_genome_out: The path to the output GTF file with scaffolds filtered.
      gtf_transcript_out: The path to the output GTF file with scaffolds
          filtered and rows without transcript ID removed.

    Raises:
      ValueError: If no overlap is found or if the GTF file is not tab delimited.
    """

    def tab_delimited(file) -> float:
    
        import statistics. # put to the top

        with open(file, "r") as f:
            data = f.read(1024)
            lines = data.split("\n")
            # most lines should have 9 tab-separated columns
            return statistics.median([line.count("\t") for line in lines])

    num_sep = tab_delimited(gtf_in)
    if num_sep != 8:
        raise ValueError("No valid tab-delimited GTF file.")

    seq_names_in_genome = extract_fasta_seq_names(fasta)
    logger.info(f"Extracted chromosome sequence names from {fasta}")
    logger.debug("All chromosome names: " + ", ".join(sorted(seq_names_in_genome)))

    with open(gtf_in) as gtf, open(gtf_in_genome_out, "w") as out, open(
        gtf_transcript_out, "w"
    ) as out2:
        line_count_all = 0
        line_count_transcript = 0
        for line in gtf:
            if (
                line.split("\t")[0] in seq_names_in_genome
                and line.count("\t") == num_sep
            ):
                out.write(line)
                line_count_all += 1
                if re.search(r'transcript_id "([^"]+)"', line):
                    out2.write(line)
                    line_count_transcript += 1
        if line_count_all == 0:
            raise ValueError("No overlapping scaffolds found.")



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
    filter_gtf(args.fasta, args.gtf, args.prefix + "_in_genome.gtf", args.prefix + "_with_transcript_ids.gtf")
