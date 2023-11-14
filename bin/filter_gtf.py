#!/usr/bin/env python
import logging
import argparse
import re
import statistics
from typing import Set

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("fasta_gtf_filter")
logger.setLevel(logging.INFO)

def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extracts the sequence names from a FASTA file."""
    with open(fasta_name) as fasta:
        return {line[1:].split(None, 1)[0] for line in fasta if line.startswith(">")}

def tab_delimited(file: str) -> float:
    """Check if file is tab-delimited and return median number of tabs."""
    with open(file, "r") as f:
        data = f.read(1024)
        return statistics.median(line.count("\t") for line in data.split("\n"))

def filter_gtf(fasta: str, gtf_in: str, gtf_in_genome_out: str, gtf_transcript_out: str) -> None:
    """Filter GTF file based on FASTA sequence names."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 8 tab-separated columns.")

    seq_names_in_genome = extract_fasta_seq_names(fasta)
    logger.info(f"Extracted chromosome sequence names from {fasta}")

    try:
        with open(gtf_in) as gtf, open(gtf_in_genome_out, "w") as out, open(gtf_transcript_out, "w") as out2:
            line_count_all, line_count_transcript = 0, 0
            for line in gtf:
                if line.split("\t")[0] in seq_names_in_genome:
                    out.write(line)
                    line_count_all += 1
                    if re.search(r'transcript_id "([^"]+)"', line):
                        out2.write(line)
                        line_count_transcript += 1
            if line_count_all == 0:
                raise ValueError("No overlapping scaffolds found.")
    except IOError as e:
        logger.error(f"File operation failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filters a GTF file based on sequence names in a FASTA file.")
    parser.add_argument("--gtf", type=str, required=True, help="GTF file")
    parser.add_argument("--fasta", type=str, required=True, help="Genome fasta file")
    parser.add_argument("--prefix", dest="prefix", default="genes", type=str, help="Prefix for output GTF files")

    args = parser.parse_args()
    filter_gtf(args.fasta, args.gtf, args.prefix + "_in_genome.gtf", args.prefix + "_with_transcript_ids.gtf")

