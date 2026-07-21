#!/usr/bin/env python

# Written by Olga Botvinnik with subsequent reworking by Jonathan Manning and Nico Trummer.

# MIT License

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import gzip
import logging
import os
import platform
import re
import statistics
from typing import Optional, Set

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("fasta_gtf_filter")
logger.setLevel(logging.INFO)

# Alias list for GTF column 9 biotype keys across different data providers
BIOTYPE_ALIASES = [
    "gene_biotype",
    "gene_type",
    "biotype",
    "transcript_biotype",
    "transcript_type",
    "so_term_name",
    "Ontology_term",
]


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string."""
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\n"
    return yaml_str


def parse_attributes(attr_str: str) -> dict:
    """Parses GTF column 9 key-value string into a clean Python dictionary."""
    attrs = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            parts = item.split(" ", 1)
            attrs[parts[0].strip()] = parts[1].strip().strip('"')
        elif "=" in item:
            parts = item.split("=", 1)
            attrs[parts[0].strip()] = parts[1].strip().strip('"')
    return attrs


def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extracts the sequence names from a FASTA file."""
    is_gz = fasta_name.endswith(".gz")
    open_fn = gzip.open if is_gz else open

    with open_fn(fasta_name) as fasta:
        sequences = set()
        for line in fasta:
            line = line.decode("utf-8") if is_gz else line
            if line.startswith(">"):
                sequences.add(line[1:].split(None, 1)[0])

        return sequences


def tab_delimited(file: str) -> float:
    """Check if file is tab-delimited and return median number of tabs."""
    with open(file) as f:
        data = f.read(102400)
        return statistics.median(line.count("\t") for line in data.split("\n"))


def filter_gtf(fasta: Optional[str], gtf_in: str, filtered_gtf_out: str, skip_transcript_id_check: bool) -> None:
    """Filter GTF file based on FASTA sequence names and propagate biotypes."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")

    seq_names_in_genome = None
    if fasta and os.path.isfile(fasta):
        seq_names_in_genome = extract_fasta_seq_names(fasta)
        logger.info(f"Extracted chromosome sequence names from {fasta}")
        logger.debug("All sequence IDs from FASTA: " + ", ".join(sorted(seq_names_in_genome)))

    is_gz = gtf_in.endswith(".gz")
    open_fn = gzip.open if is_gz else open

    # -------------------------------------------------------------
    # PASS 1: Build parent biotype mapping (gene_id -> biotype)
    # -------------------------------------------------------------
    gene_biotype_map = {}
    try:
        with open_fn(gtf_in) as gtf:
            for line in gtf:
                line = line.decode("utf-8") if is_gz else line
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                attrs = parse_attributes(fields[8])
                gene_id = attrs.get("gene_id")
                if gene_id and gene_id not in gene_biotype_map:
                    for alias in BIOTYPE_ALIASES:
                        if alias in attrs:
                            gene_biotype_map[gene_id] = attrs[alias]
                            break
    except Exception as e:
        logger.warning(f"Pre-scan for parent biotypes encountered an issue: {e}")

    # -------------------------------------------------------------
    # PASS 2: Filter lines and inject missing gene_biotype
    # -------------------------------------------------------------
    seq_names_in_gtf = set()
    try:
        with open_fn(gtf_in) as gtf, open_fn(filtered_gtf_out, "wb" if is_gz else "w") as out:
            line_count = 0
            for line in gtf:
                line = line.decode("utf-8") if is_gz else line

                if line.startswith("#"):
                    out.write(line.encode() if is_gz else line)
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    out.write(line.encode() if is_gz else line)
                    continue

                seq_name = fields[0]
                seq_names_in_gtf.add(seq_name)

                if seq_names_in_genome is None or seq_name in seq_names_in_genome:
                    if skip_transcript_id_check or re.search(r'transcript_id "([^"]+)"', line):
                        # Ensure gene_biotype exists on this line
                        attrs = parse_attributes(fields[8])
                        if "gene_biotype" not in attrs:
                            found_biotype = None

                            # 1. Check for synonymous key on current line
                            for alias in BIOTYPE_ALIASES:
                                if alias in attrs:
                                    found_biotype = attrs[alias]
                                    break

                            # 2. Inherit from parent gene_id if absent
                            if not found_biotype and attrs.get("gene_id") in gene_biotype_map:
                                found_biotype = gene_biotype_map[attrs["gene_id"]]

                            # Append gene_biotype attribute
                            if found_biotype:
                                attr_str = fields[8].strip()
                                if attr_str and not attr_str.endswith(";"):
                                    attr_str += ";"
                                fields[8] = f'{attr_str} gene_biotype "{found_biotype}";'
                                line = "\t".join(fields) + "\n"

                        out.write(line.encode() if is_gz else line)
                        line_count += 1

            if line_count == 0:
                raise ValueError("All GTF lines removed by filters")

    except OSError as e:
        logger.error(f"File operation failed: {e}")
        return

    logger.debug("All sequence IDs from GTF: " + ", ".join(sorted(seq_names_in_gtf)))
    logger.info(f"Extracted {line_count} matching sequences from {gtf_in} into {filtered_gtf_out}")


parser = argparse.ArgumentParser()
parser.add_argument("--skip_transcript_id_check", action="store_true", default=False)
parsed_args = parser.parse_args("${args}".split() if "${args}".strip() else [])

filter_gtf("${fasta}", "${gtf}", "${prefix}.${suffix}", parsed_args.skip_transcript_id_check)

# Versions

versions = {"${task.process}": {"python": platform.python_version()}}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))