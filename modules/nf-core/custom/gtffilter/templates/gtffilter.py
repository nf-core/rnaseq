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


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string."""
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extracts sequence names from a FASTA file."""
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
    """Check if file is tab-delimited and return median number of tabs (supports gzip)."""
    is_gz = file.endswith(".gz")
    open_fn = gzip.open if is_gz else open

    with open_fn(file) as f:
        data = f.read(102400)
        data = data.decode("utf-8") if is_gz else data
        lines = [line for line in data.split("\n") if line.strip()]
        return statistics.median(line.count("\t") for line in lines)


def filter_gtf(fasta: Optional[str], gtf_in: str, filtered_gtf_out: str, skip_transcript_id_check: bool) -> None:
    """Filter GTF file based on FASTA sequence names and flexible attribute rules."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")

    seq_names_in_genome = None
    if fasta and os.path.isfile(fasta):
        seq_names_in_genome = extract_fasta_seq_names(fasta)
        logger.info(f"Extracted chromosome sequence names from {fasta}")
        logger.debug("All sequence IDs from FASTA: " + ", ".join(sorted(seq_names_in_genome)))

    seq_names_in_gtf = set()
    
    # Flexible pattern matching for transcript_id across Ensembl ("..."), RefSeq/NCBI (=...), or unquoted formats
    transcript_id_pattern = re.compile(r'transcript_id[\s=]+["\']?([^"\';\s]+)')

    try:
        is_gz = gtf_in.endswith(".gz")
        open_fn = gzip.open if is_gz else open
        with open_fn(gtf_in) as gtf, open_fn(filtered_gtf_out, "wb" if is_gz else "w") as out:
            line_count = 0
            for line in gtf:
                line = line.decode("utf-8") if is_gz else line
                
                # 1. Preserve GTF header/comment lines starting with '#'
                if line.startswith("#"):
                    out.write(line.encode() if is_gz else line)
                    continue

                fields = line.split("\\t")
                if len(fields) < 9:
                    continue

                seq_name = fields[0]
                feature_type = fields[2]
                seq_names_in_gtf.add(seq_name)

                # 2. Filter by sequence name (if FASTA sequence headers provided)
                if seq_names_in_genome is None or seq_name in seq_names_in_genome:
                    has_transcript_id = bool(transcript_id_pattern.search(line))
                    is_gene_feature = (feature_type == "gene")

                    # 3. Retain feature if:
                    # - skip_transcript_id_check is set, OR
                    # - feature is a top-level 'gene' line, OR
                    # - line contains a valid transcript_id attribute
                    if skip_transcript_id_check or is_gene_feature or has_transcript_id:
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