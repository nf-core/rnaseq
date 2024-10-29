#!/usr/bin/env python3

# Written by Pranathi Vemuri, later modified by Jonathan Manning and released under the MIT license.

import logging
import os
import platform
from itertools import groupby
from typing import Iterator, Tuple


def setup_logging() -> logging.Logger:
    """Configure logging for the script.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.INFO)
    return logger


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def parse_fasta(fasta_file: str) -> Iterator[Tuple[str, str]]:
    """Parse a fasta file and yield tuples of header and sequence.

    Args:
        fasta_file (str): Path to the fasta file.

    Yields:
        Iterator[Tuple[str, str]]: Tuples of header and sequence from the fasta file.


    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence

    Fasta iterator from https://www.biostars.org/p/710/#120760
    """
    with open(fasta_file) as file_handle:
        fasta_iter = (x[1] for x in groupby(file_handle, lambda line: line[0] == ">"))
        for header in fasta_iter:
            header_str = next(header)[1:].strip()
            sequence = "".join(s.strip() for s in next(fasta_iter))
            yield (header_str, sequence)


def fasta_to_gtf(fasta: str, output_file: str, biotype: str) -> None:
    """
    Read a fasta file and create a GTF file.

    Args:
        fasta (str): Path to the fasta file.
        output_file (str): Path for the output GTF file.
        biotype (str): The biotype to use in the GTF.
    """
    fasta_iter = parse_fasta(fasta)
    lines = []

    for header, sequence in fasta_iter:
        seq_name = header.split()[0].replace(" ", "_")
        line = generate_gtf_line(seq_name, len(sequence), biotype)
        lines.append(line)

    with open(output_file, "w") as file_handle:
        file_handle.writelines(lines)


def generate_gtf_line(name: str, length: int, biotype: str) -> str:
    """Generate a single GTF line given sequence name, length, and biotype.

    Args:
        name (str): Name of the sequence.
        length (int): Length of the sequence.
        biotype (str): Biotype of the sequence.

    Returns:
        str: A formatted GTF line.
    """
    biotype_attr = f' {biotype} "transgene";' if biotype else ""
    attributes = f'exon_id "{name}.1"; exon_number "1";{biotype_attr} gene_id "{name}_gene"; gene_name "{name}_gene"; gene_source "custom"; transcript_id "{name}_gene"; transcript_name "{name}_gene";\\n'
    return f"{name}\\ttransgene\\texon\\t1\\t{length}\\t.\\t+\\t.\\t{attributes}"


def main() -> None:
    # Parse arguments using argparse (not shown for brevity)
    # Example: args = parser.parse_args()

    logger = setup_logging()
    logger.info("Starting fasta to GTF conversion.")

    # Add fasta lines to GTF
    add_name = os.path.splitext(os.path.basename("$add_fasta"))[0]
    fasta_to_gtf("$add_fasta", f"{add_name}.gtf", "$biotype")

    # Concatenate new fasta to existing fasta, and the GTF we just generated to the GTF
    genome_name = "$params.genome" if "$params.genome" != "null" else os.path.splitext(os.path.basename("$fasta"))[0]
    output_prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else f"{genome_name}_{add_name}"

    os.mkdir("out")
    os.system(f"cat $fasta $add_fasta > out/{output_prefix}.fasta")
    os.system(f"cat $gtf {add_name}.gtf > out/{output_prefix}.gtf")

    logger.info("Conversion completed successfully.")

    # Write the versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))


if __name__ == "__main__":
    main()
