#!/usr/bin/env python

# Written by Lorena Pantano with subsequent reworking by Jonathan Manning. Released under the MIT license.

import logging
import argparse
import glob
import os
import re
from collections import Counter, defaultdict, OrderedDict
from collections.abc import Set
from typing import Dict

# Configure logging
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def read_top_transcripts(quant_dir: str, file_pattern: str) -> Set[str]:
    """
    Read the top 100 transcripts from the quantification file.

    Parameters:
    quant_dir (str): Directory where quantification files are located.
    file_pattern (str): Pattern to match quantification files.

    Returns:
    set: A set containing the top 100 transcripts.
    """
    try:
        # Find the quantification file within the directory
        quant_file_path = glob.glob(os.path.join(quant_dir, "*", file_pattern))[0]
        with open(quant_file_path, "r") as file_handle:
            # Read the file and extract the top 100 transcripts
            return {line.split()[0] for i, line in enumerate(file_handle) if i > 0 and i <= 100}
    except IndexError:
        # Log an error and raise a FileNotFoundError if the quant file does not exist
        logger.error("No quantification files found.")
        raise FileNotFoundError("Quantification file not found.")


def discover_transcript_attribute(gtf_file: str, transcripts: Set[str]) -> str:
    """
    Discover the attribute in the GTF that corresponds to transcripts, prioritizing 'transcript_id'.

    Parameters:
    gtf_file (str): Path to the GTF file.
    transcripts (Set[str]): A set of transcripts to match in the GTF file.

    Returns:
    str: The attribute name that corresponds to transcripts in the GTF file.
    """

    votes = Counter()
    with open(gtf_file) as inh:
        # Read GTF file, skipping header lines
        for line in filter(lambda x: not x.startswith("#"), inh):
            cols = line.split("\t")

            # Use regular expression to correctly split the attributes string
            attributes_str = cols[8]
            attributes = dict(re.findall(r'(\S+) "(.*?)(?<!\\)";', attributes_str))

            votes.update(key for key, value in attributes.items() if value in transcripts)

    if not votes:
        # Log a warning if no matching attribute is found
        logger.warning("No attribute in GTF matching transcripts")
        return ""

    # Check if 'transcript_id' is among the attributes with the highest votes
    if "transcript_id" in votes and votes["transcript_id"] == max(votes.values()):
        logger.info("Attribute 'transcript_id' corresponds to transcripts.")
        return "transcript_id"

    # If 'transcript_id' isn't the highest, determine the most common attribute that matches the transcripts
    attribute, _ = votes.most_common(1)[0]
    logger.info(f"Attribute '{attribute}' corresponds to transcripts.")
    return attribute


def parse_attributes(attributes_text: str) -> Dict[str, str]:
    """
    Parse the attributes column of a GTF file.

    :param attributes_text: The attributes column as a string.
    :return: A dictionary of the attributes.
    """
    # Split the attributes string by semicolon and strip whitespace
    attributes = attributes_text.strip().split(";")
    attr_dict = OrderedDict()

    # Iterate over each attribute pair
    for attribute in attributes:
        # Split the attribute into key and value, ensuring there are two parts
        parts = attribute.strip().split(" ", 1)
        if len(parts) == 2:
            key, value = parts
            # Remove any double quotes from the value
            value = value.replace('"', "")
            attr_dict[key] = value

    return attr_dict


def map_transcripts_to_gene(
    quant_type: str, gtf_file: str, quant_dir: str, gene_id: str, extra_id_field: str, output_file: str
) -> bool:
    """
    Map transcripts to gene names and write the output to a file.

    Parameters:
    quant_type (str): The quantification method used (e.g., 'salmon').
    gtf_file (str): Path to the GTF file.
    quant_dir (str): Directory where quantification files are located.
    gene_id (str): The gene ID attribute in the GTF file.
    extra_id_field (str): Additional ID field in the GTF file.
    output_file (str): The output file path.

    Returns:
    bool: True if the operation was successful, False otherwise.
    """
    # Read the top transcripts based on quantification type
    transcripts = read_top_transcripts(quant_dir, "quant.sf" if quant_type == "salmon" else "abundance.tsv")
    # Discover the attribute that corresponds to transcripts in the GTF
    transcript_attribute = discover_transcript_attribute(gtf_file, transcripts)

    if not transcript_attribute:
        # If no attribute is found, return False
        return False

    # Open GTF and output file to write the mappings
    # Initialize the set to track seen combinations
    seen = set()

    with open(gtf_file) as inh, open(output_file, "w") as output_handle:
        # Parse each line of the GTF, mapping transcripts to genes
        for line in filter(lambda x: not x.startswith("#"), inh):
            cols = line.split("\t")
            attr_dict = parse_attributes(cols[8])
            if gene_id in attr_dict and transcript_attribute in attr_dict:
                # Create a unique identifier for the transcript-gene combination
                transcript_gene_pair = (attr_dict[transcript_attribute], attr_dict[gene_id])

                # Check if the combination has already been seen
                if transcript_gene_pair not in seen:
                    # If it's a new combination, write it to the output and add to the seen set
                    extra_id = attr_dict.get(extra_id_field, attr_dict[gene_id])
                    output_handle.write(f"{attr_dict[transcript_attribute]}\t{attr_dict[gene_id]}\t{extra_id}\n")
                    seen.add(transcript_gene_pair)

    return True


# Main function to parse arguments and call the mapping function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map transcripts to gene names for tximport.")
    parser.add_argument("--quant_type", type=str, help="Quantification type", default="salmon")
    parser.add_argument("--gtf", type=str, help="GTF file", required=True)
    parser.add_argument("--quants", type=str, help="Output of quantification", required=True)
    parser.add_argument("--id", type=str, help="Gene ID in the GTF file", required=True)
    parser.add_argument("--extra", type=str, help="Extra ID in the GTF file")
    parser.add_argument("-o", "--output", dest="output", default="tx2gene.tsv", type=str, help="File with output")

    args = parser.parse_args()
    if not map_transcripts_to_gene(args.quant_type, args.gtf, args.quants, args.id, args.extra, args.output):
        logger.error("Failed to map transcripts to genes.")
