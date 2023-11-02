#!/usr/bin/env python
import logging
import argparse
import glob
import os
from typing import Set, Optional, Counter, Dict, List
from collections import defaultdict

# Configure logging
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def read_top_transcripts(quant_dir: str, file_pattern: str) -> Set[str]:
    """
    Read top transcripts from quantification file.

    :param quant_dir: Directory containing quantification files.
    :param file_pattern: File pattern to match the quantification file.
    :return: A set of transcript names.
    """
    transcripts = set()
    try:
        # Find the first file that matches the pattern
        file_name = glob.glob(os.path.join(quant_dir, '*', file_pattern))[0]
    except IndexError:
        logger.error('No quantification files found.')
        raise FileNotFoundError('Quantification file not found.')

    with open(file_name, 'r') as file_handle:
        for line in file_handle:
            if line.startswith('Name'):
                continue
            transcript = line.split()[0]
            transcripts.add(transcript)
            if len(transcripts) > 100:
                break

    logger.info(f'Transcripts found: {transcripts}')
    return transcripts


def map_transcripts_to_gene(quant_type: str, gtf_file: str, quant_dir: str, gene_id: str, extra_id_field: str, output_file: str) -> Optional[bool]:
    """
    Map transcripts to gene names and write the output to a file.

    :param quant_type: Type of quantification (e.g., 'kallisto' or 'salmon').
    :param gtf_file: Path to the GTF file.
    :param quant_dir: Directory containing quantification files.
    :param gene_id: Gene identifier in the GTF file.
    :param extra_id_field: Additional identifier in the GTF file.
    :param output_file: Path to the output file.
    :return: True if mapping is successful, None otherwise.
    """
    file_pattern = 'quant.sf' if quant_type == 'salmon' else 'abundance.tsv'
    transcripts = read_top_transcripts(quant_dir, file_pattern)

    votes: Counter = Counter()
    gene_dict: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    with open(gtf_file, 'r') as file_handle:
        for line in file_handle:
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            attributes = {item.split()[0]: item.split()[1].strip('"') for item in columns[8].split(';') if item}

            if gene_id in attributes and attributes[gene_id] in transcripts:
                votes[attributes[gene_id]] += 1
                gene_dict[attributes[gene_id]].append(attributes)

    if not votes:
        logger.warning('No matching attributes found in GTF.')
        return None

    # Get the most common transcript ID
    top_transcript_id = votes.most_common(1)[0][0]
    logger.info(f'Transcript attributed to be gene: {top_transcript_id}')

    seen = set()
    with open(output_file, 'w') as output_handle:
        for gene, attrs in gene_dict.items():
            for attr in attrs:
                if top_transcript_id not in attr:
                    continue
                if (gene, attr[top_transcript_id]) not in seen:
                    seen.add((gene, attr[top_transcript_id]))
                    extra_id = attr.get(extra_id_field, gene)
                    output_handle.write(f'{attr[top_transcript_id]}\t{gene}\t{extra_id}\n')

    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Map transcripts to gene names for tximport.')
    parser.add_argument('--quant_type', type=str, help='Quantification type', default='salmon')
    parser.add_argument('--gtf', type=str, help='GTF file', required=True)
    parser.add_argument('--quants', type=str, help='Output of quantification', required=True)
    parser.add_argument('--id', type=str, help='Gene ID in the GTF file', required=True)
    parser.add_argument('--extra', type=str, help='Extra ID in the GTF file')
    parser.add_argument('-o', '--output', dest='output', default='tx2gene.tsv', type=str, help='File with output')

    args = parser.parse_args()
    map_transcripts_to_gene(args.quant_type, args.gtf, args.quants, args.id, args.extra, args.output)

