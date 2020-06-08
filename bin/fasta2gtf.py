#!/usr/bin/env python
"""
Read a custom fasta file and create a custom GTF containing each entry
"""
import argparse
from itertools import groupby
import logging

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence

    Fasta iterator from https://www.biostars.org/p/710/#120760
    """
    with open(fasta_name) as fh:
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            headerStr = header.__next__()[1:].strip()

            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.__next__())

            yield (headerStr, seq)


def fasta2gtf(fasta, output):
    fiter = fasta_iter(fasta)
    # GTF output lines
    lines = []
    attributes = \
        'gene_id "{name_sanitized}"; gene_name "{name_sanitized}";transcript_id "{name_sanitized}"; gene_biotype "{name_sanitized}"; gene_type "{name_sanitized}"\n'
    line_template = \
        "{name_sanitized}\ttransgene\texon\t1\t{length}\t.\t+\t.\t" + attributes

    for ff in fiter:
        name, seq = ff
        # Use first ID as separated by spaces as the "sequence name"
        # (equivalent to "chromosome" in other cases)
        seqname = name.split()[0]
        # Remove all spaces
        name_sanitized = seqname.replace(' ', '_')
        length = len(seq)
        line = line_template.format(
            name_sanitized=name_sanitized, length=length)
        lines.append(line)

    with open(output, 'w') as f:
        f.write(''.join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Convert a custom fasta (e.g. transgene)
        to a GTF annotation.""")
    parser.add_argument("fasta", type=str, help="Custom transgene sequence")
    parser.add_argument(
        "-o", "--output", dest='output',
        default='transgenes.gtf', type=str, help="gene annotation GTF output")
    args = parser.parse_args()
    fasta2gtf(args.fasta, args.output)
