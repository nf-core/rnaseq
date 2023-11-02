#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict, Counter
import logging
import argparse
import glob
import os

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def read_salmon_top_transcript(salmon):
    txs = set()
    fn = glob.glob(os.path.join(salmon, "*", "quant.sf"))[0]
    with open(fn) as inh:
        for line in inh:
            if line.startswith("Name"):
                continue
            txs.add(line.split()[0])
            if len(txs) > 100:
                break
    logger.info("Transcripts found in FASTA: %s" % txs)
    return txs

def read_kallisto_top_transcript(kallisto):
    txs = set()
    fn = glob.glob(os.path.join(kallisto, "*", "abundance.tsv"))[0]
    with open(fn) as inh:
        next(inh)  # Skip header line
        for line in inh:
            txs.add(line.split('\t')[0])
            if len(txs) > 100:
                break
    logger.info("Transcripts found in FASTA: %s" % txs)
    return txs

def tx2gene(quant_type, gtf, quants, gene_id, extra, out):
    if quant_type == 'kallisto':
        txs = read_kallisto_top_transcript(quants)
    else:
        txs = read_salmon_top_transcript(quants)
    
    votes = Counter()
    gene_dict = defaultdict(list)
    with open(gtf) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            attr_dict = OrderedDict()
            for gff_item in cols[8].split(";"):
                item_pair = gff_item.strip().split(" ")
                if len(item_pair) > 1:
                    value = item_pair[1].strip().replace('"', "")
                    if value in txs:
                        votes[item_pair[0].strip()] += 1

                    attr_dict[item_pair[0].strip()] = value
            try:
                gene_dict[attr_dict[gene_id]].append(attr_dict)
            except KeyError:
                continue

    if not votes:
        logger.warning("No attribute in GTF matching transcripts")
        return None

    txid = votes.most_common(1)[0][0]
    logger.info("Attributed found to be transcript: %s" % txid)
    seen = set()
    with open(out, "w") as outh:
        for gene in gene_dict:
            for row in gene_dict[gene]:
                if txid not in row:
                    continue
                if (gene, row[txid]) not in seen:
                    seen.add((gene, row[txid]))
                    if not extra in row:
                        extra_id = gene
                    else:
                        extra_id = row[extra]
                    print("%s\t%s\t%s" % (row[txid], gene, extra_id), file=outh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Get tx to gene names for tximport""")
    parser.add_argument("--quant_type", type=str, help="Quantification type", default = 'salmon')
    parser.add_argument("--gtf", type=str, help="GTF file")
    parser.add_argument("--quants", type=str, help="output of quantification")
    parser.add_argument("--id", type=str, help="gene id in the gtf file")
    parser.add_argument("--extra", type=str, help="extra id in the gtf file")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="tx2gene.tsv",
        type=str,
        help="file with output",
    )

    args = parser.parse_args()
    tx2gene(args.quant_type, args.gtf, args.quants, args.id, args.extra, args.output)
