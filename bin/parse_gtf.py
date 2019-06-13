#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict, Counter
import logging
import argparse
import glob
import os

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def read_top_transcript(salmon):
    txs = set()
    fn = glob.glob(os.path.join(salmon, "*", "quant.sf"))[1]
    with open(fn) as inh:
        for line in inh:
            if line.startswith("Name"):
                continue
            txs.add(line.split()[0])
            if len(txs) > 100:
                break
    logger.info("Transcripts found in FASTA: %s" % txs)
    return txs


def tx2gene(gtf, salmon, gene_id, extra, out):
    txs = read_top_transcript(salmon)
    votes = Counter()
    gene_dict = defaultdict(dict)
    with open(gtf) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            attr_dict = OrderedDict()
            for gff_item in cols[8].split(";"):
                item_pair = gff_item.strip().split(" ")
                if len(item_pair) > 1:
                    value = item_pair[1].strip().replace("\"", "")
                    if value in txs:
                        votes[item_pair[0].strip()] += 1

                    attr_dict[item_pair[0].strip()] = value
            gene_dict[attr_dict[gene_id]] = attr_dict

    if not votes:
        logger.warning("No attribute in GTF matching transcripts")
        return None

    txid = votes.most_common(1)[0][0]
    logger.info("Attributed found to be transcript: %s" % txid)
    with open(out, 'w') as outh:
        for gene in gene_dict:
            print("%s,%s,%s" % (gene_dict[gene][txid], gene, gene_dict[gene][extra]), file=outh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Get tx to gene names for tximport""")
    parser.add_argument("--gtf", type=str, help="GTF file")
    parser.add_argument("--salmon", type=str, help="output of salmon")
    parser.add_argument("--id", type=str, help="gene id in the gtf file")
    parser.add_argument("--extra", type=str, help="extra id in the gtf file")
    parser.add_argument("-o", "--output", dest='output', default='tx2gene.csv', type=str, help="file with output")

    args = parser.parse_args()
    tx2gene(args.gtf, args.salmon, args.id, args.extra, args.output)
