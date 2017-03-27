#!/usr/bin/env python

import argparse
import os
import re
import logging
from collections import defaultdict

def merge_featureCounts(dest_dir,out_file,input_files):
   logger = logging.getLogger(__name__)
   logger.addHandler(logging.StreamHandler())
   logger.setLevel(logging.INFO)
   
   table_dict=defaultdict(dict)
   sample_names=[]
   genes=set()
   for input_file in input_files:
       logger.info("Reading from {}".format(input_file))
       sample_name=os.path.basename(input_file)
       sample_name= sample_name.replace("Aligned.sortedByCoord.out_gene.featureCounts.txt","")
       sample_names.append(sample_name)
       table_dict[sample_name]=dict()
     
       with open(input_file, 'r') as f:
           f.readline()
           f.readline()
           for line in f:
               #save the genes to a list for the first file
               line_info=line.split('\t')
               gene=line_info[0]
               genes.add(gene)
               try:
                   gene_count = int(line_info[-1].rstrip())
                   table_dict[sample_name][gene] = gene_count
               except TypeError:
                       logger.warning("Detected discrepancy in {}  line {}".format(input_file, line))
                
               table_dict[sample_name][gene]=gene_count

   #write Output
   logger.info("Writing to file {}".format(out_file))
   with open(out_file, 'w') as f:
       #Generate header
       line_to_write="ENSEMBL_ID"
       sample_names.sort()
       for sample_name in sample_names:
           line_to_write+=('\t{}'.format(sample_name))
       line_to_write += "\n"
       f.write(line_to_write)    
       #Write the rest of the lines
       gene_list=list(genes)
       gene_list.sort()
       for gene in gene_list:
           line_to_write=gene
           for sample_name in sample_names:
               try:
                   line_to_write += ('\t{}'.format(table_dict[sample_name][gene]))
               except KeyError:
                    #Missing gene in one of the filess.
                    line_to_write += ('\tNa')
           line_to_write += "\n"
           f.write(line_to_write)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples in a project
    """)
    parser.add_argument("-d", "--dest_dir", dest='dest_dir', default='.',
                                   help="Path to output.")
    parser.add_argument("-o", "--results_file_name",dest='out_file', default='all_counts.txt',
                                   help= "Name of the output file that will be created")
    parser.add_argument("-i", "--input_files", metavar='<input_files>', nargs='+', default='*.featureCounts.txt',
                                   help="Path to the outputfiles from FeatureCounts. ")
    args = parser.parse_args()
    merge_featureCounts(args.dest_dir, args.out_file, args.input_files)
