#!/usr/bin/env/python

import argparse
import os
import re
import logging
from collections import defaultdict


def merge_featureCounts(dest_dir,out_file,input_files):
   table_dict=defaultdict(dict)
   sample_names=[]
   gene_list=[]
   for input_file in input_files:
       if not gene_list:
           first_file=True
       print"Reading from {}".format(input_file)
       sample_name=os.path.basename(input_file)
       sample_name= sample_name.replace("Aligned.sortedByCoord.out_gene.featureCounts.txt","")
       sample_names.append(sample_name)
       table_dict[sample_name]=dict()
     
       with open(input_file, 'r') as f:
           for line in f:
               if not line.startswith('E'):
                   continue
               #save the genes to a list for the first file
               if first_file:
                   gene_list.append(line.split('\t')[0])
               line_info=line.split('\t')
               gene=line_info[0]
               try:
                   gene_count = int(line_info[-1].rstrip())
                   table_dict[sample_name][gene] = gene_count
               except TypeError:
                       logging.info("Detected discrepancy in {}  line {}".format(input_file, line))
                
               table_dict[sample_name][gene]=gene_count
               #Done with first file
               first_file=False

   #write Output
   logging.info("Writing to file {}".format(out_file))
   with open(out_file, 'w') as f:
       #Generate header
       line_to_write="ENSEMBL_ID"
       sample_names.sort()
       for sample_name in sample_names:
           line_to_write+=('\t'+ sample_name)
       line_to_write += "\n"
       f.write(line_to_write)    
       #Write the rest of the lines
       gene_list.sort()
       for gene in gene_list:
           line_to_write=gene
           for sample_name in sample_names:
               line_to_write += ('\t'+table_dict[sample_name][gene])
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
