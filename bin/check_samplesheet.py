#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context='Line', context_str=''):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != '' and context_str != '':
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(error, context.strip(), context_str.strip())
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    group,replicate,fastq_1,fastq_2
    WT,1,WT_LIB1_REP1_1.fastq.gz,WT_LIB1_REP1_2.fastq.gz
    WT,1,WT_LIB2_REP1_1.fastq.gz,WT_LIB2_REP1_2.fastq.gz
    WT,2,WT_LIB1_REP2_1.fastq.gz,WT_LIB1_REP2_2.fastq.gz
    KO,1,KO_LIB1_REP1_1.fastq.gz,KO_LIB1_REP1_2.fastq.gz
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2']
        header = fin.readline().strip().split(",")
        if header[:len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip() for x in line.strip().split(",")][:len(HEADER)]
            sample, replicate, fastq_1, fastq_2 = lspl

            ## Check valid number of columns per row
            if len(lspl) != len(HEADER):
                print_error("Invalid number of columns (minimum = {})!".format(len(HEADER)), 'Line', line)

            num_cols = len([x for x in lspl if x])
            if num_cols < 3:
                print_error("Invalid number of populated columns (minimum = 3)!", 'Line', line)

            ## Check sample name entries
            if sample:
                if sample.find(" ") != -1:
                    print_error("Group entry contains spaces!", 'Line', line)
            else:
                print_error("Group entry has not been specified!", 'Line', line)

            ## Check replicate entry is integer
            if not replicate.isdigit():
                print_error("Replicate id not an integer!", 'Line', line)
            replicate = int(replicate)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", 'Line', line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error("FastQ file does not have extension '.fastq.gz' or '.fq.gz'!", 'Line', line)

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", 'Line', line)

            ## Create sample mapping dictionary = {sample: {replicate : [ single_end, fastq_1, fastq_2 ]}}
            if sample not in sample_run_dict:
                sample_run_dict[sample] = {}
            if replicate not in sample_run_dict[sample]:
                sample_run_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_run_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", 'Line', line)
                else:
                    sample_run_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:

            fout.write(",".join(['sample', 'single_end', 'fastq_1', 'fastq_2']) + "\n")
            for sample in sorted(sample_run_dict.keys()):

                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_run_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error("Replicate IDs must start with 1..<num_replicates>!", 'Group', sample)

                for replicate in sorted(sample_run_dict[sample].keys()):

                    ## Check that multiple runs of the same sample are of the same datatype
                    if not all(x[0] == sample_run_dict[sample][replicate][0][0] for x in sample_run_dict[sample][replicate]):
                        print_error("Multiple runs of a sample must be of the same datatype!", 'Group', sample)

                    ## Write to file
                    for idx, sample_info in enumerate(sample_run_dict[sample][replicate]):
                        sample_id = "{}_R{}_T{}".format(sample,replicate,idx+1)
                        fout.write(','.join([sample_id] + sample_info) + '\n')


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
