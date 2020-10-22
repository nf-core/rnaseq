#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Create samplesheet with FTP download links and md5ums from sample information obtained via 'sra_ids_to_runinfo.py' script."
    Epilog = 'Example usage: python sra_runinfo_to_ftp.py <FILES_IN> <FILE_OUT>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILES_IN', help="Comma-separated list of metadata file created from 'sra_ids_to_runinfo.py' script.")
    parser.add_argument('FILE_OUT', help="Output file containing paths to download FastQ files along with their associated md5sums.")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def parse_sra_runinfo(file_in):
    runinfo_dict = {}
    with open(file_in, "r") as fin:
        header = fin.readline().strip().split('\t')
        for line in fin:
            line_dict   = dict(zip(header,line.strip().split('\t')))
            run_id      = line_dict['run_accession']
            exp_id      = line_dict['experiment_accession']
            library     = line_dict['library_layout']
            fastq_files = line_dict['fastq_ftp']
            fastq_md5   = line_dict['fastq_md5']

            db_id = exp_id
            sample_dict = {}
            if library == 'SINGLE':
                sample_dict = {'single_end':True, 'is_ftp':False, 'fastq':[], 'md5':[]}
                if fastq_files:
                    sample_dict['is_ftp'] = True
                    sample_dict['fastq']  = [fastq_files, '']
                    sample_dict['md5']    = [fastq_md5, '']
                else:
                    ## In some instances fastq files don't exist for an entry and have to be downloaded via parallel-fastq-dump
                    db_id = run_id
            
            elif library == 'PAIRED':
                sample_dict = {'single_end':False, 'is_ftp':False, 'fastq':[], 'md5':[]}
                if fastq_files:
                    fq_files = fastq_files.split(';')[-2:]
                    if fq_files[0].find('_1.fastq.gz') != -1 and fq_files[1].find('_2.fastq.gz') != -1:
                        sample_dict['is_ftp'] = True
                        sample_dict['fastq']  = fq_files
                        sample_dict['md5']    = fastq_md5.split(';')[-2:]
                    else:
                        print("Invalid FastQ files found for database id:'{}'!.".format(run_id))
                else:
                    db_id = run_id

            if sample_dict:
                sample_info  = [str(int(sample_dict['single_end'])), str(int(sample_dict['is_ftp']))]
                sample_info += sample_dict['fastq']
                sample_info += sample_dict['md5']
                if db_id not in runinfo_dict:
                    runinfo_dict[db_id] = [sample_info]
                else:
                    if sample_info in runinfo_dict[db_id]:
                        print("Input run info file contains duplicate rows!\nLine: '{}'".format(line))
                    else:
                        runinfo_dict[db_id].append(sample_info)
    return runinfo_dict


def sra_runinfo_to_ftp(files_in,file_out):
    samplesheet_dict = {}
    for file_in in files_in:
        runinfo_dict = parse_sra_runinfo(file_in)        
        for db_id in runinfo_dict.keys():
            if db_id not in samplesheet_dict:
                samplesheet_dict[db_id] = runinfo_dict[db_id]
            else:
                print("Duplicate sample identifier found!\nID: '{}'".format(db_id))
    
    ## Write samplesheet with paths to FastQ files and md5 sums
    if len(samplesheet_dict) != 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(['id', 'single_end', 'is_ftp', 'fastq_1', 'fastq_2', 'md5_1', 'md5_2']) + "\n")
            for db_id in sorted(samplesheet_dict.keys()):
                for idx,val in enumerate(samplesheet_dict[db_id]):
                    fout.write(','.join(["{}_T{}".format(db_id,idx+1)] + val) + '\n')


def main(args=None):
    args = parse_args(args)
    sra_runinfo_to_ftp([x.strip() for x in args.FILES_IN.split(',')], args.FILE_OUT)
    

if __name__ == '__main__':
    sys.exit(main())
