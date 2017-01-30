#!/usr/bin/env/python
#Merges all fastq_files from all samples into one file. 

import re
import os
import glob
import shutil
def merge_files(dest_dir, fastq_files):
    #sample_pattern=re.compile("^(.+)_S[0-9]+_.+_R([1-2])_")
    #Below is used for BCBP above for Irma style projects.
    sample_pattern=re.compile("^[0-9]_[0-9]+_[A-Z0-9]+_(P[0-9]+_[0-9]+)_([1-2])")
    while fastq_files:
        tomerge=[]
        first=fastq_files[0]
        fq_bn=os.path.basename(first)
        sample_name=sample_pattern.match(fq_bn).group(1)
        read_nb=sample_pattern.match(fq_bn).group(2)
        for fq in fastq_files:
            if sample_name in os.path.basename(fq) and "_{}.".format(read_nb) in os.path.basename(fq):
                tomerge.append(fq)
        for fq in tomerge:
            fastq_files.remove(fq)

        outfile=os.path.join(dest_dir, "{}_R{}.fastq.gz".format(sample_name, read_nb))
        with open(outfile, 'wb') as wfp:
            for fn in tomerge:
                with open(fn, 'rb') as rfp:
                    shutil.copyfileobj(rfp, wfp)

def main():
    destination_dir = str(raw_input("Destiantion dir:"))
    #all fastq.gz files in CWD 
    fastq_files=glob.glob('*.fastq.gz')
    merge_files(destination_dir,fastq_files)

if __name__ == "__main__":
    main()
