#!/usr/bin/env python

import os
import re
import sys
import csv
import errno
import requests
import argparse


## Example ids supported by this script
SRA_IDS = ['PRJNA63463', 'SAMN00765663', 'SRA023522', 'SRP003255', 'SRR390278', 'SRS282569', 'SRX111814']
ENA_IDS = ['ERA2421642', 'ERP120836', 'ERR674736', 'ERS4399631', 'ERX629702', 'PRJEB7743', 'SAMEA3121481']
GEO_IDS = ['GSE18729', 'GSM465244']
ID_REGEX = r'^[A-Z]+'
PREFIX_LIST = sorted(list(set([re.search(ID_REGEX,x).group() for x in SRA_IDS + ENA_IDS + GEO_IDS])))


def parse_args(args=None):
    Description = 'Download and create a run information metadata file from SRA/ENA/GEO identifiers.'
    Epilog = 'Example usage: python fetch_sra_runinfo.py <FILE_IN> <FILE_OUT>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="File containing database identifiers, one per line.")
    parser.add_argument('FILE_OUT', help="Output file in tab-delimited format.")
    parser.add_argument('-pl', '--platform', type=str, dest="PLATFORM", default='', help="Comma-separated list of platforms to use for filtering. Accepted values = 'ILLUMINA', 'OXFORD_NANOPORE' (default: '').")
    parser.add_argument('-ll', '--library_layout', type=str, dest="LIBRARY_LAYOUT", default='', help="Comma-separated list of library layouts to use for filtering. Accepted values = 'SINGLE', 'PAIRED' (default: '').")
    return parser.parse_args(args)


def validate_csv_param(param,validVals,param_desc):
    validList = []
    if param:
        userVals = param.split(',')
        intersect = list(set(userVals) & set(validVals))
        if len(intersect) == len(userVals):
            validList = intersect
        else:
            print("ERROR: Please provide a valid {} parameter!\nProvided values = {}\nAccepted values = {}".format(param_desc,param,','.join(validVals)))
            sys.exit(1)
    return validList


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def fetch_url(url,encoding='utf-8'):
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    if r.status_code != 200:
        print("ERROR: Connection failed\nError code '{}'".format(r.status_code))
        sys.exit(1)
    return r.content.decode(encoding).splitlines()


def id_to_srx(db_id):
    ids = []
    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={}'.format(db_id)
    for row in csv.DictReader(fetch_url(url), delimiter=','):
        ids.append(row['Experiment'])
    return ids


def id_to_erx(db_id):
    ids = []
    fields = ['run_accession', 'experiment_accession']
    url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run&fields={}'.format(db_id,','.join(fields))
    for row in csv.DictReader(fetch_url(url), delimiter='\t'):
        ids.append(row['experiment_accession'])
    return ids


def gse_to_srx(db_id):
    ids = []
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=gsm&view=data&form=text'.format(db_id)
    gsm_ids = [x.split('=')[1].strip() for x in fetch_url(url) if x.find('GSM') != -1]
    for gsm_id in gsm_ids:
        ids += id_to_srx(gsm_id)
    return ids


def get_ena_fields():
    fields = []
    url = 'https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&format=tsv&result=read_run'
    for row in csv.DictReader(fetch_url(url), delimiter='\t'):
        fields.append(row['columnId'])
    return fields


def fetch_sra_runinfo(FileIn,FileOut,platformList=[],libraryLayoutList=[]):
    total_out = 0
    seen_ids = []; run_ids = []
    header = []
    make_dir(os.path.dirname(FileOut))
    ena_fields = get_ena_fields()
    fin = open(FileIn,'r')
    fout = open(FileOut,'w')
    while True:
        line = fin.readline()
        if line:
            db_id = line.strip()
            match = re.search(ID_REGEX, db_id)
            if match:
                prefix = match.group()
                if prefix in PREFIX_LIST:
                    if not db_id in seen_ids:

                        ids = [db_id]
                        ## Resolve/expand these ids against GEO URL
                        if prefix in ['GSE']:
                            ids = gse_to_srx(db_id)

                        ## Resolve/expand these ids against SRA URL
                        elif prefix in ['GSM', 'PRJNA', 'SAMN', 'SRR']:
                            ids = id_to_srx(db_id)

                        ## Resolve/expand these ids against ENA URL
                        elif prefix in ['ERR']:
                            ids = id_to_erx(db_id)

                        ## Resolve/expand to get run identifier from ENA and write to file
                        for id in ids:
                            url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run&fields={}'.format(id,','.join(ena_fields))
                            csv_dict = csv.DictReader(fetch_url(url), delimiter='\t')
                            for row in csv_dict:
                                run_id = row['run_accession']
                                if not run_id in run_ids:

                                    writeID = True
                                    if platformList:
                                        if row['instrument_platform'] not in platformList:
                                            writeID = False
                                    if libraryLayoutList:
                                        if row['library_layout'] not in libraryLayoutList:
                                            writeID = False

                                    if writeID:
                                        if total_out == 0:
                                            header = sorted(row.keys())
                                            fout.write('{}\n'.format('\t'.join(sorted(header))))
                                        else:
                                            if header != sorted(row.keys()):
                                                print("ERROR: Metadata columns do not match for id {}!\nLine: '{}'".format(run_id,line.strip()))
                                                sys.exit(1)
                                        fout.write('{}\n'.format('\t'.join([row[x] for x in header])))
                                        total_out += 1
                                    run_ids.append(run_id)
                        seen_ids.append(db_id)
                else:
                    id_str = ', '.join([x + "*" for x in PREFIX_LIST])
                    print("ERROR: Please provide a valid database id starting with {}!\nLine: '{}'".format(id_str,line.strip()))
                    sys.exit(1)
            else:
                id_str = ', '.join([x + "*" for x in PREFIX_LIST])
                print("ERROR: Please provide a valid database id starting with {}!\nLine: '{}'".format(id_str,line.strip()))
                sys.exit(1)
        else:
            break
            fin.close()
            fout.close()


def main(args=None):
    args = parse_args(args)
    platformList = validate_csv_param(args.PLATFORM,validVals=['ILLUMINA'],param_desc='--platform')
    libraryLayoutList = validate_csv_param(args.LIBRARY_LAYOUT,validVals=['SINGLE', 'PAIRED'],param_desc='--library_layout')
    fetch_sra_runinfo(args.FILE_IN,args.FILE_OUT,platformList,libraryLayoutList)


if __name__ == '__main__':
    sys.exit(main())
