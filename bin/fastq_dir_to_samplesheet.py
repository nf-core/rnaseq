#!/usr/bin/env python

import os
import sys
import glob
import argparse


def parse_args(args=None):
    Description = "Generate nf-core/rnaseq samplesheet from a directory of FastQ files."
    Epilog = "Example usage: python fastq_dir_to_samplesheet.py <FASTQ_DIR> <SAMPLESHEET_FILE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FASTQ_DIR", help="Folder containing raw FastQ files.")
    parser.add_argument("SAMPLESHEET_FILE", help="Output samplesheet file.")
    parser.add_argument(
        "-st",
        "--strandedness",
        type=str,
        dest="STRANDEDNESS",
        default="unstranded",
        help="Value for 'strandedness' in samplesheet. Must be one of 'unstranded', 'forward', 'reverse'.",
    )
    parser.add_argument(
        "-r1",
        "--read1_pattern",
        type=str,
        dest="READ1_PATTERN",
        default="*_R1_001.fastq.gz",
        help="File pattern for read 1.",
    )
    parser.add_argument(
        "-r2",
        "--read2_pattern",
        type=str,
        dest="READ2_PATTERN",
        default="*_R2_001.fastq.gz",
        help="File pattern for read 2.",
    )
    parser.add_argument(
        "-se",
        "--single_end",
        dest="SINGLE_END",
        action="store_true",
        help="Single-end information will be auto-detected but this option forces paired-end FastQ files to be treated as single-end so only read 1 information is included in the samplesheet.",
    )
    parser.add_argument(
        "-sn",
        "--sanitise_name",
        dest="SANITISE_NAME",
        action="store_true",
        help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with --sanitise_name_delimiter and --sanitise_name_index.",
    )
    parser.add_argument(
        "-sd",
        "--sanitise_name_delimiter",
        type=str,
        dest="SANITISE_NAME_DELIMITER",
        default="_",
        help="Delimiter to use to sanitise sample name.",
    )
    parser.add_argument(
        "-si",
        "--sanitise_name_index",
        type=int,
        dest="SANITISE_NAME_INDEX",
        default=1,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    return parser.parse_args(args)


def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    strandedness="unstranded",
    read1_pattern="*_R1_001.fastq.gz",
    read2_pattern="*_R2_001.fastq.gz",
    single_end=False,
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
):
    def sanitize_sample(path):
        sample = os.path.splitext(os.path.basename(path))[0]
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                    :sanitise_name_index
                ]
            )
        return sample

    ## Get read 1 files
    read_dict = {}
    read1_files = glob.glob(os.path.join(fastq_dir, read1_pattern), recursive=False)
    for read1_file in read1_files:
        sample = sanitize_sample(read1_file)
        read_dict[sample] = [read1_file]

    ## Get read 2 files
    if not single_end:
        read2_files = glob.glob(os.path.join(fastq_dir, read2_pattern), recursive=False)
        for read2_file in read2_files:
            sample = sanitize_sample(read2_file)
            read_dict[sample] += [read2_file]

    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fastq_1", "fastq_2", "strandedness"]
            fout.write(",".join(header) + "\n")
            for sample, reads in read_dict.items():
                sample_info = ",".join([sample] + reads + [strandedness])
                if len(reads) == 1:
                    sample_info += ","
                fout.write(f"{sample_info}\n")
    else:
        error_str = (
            "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        )
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_pattern' parameter\n"
        error_str += "  - '--read2_pattern' parameter\n"
        print(error_str)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)

    strandedness = "unstranded"
    if args.STRANDEDNESS in ["unstranded", "forward", "reverse"]:
        strandedness = args.STRANDEDNESS

    fastq_dir_to_samplesheet(
        fastq_dir=args.FASTQ_DIR,
        samplesheet_file=args.SAMPLESHEET_FILE,
        strandedness=strandedness,
        read1_pattern=args.READ1_PATTERN,
        read2_pattern=args.READ2_PATTERN,
        single_end=args.SINGLE_END,
        sanitise_name=args.SANITISE_NAME,
        sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
        sanitise_name_index=args.SANITISE_NAME_INDEX,
    )


if __name__ == "__main__":
    sys.exit(main())
