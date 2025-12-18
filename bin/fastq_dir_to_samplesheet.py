#!/usr/bin/env python3

# Written by Harshil Patel and released under the MIT license.

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
        default="auto",
        help="Value for 'strandedness' in samplesheet. Must be one of 'unstranded', 'forward', 'reverse', 'auto'.",
    )
    parser.add_argument(
        "-r1",
        "--read1_extension",
        type=str,
        dest="READ1_EXTENSION",
        default="_R1_001.fastq.gz",
        help="File extension for read 1.",
    )
    parser.add_argument(
        "-r2",
        "--read2_extension",
        type=str,
        dest="READ2_EXTENSION",
        default="_R2_001.fastq.gz",
        help="File extension for read 2.",
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
    parser.add_argument(
        "-re",
        "--recursive",
        dest="RECURSIVE",
        action="store_true",
        help="Whether or not to search for FastQ files recursively in <FASTQ_DIR>.",
    )
    return parser.parse_args(args)


def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    strandedness="auto",
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    single_end=False,
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
    recursive=False,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[:sanitise_name_index]
            )
        return sample

    def get_fastqs(extension, recursive=False):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        search_path = f"*{extension}"
        if recursive:
            search_path = f"**/*{extension}"
        return sorted(glob.glob(os.path.join(fastq_dir, search_path), recursive=recursive))

    read_dict = {}

    ## Get read 1 files
    for read1_file in get_fastqs(read1_extension, recursive):
        sample = sanitize_sample(read1_file, read1_extension)
        if sample not in read_dict:
            read_dict[sample] = {"R1": [], "R2": []}
        read_dict[sample]["R1"].append(read1_file)

    ## Get read 2 files
    if not single_end:
        for read2_file in get_fastqs(read2_extension, recursive):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fastq_1", "fastq_2", "strandedness"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, read_1 in enumerate(reads["R1"]):
                    read_2 = ""
                    if idx < len(reads["R2"]):
                        read_2 = reads["R2"][idx]
                    sample_info = ",".join([sample, read_1, read_2, strandedness])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)

    strandedness = "auto"
    if args.STRANDEDNESS in ["unstranded", "forward", "reverse", "auto"]:
        strandedness = args.STRANDEDNESS

    fastq_dir_to_samplesheet(
        fastq_dir=args.FASTQ_DIR,
        samplesheet_file=args.SAMPLESHEET_FILE,
        strandedness=strandedness,
        read1_extension=args.READ1_EXTENSION,
        read2_extension=args.READ2_EXTENSION,
        single_end=args.SINGLE_END,
        sanitise_name=args.SANITISE_NAME,
        sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
        sanitise_name_index=args.SANITISE_NAME_INDEX,
        recursive=args.RECURSIVE,
    )


if __name__ == "__main__":
    sys.exit(main())
