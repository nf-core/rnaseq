#!/usr/bin/env python

## This script is used to generate scan *.nf.test files for function/process/workflow name and return as a JSON list
# It is functionally similar to nf-test list but fills a gap until feature https://github.com/askimed/nf-test/issues/196 is added

import argparse
import json
import logging
import re

from itertools import chain
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments and return an ArgumentParser object.

    Returns:
        argparse.ArgumentParser: The ArgumentParser object with the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Scan *.nf.test files for function/process/workflow name and return as a JSON list"
    )
    parser.add_argument(
        "-p",
        "--paths",
        nargs="+",
        default=["."],
        help="List of directories or files to scan",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level",
    )
    parser.add_argument(
        "-t",
        "--types",
        nargs="+",
        choices=["function", "process", "workflow", "pipeline"],
        default=["function", "process", "workflow", "pipeline"],
        help="Types of tests to include.",
    )
    return parser.parse_args()


def find_files(paths: list[str]) -> list[Path]:
    """
    Find all files matching pattern *.nf.test recursively from a list of paths.

    Args:
        paths (list): List of directories or files to scan.

    Returns:
        list: List of files matching the pattern *.nf.test.
    """
    # this is a bit clunky
    result = []
    for path in paths:
        path_obj = Path(path)
        # If Path is the exact nf-test file add to list:
        if path_obj.match("*.nf.test"):
            result.append(path_obj)
        # Else recursively search for nf-test files:
        else:
            for file in path_obj.rglob("*.nf.test"):
                result.append(file)
    return result


def process_files(files: list[Path]) -> list[str]:
    """
    Process the files and return lines that begin with 'workflow', 'process', or 'function' and have a single string afterwards.

    Args:
        files (list): List of files to process.

    Returns:
        list: List of lines that match the criteria.
    """
    result = []
    for file in files:
        with open(file, "r") as f:
            is_pipeline_test = True
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith(("workflow", "process", "function")):
                    words = line.split()
                    if len(words) == 2 and re.match(r'^".*"$', words[1]):
                        result.append(line)
                        is_pipeline_test = False

            # If no results included workflow, process or function
            # Add a dummy result to fill the 'pipeline' category
            if is_pipeline_test:
                result.append("pipeline 'PIPELINE'")

    return result


def generate(
    lines: list[str], types: list[str] = ["function", "process", "workflow", "pipeline"]
) -> dict[str, list[str]]:
    """
    Generate a dictionary of function, process and workflow lists from the lines.

    Args:
        lines (list): List of lines to process.
        types (list): List of types to include.

    Returns:
        dict: Dictionary with function, process and workflow lists.
    """
    result: dict[str, list[str]] = {
        "function": [],
        "process": [],
        "workflow": [],
        "pipeline": [],
    }
    for line in lines:
        words = line.split()
        if len(words) == 2:
            keyword = words[0]
            name = words[1].strip("'\"")  # Strip both single and double quotes
            if keyword in types:
                result[keyword].append(name)
    return result


if __name__ == "__main__":

    # Utility stuff
    args = parse_args()
    logging.basicConfig(level=args.log_level)

    # Parse nf-test files for targets of tests
    files = find_files(args.paths)
    lines = process_files(files)
    result = generate(lines)

    # Get only relevant results (specified by -t)
    # Unique using a set
    target_results = list(
        {item for sublist in map(result.get, args.types) for item in sublist}
    )

    # Print to stdout
    print(json.dumps(target_results))
