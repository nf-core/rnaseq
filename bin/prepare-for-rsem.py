#!/usr/bin/env python3

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Credits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script is a clone of the "prepare-for-rsem.py" script written by
Ian Sudbury, Tom Smith and other contributors to the UMI-tools package:
https://github.com/CGATOxford/UMI-tools

It has been included here to address problems encountered with
Salmon quant and RSEM as discussed in the issue below:
https://github.com/CGATOxford/UMI-tools/issues/465

When the "umi_tools prepare-for-rsem" command becomes available in an official
UMI-tools release this script will be replaced and deprecated.

Commit:
https://github.com/CGATOxford/UMI-tools/blob/bf8608d6a172c5ca0dcf33c126b4e23429177a72/umi_tools/prepare-for-rsem.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prepare_for_rsem - make the output from dedup or group compatible with RSEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The SAM format specification states that the mnext and mpos fields should point
to the primary alignment of a read's mate. However, not all aligners adhere to
this standard. In addition, the RSEM software requires that the mate of a read1
appears directly after it in its input BAM. This requires that there is exactly
one read1 alignment for every read2 and vice versa.

In general (except in a few edge cases) UMI tools outputs only the read2 to that
corresponds to the read specified in the mnext and mpos positions of a selected
read1, and only outputs this read once, even if multiple read1s point to it.
This makes UMI-tools outputs incompatible with RSEM. This script takes the output
from dedup or groups and ensures that each read1 has exactly one read2 (and vice
versa), that read2 always appears directly after read1,and that pairs point to
each other (note this is technically not valid SAM format). Copy any specified
tags from read1 to read2 if they are present (by default, UG and BX, the unique
group and correct UMI tags added by _group_)

Input must to name sorted.


https://raw.githubusercontent.com/CGATOxford/UMI-tools/master/LICENSE

"""

from umi_tools import Utilities as U
from collections import defaultdict, Counter
import pysam
import sys

usage = """
prepare_for_rsem - make output from dedup or group compatible with RSEM

Usage: umi_tools prepare_for_rsem [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is omited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr """


def chunk_bam(bamfile):
    """Take in a iterator of pysam.AlignmentSegment entries and yield
    lists of reads that all share the same name"""

    last_query_name = None
    output_buffer = list()

    for read in bamfile:
        if last_query_name is not None and last_query_name != read.query_name:
            yield (output_buffer)
            output_buffer = list()

        last_query_name = read.query_name
        output_buffer.append(read)

    yield (output_buffer)


def copy_tags(tags, read1, read2):
    """Given a  list of tags, copies the values of these tags from read1
    to read2, if the tag is set"""

    for tag in tags:
        try:
            read1_tag = read1.get_tag(tag, with_value_type=True)
            read2.set_tag(tag, value=read1_tag[0], value_type=read1_tag[1])
        except KeyError:
            pass

    return read2


def pick_mate(read, template_dict, mate_key):
    """Find the mate of read in the template dict using key. It will retrieve
    all reads at that key, and then scan to pick the one that refers to _read_
    as it's mate. If there is no such read, it picks a first one it comes to"""

    mate = None

    # get a list of secondary reads at the correct alignment position
    potential_mates = template_dict[not read.is_read1][mate_key]

    # search through one at a time to find a read that points to the current read
    # as its mate.
    for candidate_mate in potential_mates:
        if (
            candidate_mate.next_reference_name == read.reference_name
            and candidate_mate.next_reference_start == read.pos
        ):
            mate = candidate_mate

    # if no such read is found, then pick any old secondary alignment at that position
    # note: this happens when UMI-tools outputs the wrong read as something's pair.
    if mate is None and len(potential_mates) > 0:
        mate = potential_mates[0]

    return mate


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$", usage=usage, description=globals()["__doc__"])
    group = U.OptionGroup(parser, "RSEM preparation specific options")

    group.add_option(
        "--tags",
        dest="tags",
        type="string",
        default="UG,BX",
        help="Comma-separated list of tags to transfer from read1 to read2",
    )
    group.add_option(
        "--sam", dest="sam", action="store_true", default=False, help="input and output SAM rather than BAM"
    )

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(
        parser, argv=argv, add_group_dedup_options=False, add_umi_grouping_options=False, add_sam_options=False
    )

    skipped_stats = Counter()

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        in_name = "-"

    if options.sam:
        mode = ""
    else:
        mode = "b"

    inbam = pysam.AlignmentFile(in_name, "r" + mode)

    if options.stdout != sys.stdout:
        out_name = options.stdout.name
        options.stdout.close()
    else:
        out_name = "-"

    outbam = pysam.AlignmentFile(out_name, "w" + mode, template=inbam)

    options.tags = options.tags.split(",")

    for template in chunk_bam(inbam):
        assert len(set(r.query_name for r in template)) == 1
        current_template = {True: defaultdict(list), False: defaultdict(list)}

        for read in template:
            key = (read.reference_name, read.pos, not read.is_secondary)
            current_template[read.is_read1][key].append(read)

        output = set()

        for read in template:
            mate = None

            # if this read is a non_primary alignment, we first want to check if it has a mate
            # with the non-primary alignment flag set.

            mate_key_primary = True
            mate_key_secondary = (read.next_reference_name, read.next_reference_start, False)

            # First look for a read that has the same primary/secondary status
            # as read (i.e. secondary mate for secondary read, and primary mate
            # for primary read)
            mate_key = (read.next_reference_name, read.next_reference_start, read.is_secondary)
            mate = pick_mate(read, current_template, mate_key)

            # If none was found then look for the opposite (primary mate of secondary
            # read or seconadary mate of primary read)
            if mate is None:
                mate_key = (read.next_reference_name, read.next_reference_start, not read.is_secondary)
                mate = pick_mate(read, current_template, mate_key)

            # If we still don't have a mate, then their can't be one?
            if mate is None:
                skipped_stats["no_mate"] += 1
                U.warn(
                    "Alignment {} has no mate -- skipped".format(
                        "\t".join(map(str, [read.query_name, read.flag, read.reference_name, int(read.pos)]))
                    )
                )
                continue

            # because we might want to make changes to the read, but not have those changes reflected
            # if we need the read again,we copy the read. This is only way I can find to do this.
            read = pysam.AlignedSegment().from_dict(read.to_dict(), read.header)
            mate = pysam.AlignedSegment().from_dict(mate.to_dict(), read.header)

            # Make it so that if our read is secondary, the mate is also secondary. We don't make the
            # mate primary if the read is primary because we would otherwise end up with mulitple
            # primary alignments.
            if read.is_secondary:
                mate.is_secondary = True

            # In a situation where there is already one mate for each read, then we will come across
            # each pair twice - once when we scan read1 and once when we scan read2. Thus we need
            # to make sure we don't output something already output.
            if read.is_read1:
                mate = copy_tags(options.tags, read, mate)
                output_key = str(read) + str(mate)

                if output_key not in output:
                    output.add(output_key)
                    outbam.write(read)
                    outbam.write(mate)
                    skipped_stats["pairs_output"] += 1

            elif read.is_read2:
                read = copy_tags(options.tags, mate, read)
                output_key = str(mate) + str(read)

                if output_key not in output:
                    output.add(output_key)
                    outbam.write(mate)
                    outbam.write(read)
                    skipped_stats["pairs_output"] += 1

            else:
                skipped_stats["skipped_not_read12"] += 1
                U.warn(
                    "Alignment {} is neither read1 nor read2 -- skipped".format(
                        "\t".join(map(str, [read.query_name, read.flag, read.reference_name, int(read.pos)]))
                    )
                )
                continue

    if not out_name == "-":
        outbam.close()

    U.info(
        "Total pairs output: {}, Pairs skipped - no mates: {},"
        " Pairs skipped - not read1 or 2: {}".format(
            skipped_stats["pairs_output"], skipped_stats["no_mate"], skipped_stats["skipped_not_read12"]
        )
    )
    U.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
