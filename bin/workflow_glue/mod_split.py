#!/usr/bin/env python
"""Split bedMethyl into subfiles."""

import gzip
import os

from ezcharts.components.common import MOD_CONVERT

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    # Check that the input files exist
    if not os.path.exists(args.bedmethyl):
        raise FileNotFoundError(f"File {args.bedmethyl} not found.")

    # Load other input files
    logger = get_named_logger("mod_split")
    logger.info(f'Splitting: {args.bedmethyl}')

    # Define file root name and file reader
    file_read = open
    f_root_name = args.bedmethyl
    if args.bedmethyl.endswith('.gz'):
        file_read = gzip.open
        f_root_name = args.bedmethyl.strip('.gz')

    # Dict of output files
    output_files = {}
    # Start processing the files
    for line in file_read(args.bedmethyl):
        # Get change type
        line = line.decode()
        change = MOD_CONVERT.get(line.split()[3], line.split()[3])
        # Get strandedness
        strand = line.split()[5]
        # Define a prefix (which is also the key)
        if strand == '.':
            prefix = change
        else:
            prefix = f"{change}_{strand}"
        # Create the file, if it doesn't exits
        if prefix not in output_files:
            output_files[prefix] = open(f"{prefix}.{f_root_name}", 'w')
        # Write the line
        output_files[prefix].write(line)

    # Close everything
    for change, filename in output_files.items():
        filename.close()

    # At-a-glance report
    logger.info("All done.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("mod_split")
    parser.add_argument("bedmethyl", help="Input bedMethyl file")
    return parser
