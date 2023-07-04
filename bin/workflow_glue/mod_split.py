#!/usr/bin/env python
"""Split bedMethyl into subfiles."""

import os

from ezcharts.components.common import CATEGORICAL, MOD_CONVERT
from ezcharts.components.modkit import load_bedmethyl
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


relevant_stats_cols_dtypes = {
    "chrom": CATEGORICAL,
    "start": int,
    "end": int,
    "mod": CATEGORICAL,
    "score": int,
    "strand": CATEGORICAL,
    "startp": int,
    "endp": int,
    "colour": CATEGORICAL,
    "Nvalid": int,
    "fraction": float,
    "Nmod": int,
    "Ncanonical": int,
    "Nother_mod": int,
    "Ndelete": int,
    "Nfail": int,
    "Ndiff": int,
    "Nnocall": int,
}


def main(args):
    """Run the entry point."""
    # Check that the input files exist
    if not os.path.exists(args.bedmethyl):
        raise FileNotFoundError(f"File {args.bedmethyl} not found.")

    # Load other input files
    logger = get_named_logger("mod_split")
    logger.info(f'Load: {args.bedmethyl}')
    bedmethyl = load_bedmethyl(args.bedmethyl)

    # If bedmethyl is empty, then exit with error:
    if bedmethyl.empty:
        raise pd.errors.EmptyDataError(
            f'Input {args.bedmethyl} is empty. Check that the BAM file ' +
            'has the MM/ML flags required by modkit.')

    # Add missing columns
    bedmethyl.insert(
        loc=6, column='endp', value=bedmethyl['end'].values)
    bedmethyl.insert(
        loc=6, column='startp', value=bedmethyl['start'].values)
    # ... and drop the unnecessary ones
    bedmethyl = bedmethyl.drop(columns=['total_mean_pos', 'filename'])

    # Split bedmethyl based on change type
    for change in bedmethyl['mod'].unique():
        for strand in bedmethyl['strand'].unique():
            logger.info(f'Processing: {change}; Strand: {strand}')
            subdf = bedmethyl[
                (bedmethyl['mod'] == change) &
                (bedmethyl['strand'] == strand)
                ]
            # Get modification code
            mod = MOD_CONVERT.get(change, change)
            # If strand-aware calling, then make the prefix contain the info
            # to keep the data separate.
            prefix = f'{mod}_{strand}' if strand != '.' else mod
            subdf.to_csv(
                f'{prefix}.{args.bedmethyl}',
                sep='\t',
                header=False,
                index=False)

    # At-a-glance report
    logger.info("All done.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("mod_split")
    parser.add_argument("bedmethyl", help="Input bedMethyl file")
    return parser
