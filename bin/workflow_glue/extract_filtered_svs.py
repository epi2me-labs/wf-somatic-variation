#!/usr/bin/env python
"""check_seq_ref.

Compare a FASTA reference to BAM/CRAM SQ lines and determine
if realignment is required.
"""

import csv
import os
import sys

import pysam
from .util import wf_parser  # noqa: ABS101

# NOTE Both OK and DATAERR are permissible exits from this script so
#      if you want to raise an error to stop Nextflow, it better be
#      something else!


def argparser():
    """Create argument parser."""
    parser = wf_parser("extract_filtered_svs")

    parser.add_argument(
        "--in_vcf",
        required=True,
    )
    parser.add_argument(
        "--filtered",
        required=True,
    )
    parser.add_argument(
        "--out_vcf",
        required=True,
    )
    return parser


# Fast dictionary fetcher
def fetch_value(mydict, mykey, default=None):
    """Quick value fetching from a dictionary."""
    try:
        return mydict[mykey]
    except KeyError:
        return default


def main(args):
    """Run entry point."""
    desc = 'Nanomonsv repeat filtering'
    try:
        i_vcf = pysam.VariantFile(args.in_vcf)
        insites = csv.DictReader(open(args.filtered, 'r'), delimiter='\t')
        sites_id = {}
        for site in insites:
            if site["Chr_1"] != site["Chr_2"]:
                sites_id[site["SV_ID"] + "_0"] = site["Is_Filter"]
            else:
                sites_id[site["SV_ID"]] = site["Is_Filter"]
            if site["Is_Filter"] not in list(i_vcf.header.filters):
                i_vcf.header.filters.add(site["Is_Filter"], None, None, desc)
    except ValueError:
        sys.stderr.write(
            "[FAIL] One (or both) of the input files could not be"
            " prepared. Are they the right format?"
        )
        sys.exit(os.EX_NOINPUT)

    try:
        o_vcf = pysam.VariantFile(args.out_vcf, 'w', header=i_vcf.header)
    except ValueError:
        sys.stderr.write(
            "[FAIL] Impossible to create the output files."
        )
        sys.exit()
    for rec in i_vcf:
        filt = fetch_value(sites_id, rec.id, default='PASS')
        if filt not in rec.filter.keys() or filt == 'PASS':
            rec.filter.add(fetch_value(sites_id, rec.id, default='PASS'))
        o_vcf.write(rec)
    o_vcf.close()
