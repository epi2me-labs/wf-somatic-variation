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
    parser = wf_parser("classify_vcf_svs")

    parser.add_argument(
        "--in_vcf",
        required=True,
    )
    parser.add_argument(
        "--annotated",
        required=True,
    )
    parser.add_argument(
        "--original",
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


# Get annotation columns
def fetch_annotation(annotated, header):
    """Extract annotation from a second file header."""
    for i in annotated.readline().strip().split():
        if i not in header:
            yield i


def main(args):
    """Run entry point."""
    try:
        filter_description = 'Transposable elements inferred by'
        filter_description += ' nanomonsv insert_classify with structure'
        i_vcf = pysam.VariantFile(args.in_vcf)
        with open(args.original) as preclass:
            header = preclass.readline().strip().split()
        with open(args.annotated) as annotated:
            annotations = [i for i in fetch_annotation(annotated, header)]
            filter_description += ' ' + '|'.join(annotations)
        i_vcf.header.add_meta('INFO', items=[('ID', "REPCLASS"),
                                             ('Number', 1),
                                             ('Type', 'String'),
                                             ('Description', filter_description)])
        insites = csv.DictReader(open(args.annotated, 'r'), delimiter='\t')
        site_annots = {}
        for site in insites:
            for k in annotations:
                for orig, new in zip([',', ';', '(', ')'], ['-', '-', '', '']):
                    site[k] = site[k].replace(orig, new)
            if site["Chr_1"] != site["Chr_2"]:
                ann = '|'.join([site[k] for k in annotations])
                site_annots[site["SV_ID"] + "_0"] = ann
                site_annots[site["SV_ID"] + "_1"] = ann
            else:
                site_annots[site["SV_ID"]] = '|'.join([site[k] for k in annotations])
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
        rec.info.__setitem__('REPCLASS', fetch_value(site_annots, rec.id, default='.'))
        o_vcf.write(rec)
    o_vcf.close()
