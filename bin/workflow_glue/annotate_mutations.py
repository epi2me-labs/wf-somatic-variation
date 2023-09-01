"""Annotate mutation type."""

import itertools
import json
import re

import pandas as pd
import pysam

from .util import wf_parser  # noqa: ABS101


def reverse(seq):
    """Reverse-complement sequence."""
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def main(args):
    """Run the entry point."""
    # Define input files and prepare the required datasets
    i_vcf = pysam.VariantFile(args.i_vcf, 'r')
    # Get sample IDs
    if len(i_vcf.header.samples) > 1:
        raise Exception("VCF files with more than one sample are not supported.")
    sample_id = i_vcf.header.samples[0]
    # Add mutation type header
    i_vcf.header.add_meta(
        'INFO',
        items=[
                ('ID', 'mutation_type'),
                ('Number', 'A'),
                ('Type', 'Character'),
                ('Description', f"{args.k}-mer mutation type")
            ]
    )
    fasta = pysam.FastaFile(args.genome)

    # Check given K-mer size is odd to ensure equal sized flanking
    if args.k % 2 != 1:
        raise ValueError('K-mer size has to be an odd number.')
    # Define mutation types
    # The 6 changes have been chosen to match the MutationalPatterns package
    fsize = int((args.k - 1) / 2)
    flanks = list(map(''.join, itertools.product('ATCG', repeat=fsize)))
    muts = ['[C>A]', '[T>A]', '[C>T]', '[T>G]', '[C>G]', '[T>C]']
    mut_count = dict.fromkeys(map(''.join, itertools.product(flanks, muts, flanks)), 0)
    # Define output file and process inputs
    with pysam.VariantFile(args.o_vcf, 'w', header=i_vcf.header) as o_vcf:
        for rec in i_vcf:
            # If it's a no-alt site (i.e. from genotyping), save
            # the record as-is and continue.
            if not rec.alts:
                o_vcf.write(rec)
                continue
            # First, we get the upper case reference K-mer
            # It now accounts for the size of the flanks when selecting
            # the region (pos is 1-based)
            ref_kmer = fasta.fetch(
                reference=rec.chrom,
                start=rec.pos - 1 - fsize,
                end=rec.pos + fsize
                ).upper()
            # Get filters
            filters = [f for f in rec.filter.keys() if f not in ['PASS', '.']]
            # Get SNVs
            is_snv = len(rec.ref) == len(rec.alts[0]) == 1
            # Check that there are valid nucleotides in the K-mers, that it is a SNP
            # and that it is sorrounded by actual flanks
            if re.match(r'[ACTGactg]+$', ref_kmer) and is_snv and rec.pos > fsize:
                # If it is an snv and all Nts are ACTG, check if the central base == ref
                if ref_kmer[fsize] != rec.ref:
                    raise ValueError('Reference allele does not match fasta sequence.')
                # If so, create the non-ref state using the same flanking
                alt_kmer = ref_kmer[0:fsize] + rec.alts[0] + ref_kmer[fsize+1:]
                # And then prepare the mutation.
                mut = f"{ref_kmer[0:fsize]}[{ref_kmer[fsize]}>" +\
                    f"{alt_kmer[fsize]}]{ref_kmer[fsize+1:]}"
                # Check if mutation is among the mutation types
                if mut not in mut_count.keys():
                    ref_kmer = reverse(ref_kmer)
                    alt_kmer = reverse(alt_kmer)
                    mut = f"{ref_kmer[0:fsize]}[{ref_kmer[fsize]}>" +\
                        f"{alt_kmer[fsize]}]{ref_kmer[fsize+1:]}"
                    # Check if the rev-comp is among mutation types
                    if mut not in mut_count.keys():
                        raise ValueError(f"Change {mut} is not valid.")
                # If so, set mutation type
                rec.info['mutation_type'] = mut
                # Count only if PASS
                if len(filters) == 0:
                    mut_count[mut] += 1
            # If not, save as missing mutation_type
            else:
                rec.info.__setitem__('mutation_type', '.')
            o_vcf.write(rec)

    # Save output matrix of counts
    with open(f'{sample_id}_changes.csv', 'w') as o_file:
        o_file.write(f'Type,{sample_id}\n')
        for key in sorted(mut_count.keys()):
            o_file.write(f'{key},{mut_count[key]}\n')

    # Save json for sankey
    if args.json:
        outjson = {}
        min_rank = -1 * int(args.k/2)
        regexp = r'(?P<Before>[\s\S]+)\[(?P<Change>[\s\S\>]+)\](?P<After>[\s\S]+)'  # noqa: E501
        df = pd.read_csv(f'{sample_id}_changes.csv')
        df[['Before', 'Change', 'After']] = df['Type'].str.extract(regexp, expand=True)
        df = df[['Before', 'Change', 'After', sample_id]]

        # Add before-change site to the json
        before = df[['Before', sample_id]].groupby('Before').sum().reset_index()
        for _, row in before.iterrows():
            b, t = row[['Before', sample_id]]
            outjson[b] = {
                'rank': str(min_rank),
                'children': {},
                'count': str(t)
                }

        # Then add the changes to the sublevel
        change = df[['Before', 'Change', sample_id]]\
            .groupby(['Before', 'Change'])\
            .sum()\
            .reset_index()
        for _, row in change.iterrows():
            b, c, t = row[['Before', 'Change', sample_id]]
            outjson[b]['children'][c] = {
                'rank': str(min_rank + 1),
                'children': {},
                'count': str(t)
                }

        # Finally, the individual level change
        for _, (b, c, a, t) in df.iterrows():
            outjson[b]['children'][c]['children'][a] = {
                'rank': str(min_rank + 2),
                'children': {},
                'count': str(t)
                }

        # Finally save the json
        with open(f'{sample_id}_changes.json', 'w') as j_file:
            json.dump(outjson, j_file)
    return 0


def argparser():
    """Create argument parser."""
    parser = wf_parser("annotate_mutations")

    parser.add_argument("i_vcf", help="Input vcf file")
    parser.add_argument("o_vcf", help="Input vcf file")
    parser.add_argument(
        "-k", default=3, type=int,
        help="K-mer size (has to be an odd number)")
    parser.add_argument(
        "--genome", required=True,
        help="Input fasta file")
    parser.add_argument(
        "--json", action="store_true",
        help="Save the spectrum as json for the sankey plot")

    return parser
