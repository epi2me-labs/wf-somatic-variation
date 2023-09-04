#!/usr/bin/env python
"""Convert two-sample VCF to single-sample VCF."""

import argparse
import gzip


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--vcf', required=True,
        help="Output from samtools faidx")
    parser.add_argument(
        '--sample_id', required=False, default='SAMPLE',
        help="Sample ID")
    parser.add_argument(
        '--control_id', required=False, default="CONTROL",
        help="Control sample ID")
    parser.add_argument(
        '--tumor_id', required=False, default="TUMOR",
        help="Tumor sample ID")
    parser.add_argument(
        '-o', '--output', required=True,
        help="Output genome")
    args = parser.parse_args()

    # Prepare additional format details.
    # These are consistent with some in ClairS outputs.
    gtfield = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    naffield = '##FORMAT=<ID=NAF,Number=1,Type=Float,Description=' +\
        '"Estimated allele frequency in the normal BAM">\n'
    nadfield = '##FORMAT=<ID=NAD,Number=1,Type=Float,Description=' +\
        '"Alternate allele depth in the normal BAM">\n'
    ndpfield = '##FORMAT=<ID=NDP,Number=1,Type=Float,Description=' +\
        '"Read depth in the tumor BAM">\n'
    vaffield = '##FORMAT=<ID=AF,Number=1,Type=Float,Description=' +\
        '"Estimated allele frequency in the tumor BAM">\n'
    vadfield = '##FORMAT=<ID=AD,Number=1,Type=Float,Description=' +\
        '"Alternate allele depth in the tumor BAM">\n'
    vdpfield = '##FORMAT=<ID=DP,Number=1,Type=Float,Description=' +\
        '"Read depth in the tumor BAM">\n'

    # Prepare output file
    output_file = open(args.output, 'w')

    # Define correct loader
    if args.vcf.endswith('.gz'):
        loader = gzip.open(args.vcf)
    else:
        loader = open(args.vcf)

    # explode on bad genome
    for line in loader:
        try:
            line = line.decode()
        except AttributeError:
            pass

        if '##' in line:
            output_file.write(line)
            continue
        if line.startswith('#CHROM'):
            # Add missing header lines
            output_file.write(gtfield)
            output_file.write(vaffield)
            output_file.write(vadfield)
            output_file.write(vdpfield)
            output_file.write(naffield)
            output_file.write(nadfield)
            output_file.write(ndpfield)
            header = line.strip().split()
            # Get field index in each record
            tumor_idx = header.index(args.tumor_id)
            control_idx = header.index(args.control_id)
            # Replace TUMOR/CONTROL with provided sample ID
            header = [
                f for n, f in enumerate(header) if n not in [tumor_idx, control_idx]] +\
                [args.sample_id]
            # Write header line
            output_file.write('\t'.join(header) + '\n')
            continue
        # Process the records
        record = line.strip().split()
        # Replace the format field
        record[header.index('FORMAT')] = 'GT:DP:AF:AD:NDP:NAF:NAD'
        # Prepare the allele frequencies, alleles depths and total depths
        # for normal and tumor samples.
        ndp, nad = record[control_idx].split(':')
        tdp, tad = record[tumor_idx].split(':')
        naf = round(float(nad)/float(ndp), 4)
        taf = round(float(tad)/float(tdp), 4)
        # Remove old sample fields
        record = [f for n, f in enumerate(record) if n not in [tumor_idx, control_idx]]
        # Add the new sample
        record.append(f'0/1:{tdp}:{taf}:{tad}:{ndp}:{naf}:{nad}')
        # Save the record
        record = "\t".join(record)
        output_file.write(f'{record}\n')


if __name__ == '__main__':
    main()
