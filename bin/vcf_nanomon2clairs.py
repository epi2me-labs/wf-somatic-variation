#!/usr/bin/env python
"""Convert two-sample VCF to single-sample VCF."""

import argparse
import gzip

from workflow_glue.report_utils.utils import PRECISION  # noqa: ABS101


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
        '--tumor_only', required=False, action="store_true",
        help="The run was tumor-only")
    parser.add_argument(
        '--normal_id', required=False, default="CONTROL",
        help="Normal sample ID")
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
            if not args.tumor_only:
                normal_idx = header.index(args.normal_id)
            # Replace TUMOR/CONTROL with provided sample ID
            header = header[0:9] + [args.sample_id]
            # Write header line
            output_file.write('\t'.join(header) + '\n')
            continue
        # Process the records
        record = line.strip().split()
        # Replace the format field
        record[header.index('FORMAT')] = 'GT:DP:AF:AD:NDP:NAF:NAD'
        # Prepare the allele frequencies, alleles depths and total depths
        # for normal and tumor samples.
        tumor_dp, tumor_ad = record[tumor_idx].split(':')
        tumor_af = round(float(tumor_ad)/float(tumor_dp), PRECISION)
        # If not tumor-only, extract values; otherwise set to NA
        if args.tumor_only:
            normal_af = normal_dp = normal_ad = 0
        else:
            normal_dp, normal_ad = record[normal_idx].split(':')
            normal_af = round(float(normal_ad)/float(normal_dp), PRECISION)
        # Remove old sample fields
        record = record[0:9]
        # Add the new sample
        record.append(
            f'0/1:{tumor_dp}:{tumor_af}:{tumor_ad}' +
            f':{normal_dp}:{normal_af}:{normal_ad}'
        )
        # Save the record
        record = "\t".join(record)
        output_file.write(f'{record}\n')


if __name__ == '__main__':
    main()
