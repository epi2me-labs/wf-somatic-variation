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
        '--min_ref_support', required=False, default=3, type=int,
        help="Number of REF reads required to call a heterozygote site")
    parser.add_argument(
        '--genotype', action="store_true",
        help="Add a genotype field to the output record")
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
            normal_idx = header.index(args.normal_id)
            # Replace TUMOR/CONTROL with provided sample ID
            header = header[0:9] + [args.sample_id]
            # Write header line
            output_file.write('\t'.join(header) + '\n')
            continue
        # Process the records
        record = line.strip().split()
        # Replace the format field
        if args.genotype:
            format_string = 'GT:DP:AF:AD:NDP:NAF:NAD'
        else:
            format_string = 'DP:AF:AD:NDP:NAF:NAD'
        record[header.index('FORMAT')] = format_string
        # Prepare the allele frequencies, alleles depths and total depths
        # for normal and tumor samples.
        tumor_dp, tumor_ad = record[tumor_idx].split(':')
        tumor_af = round(float(tumor_ad)/float(tumor_dp), PRECISION)
        # Define the genotype based on the number of reads supporting the
        # reference allele. If there are < args.min_ref_support, then call
        # the site as homozygote for the alternative allele (1/1). Otherwise,
        # call it as heterozygote (0/1).
        genotype = (
            "0/1" if (int(tumor_dp) - int(tumor_ad)) >= args.min_ref_support else '1/1')
        # Extract values
        normal_dp, normal_ad = record[normal_idx].split(':')
        normal_af = round(float(normal_ad)/float(normal_dp), PRECISION)
        # Remove old sample fields
        record = record[0:9]
        # Create the GT field for the VCF, if requested.
        if args.genotype:
            gt_field = f"{genotype}:"
        else:
            gt_field = ""
        # Save record.
        record.append(
            f'{gt_field}{tumor_dp}:{tumor_af}:{tumor_ad}' +
            f':{normal_dp}:{normal_af}:{normal_ad}'
        )

        # Save the record
        record = "\t".join(record)
        output_file.write(f'{record}\n')


if __name__ == '__main__':
    main()
