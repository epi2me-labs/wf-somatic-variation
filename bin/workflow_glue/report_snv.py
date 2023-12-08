#!/usr/bin/env python
"""Create SNV report."""

import os

from dominate.tags import a, h6, p
from ezcharts.components.bcfstats import load_bcfstats
from ezcharts.components.clinvar import load_clinvar_vcf
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
import numpy as np
import pandas as pd
import pysam

from .report_utils.utils import COLORS, display_alert, PRECISION  # noqa: ABS101
from .report_utils.visualizations import plot_profile  # noqa: ABS101
from .report_utils.visualizations import plot_spectra  # noqa: ABS101
from .report_utils.visualizations import scatter_plot  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


# Mutation profile palette to match the COSMIC
# patterns.
cmap = {
    'C>A': COLORS.cerulean,
    'C>G': COLORS.black,
    'C>T': COLORS.cinnabar,
    'T>A': COLORS.grey70,
    'T>C': COLORS.medium_spring_bud,
    'T>G': COLORS.fandango,
}


def vcf_parse(args):
    """Parse the vcf file."""
    vcf_df = pysam.VariantFile(args.vcf, 'r')
    # Sample name
    samplename = vcf_df.header.samples[0]
    # Output lists
    flt = []
    vaf = []
    naf = []
    vtype = []
    for rec in vcf_df:
        if not rec.alts:
            continue
        flt.append(rec.filter.keys()[0])
        vaf.append(rec.samples[samplename]['AF'])
        naf.append(rec.samples[samplename]['NAF'])
        if len(rec.alts[0]) == len(rec.ref) and len(rec.ref) == 1:
            vtype.append('SNV')
        else:
            vtype.append('Indel')
    return samplename, flt, vaf, naf, vtype


def filt_stats(filters, vaf, naf, noalt=0, thresholds=[0.2, 0.1, 0.05]):
    """Plot the filtering stats."""
    df = pd.DataFrame({'Filter': filters, 'Tumor_AF': vaf, 'Normal_AF': naf})
    # Filtering table
    summary_table = {
        'N sites': [df.count()['Filter'],
                    ""],
        'No ALTS': [noalt,
                    ""],
        'PASS': [df[df['Filter'] == 'PASS'].count()['Filter'],
                 ""],
        'not PASS': [df[df['Filter'] != 'PASS'].count()['Filter'],
                     ""],
        'Type': ['Tumor', 'Normal']}
    for threshold in thresholds:
        summary_table.update({'VAF > {:.2f}'.format(threshold): [
            df.query(f'Tumor_AF >= {threshold}').shape[0],
            df.query(f'Normal_AF >= {threshold}').shape[0]]})

    # Plot the histogram of the allele frequencies in the tumor and normal
    return pd.DataFrame(summary_table)


def process_spectra(inspectra):
    """Process the mutation spectra."""
    # Mutation order
    try:
        spectra = pd.read_csv(inspectra)
        # Define middle point
        midpoint = int(np.floor(len(spectra['Type'].iloc[0])/2))
        # Define change type and flankings
        spectra["Change"] = spectra['Type'].str.slice(start=midpoint-1, stop=midpoint+2)
        spectra["Flanks"] = spectra['Type'].str.slice(start=0, stop=midpoint-2) + \
            '-' + spectra['Type'].str.slice(start=midpoint+3)
        # Re-arrange the columns to have the Change first and then the flanks
        change_col = spectra.pop("Change")
        flank_col = spectra.pop("Flanks")
        spectra.insert(0, "Flanks", flank_col)
        spectra.insert(0, "Change", change_col)
        # Sort the data by change x mutation
        spectra = spectra.sort_values(by=['Change', 'Type'])
    except pd.errors.EmptyDataError:
        spectra = pd.DataFrame()
    return spectra


def main(args):
    """Run the entry point."""
    logger = get_named_logger("report_snp")
    # Check that the input files exist
    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"File {args.vcf} not found.")
    if not os.path.exists(args.mut_spectra):
        raise FileNotFoundError(f"File {args.mut_spectra} not found.")
    if not os.path.exists(args.vcf_stats):
        raise FileNotFoundError(f"File {args.vcf_stats} not found.")

    # Import var. alleles thresholds
    vaf_thresholds = [float(i) for i in args.vaf_thresholds.split(',') if float(i) < 1]
    if not vaf_thresholds:
        raise Exception("Invalid range of variant allele frequencies thresholds.")

    # Load the data
    sample_id, filters, var_af, nvar_af, vtype = vcf_parse(args)
    filtstats = filt_stats(filters, var_af, nvar_af, thresholds=vaf_thresholds)
    spectra = process_spectra(args.mut_spectra)
    try:
        bcfstats = load_bcfstats(
            args.vcf_stats,
            sample_names=[sample_id])
    except IndexError:
        bcfstats = {'SN': pd.DataFrame(), 'TSTV': pd.DataFrame()}
    # Instantiate the report
    report = LabsReport(
        f"{sample_id} | Small variation statistics", "wf-somatic-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    # VCF At-a-glance report

    with report.add_section('At a glance', 'Summary'):
        if args.no_germline:
            display_alert(
                "This report was produced in somatic-only mode,"
                " no germline calling was performed."
            )
        if args.normal_vcf:
            display_alert(
                "Germline calls for the normal sample were provided by the user",
                f" in file {args.normal_vcf}"
            )
        if args.hybrid_mode_vcf_fn:
            display_alert(
                "Hybrid calling from VCF file provided by the user",
                f": {args.hybrid_mode_vcf_fn}"
            )
        if args.genotyping_mode_vcf_fn:
            display_alert(
                "Genotyping sites in VCF file provided by the user",
                f": {args.genotyping_mode_vcf_fn}"
            )
        tabs = Tabs()
        with tabs.add_tab(sample_id):
            if bcfstats['SN'].empty:
                p('The bcftools stats file is empty.')
            else:
                bcfstats['SN'].columns = [
                    i.replace('SNP', 'SNV').replace('MNP', 'MNV')
                    for i in bcfstats['SN'].columns
                ]
                # Extract values to avoid displaying as list in Stats
                titv = bcfstats['TSTV']['ts/tv'].values[0]
                nsites = bcfstats['SN']['records'].values[0]
                nsnvs = bcfstats['SN']['SNVs'].values[0]
                nindels = bcfstats['SN']['indels'].values[0]
                nvars = int(bcfstats['SN']['records'].values[0]) -\
                    int(bcfstats['SN']['no-ALTs'].values[0])
                Stats(
                    columns=4,
                    items=[(f'{"{:,} ({:,})".format(int(nsites), nvars)}',
                            'Sites (of which variants)'),
                           (f'{"{:,}".format(int(nsnvs))}', 'SNVs'),
                           (f'{"{:,}".format(int(nindels))}', 'Indels'),
                           (f'{titv}', 'Ti/Tv')])

    # Base statistics
    with report.add_section('Statistics', 'Stats'):
        tabs = Tabs()
        with tabs.add_tab(sample_id):
            if bcfstats['SN'].empty:
                p('The bcftools stats file is empty.')
            else:
                DataTable.from_pandas(
                    bcfstats['SN'].drop(columns=['id']),
                    use_index=False)
                DataTable.from_pandas(
                    bcfstats['TSTV'].drop(columns='id'),
                    use_index=False)

    # Plot tumor/normal AF if provided
    with report.add_section('Variant allele frequencies', 'VAF'):
        p(
            'Variant allele frequency (VAF) comparison between normal (x-axis)'
            ' and tumor (y-axis). The tooltips display the tumor VAF of each site.')
        tabs = Tabs()
        with tabs.add_tab(sample_id):
            if len(filters) == 0:
                p('No variants are in the VCF file.')
            else:
                DataTable.from_pandas(filtstats, use_index=False)
                vcf_df = pd.DataFrame(
                    data={
                        'Filters': filters, 'VAF': var_af,
                        'NVAF': nvar_af, 'TYPE': vtype})\
                    .query('Filters=="PASS"')
                with Grid():
                    for vt in ('SNV', 'Indel'):
                        sub_df = vcf_df.round(PRECISION).query(f'TYPE=="{vt}"')
                        if sub_df.empty:
                            p(f'No {vt}s to display.')
                            continue
                        plt = scatter_plot(
                            sub_df, 'NVAF', 'VAF', None,
                            f'Filtered Tumor vs Normal VAF ({vt})',
                            xaxis='Normal VAF', yaxis='Tumor VAF',
                            min_x=0, max_x=1, min_y=0, max_y=1)
                        plt.color = [
                            COLORS.cerulean if vt == 'SNV' else COLORS.cinnabar
                        ]
                        for s in plt.series:
                            s.symbolSize = 3
                            s.encode = {
                                'x': 'NVAF', 'y': 'VAF',
                                'itemName': 'NVAF', 'tooltip': ['VAF']}
                        # Add tooltip to facilitate reading
                        plt.tooltip = dict({'trigger': 'item'})
                        EZChart(plt, 'epi2melabs')

    # Plot mutation spectra if provided
    with report.add_section('Mutational characterisation', 'Changes'):
        p('This section shows the type of changes found and their '
            'distribution.')
        tabs = Tabs()
        with tabs.add_tab(sample_id):
            if spectra.empty:
                p('Mutation counts file is empty.')
            else:
                h6('Change distribution')
                spect = plot_spectra(spectra, sample_id, cmap=cmap)
                EZChart(spect, 'epi2melabs')
                h6('96 mutational profile')
                prof = plot_profile(spectra, sample_id, cmap=cmap)
                EZChart(prof, 'epi2melabs')

    # ClinVar variants
    if args.clinvar_vcf is not None:
        if os.path.exists(args.clinvar_vcf):
            with report.add_section('ClinVar variant annotations', 'ClinVar'):
                clinvar_docs_url = "https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/"
                p(
                    "The ",
                    a("SnpEff", href="https://pcingola.github.io/SnpEff/"),
                    " annotation tool has been used to annotate with",
                    a("ClinVar", href="https://www.ncbi.nlm.nih.gov/clinvar/"), '.'
                    " If any variants have ClinVar annotations, they will appear in a ",
                    "table below. Please note, this table excludes variants with",
                    " 'Benign' or 'Likely benign' significance, however these ",
                    "variants will appear in the VCF output by the workflow. For ",
                    "further details on the terms in the 'Significance' column,",
                    "  please visit ",
                    a("this page", href=clinvar_docs_url),
                    '.')

                # check if there are any ClinVar sites to report
                clinvar_for_report = load_clinvar_vcf(args.clinvar_vcf)
                if clinvar_for_report.empty:
                    h6('No ClinVar sites to report.')
                else:
                    DataTable.from_pandas(
                        clinvar_for_report, export=True, use_index=False)
    else:
        # Annotations were skipped
        with report.add_section('ClinVar variant annotations', 'ClinVar'):
            p(
                "This report was generated without annotations. To see"
                " them, re-run the workflow with --annotation true.")

    # write report
    report.write(args.report)
    logger.info(f"Written report to '{args.report}'.")


def argparser():
    """Create argument parser."""
    parser = wf_parser("report")

    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--read_stats", default='unknown',
        help="read statistics output from bamstats")
    parser.add_argument(
        "--vcf", default='unknown',
        help="input vcf file")
    parser.add_argument(
        "--clinvar_vcf", required=False,
        help="VCF file of variants annotated in ClinVar")
    parser.add_argument(
        "--vcf_stats", default='unknown',
        help="final vcf stats file")
    parser.add_argument(
        "--mut_spectra", default='unknown',
        help="final vcf stats file")
    parser.add_argument(
        "--read_depth", default="unknown",
        help="read coverage output from mosdepth")
    parser.add_argument(
        "--vaf_thresholds", default="0.2,0.1,0.05",
        help="read coverage output from mosdepth")
    parser.add_argument(
        "--no_germline", action="store_true",
        help="workflow run without germline call")
    parser.add_argument(
        "--normal_vcf", type=str,
        help="workflow run with germline VCF for normal sample")
    parser.add_argument(
        "--hybrid_mode_vcf_fn", type=str,
        help="workflow run hybrid typing with the given VCF")
    parser.add_argument(
        "--genotyping_mode_vcf_fn", type=str,
        help="workflow run genotyping with the given VCF")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="directory containing CSVs containing name,version.")

    return parser
