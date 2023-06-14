#!/usr/bin/env python
"""Create SNV report."""

import os

from aplanat.parsers.bcfstats import parse_bcftools_stats_multi
from dominate.tags import h6, p
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
from ezcharts.plots import util
from ezcharts.plots.categorical import barplot
from ezcharts.plots.distribution import histplot
import numpy as np
import pandas as pd
import pysam

from .report_sv import scatter_plot  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101

# Global variables
Colors = util.Colors
# Number of digits
PRECISION = 4

# Mutation profile palette to match the COSMIC
# patterns.
cmap = {
    'C>A': Colors.cerulean,
    'C>G': Colors.black,
    'C>T': Colors.cinnabar,
    'T>A': Colors.grey70,
    'T>C': Colors.medium_spring_bud,
    'T>G': Colors.fandango,
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
        flt.append(rec.filter.keys()[0])
        vaf.append(rec.samples[samplename]['AF'])
        naf.append(rec.samples[samplename]['NAF'])
        if len(rec.alts[0]) == len(rec.ref) and len(rec.ref) == 1:
            vtype.append('SNV')
        else:
            vtype.append('Indel')
    return samplename, flt, vaf, naf, vtype


def vaf_plot(vaf, title, xaxis="", yaxis="", color=None):
    """Plot the allele frequencies for tumor sample."""
    # Plot the histogram of the allele frequencies in the tumor and normal
    plt = histplot(
        vaf, binrange=(0, 1), binwidth=0.01,
        stat='proportion', color=color)
    plt.title = dict(text=title)
    # Update axis parameters
    plt.xAxis.name = xaxis
    # Define y-axis
    plt.yAxis = dict(
        name=yaxis,
        max=1.0,
        min=0.0,
        type='value',
        axisTick=dict(alignWithLabel=True),
        data=[dict(value=i) for i in plt.dataset[0].source[:, 2].round(1)])
    return plt


def filt_stats(filters, vaf, naf, thresholds=[0.2, 0.1, 0.05]):
    """Plot the filtering stats."""
    df = pd.DataFrame({'Filter': filters, 'Tumor_AF': vaf, 'Normal_AF': naf})
    # Filtering table
    summary_table = {
        'N sites': [df.count()['Filter'],
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


def plot_spectra(spectra, sample):
    """Plot the mutation spectra."""
    # Plot change spectrum
    df = spectra[['Change', sample]]
    # Count by change type
    tots = df.groupby('Change').sum().reset_index()
    tots.columns = ['Change', sample]
    # Create bar plot
    plot = barplot(data=tots, x='Change', y=sample)
    plot.color = list(cmap.values())
    plot.xAxis.name = 'Change'
    plot.yAxis.name = "Counts"
    for s in plot.series:
        s.colorBy = 'data'
    plot.tooltip = dict({'trigger': 'item'})
    return plot


def plot_profile(df, sample):
    """Plot the mutation profile."""
    # Ensure sorting
    df = df.sort_values(["Change", "Flanks"])

    # Create base barplot
    plt = barplot(
        data=df.assign(order=range(df.shape[0])),
        x="order", y=sample,
        hue='Change', dodge=False)

    # Replace NaN with zeros to avoid ValueError('invalid input Character)
    plt.dataset[0].source = np.nan_to_num(plt.dataset[0].source, nan=0)

    # Define appropriate palette
    plt.color = list(cmap.values())

    # Definition of the x-axis
    # in distinct axes
    flanks_axis = dict(
        data=[dict(value=x) for x in df["Flanks"]],
        axisLabel=dict(
            rotate=90,
            fontSize=10,
            interval=0,  # this forces echarts to show all tick labels
        ),
        axisTick=dict(alignWithLabel=True),
    )
    change_axis = dict(
        position="bottom",
        # make sure the df is sorted by `change` (otherwise the labels will be wrong)
        data=[dict(value=x) for x in df["Change"].unique()],
        # tick length and label margin below might need tweaking
        axisTick=dict(length=50),
        axisLabel=dict(margin=40),
        splitLine=dict(show=True),
    )

    # replace the original axis with the new ones
    plt.xAxis = [flanks_axis, change_axis]
    plt.yAxis.name = 'Count'

    # Add tooltip to facilitate reading
    tmp_df = plt.dataset[0].source
    tmp_sum = np.array([tmp_df[:, 0], np.sum(tmp_df[:, 1:], axis=1)]).T
    plt.add_series(
        dict(
            type="scatter",
            name="Mean",
            data=[dict(value=tmp_sum[i, :]) for i in range(0, tmp_sum.shape[0])],
            symbolSize=0,
            tooltip=dict(formatter='{c}')
        ))
    plt.tooltip = dict({'trigger': 'item'})
    return plt


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
        bcfstats = parse_bcftools_stats_multi(
            [args.vcf_stats],
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
        tabs = Tabs()
        with tabs.add_tab(sample_id):
            if bcfstats['SN'].empty:
                p('The bcftools stats file is empty.')
            else:
                bcfstats['SN'].columns = [
                    i.replace('SNP', 'SNV').replace('MNP', 'MNV')
                    for i in bcfstats['SN'].columns
                ]
                titv = bcfstats['TSTV']['ts/tv'].values[0]
                nsites = bcfstats['SN']['records'].values[0]
                nsnvs = bcfstats['SN']['SNVs'].values[0]
                nindels = bcfstats['SN']['indels'].values[0]
                Stats(
                    columns=4,
                    items=[(f'{"{:,}".format(int(nsites))}', 'Variants'),
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
                    bcfstats['SN'].drop(columns=['sample', 'samples']),
                    use_index=False)
                DataTable.from_pandas(
                    bcfstats['TSTV'].drop(columns='sample'),
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
                        plt.color = [Colors.cerulean if vt == 'SNV' else Colors.green]
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
                spect = plot_spectra(spectra, sample_id)
                EZChart(spect, 'epi2melabs')
                h6('96 mutational profile')
                prof = plot_profile(spectra, sample_id)
                EZChart(prof, 'epi2melabs')

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
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="directory containing CSVs containing name,version.")

    return parser
