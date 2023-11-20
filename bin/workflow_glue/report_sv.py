#!/usr/bin/env python
"""Create workflow report."""
from dominate.tags import p
from ezcharts.components.common import CATEGORICAL
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
from ezcharts.plots.ideogram import ideogram
import numpy as np
import pandas as pd
import pysam

from .report_utils.utils import compare_max_axes  # noqa: ABS101
from .report_utils.utils import COLORS, PRECISION  # noqa: ABS101
from .report_utils.visualizations import hist_plot  # noqa: ABS101
from .report_utils.visualizations import scatter_plot  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101

chroms_37 = [str(x) for x in range(1, 23)] + ['X', 'Y']


def read_vcf(fname):
    """Read input VCF as pandas dataframe."""
    vcf = pysam.VariantFile(fname)
    sample_name = vcf.header.samples[0]
    cols = {
        'CHROM': CATEGORICAL,
        'POS': int,
        'ID': str,
        'REF': str,
        'ALT': str,
        'FILTER': str,
        'VAF': float,
        'NVAF': float,
        'SVLEN': int,
        'SVTYPE': CATEGORICAL,
        }
    df = pd.DataFrame(columns=cols).astype(cols)
    try:
        dropped = 0
        idx = 0

        # Process one entry at time using pysam
        for n, rec in enumerate(vcf):
            # Ignore non-chromosomal SVs
            if rec.chrom.replace('chr', '') not in chroms_37:
                dropped += 1
                continue
            # Ignore filtered SVs
            if rec.filter.keys()[0] != 'PASS':
                dropped += 1
                continue
            vaf = float(rec.samples[sample_name]['AF'])
            nvaf = float(rec.samples[sample_name]['NAF'])
            sv_len = rec.info.get('SVLEN')
            if not sv_len:
                sv_len = rec.info.get('SVINSLEN')
            if not sv_len:
                sv_len = np.nan
            # Deal with multiple allele lines by treating them as
            # independent sites, changing the site ID
            for m, alt in enumerate(rec.alts):
                if alt == '<DEL>':
                    sv_type = 'DEL'
                elif alt == '<INS>':
                    sv_type = 'INS'
                else:
                    sv_type = 'BND'
                rec_id = rec.id if len(rec.alts) > 1 else f"{rec.id}_{m}"
                new_df = pd.DataFrame(
                    data={
                        'CHROM': rec.chrom.replace('chr', ''),
                        'POS': rec.pos,
                        'ID': rec_id,
                        'REF': rec.ref,
                        'ALT': alt,
                        'FILTER': ','.join(rec.filter.keys()),
                        'VAF': vaf,
                        'NVAF': nvaf,
                        'SVLEN': sv_len,
                        'SVTYPE': sv_type,
                    }, index=[idx]
                )
                idx += 1
            df = pd.concat([df, new_df])
        # Use appropriate typing
        df.reset_index(drop=True)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    return df, dropped


def get_sv_summary_table(vcf_df):
    """Aggregate summary info for SV calls per type."""
    return vcf_df.groupby('SVTYPE').agg(**{
        'Count': ('POS', 'count'),
        'Min. Length': ('SVLEN', lambda x: np.min(x.abs())),
        'Ave. Length': ('SVLEN', lambda x: np.median(x.abs())),
        'Max. Length': ('SVLEN', lambda x: np.max(x.abs()))}).transpose()


def sv_stats(vcf_data):
    """Display base stats."""
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                inserts = vcf_df.loc[vcf_df['SVTYPE'] == 'INS']
                delets = vcf_df.loc[vcf_df['SVTYPE'] == 'DEL']
                breakends = vcf_df.loc[vcf_df['SVTYPE'] == 'BND']
                Stats(
                    columns=4,
                    items=[
                        (f'{f"{inserts.shape[0]}"}',
                         'Number of Insertions'),
                        (f'{f"{delets.shape[0]}"}',
                         'Number of Deletions'),
                        (f'{f"{breakends.shape[0]}"}',
                         'Other SVs'),
                        (f'{f"{len(vcf_df.CHROM.unique())}"}',
                         'Chromosomes with SVs'),
                    ])


def sv_size_plots(vcf_data):
    """Plot size distributions of SV calls per type."""
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                with Grid():
                    inserts = vcf_df[vcf_df['SVTYPE'] == 'INS']
                    delets = vcf_df[vcf_df['SVTYPE'] == 'DEL']
                    # Force deletion's length to be negative
                    delets = delets.eval('SVLEN=abs(SVLEN)')
                    # Define a reasonable bin number and size based on the number
                    # of variants in the dataset
                    n_bins = int(np.ceil(2*np.cbrt(vcf_df.shape[0])))
                    max_x = int(vcf_df.eval('SVLEN=abs(SVLEN)').SVLEN.max())
                    binwidth = int(vcf_df.eval('SVLEN=abs(SVLEN)').SVLEN.max() / n_bins)
                    binrange = [0, 0]
                    while binrange[-1] < max_x:
                        binrange[-1] += binwidth
                    # Compute max Y axis
                    if inserts.empty or delets.empty:
                        max_y = None
                    else:
                        max_y = compare_max_axes(
                            inserts, delets, 'SVLEN', binwidth=binwidth,
                            ptype='hist', buffer=1.2, precision=1)
                    if inserts.shape[0] > 0:
                        plt = hist_plot(
                            inserts, 'SVLEN', 'Insertion lengths', no_stats=True,
                            xaxis='abs. Length', yaxis='Count', rounding=0,
                            color=COLORS.cinnabar, binwidth=binwidth,
                            binrange=binrange, max_y=max_y)
                        EZChart(plt, 'epi2melabs')
                    else:
                        p('No insertions to show.')
                    if delets.shape[0] > 0:
                        plt = hist_plot(
                            delets, 'SVLEN', 'Deletion lengths', no_stats=True,
                            xaxis='abs. Length', yaxis='Count', rounding=0,
                            color=COLORS.cerulean, binwidth=binwidth,
                            binrange=binrange, max_y=max_y)
                        EZChart(plt, 'epi2melabs')
                    else:
                        p('No deletions to show.')


def karyoplot(vcf_data, args):
    """Karyogram plot."""
    p("Chromosomal hotspots of structural variation.")
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                colors = {
                    "INS": COLORS.cinnabar,
                    "DEL": COLORS.cerulean
                }
                # Extract relevant columns
                df = vcf_df[vcf_df['SVTYPE'] != 'BND']
                df = df[['CHROM', 'POS', 'ID', 'SVLEN', 'SVTYPE']]
                # set SV length to absolute value
                df['SVLEN'] = df['SVLEN'].abs()
                df.columns = ['chr', 'start', 'name', 'length', 'color']
                # Define end point to start+length...
                # ... unless its too small, in which case set it
                # to at least 100Kb to allow visualization of the hotspots.
                df['end'] = df.apply(lambda x: x.start + max(x.length, 1e5), axis=1)
                # Define color mappings
                df['color'] = df['color'].map(colors).fillna(COLORS.black)
                # df['color'] = df['color'].cat.rename_categories(colors)
                df = df[['chr', 'start', 'end', 'name', 'color']]
                # Prepare the ideogram
                plt = ideogram(blocks=df, genome=args.genome)
                EZChart(plt, height='600px', width='90%', theme='epi2melabs')
        p("""Red: Insertion""")
        p("""Blue: Deletion""")


def main(args):
    """Run the entry point."""
    # Input all VCFs
    vcf_data = []
    for index, sample_vcf in enumerate(args.vcf):
        vcf_df, dropped_svs = read_vcf(sample_vcf)
        vcf_data.append((index, sample_vcf.split('.')[0], vcf_df))

    # Create report file
    report = LabsReport(
        f"{args.vcf[0].split('.')[0]} | Structural variants statistics",
        "wf-somatic-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    with report.add_section('At a glance', 'Summary'):
        p(
            "This section displays a description"
            " of the variant calls made by nanomonsv.")
        sv_stats(vcf_data)
        if dropped_svs > 0:
            p(
                f"A total of {dropped_svs} SVs were not considered because"
                " either soft filtered or on small contigs.")

    with report.add_section('Variant calling results', 'Variants'):
        p(
            "This section displays summary statistics"
            " of the variant calls made by nanomonsv.")
        tabs = Tabs()
        for (index, sample_name, vcf_df) in vcf_data:
            with tabs.add_tab(sample_name):
                if vcf_df.empty:
                    p("The workflow found no structural variants to report.")
                else:
                    DataTable.from_pandas(get_sv_summary_table(vcf_df))

    with report.add_section('Karyogram', 'Karyogram'):
        karyoplot(vcf_data, args)

    with report.add_section('Size distribution', 'Size'):
        sv_size_plots(vcf_data)

    with report.add_section('Variant allele frequencies', 'VAF'):
        tabs = Tabs()
        for (index, sample_name, vcf_df) in vcf_data:
            with tabs.add_tab(sample_name):
                if vcf_df.empty:
                    p('The workflow found no structural variants to report.')
                else:
                    plt = scatter_plot(
                        vcf_df.round(PRECISION), 'NVAF', 'VAF', None,
                        'Tumor vs Normal variant allele frequency (VAF)',
                        xaxis='Normal VAF', yaxis='Tumor VAF',
                        min_x=0, max_x=1, min_y=0, max_y=1)
                    for s in plt.series:
                        s.symbolSize = 3
                    EZChart(plt, 'epi2melabs')

    #
    # write report
    #
    report.write(args.report)


def argparser():
    """Create argument parser."""
    parser = wf_parser("report_sv")
    parser.add_argument(
        "report",
        help="Report output file.")
    parser.add_argument(
        "--vcf",
        nargs='+',
        required=True)
    parser.add_argument(
        "--genome",
        default='hg38',
        required=False)
    parser.add_argument(
        "--eval_results",
        nargs='+',
        required=False)
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--params-hidden", default="",
        help="Comma delimited list of keys to hide from parameters table")
    parser.add_argument(
        "--versions", required=True,
        help="directory contained CSVs containing name,version.")

    return parser
