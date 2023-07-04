#!/usr/bin/env python
"""Create mod base report."""

import os

from dominate.tags import a, p
from ezcharts.components.common import add_missing_windows, CATEGORICAL
from ezcharts.components.common import fasta_idx, HSA_CHROMOSOME_ORDER
from ezcharts.components.common import MOD_ORDER
from ezcharts.components.dss import load_dml, load_dmr
from ezcharts.components.ezchart import EZChart
from ezcharts.components.modkit import load_modkit_summary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Stats, Tabs
from ezcharts.plots.karyomap import karyomap
import numpy as np
import pandas as pd
from pkg_resources import resource_filename
import sigfig as sg

# Ignoring F401 for CHROMOSOMES since used within eval expressions
from .report_utils.utils import CHROM_RENAME, CHROMOSOMES  # noqa: ABS101, F401
from .report_utils.visualizations import line_plot  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101

dss_url = 'https://bioconductor.org/packages/release' +\
    '/bioc/vignettes/DSS/inst/doc/DSS.html'


def central_value(df, x, col='values'):
    """Return the average of a given column in a region."""
    return df[(df['chrom'] == x.chr) &
              (df['pos'] >= x.start) &
              (df['pos'] < x.end)].median(numeric_only=True)[col]


# Load diff. modified regions
def intervals_to_points(dmr, faidx, total_ref_starts, colname, logger):
    """Verticalize dmr intervals."""
    # output columns
    outcols = {'chrom': CATEGORICAL, 'pos': int, 'value': float}
    # Create vertical version of the DMRs
    outdmr = pd.DataFrame(columns=outcols).astype(outcols)
    # Ensure it is sorted
    dmr = dmr.sort_values(['chrom', 'mean_pos'])
    # Add missing windows to depth intervals
    logger.info('Add missing windows...')
    dmr = add_missing_windows(dmr, faidx, value='areaStat')
    # Compute progressive position
    dmr["start"] = dmr.apply(
        lambda x: x.start + total_ref_starts[x.chrom], axis=1)
    dmr["end"] = dmr.apply(
        lambda x: x.end + total_ref_starts[x.chrom], axis=1)
    # Re-arrange the data in vertical format. Given a window:
    # chrom start end value
    # chr1 0 25000 124
    # chr1 25000 50000 12
    # It will generate a dataframe with the format:
    # chrom pos variable value
    # chr1 0 "start" 124
    # chr1 25000 "end" 124
    # chr1 50000 "start" 12
    # chr1 50000 "end" 12
    # This will allow to generate a line chart
    # representing the increase and decrease of
    # modified scores, while keeping the window size.
    logger.info('Verticalize...')
    outdmr = pd.melt(
        dmr[['chrom', 'start', 'end', colname]],
        id_vars=['chrom', colname],
        value_name='pos')
    # Reordering by chromosome, position and variable, obtaining
    # where variables can be start/end as mentioned above.
    logger.info('Re-order the data...')
    outdmr['chrom'] = outdmr['chrom'].cat.rename_categories(CHROM_RENAME)
    outdmr = outdmr\
        .eval("chr_id = chrom.map(@CHROMOSOMES)") \
        .sort_values(["chr_id", "pos", "variable"]) \
        .drop(columns=['chr_id'])
    outdmr['chrom'] = \
        outdmr.chrom.cat.remove_unused_categories()
    outdmr['chrom'] = \
        outdmr.chrom.cat.reorder_categories(
        [i for i in outdmr.chrom.unique()])
    # Drop unused cols and fix column IDs
    outdmr = outdmr\
        .drop(columns=['variable'])\
        .rename(columns={colname: 'value'})
    outdmr['value'] = outdmr['value'].fillna(0)
    return outdmr


def summary_values(df, n=100, columns=None):
    """Extract top and bottom."""
    if isinstance(columns, str):
        columns = [columns]
    dfs = []
    for column in columns:
        dfs.append(df.sort_values(
            column, ascending=False)[0:n])
        dfs.append(df.sort_values(
            column, ascending=True)[0:n])
    out_df = pd.concat(dfs).drop_duplicates()
    return out_df


def main(args):
    """Run the entry point."""
    logger = get_named_logger("report_mod")

    # Create ref lengths from faidx
    logger.info(f'Load: {args.reference_fai}')
    faidx = fasta_idx(args.reference_fai, rename=CHROM_RENAME)
    # Plot all chromosomes in the heatmap if --genome is provided
    if args.genome in ['hg38', 'hg19']:
        src = resource_filename(
            "ezcharts", f"data/reference/{args.genome}/cytoBand.txt.gz")
        names = ['chrom', 'start', 'end', 'name', 'type']
        bands_df = pd.read_csv(src, sep="\t", header=None, names=names)
        ref_lengths = bands_df\
            .eval("chrom = chrom.map(@CHROM_RENAME)") \
            .eval("chr_id = chrom.map(@CHROMOSOMES)") \
            .sort_values(["chr_id", "start"]) \
            .drop(columns=['chr_id']) \
            .groupby('chrom', observed=True, sort=False)["start"] \
            .last()
    # Otherwise, use faidx
    else:
        # Summarise lengths
        ref_lengths = faidx \
            .eval("chr_id = chrom.map(@CHROMOSOMES)") \
            .sort_values(["chr_id", "length"]) \
            .drop(columns=['chr_id']) \
            .groupby('chrom', observed=True, sort=False)["length"] \
            .last()
    # Cumulative length, used in visualization of line chart
    tot_ref_starts = ref_lengths \
        .cumsum() \
        .shift(1, fill_value=0)

    # Load other input files
    logger.info(f'Load: {args.tumor_summary}')
    t_df = load_modkit_summary(args.tumor_summary)
    t_df[['sample', 'type']] = t_df.filename.str.split('.', expand=True)[[0, 1]]
    t_df = t_df.replace('modC', '-')

    logger.info(f'Load: {args.normal_summary}')
    n_df = load_modkit_summary(args.normal_summary)
    n_df[['sample', 'type']] = n_df.filename.str.split('.', expand=True)[[0, 1]]
    n_df = n_df.replace('modC', '-')

    logger.info(f'Load: {args.dml}')
    if os.path.isfile(args.dml):
        mod = args.dml.split('.')[1]
        dmls = {mod: load_dml(args.dml, faidx=faidx)}
    else:
        dmls = {}
        for fname in os.listdir(args.dml):
            mod = fname.split('.')[1]
            dmls.update({mod: load_dml(f"{args.dml}/{fname}", faidx=faidx)})

    logger.info(f'Load: {args.dmr}')
    if os.path.isfile(args.dmr):
        mod = args.dmr.split('.')[1]
        dmrs = {mod: load_dmr(args.dmr, faidx=faidx)}
    else:
        dmrs = {}
        for fname in os.listdir(args.dmr):
            mod = fname.split('.')[1]
            dmrs.update({mod: load_dmr(f"{args.dmr}/{fname}", faidx=faidx)})

    # Instantiate the report
    report = LabsReport(
        f"{args.sample_name} | Modified bases report", "wf-somatic-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    # At-a-glance report
    logger.info('Report summary')
    with report.add_section('At a glance', 'Summary'):
        tabs = Tabs()
        with tabs.add_tab(f'{args.sample_name}'):
            t_mod_rate = t_df[
                (t_df['sample'] == args.sample_name) & (t_df['mod'] != '-')
                ].all_frac.sum()
            n_mod_rate = n_df[
                (n_df['sample'] == args.sample_name) & (n_df['mod'] != '-')
                ].all_frac.sum()
            n_dml = pd.concat(dmls)[['chrom', 'pos']]\
                .drop_duplicates().shape[0]
            n_dmr = pd.concat(dmrs)[['chrom', 'start']]\
                .drop_duplicates().shape[0]
            Stats(
                columns=4,
                items=[
                    (f'{"{:,}".format(round(t_mod_rate, 4))}',
                     'Modified rate (tumor)'),
                    (f'{"{:,}".format(round(n_mod_rate, 4))}',
                     'Modified rate (normal)'),
                    (f'{"{:,}".format(int(n_dml))}', 'DML'),
                    (f'{"{:,}".format(int(n_dmr))}', 'DMR')
                    ])
        p(
            'DML: differentially modified loci. ',
            'DMR: differentially modified regions.'
            )

    # Stat table
    logger.info('Report base stats')
    kept_columns = [
        'type', 'base', 'mod', 'threshold', 'pass_count',
        'pass_frac', 'all_count', 'all_frac', 'filename']
    with report.add_section('Modified base stats', 'Stats'):
        p(
            'Modified base ',
            a(
                'summary table',
                href='https://nanoporetech.github.io/modkit/intro_summary.html'),
            ' computed by ',
            a('modkit', href='https://nanoporetech.github.io/modkit/'),
            ' (subsample of up to 10,000 reads).')

        tabs = Tabs()
        with tabs.add_tab(f'{args.sample_name}'):
            stat_table = pd.concat([n_df, t_df])\
                .drop(columns=['sample'])[kept_columns]
            DataTable.from_pandas(
                stat_table,
                use_index=False,
                export=True,
                file_name=(f'{args.sample_name}-wf-somatic-methyl-stats'))

    # DML plot
    logger.info('Report DML')
    # Define color palette
    with report.add_section('Differentially modified loci', 'DML'):
        p(
            'Differentially modified loci (DMLs) distribution in the different',
            ' chromosomes computed by',
            a(
                'DSS',
                href=dss_url),
            '. Values shown are the median -log10(P-value) ',
            'in {} windows.'.format(sg.round(args.window_size, prefix=True, sigfigs=3)))
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for mod in MOD_ORDER:
                if mod not in dmls.keys():
                    continue
                else:
                    dml = dmls[mod]
                if '_' in mod:
                    m, s = mod.split('_')
                    tab_name = f'{args.sample_name} - {m} ({s})'
                else:
                    tab_name = f'{args.sample_name} - {mod}'
                with tabs.add_dropdown_tab(tab_name):
                    if dml.empty:
                        p('No differentially modified loci to show.')
                        continue
                    # Create intervals for the density plot
                    plt = karyomap(
                        dml, 'chrom', 'pos', 'neg_log10_p', stats='median',
                        order=HSA_CHROMOSOME_ORDER, ref_lengths=faidx)
                    # Prepare the plot.
                    EZChart(plt, theme='epi2melabs')

    # DMR plot
    logger.info('Save DMR')
    # Rename column in a readable form
    dmr_colnames = {
        'chrom': 'Chromosome',
        'start': 'Start',
        'end': 'End',
        'length': 'Length',
        'meanMethy1': 'Avg. methylation (T)',
        'meanMethy2': 'Avg. methylation (N)',
        'diff.Methy': 'Diff. methylation',
        'abs_diff': '|Diff. methylation|'
        }
    # Columns to keep
    kept_columns = [
        'Chromosome',
        'Start',
        'End',
        'Length',
        'nCG',
        'Avg. methylation (T)',
        'Avg. methylation (N)',
        'Diff. methylation',
        '|Diff. methylation|',
        'areaStat'
        ]
    # Generate the report, multiple tabs per sample, one option per mod
    # for each sample in a drop-down menu.
    with report.add_section('Differentially modified regions', 'DMR'):
        p(
            'Differentially modified regions computed by ',
            a('DSS', href=dss_url), '.')
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for mod in MOD_ORDER:
                if mod not in dmrs.keys():
                    continue
                else:
                    dmr = dmrs[mod]
                if '_' in mod:
                    m, s = mod.split('_')
                    tab_name = f'{args.sample_name} - {m} ({s})'
                else:
                    tab_name = f'{args.sample_name} - {mod}'
                with tabs.add_dropdown_tab(tab_name):
                    if dmr.empty:
                        p('No differentially modified regions to show.')
                        continue
                    # Plot DMRs as lines to show dips and constant regions
                    logger.info('Convert to vertical DMR')
                    # Turn each interval into two single consecutive points.
                    area = intervals_to_points(
                        dmr, faidx, tot_ref_starts, 'areaStat', logger)
                    logger.info('Save line plots')
                    plt = line_plot(
                        area.round(4), 'pos', 'value',
                        'chrom', f'DMR (areaStat for modification: {mod})',
                        xaxis='cumulative position', yaxis='Area Statistic')
                    for s in plt.series:
                        s.symbolSize = 3
                    p('Differentially modified regions along the genome.')
                    EZChart(plt, theme='epi2melabs')

                    # Create summary table
                    dmr['abs_diff'] = np.abs(dmr['diff.Methy'])
                    # Save top- and bottom-N values by multiple cols
                    datatable = summary_values(
                        dmr,
                        n=100,
                        columns=['length', 'areaStat', 'diff.Methy'])\
                        .rename(columns=dmr_colnames).drop(
                            columns=['mean_pos'])[kept_columns]
                    p(
                        'Top and bottom 100 differentially modified regions for',
                        ' size, areaStat and diff.Methy.'
                        )
                    DataTable.from_pandas(
                        datatable,
                        use_index=False)

    # write report
    report.write(args.report)
    logger.info(f"Written report to '{args.report}'.")


def path_check(infile):
    """Check that input file exists and return str type."""
    if not os.path.exists(infile):
        raise FileNotFoundError(f"File {infile} not found.")
    return infile


def argparser():
    """Create argument parser."""
    parser = wf_parser("report_mod")

    parser.add_argument("report", help="Report mod file")
    parser.add_argument(
        "--tumor_summary", default=None,
        help="input tumor modification ratio",
        type=lambda x: path_check(x))
    parser.add_argument(
        "--normal_summary", default=None,
        help="input normal modification ratio",
        type=lambda x: path_check(x))
    parser.add_argument(
        "--dml", default=None,
        help="DSS differentially modified loci",
        type=lambda x: path_check(x))
    parser.add_argument(
        "--dmr", default=None,
        help="DSS differentially modified regions",
        type=lambda x: path_check(x))
    parser.add_argument(
        "--reference_fai", default=None,
        help="Reference fai file.",
        type=lambda x: path_check(x))
    parser.add_argument(
        "--genome", default=None, required=False,
        help="Reference genome id.")
    parser.add_argument(
        "--window_size", default=1000000, type=int,
        help="Window size for visualization DML heatmap.")
    parser.add_argument(
        "--sample_name", default="SAMPLE",
        help="Sample name.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="directory containing CSVs containing name,version.")

    return parser
