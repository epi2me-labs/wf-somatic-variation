#!/usr/bin/env python
"""Plot QC metrics."""

from dominate.tags import a, p
from ezcharts import lineplot
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
from ezcharts.plots import util
from ezcharts.plots.distribution import histplot
import numpy as np
import pandas as pd
from pandas.api import types as pd_types
from seaborn._statistics import Histogram

from .util import get_named_logger, wf_parser  # noqa: ABS101

# Global variables
COLORS = util.Colors
CATEGORICAL = pd_types.CategoricalDtype(ordered=True)
CHROMOSOMES = {str(i): i for i in range(1, 23)}
CHROMOSOMES.update({f"chr{i}": int(i) for i in CHROMOSOMES})
CHROMOSOMES.update({'X': 23, 'Y': 24, 'chrX': 23, 'chrY': 24})


# File loaders
def fasta_idx(faidx):
    """Read faidx for the reference fasta."""
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "length": int,
        "offset1": int,
        "offset2": int,
        "offset3": int,
    }
    try:
        df = pd.read_csv(
            faidx,
            sep="\t",
            names=relevant_stats_cols_dtypes.keys(),
            dtype=relevant_stats_cols_dtypes,
        ).drop(columns=['offset1', 'offset2', 'offset3'])
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=relevant_stats_cols_dtypes)
    return df


def process_fastcat(fastcat_file):
    """Load fastcat results into dataframe."""
    relevant_stats_cols_dtypes = {
        "name": str,
        "sample_name": CATEGORICAL,
        "ref": CATEGORICAL,
        "coverage": float,
        "ref_coverage": float,
        "read_length": int,
        "mean_quality": float,
        "acc": float,
    }
    try:
        d = pd.read_csv(
            fastcat_file,
            sep="\t",
            header=0,
            usecols=relevant_stats_cols_dtypes,
            dtype=relevant_stats_cols_dtypes
            )
    except pd.errors.EmptyDataError:
        d = pd.DataFrame(columns=relevant_stats_cols_dtypes)
    return d


def make_breaks(minval, maxval, winsize):
    """Create intervals inclusive of last value."""
    # Create the intervals
    breaks = list(range(minval, maxval, winsize))
    # Check that the max value is the chr length
    if breaks[-1] > maxval:
        breaks[-1] = maxval
    if breaks[-1] < maxval:
        breaks += [maxval]
    return breaks


def add_missing_windows(intervals, faidx, winsize=25000):
    """Add missing windows to depth dataframe.

    If the depth refers to a narrow region of the chromosome, this will
    cause to generate a dataframe with either a single-entry, or only
    few entries. This causes the production of poor quality or impossible to
    understand plots.
    This function fixes this by adding the missing leading or trailing
    intervals, matching the full length of the chromosome. Moreover, it will
    check for the presence of gaps in the region, adding back the missing
    points. All added intervals are set to 0-coverage values.
    """
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "start": int,
        "end": int,
        "depth": float,
    }
    # Get unique chromosomes
    chrs = intervals.chrom.unique().tolist()
    # New intervals
    final_intervals = []
    # Add missing windows for each chromosome
    for chr_id in chrs:
        # Get chromosome entries
        chr_entry = intervals.loc[intervals['chrom'] == chr_id].reset_index(drop=True)
        # First window start
        minval = chr_entry.start.min()
        # Final window end
        maxval = chr_entry.end.max()
        # Chromosome length
        chr_len = faidx[faidx['chrom'] == chr_id].length.max()
        # If the minimum value is != 0, add intermediate windows
        if minval != 0:
            # Create the breaks for the intervals
            breaks = make_breaks(0, minval, winsize)
            # Add leading intervals
            final_intervals.append(
                pd.DataFrame(data={
                    'chrom': chr_id,
                    'start': breaks[0:-1],
                    'end': breaks[1:],
                    'depth': 0}))
        # Append precomputed intervals, checking for breaks
        # in between regions.
        for idx, region in chr_entry.iterrows():
            # To DF
            region = region.to_frame().T
            # Add the first window as default
            if idx == 0:
                final_intervals.append(region)
                continue
            # If the new window start is not the end of the previous,
            # then create new regions
            if region.start.min() != final_intervals[-1].end.max():
                # Create intervals
                breaks = make_breaks(
                    final_intervals[-1].end.max(),
                    region.start.min(),
                    winsize)
                # Add heading    intervals
                final_intervals.append(
                    pd.DataFrame(data={
                        'chrom': chr_id,
                        'start': breaks[0:-1],
                        'end': breaks[1:],
                        'depth': 0}))
            final_intervals.append(region)
        # If the max value is less than the chromosome length, add these too
        if maxval < chr_len:
            # Create the breaks for the intervals
            breaks = make_breaks(maxval, chr_len, winsize)
            # Create vectors to populate the DF
            chr_vec = [chr_id] * (len(breaks) - 1)
            depths = [0] * (len(breaks) - 1)
            # Add tailing intervals
            final_intervals.append(
                pd.DataFrame(data={
                    'chrom': chr_vec,
                    'start': breaks[0:-1],
                    'end': breaks[1:],
                    'depth': depths}))
    return pd.concat(final_intervals).astype(
        relevant_stats_cols_dtypes).reset_index(drop=True)


def load_mosdepth(mosdepth_file, faidx, winsize):
    """Load mosdepth results into dataframe."""
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "start": int,
        "end": int,
        "depth": float,
    }
    try:
        df = pd.read_csv(
            mosdepth_file,
            sep="\t",
            names=relevant_stats_cols_dtypes.keys(),
            dtype=relevant_stats_cols_dtypes
        )
        df = add_missing_windows(df, faidx, winsize)
        df = df.loc[df['chrom'].isin(CHROMOSOMES)]
        df = (
            df.eval("mean_pos = (start + end) / 2")
            .eval("step = end - start")
            .reset_index()
            )
        # Sort column by chromosome number
        df = df.eval("chr_id = chrom.map(@CHROMOSOMES)") \
            .sort_values(["chr_id", "start"]) \
            .drop(columns="chr_id")
        # Reordering
        df['chrom'] = \
            df.chrom.cat.remove_unused_categories()
        df['chrom'] = \
            df.chrom.cat.reorder_categories(
            [i for i in df.chrom.unique()])
        # Extract ref lengths
        ref_lengths = df.groupby(
            "chrom", observed=True, sort=False)["end"].last()
        total_ref_starts = ref_lengths.cumsum().shift(1, fill_value=0)
        # Add cumulative depth
        df["total_mean_pos"] = df.apply(
            lambda x: x.mean_pos + total_ref_starts[x.chrom], axis=1)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=relevant_stats_cols_dtypes) \
            .astype({"chrom": CATEGORICAL})
    return df


def load_mosdepth_summary(summary_file):
    """Load mosdepth results into dataframe."""
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "length": int,
        "bases": int,
        "mean": float,
        "min": int,
        "max": int
    }
    try:
        df = pd.read_csv(
            summary_file,
            sep="\t",
            usecols=relevant_stats_cols_dtypes.keys(),
            dtype=relevant_stats_cols_dtypes
            )
        df_tot = df[~df['chrom'].str.contains("_region")]
        df_reg = df[df['chrom'].str.contains("_region")]
    except pd.errors.EmptyDataError:
        df_tot = pd.DataFrame(columns=relevant_stats_cols_dtypes)
        df_reg = pd.DataFrame(columns=relevant_stats_cols_dtypes)
    return df_tot, df_reg


def load_flagstat(flagstat_file):
    """Load and process flagstat."""
    relevant_stats_cols_dtypes = {
        "ref": CATEGORICAL,
        "sample_name": CATEGORICAL,
        "total": int,
        "primary": int,
        "secondary": int,
        "supplementary": int,
        "unmapped": int,
        "qcfail": int,
        "duplicate": int,
    }
    try:
        df = pd.read_csv(
            flagstat_file, sep="\t",
            usecols=relevant_stats_cols_dtypes,
            dtype=relevant_stats_cols_dtypes)
        conditions = [(df['ref'] == '*'), (df['ref'] != '*')]
        choices = ['Unmapped', 'Mapped']
        df['Status'] = np.select(conditions, choices, default='Unmapped')
        df.pop('ref')
        movecol(df, "Status", 0)
        movecol(df, "sample_name", 0)
        df = df.groupby(['sample_name', 'Status']).sum().reset_index()
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    return df


# Data manipulation
def movecol(df, colname, pos):
    """Move a column in a dataframe."""
    df.insert(pos, colname, df.pop(colname))


def process_chr_sizes(fai):
    """Load the chromosome sizes."""
    sizes = pd.read_csv(
        fai,
        sep="\t",
        header=None, names=['chr', 'length', 'offset', 'b1', 'b2'])
    # Add percentage and return
    sizes['percent'] = sizes['length']/sizes['length'].sum()*100
    return sizes


def compute_n50(lengths):
    """Compute read N50."""
    # Sort the read lengths
    sorted_l = np.sort(lengths)[::-1]
    # Generate cumsum
    cumsum = np.cumsum(sorted_l)
    # Get lowest cumulative value >= (total_length/2)
    n50 = sorted_l[cumsum >= cumsum[-1]/2][0]
    return n50


def hist_max(variable_data, binwidth=None):
    """Compute max value to set in a plot."""
    estimate_kws = dict(
        stat='count',
        bins='auto',
        binwidth=binwidth,
        binrange=None,
        discrete=None,
        cumulative=False,
    )
    estimator = Histogram(**estimate_kws)
    heights, edges = estimator(variable_data, weights=None)
    return max(heights)


def compare_max_axes(df1, df2, col, ptype='val', binwidth=None):
    """Compute max value to set in a plot."""
    if ptype == 'hist':
        v1max = hist_max(df1[col].dropna(), binwidth=binwidth)
        v2max = hist_max(df2[col].dropna(), binwidth=binwidth)
    else:
        v1max = df1[col].max()
        v2max = df2[col].max()
    return np.ceil(max(v1max, v2max) * 1.1)


def add_cumulative(df):
    """Compute cumulative length."""
    ref_lengths = df.groupby("chrom", observed=True)["stop"].last()
    total_ref_starts = ref_lengths.cumsum().shift(1, fill_value=0)
    df["total_mean_pos"] = df.groupby(
        "chrom", observed=True, group_keys=False
    )["mean_pos"].apply(lambda s: s + total_ref_starts[s.name])
    return df


# Histogram
def hist_plot(
        df, col, title, xaxis='', yaxis='', rounding=None,
        n50=None, color=None, binwidth=None, binrange=None, bins='auto',
        max_y=None, max_x=None, min_x=None, min_y=None):
    """Make a histogram of given parameter."""
    histogram_data = df[col].values

    plt = histplot(
        data=histogram_data,
        bins=bins,
        binwidth=binwidth,
        binrange=binrange,
        color=color)
    if isinstance(rounding, int):
        meanv = df[col].mean().round(rounding)
        medianv = df[col].median().round(rounding)
    else:
        meanv = df[col].mean()
        medianv = df[col].median()

    if n50 is None:
        plt.title = dict(
            text=title,
            subtext=(
                f"Mean: {meanv}. "
                f"Median: {medianv}. "
            ),
        )
    else:
        plt.title = dict(
            text=title,
            subtext=(
                f"Mean: {meanv}. "
                f"Median: {medianv}. "
                f"N50: {n50}. "
            ),
        )

    # Add mean and median values (Thanks Julian!)
    plt.add_series(
        dict(
            type="line",
            name="Mean",
            data=[dict(value=[meanv, 0]), dict(value=[meanv, max_y])],
            itemStyle=(dict(color=COLORS.sandstorm)),
            symbolSize=0
        )
    )
    plt.add_series(
        dict(
            type="line",
            name="Median",
            data=[dict(value=[medianv, 0]), dict(value=[medianv, max_y])],
            itemStyle=(dict(color=COLORS.fandango)),
            symbolSize=0
        )
    )
    # Add N50 if required
    if n50 is not None:
        plt.add_series(
            dict(
                type="line",
                name="N50",
                data=[dict(value=[n50, 0]), dict(value=[n50, max_y])],
                itemStyle=(dict(color=COLORS.cinnabar)),
                symbolSize=0
            )
        )
    # Change color if requested
    if color is not None:
        plt.color = [color]
    # Customize X-axis
    plt.xAxis.name = xaxis
    plt.yAxis.name = yaxis
    plt.yAxis.nameGap = 0
    if max_x is not None:
        plt.xAxis.max = max_x
    if min_x is not None:
        plt.xAxis.min = min_x
    if max_y is not None:
        plt.yAxis.max = max_y
    if min_y is not None:
        plt.yAxis.min = min_y
    return plt


# Line plot
def line_plot(
        df, x, y, hue, title, add_mean=None,
        xaxis='', yaxis='', max_y=None):
    """Make a histogram of given parameter."""
    plt = lineplot(
        data=df,
        x=x,
        y=y,
        hue=hue
    )
    # Change axes names
    plt.xAxis.name = xaxis
    plt.yAxis.name = yaxis
    plt.yAxis.nameGap = 0
    # Set y limit
    if max_y is not None:
        plt.yAxis.max = max_y
    # Remove markers
    for s in plt.series:
        s.showSymbol = False
    # Change title and add mean as subtitle and horiz. line, if provided
    if add_mean is not None:
        plt.title = {
            "text": title,
            "subtext": f"Mean coverage: {round(add_mean, 2)}."
            }
        max_x_val = df.max()[x]
        plt.add_series(
            dict(
                type="line",
                name="Mean coverage",
                data=[dict(value=[0, add_mean]), dict(value=[max_x_val, add_mean])],
                itemStyle=(dict(color=COLORS.black)),
                lineStyle=(dict(type='dashed')),
                symbolSize=0
            )
        )
    else:
        plt.title = {"text": title}
    return plt


# Reporting function
def populate_report(report, args, **kwargs):
    """Populate the report with the different sections."""
    # Extract data used
    stats_df_t = kwargs['stats_df_t']
    stats_df_n = kwargs['stats_df_n']
    flags_df_t = kwargs['flags_df_t']
    flags_df_n = kwargs['flags_df_n']
    depth_su_t = kwargs['depth_su_t']
    depth_su_n = kwargs['depth_su_n']
    depth_df_t = kwargs['depth_df_t']
    depth_df_n = kwargs['depth_df_n']
    logger = kwargs['logger']

    # Average total coverage
    logger.info('Compute average cvg...')
    t_cov = depth_su_t.loc[depth_su_t['chrom'] == 'total', 'mean'].values[0]
    n_cov = depth_su_n.loc[depth_su_n['chrom'] == 'total', 'mean'].values[0]

    # Define collected pandas coverage dataframe
    logger.info('Collect cvg dataframes...')
    t_pass = 'PASS' if t_cov > args.tumor_cov_threshold else 'FAIL'
    n_pass = 'PASS' if n_cov > args.normal_cov_threshold else 'FAIL'

    # Compute N50s for later use
    logger.info('Compute N50...')
    t_n50 = compute_n50(stats_df_t['read_length'].values)
    n_n50 = compute_n50(stats_df_n['read_length'].values)

    # Save coverages
    t_cov = depth_su_t.loc[depth_su_t['chrom'] == 'total', 'mean'].values[0]
    n_cov = depth_su_n.loc[depth_su_n['chrom'] == 'total', 'mean'].values[0]

    # Start structuring the report
    logger.info('Start reporting.')
    logger.info('Saving summary section...')
    with report.add_section('At a glance', 'Description'):
        p(
            """
            This report contains visualisations of alignment statistics for paired
            tumor/normal samples that can help in understanding the results from the
            wf-somatic-variation. Each section contains different plots or tables,
            and in general the results are broken down by sample or the reference file
            to which alignments were made. You can quickly jump to an individual
            section with the links in the header bar.
            """
            )
        # Create tabs for each sample
        tabs = Tabs()
        with tabs.add_tab(args.sample_id):
            Stats(
                columns=2,
                items=[
                    (f'{"{:,}".format(len(stats_df_t.index))}', 'Tumor Total Reads'),
                    (f'{"{:,}".format(len(stats_df_n.index))}', 'Normal Total Reads'),
                    (
                        f'{"{:,}".format(t_n50)} bp',
                        'Tumor Read N50'),
                    (
                        f'{"{:,}".format(n_n50)} bp',
                        'Normal Read N50'),
                    (
                        f"{t_cov}x",
                        'Tumor mean coverage'),
                    (
                        f"{n_cov}x",
                        'Normal mean coverage')
                ])

    # Add coverage thresholds passing/failing
    logger.info('Saving filtering section...')
    with report.add_section('Base statistics', 'Stats'):
        p(
            f"""
            The aligned bam files have been tested for coverage thresholds of
            {args.tumor_cov_threshold}x for the tumor and {args.normal_cov_threshold}x
            for the normal sequences.
            """
            )
        # Create tabs for each sample
        df = pd.DataFrame({
            'Sample': [args.sample_id, args.sample_id],
            'Type': ['Tumor', 'Normal'],
            'Total Reads': [len(stats_df_t.index), len(stats_df_n.index)],
            'Median Read Length': [int(stats_df_t['read_length'].median()),
                                   int(stats_df_n['read_length'].median())],
            'Read N50': [t_n50, n_n50],
            'Min chrom. coverage': [depth_su_t['mean'].min(), depth_su_n['mean'].min()],
            'Mean chrom. coverage': [t_cov, n_cov],
            'Max chrom. coverage': [depth_su_t['mean'].max(), depth_su_n['mean'].max()],
            'Pass threshold': [t_pass, n_pass],
        })
        DataTable.from_pandas(df, use_index=False)

    # Add comparison of read lengths
    logger.info('Saving read distribution...')
    with report.add_section('Read length distribution', 'Read Length'):
        tabs = Tabs()
        with tabs.add_tab(args.sample_id):
            with Grid():
                max_y = compare_max_axes(
                    stats_df_t, stats_df_n, 'read_length',
                    ptype='hist', binwidth=1000)
                inputs = (
                    (stats_df_t, t_n50, "Tumor Read Length", COLORS.cerulean),
                    (stats_df_n, n_n50, "Normal Read Length", COLORS.green)
                )
                for (df, n50, header, color) in inputs:
                    plt = hist_plot(
                        df, 'read_length', header, xaxis='Read length', color=color,
                        yaxis='Number of reads', n50=n50, binwidth=1000, max_y=max_y,
                        rounding=0)
                    EZChart(plt, 'epi2melabs')
            p("""Red: read N50; Yellow: mean length; Purple: median length.""")

    # Add comparison of read quality
    logger.info('Saving mean read quality...')
    with report.add_section('Mean read quality', 'Read Quality'):
        tabs = Tabs()
        with tabs.add_tab(args.sample_id):
            with Grid():
                max_y = compare_max_axes(
                    stats_df_t, stats_df_n, 'mean_quality',
                    ptype='hist', binwidth=0.5)
                inputs = (
                    (stats_df_t, "Tumor Mean Read Quality", COLORS.cerulean),
                    (stats_df_n, "Normal Mean Read Quality", COLORS.green)
                )
                for (df, header, color) in inputs:
                    plt = hist_plot(
                        df, 'mean_quality', header, xaxis='Mean read quality',
                        color=color, yaxis='Number of reads', binwidth=0.5, max_y=max_y,
                        rounding=1)
                    EZChart(plt, 'epi2melabs')
            p("""Yellow: mean length; Purple: median length.""")

    # Create the alignment stats
    logger.info('Saving alignment statistics...')
    with report.add_section('Alignment statistics', 'Alignments'):
        tabs = Tabs()
        with tabs.add_tab(args.sample_id):
            flags_df_t.insert(0, 'Type', 'Tumor')
            flags_df_n.insert(0, 'Type', 'Normal')
            data_table = pd.concat((flags_df_t, flags_df_n)).drop(columns=['unmapped'])
            DataTable.from_pandas(data_table, use_index=False)

            # Add accuracy plot
            with Grid():
                max_y = compare_max_axes(
                    stats_df_t, stats_df_n, 'acc',
                    ptype='hist', binwidth=0.1)
                inputs = (
                    (stats_df_t, "Tumor Alignment Accuracy", COLORS.cerulean),
                    (stats_df_n, "Normal Alignment Accuracy", COLORS.green)
                )
                for (df, header, color) in inputs:
                    plt = hist_plot(
                        df, 'acc', header, xaxis='Accuracy [%]',
                        rounding=1, color=color, yaxis='Number of reads',
                        binwidth=0.1, max_y=max_y, min_x=80, max_x=100)
                    EZChart(plt, 'epi2melabs')
            p("""Yellow: mean length; Purple: median length.""")
            p("""The distribution is truncated to show the range between 80-100%.""")

    with report.add_section('Coverage', 'Coverage'):
        # Predispose to tabs
        tabs = Tabs()
        # Create plots
        with tabs.add_tab(args.sample_id):
            # depth vs genomic coordinate plot on the left and cumulative depth
            # plot on the right
            with Grid():
                max_y = compare_max_axes(
                    depth_df_t, depth_df_n, 'depth',
                    ptype='val')
                inputs = (
                    (depth_df_t, t_cov, "Tumor coverage along reference"),
                    (depth_df_n, n_cov, "Normal coverage along reference")
                )
                for (df, mean_cov, header) in inputs:
                    # Consider replacing mean_pos with total_mean_pos to have a
                    # single line, instead of one line per chromosome.
                    plt = line_plot(
                        df, 'mean_pos', 'depth', 'chrom', header, add_mean=mean_cov,
                        xaxis='Position along reference', yaxis='Sequencing depth',
                        max_y=max_y)
                    EZChart(plt, 'epi2melabs')
        p(
            'Depth of coverage computed by',
            a("mosdepth", href="https://github.com/brentp/mosdepth"),
            '. The dashed line shows the total mean coverage.'
            )


# Finally, main
def main(args):
    """Run entry point."""
    logger = get_named_logger("report_qc")

    # Instantiate the report.
    report = LabsReport(
        f"{args.sample_id} | Read alignment statistics", "wf-somatic-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    # Load data separately for debugging
    logger.info(f'Loading {args.reference_fai}')
    faidx = fasta_idx(args.reference_fai)

    logger.info(f'Loading {args.read_stats_tumor}')
    stats_df_t = process_fastcat(args.read_stats_tumor)
    logger.info(f'Loading {args.read_stats_normal}')
    stats_df_n = process_fastcat(args.read_stats_normal)

    logger.info(f'Loading {args.flagstat_tumor}')
    flags_df_t = load_flagstat(args.flagstat_tumor)
    logger.info(f'Loading {args.flagstat_normal}')
    flags_df_n = load_flagstat(args.flagstat_normal)

    logger.info(f'Loading {args.mosdepth_summary_tumor}')
    depth_su_t = load_mosdepth_summary(args.mosdepth_summary_tumor)
    logger.info(f'Loading {args.mosdepth_summary_normal}')
    depth_su_n = load_mosdepth_summary(args.mosdepth_summary_normal)

    logger.info(f'Loading {args.depth_tumor}')
    depth_df_t = load_mosdepth(args.depth_tumor, faidx, args.window_size)
    logger.info(f'Loading {args.depth_normal}')
    depth_df_n = load_mosdepth(args.depth_normal, faidx, args.window_size)

    # Populate report
    populate_report(
        report=report,
        args=args,
        stats_df_t=stats_df_t,
        stats_df_n=stats_df_n,
        flags_df_t=flags_df_t,
        flags_df_n=flags_df_n,
        depth_su_t=depth_su_t[0],
        depth_su_n=depth_su_n[0],
        depth_df_t=depth_df_t,
        depth_df_n=depth_df_n,
        logger=logger)

    # Save report
    report_fname = f"{args.name}-report.html"
    report.write(report_fname)
    logger.info(f"Written report to '{report_fname}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--name",
        help="report name",
    )
    parser.add_argument(
        "--sample_id",
        help="Sample identifier",
    )
    parser.add_argument(
        "--read_stats_tumor",
        help="`bamstats` file of per-read stats for the tumor sample",
    )
    parser.add_argument(
        "--read_stats_normal",
        help="`bamstats` file of per-read stats for the normal sample",
    )
    parser.add_argument(
        "--flagstat_tumor",
        help="`bamstats` flagstats for the tumor sample",
    )
    parser.add_argument(
        "--flagstat_normal",
        help="`bamstats` flagstats for the normal sample",
    )
    parser.add_argument(
        "--mosdepth_summary_tumor",
        required=False,
        help="`mosdepth` summary for the tumor sample",
    )
    parser.add_argument(
        "--mosdepth_summary_normal",
        required=False,
        help="`mosdepth` summary for the normal sample",
    )
    parser.add_argument(
        "--depth_tumor",
        required=False,
        help="`mosdepth` depth files for the tumor sample",
    )
    parser.add_argument(
        "--depth_normal",
        required=False,
        help="`mosdepth` depth files for the normal sample",
    )
    parser.add_argument(
        "--tumor_cov_threshold",
        required=True,
        type=float,
        help="Coverage threshold for the tumor sample",
    )
    parser.add_argument(
        "--normal_cov_threshold",
        required=True,
        type=float,
        help="Coverage threshold for the tumor sample",
    )
    parser.add_argument(
        "--reference_fai",
        required=True,
        help="Fai index for the reference genome",
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=50000,
        help="Fai index for the reference genome",
    )
    parser.add_argument(
        "--params",
        default=None,
        help="CSV file with workflow parameters",
    )
    parser.add_argument(
        "--versions",
        help="CSV file with software versions",
    )
    return parser
