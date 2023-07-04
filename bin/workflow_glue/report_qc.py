#!/usr/bin/env python
"""Plot QC metrics."""

from dominate.tags import a, p
from ezcharts.components.common import fasta_idx
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import load_bamstats_flagstat, load_stats
from ezcharts.components.mosdepth import load_mosdepth_regions, load_mosdepth_summary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
import pandas as pd

from .report_utils.utils import COLORS, compute_n50  # noqa: ABS101
from .report_utils.utils import compare_max_axes  # noqa: ABS101
from .report_utils.visualizations import hist_plot, line_plot  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


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
                        yaxis='Number of reads', extra_metric={'N50': n50},
                        binwidth=1000, max_y=max_y, rounding=0)
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
                        df, 'total_mean_pos', 'depth', 'chrom', header,
                        xaxis='Position along reference', yaxis='Sequencing depth',
                        max_y=max_y, add_mean=mean_cov)
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
    stats_df_t = load_stats(args.read_stats_tumor, format='bamstats')
    logger.info(f'Loading {args.read_stats_normal}')
    stats_df_n = load_stats(args.read_stats_normal, format='bamstats')

    logger.info(f'Loading {args.flagstat_tumor}')
    flags_df_t = load_bamstats_flagstat(args.flagstat_tumor)
    logger.info(f'Loading {args.flagstat_normal}')
    flags_df_n = load_bamstats_flagstat(args.flagstat_normal)

    logger.info(f'Loading {args.mosdepth_summary_tumor}')
    depth_su_t = load_mosdepth_summary(args.mosdepth_summary_tumor)
    logger.info(f'Loading {args.mosdepth_summary_normal}')
    depth_su_n = load_mosdepth_summary(args.mosdepth_summary_normal)

    logger.info(f'Loading {args.depth_tumor}')
    depth_df_t = load_mosdepth_regions(
        args.depth_tumor, faidx=faidx, winsize=args.window_size, min_size=10000000)
    logger.info(f'Loading {args.depth_normal}')
    depth_df_n = load_mosdepth_regions(
        args.depth_normal, faidx=faidx, winsize=args.window_size, min_size=10000000)

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
