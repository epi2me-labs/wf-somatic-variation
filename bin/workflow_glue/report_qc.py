#!/usr/bin/env python
"""Plot QC metrics."""

from dominate.tags import a, p
from ezcharts.components import fastcat
from ezcharts.components.common import fasta_idx
from ezcharts.components.ezchart import EZChart
from ezcharts.components.mosdepth import load_mosdepth_regions, load_mosdepth_summary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
import pandas as pd

from .report_utils.utils import compute_n50  # noqa: ABS101
from .report_utils.utils import compare_max_axes  # noqa: ABS101
from .report_utils.utils import display_alert  # noqa: ABS101
from .report_utils import visualizations  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


# Reporting function
def populate_report(report, args, **kwargs):
    """Populate the report with the different sections."""
    # Extract data used
    tumor_hist_df = kwargs["tumor_hist_df"]
    normal_hist_df = kwargs["normal_hist_df"]
    tumor_depth_su = kwargs["tumor_depth_su"]
    normal_depth_su = kwargs["normal_depth_su"]
    tumor_depth_df = kwargs["tumor_depth_df"]
    normal_depth_df = kwargs["normal_depth_df"]
    logger = kwargs["logger"]

    # Define depth statistics and variables for tumor.
    logger.info("Compute average cvg...")
    tumor_cov = tumor_depth_su.loc[tumor_depth_su["chrom"] == "total", "mean"].values[0]
    tumor_pass = "PASS" if tumor_cov > args.tumor_cov_threshold else "FAIL"
    tumor_cov_str = "{}x".format(tumor_cov)
    tumor_min_chrom_cov = tumor_depth_su["mean"].min()
    tumor_max_chrom_cov = tumor_depth_su["mean"].max()

    # Define depth statistics and variables for normal, if available.
    if normal_depth_su.empty:
        normal_min_chrom_cov = normal_max_chrom_cov = "NA"
        normal_pass = normal_cov_str = normal_cov = "NA"
    else:
        normal_cov = normal_depth_su.loc[
            normal_depth_su["chrom"] == "total", "mean"
        ].values[0]
        normal_pass = "PASS" if normal_cov > args.normal_cov_threshold else "FAIL"
        normal_cov_str = "{}x".format(normal_cov)
        normal_min_chrom_cov = normal_depth_su["mean"].min()
        normal_max_chrom_cov = normal_depth_su["mean"].max()

    # Define read stats for tumor.
    logger.info("Compute N50...")
    tumor_read_number_str = str(tumor_hist_df["count"].sum())
    tumor_n50 = compute_n50(tumor_hist_df)
    tumor_n50_str = str(tumor_n50) + " bp"
    tumor_median_read_length = int(fastcat.histogram_median(tumor_hist_df))

    # Define read stats for normal, if available.
    if normal_hist_df.empty:
        normal_median_read_length = normal_read_number_str = "NA"
        normal_n50_str = normal_n50 = "NA"
    else:
        normal_read_number_str = str(normal_hist_df["count"].sum())
        normal_n50 = compute_n50(normal_hist_df)
        normal_n50_str = str(normal_n50) + " bp"
        normal_median_read_length = int(fastcat.histogram_median(normal_hist_df))

    # Start structuring the report
    logger.info("Start reporting.")
    logger.info("Saving summary section...")
    with report.add_section("At a glance", "Description"):
        if normal_read_number_str == "NA":
            display_alert(
                "Workflow run in tumor-only mode; statistics for only the",
                " tumor sample will be reported.",
            )
            sample_string = "a tumor sample"
        else:
            sample_string = "paired tumor/normal samples"
        p(
            f"""
            This report contains visualisations of alignment statistics for
            {sample_string} that can help in understanding the results from
            wf-somatic-variation. Each section contains different plots or tables.
            You can quickly jump to an individual section with the links in the
            header bar.
            """
        )

        # Create tabs for each sample
        tabs = Tabs()
        with tabs.add_tab(args.sample_id):
            Stats(
                columns=2,
                items=[
                    (tumor_read_number_str, "Tumor Total Reads"),
                    (normal_read_number_str, "Normal Total Reads"),
                    (tumor_n50_str, "Tumor Read N50"),
                    (normal_n50_str, "Normal Read N50"),
                    (tumor_cov_str, "Tumor mean coverage"),
                    (normal_cov_str, "Normal mean coverage"),
                ],
            )

    # Add coverage thresholds passing/failing
    logger.info("Saving filtering section...")
    with report.add_section("Base statistics", "Stats"):
        normal_threshold_string = ""
        if args.normal_cov_threshold:
            normal_threshold_string = (
                f" and {args.normal_cov_threshold}x for the normal sequences"
            )
        p(
            f"""
            The aligned bam files have been tested for coverage thresholds of
            {args.tumor_cov_threshold}x for the tumor{normal_threshold_string}.
            """
        )

        # Create tabs for each sample
        df = pd.DataFrame(
            {
                "Sample": [args.sample_id, args.sample_id],
                "Type": ["Tumor", "Normal"],
                "Total Reads": [tumor_read_number_str, normal_read_number_str],
                "Median Read Length": [
                    tumor_median_read_length,
                    normal_median_read_length,
                ],
                "Read N50": [tumor_n50, normal_n50],
                "Min chrom. coverage": [tumor_min_chrom_cov, normal_min_chrom_cov],
                "Mean chrom. coverage": [tumor_cov, normal_cov],
                "Max chrom. coverage": [tumor_max_chrom_cov, normal_max_chrom_cov],
                "Pass threshold": [tumor_pass, normal_pass],
            }
        )
        DataTable.from_pandas(df, use_index=False)

    # Add comparison of read lengths
    logger.info("Saving read distribution...")
    with report.add_section("Read statistics", "Read Stats"):
        # Show only one plot if no normal is provided.
        if normal_hist_df.empty:
            fastcat.SeqSummary(
                args.hists_tumor,
                flagstat=args.flagstat_tumor,
                sample_names=f"{args.sample_id} (tumor)"
            )
        else:
            fastcat.SeqCompare(
                (args.hists_tumor, args.hists_normal),
                flagstat=(args.flagstat_tumor, args.flagstat_normal),
                sample_names=(f"{args.sample_id} (tumor)", f"{args.sample_id} (normal)")
            )

    with report.add_section("Coverage", "Coverage"):
        # Predispose to tabs
        tabs = Tabs()
        # Create plots
        with tabs.add_tab(args.sample_id):
            # depth vs genomic coordinate plot on the left and cumulative depth
            # plot on the right
            if normal_depth_df.empty:
                max_y = tumor_depth_df["depth"].max()
                # Consider replacing mean_pos with total_mean_pos to have a
                # single line, instead of one line per chromosome.
                plt = visualizations.line_plot(
                    tumor_depth_df,
                    "total_mean_pos",
                    "depth",
                    "chrom",
                    "Tumor coverage along reference",
                    xaxis="Position along reference",
                    yaxis="Sequencing depth",
                    max_y=max_y,
                    add_mean=tumor_cov,
                )
                EZChart(plt, "epi2melabs")
            else:
                with Grid():
                    max_y = compare_max_axes(
                        tumor_depth_df, normal_depth_df, "depth", ptype="val"
                    )
                    inputs = (
                        (tumor_depth_df, tumor_cov, "Tumor coverage along reference"),
                        (
                            normal_depth_df,
                            normal_cov,
                            "Normal coverage along reference",
                        ),
                    )
                    for df, mean_cov, header in inputs:
                        # Consider replacing mean_pos with total_mean_pos to have a
                        # single line, instead of one line per chromosome.
                        plt = visualizations.line_plot(
                            df,
                            "total_mean_pos",
                            "depth",
                            "chrom",
                            header,
                            xaxis="Position along reference",
                            yaxis="Sequencing depth",
                            max_y=max_y,
                            add_mean=mean_cov,
                        )
                        EZChart(plt, "epi2melabs")
        p(
            "Depth of coverage computed by",
            a("mosdepth", href="https://github.com/brentp/mosdepth"),
            ". The dashed line shows the total mean coverage.",
        )


def read_stats(
    logger,
    hists_dir=None,
    depth=None,
    mosdepth_summary=None,
    faidx=None,
    winsize=50000,
):
    """Load all stats files for a sample."""
    if hists_dir:
        logger.info(f"Loading length hists from: {hists_dir}")
        hist_df = fastcat.sum_hists((
            fastcat.load_histogram(hists_dir, "length"),
            fastcat.load_histogram(hists_dir, "length.unmap")
        ))
    else:
        hist_df = pd.DataFrame()

    if mosdepth_summary:
        logger.info(f"Loading {mosdepth_summary}")
        depth_su = load_mosdepth_summary(mosdepth_summary)
    else:
        depth_su = pd.DataFrame()

    if depth:
        logger.info(f"Loading {depth}")
        depth_df = load_mosdepth_regions(
            depth, faidx=faidx, winsize=winsize, min_size=10000000
        )
        if not depth_df.empty:
            depth_df = depth_df[["chrom", "total_mean_pos", "depth"]]
    else:
        depth_df = pd.DataFrame()

    return hist_df, depth_su, depth_df


# Finally, main
def main(args):
    """Run entry point."""
    logger = get_named_logger("report_qc")

    # Instantiate the report.
    report = LabsReport(
        f"{args.sample_id} | Read alignment statistics",
        "wf-somatic-variation",
        args.params,
        args.versions,
        args.workflow_version,
        head_resources=[*LAB_head_resources],
    )

    # Load data separately for debugging
    logger.info(f"Loading {args.reference_fai}")
    faidx = fasta_idx(args.reference_fai)

    logger.info("Loading tumor stats")
    tumor_hist_df, tumor_depth_su, tumor_depth_df = read_stats(
        logger,
        hists_dir=args.hists_tumor,
        depth=args.depth_tumor,
        mosdepth_summary=args.mosdepth_summary_tumor,
        faidx=faidx,
        winsize=args.window_size,
    )

    if args.hists_normal:
        logger.info("Loading normal stats")
        normal_hist_df, normal_depth_su, normal_depth_df = read_stats(
            logger,
            hists_dir=args.hists_normal,
            depth=args.depth_normal,
            mosdepth_summary=args.mosdepth_summary_normal,
            faidx=faidx,
            winsize=args.window_size,
        )
    else:
        normal_hist_df = pd.DataFrame()
        normal_depth_su = [pd.DataFrame()]
        normal_depth_df = pd.DataFrame()

    # Populate report
    populate_report(
        report=report,
        args=args,
        tumor_hist_df=tumor_hist_df,
        normal_hist_df=normal_hist_df,
        tumor_depth_su=tumor_depth_su[0],
        normal_depth_su=normal_depth_su[0],
        tumor_depth_df=tumor_depth_df,
        normal_depth_df=normal_depth_df,
        logger=logger,
    )

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
        "--hists_tumor",
        help="`bamstats` histograms folder for the tumor sample",
    )
    parser.add_argument(
        "--hists_normal",
        help="`bamstats` histograms folder for the normal sample",
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
        required=False,
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
    parser.add_argument(
        "--workflow_version",
        required=True,
        help="Workflow version",
    )
    return parser
