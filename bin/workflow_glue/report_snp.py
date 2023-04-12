#!/usr/bin/env python
"""Create workflow report."""

from aplanat import annot, hist, report
from aplanat.components import bcfstats, depthcoverage
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import gridplot
from bokeh.models import FactorRange
from bokeh.plotting import figure
import numpy as np
import pandas as pd
import pysam

from .util import wf_parser  # noqa: ABS101


def vcf_parse(args):
    """Parse the vcf file."""
    vcf_df = pysam.VariantFile(args.vcf, 'r')

    if len(vcf_df.header.samples) > 1:
        raise Exception("VCF files with more than one sample are not supported.")
    samplename, = vcf_df.header.samples
    flt = []
    vaf = []
    naf = []
    for rec in vcf_df:
        flt.append(rec.filter.keys()[0])
        vaf.append(rec.samples[samplename]['AF'])
        naf.append(rec.samples[samplename]['NAF'])
    return flt, vaf, naf


def vaf_plot(vaf, naf, section):
    """Plot the allele frequencies for cancer sample."""
    # Plot the histogram of the allele frequencies in the tumor and normal
    plot_af = hist.histogram(
        [vaf],
        bins=20,
        xlim=(0, 1.0),
        title='Tumor allele frequency'
        )
    plot_nf = hist.histogram(
        [naf],
        bins=20,
        xlim=(0, 1.0),
        title='Control allele frequency'
        )

    section.markdown("""
### Allele frequencies
This section shows a histogram of the variant allele frequencies in the tumor
vs normal samples.
""")
    section.plot(gridplot(
        [plot_af, plot_nf],
        ncols=2,
        width=500)
    )


def filt_stats(filters, vaf, naf, af_filters, section):
    """Plot the filtering stats."""
    df = pd.DataFrame({'Filter': filters, 'Tumor_AF': vaf, 'Normal_AF': naf})
    # Pivot the values
    highest = sorted(map(float, af_filters))[-1]
    lowest = sorted(map(float, af_filters))[0]
    passing = (df['Filter'] == 'PASS') & \
        (df['Tumor_AF'] >= highest) & \
        (df['Normal_AF'] < lowest)
    n_passing = df.loc[passing].shape[0]
    stats = {
        'N sites': [df.count()['Filter']],
        'PASS': [df.query('Filter == "PASS"').shape[0]],
        'not PASS': [df.query('Filter != "PASS"').shape[0]]
        }
    # Add filtered variant counts for the cancer
    for c_name in ('Tumor_AF', 'Normal_AF'):
        for filt in af_filters:
            colname = "{} >= {}".format(c_name.replace('_', ' '), filt)
            stats[colname] = [df.query(f'{c_name} >= {filt}').shape[0]]
    summary_table = pd.DataFrame(stats)
    # Plot the histogram of the allele frequencies in the tumor and normal
    section.markdown(f"""
### Filtering
This section shows the count for quality filtering in the dataset.
There are {n_passing} variants labelled as 'PASS', with tumor VAF >= {highest} and
normal VAF < {lowest}.
""")
    section.table(summary_table, index=False)


def plot_spectra(args, section):
    """Plot the mutation spectra."""
    # Load and combine the spectra
    cnt = pd.read_csv(args.mut_spectra)
    sample = cnt.columns[1]

    # Mutation order
    midpoint = int(np.floor(len(cnt['Type'].iloc[0])/2))
    cnt["Change"] = cnt['Type'].str.slice(start=midpoint-1, stop=midpoint+2)
    cnt["Flanks"] = cnt['Type'].str.slice(start=0, stop=midpoint-2) + \
        '-' + cnt['Type'].str.slice(start=midpoint+3)
    change_col = cnt.pop("Change")
    flank_col = cnt.pop("Flanks")
    cnt.insert(0, "Flanks", flank_col)
    cnt.insert(0, "Change", change_col)
    cnt = cnt.sort_values(by=['Change', 'Type'])

    # Plot profile
    factors = list(zip(cnt['Change'], cnt['Flanks']))
    x = cnt[sample].tolist()
    pr_plot2 = figure(
        x_range=FactorRange(*factors),
        height=250,
        toolbar_location=None,
        tools=""
        )

    pr_plot2.vbar(
        x=factors,
        top=x,
        width=0.9,
        color=Colors.cerulean,
        )
    pr_plot2.y_range.start = 0
    pr_plot2.xaxis.major_label_orientation = np.pi/2
    pr_plot2.xgrid.grid_line_color = None

    # Generate the plots in the report
    section.markdown("""
### Mutational profile
This section displays the mutational profile as counts of changes in their context.
""")
    section.plot(gridplot(
        [pr_plot2],
        ncols=1,
        width=1000)
    )


def plot_qc_stats(read_stats):
    """Plot read QC stats."""
    # read length
    tp = pd.read_csv(
        read_stats, sep='\t', chunksize=1000, iterator=True)
    seq_summary = pd.concat(tp, ignore_index=True)
    total_bases = seq_summary['read_length'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = seq_summary['read_length'].median()
    datas = [seq_summary['read_length']]
    length_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, None))
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))

    # read quality
    datas = [seq_summary['acc']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[Colors.cerulean], bins=100,
        title="Read quality (wrt reference sequence)",
        x_axis_label="Read Quality",
        y_axis_label="Number of reads",
        xlim=(85, 100))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))
    return gridplot([[length_hist, q_hist]])


def depth_plot_section(read_depth, section):
    """Add read depth plots to a given section."""
    rd_plot = depthcoverage.cumulative_depth_from_dist(read_depth)
    section.markdown("""
### Genome coverage
This section displays basic metrics relating to genome coverage.
""")
    section.plot(gridplot(
        [rd_plot],
        ncols=2)
    )


def main(args):
    """Run the entry point."""
    report_doc = report.HTMLReport(
        "Somatic variant calling Summary Report",
        ("Results generated through the wf-somatic-variation nextflow "
            "workflow provided by Oxford Nanopore Technologies"))

    report_doc = WFReport(
        "Workflow Somatic SNP", "wf-somatic-variation",
        revision=args.revision, commit=args.commit)
    if args.read_stats:
        section = report_doc.add_section()
        section.markdown("""
    ### Read Quality control
    This section displays basic QC metrics indicating read data quality.
    """)
        section.plot(plot_qc_stats(args.read_stats))
    # canned VCF stats report component
    if args.vcf_stats:
        section = report_doc.add_section()
        bcfstats.full_report(args.vcf_stats, report=section)

    if args.read_depth:
        section = report_doc.add_section()
        depth_plot_section(args.read_depth, section)

    # Plot tumor/normal AF if provided
    if args.vcf:
        section = report_doc.add_section()
        filters, var_af, nvar_af = vcf_parse(args)
        af_filters = [i for i in args.af_filter_intervals.split(',') if float(i) < 1]
        if not af_filters:
            raise Exception("Invalid range of allele frequencies.")
        vaf_plot(var_af, nvar_af, section)
        filt_stats(filters, var_af, nvar_af, af_filters, section)

    # Plot mutation spectra if provided
    if args.mut_spectra:
        section = report_doc.add_section()
        plot_spectra(args, section)

    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    # write report
    report_doc.write(args.report)


def argparser():
    """Create argument parser."""
    parser = wf_parser("report_snp")

    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--read_stats",
        help="read statistics output from bamstats")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--vcf",
        help="input vcf file")
    parser.add_argument(
        "--vcf_stats",
        help="final vcf stats file")
    parser.add_argument(
        "--mut_spectra",
        help="final vcf stats file")
    parser.add_argument(
        "--read_depth",
        help="read coverage output from mosdepth")
    parser.add_argument(
        "--af_filter_intervals", default="0.05,0.10,0.20",
        help="Allele frequency intervals to use to characterize the variants")

    return parser
