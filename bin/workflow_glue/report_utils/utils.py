"""Global variables commonly used in other scripts."""
from dominate.tags import div, p
from ezcharts.plots import util
import numpy as np
import pandas as pd
from seaborn._statistics import Histogram


# Global variables
COLORS = util.Colors
# Create chromosome order for sorting
CHROMOSOMES = {f"chr{i}": int(i) for i in range(0, 23)}
CHROMOSOMES.update({'chrX': 23, 'chrY': 24})
# Fix chromosome naming to ensure all match chrN
CHROM_RENAME = {str(i): f"chr{i}" for i in range(1, 23)}
CHROM_RENAME.update({'X': "chrX", 'Y': "chrY", 'M': "chrM", 'MT': "chrMT"})
CHROM_RENAME.update({val: val for key, val in CHROM_RENAME.items()})
# Digit precision
PRECISION = 4


# Utility functions
def display_alert(*argv):
    """Display alert in report."""
    with div(cls="alert alert-warning"):
        p(argv)


def hist_max(variable_data, binwidth=None, bins='auto'):
    """Compute max value to set in a plot."""
    estimate_kws = dict(
        stat='count',
        bins=bins,
        binwidth=binwidth,
        binrange=None,
        discrete=None,
        cumulative=False,
    )
    estimator = Histogram(**estimate_kws)
    heights, edges = estimator(variable_data, weights=None)
    return max(heights)


def compare_max_axes(
        df1, df2, col, ptype='val',
        bins='auto', bins_2=None,
        binwidth=None, binwidth_2=None,
        buffer=1.1, precision=0):
    """Compute max value to set in a plot."""
    # If not specified, consider the two hists as having the same binning
    if not bins_2:
        bins_2 = bins
    if not binwidth_2:
        binwidth_2 = binwidth
    if ptype == 'hist':
        v1max = hist_max(df1[col].dropna(), bins=bins, binwidth=binwidth)
        v2max = hist_max(df2[col].dropna(), bins=bins_2, binwidth=binwidth_2)
    else:
        v1max = df1[col].max()
        v2max = df2[col].max()
    return np.ceil(max(v1max, v2max) * buffer * (10**precision))/(10**precision)


def compute_n50(data, x='start', y='count'):
    """Automatic detection of the data type for N50."""
    if isinstance(data, pd.DataFrame):
        n50_value = n50_hist(data, x=x, y=y)
    elif isinstance(data, np.ndarray):
        n50_value = n50_array(data)
    elif isinstance(data, list) or isinstance(data, tuple):
        n50_value = n50_array(np.array(data))
    else:
        raise ValueError(f"Unsupported data type {type(data)}")
    return n50_value


def n50_array(lengths):
    """Compute read N50.

    :param length: numpy vector of lengths
    """
    sorted_l = np.sort(lengths)[::-1]
    cumsum = np.cumsum(sorted_l)
    mid = cumsum[-1] / 2
    n50_idx = np.searchsorted(cumsum, mid)
    return sorted_l[n50_idx]


def n50_hist(length_hist, x='start', y='count'):
    """Compute read N50 from histogram data."""
    cumsum = np.cumsum(length_hist[y].to_numpy() * length_hist[x].to_numpy())
    mid = cumsum[-1] / 2
    n50_idx = np.searchsorted(cumsum, mid)
    return length_hist.iloc[n50_idx].start


def hist_binning(values):
    """Define histogram bin number and size."""
    n_bins = int(np.ceil(2*np.cbrt(values.shape[0])))
    max_x = int(values.max())
    binwidth = int(values.max() / n_bins)
    binrange = [0, 0]
    while binrange[-1] < max_x:
        binrange[-1] += binwidth
    return binrange, binwidth
