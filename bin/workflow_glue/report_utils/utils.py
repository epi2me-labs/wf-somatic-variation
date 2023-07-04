"""Global variables commonly used in other scripts."""
from ezcharts.plots import util
import numpy as np
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
        bins='auto', binwidth=None,
        buffer=1.1, precision=0):
    """Compute max value to set in a plot."""
    if ptype == 'hist':
        v1max = hist_max(df1[col].dropna(), bins=bins, binwidth=binwidth)
        v2max = hist_max(df2[col].dropna(), bins=bins, binwidth=binwidth)
    else:
        v1max = df1[col].max()
        v2max = df2[col].max()
    return np.ceil(max(v1max, v2max) * buffer * (10**precision))/(10**precision)


def compute_n50(lengths):
    """Compute read N50."""
    # Sort the read lengths
    sorted_l = np.sort(lengths)[::-1]
    # Generate cumsum
    cumsum = np.cumsum(sorted_l)
    # Get lowest cumulative value >= (total_length/2)
    n50 = sorted_l[cumsum >= cumsum[-1]/2][0]
    return n50
