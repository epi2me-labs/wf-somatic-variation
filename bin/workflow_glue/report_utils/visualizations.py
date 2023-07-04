"""Commonly used code in visualizations."""
from ezcharts import barplot, histplot
from ezcharts import lineplot, scatterplot
import numpy as np
from .utils import COLORS  # noqa: ABS101


# Scatter plot
def scatter_plot(
        df, x, y, hue, title, add_mean=None, xaxis='', yaxis='',
        min_x=None, max_x=None, min_y=None, max_y=None, tooltip_label=None):
    """Make a scatterplot."""
    plt = scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue
    )
    # Change axes names
    plt.xAxis.name = xaxis
    plt.yAxis.name = yaxis
    # Set y limit
    if min_y is not None:
        plt.yAxis.min = min_y
    if max_y is not None:
        plt.yAxis.max = max_y
    if min_x is not None:
        plt.xAxis.min = min_x
    if max_x is not None:
        plt.xAxis.max = max_x
    # Remove markers
    for s in plt.series:
        s.showSymbol = True
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


# Histogram
def hist_plot(
        df, col, title, xaxis='', yaxis='', rounding=None,
        extra_metric=None, color=None, binwidth=None, binrange=None, bins='auto',
        max_y=None, max_x=None, min_x=None, min_y=None, no_stats=False):
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

    if isinstance(extra_metric, dict):
        if len(extra_metric) > 1:
            raise ValueError('Too many input extra metrics')
        key = list(extra_metric.keys())[0]
        val = list(extra_metric.values())[0]
        plt.title = dict(
            text=title,
            subtext=(
                f"Mean: {meanv}. "
                f"Median: {medianv}. "
                f"{key}: {val}."
            ),
        )
    else:

        plt.title = dict(
            text=title,
            subtext=(
                f"Mean: {meanv}. "
                f"Median: {medianv}. "
            )
        )

    # Add mean and median values (Thanks Julian!)
    if not no_stats:
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
    if isinstance(extra_metric, dict):
        for key, val in extra_metric.items():
            plt.add_series(
                dict(
                    type="line",
                    name=f"{key}",
                    data=[dict(value=[val, 0]), dict(value=[val, max_y])],
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


def plot_spectra(spectra, sample, cmap=None):
    """Plot the mutation spectra."""
    # Plot change spectrum
    df = spectra[['Change', sample]]
    # Count by change type
    tots = df.groupby('Change').sum().reset_index()
    tots.columns = ['Change', sample]
    # Create bar plot
    plot = barplot(data=tots, x='Change', y=sample)
    if cmap:
        plot.color = list(cmap.values())
    plot.xAxis.name = 'Change'
    plot.yAxis.name = "Counts"
    for s in plot.series:
        s.colorBy = 'data'
    plot.tooltip = dict({'trigger': 'item'})
    return plot


def plot_profile(df, sample, cmap=None):
    """Plot the mutation profile."""
    # Ensure sorting
    df = df.sort_values(["Change", "Flanks"])

    # Create base barplot
    plt = barplot(
        data=df.assign(order=range(df.shape[0])),
        x="order", y=sample,
        hue='Change', dodge=False)

    # Define appropriate palette
    if cmap:
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
            name="",
            data=[dict(value=tmp_sum[i, :]) for i in range(0, tmp_sum.shape[0])],
            symbolSize=0,
            tooltip=dict(formatter='{c}')
        ))
    plt.tooltip = dict({'trigger': 'item'})
    return plt
