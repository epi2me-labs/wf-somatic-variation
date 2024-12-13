"""Commonly used code in visualizations."""
from bokeh.models import HoverTool, Title
from ezcharts import barplot, histplot
from ezcharts import lineplot, scatterplot


# Scatter plot
def scatter_plot(
        df, x, y, hue, title, add_mean=None, xaxis='', yaxis='', color=None,
        min_x=None, max_x=None, min_y=None, max_y=None, tooltip_label=None):
    """Make a scatterplot."""
    plt = scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        color=color
    )
    # Change axes names
    plt._fig.xaxis.axis_label = xaxis
    plt._fig.yaxis.axis_label = yaxis
    # Customize X-axis
    if max_x is not None:
        plt._fig.x_range.end = max_x
    if min_x is not None:
        plt._fig.x_range.start = min_x
    if max_y is not None:
        plt._fig.y_range.end = max_y
    if min_y is not None:
        plt._fig.y_range.start = min_y
    if add_mean:
        plt._fig.add_layout(
            Title(
                text=f"Mean: {round(add_mean, 2)}",
                text_font_size="0.8em"),
            'above'
        )
    plt._fig.add_layout(
        Title(text=title, text_font_size="1.5em"),
        'above'
    )
    return plt


# Histogram
def hist_plot(
        df, col, title, xaxis='', yaxis='', rounding=None,
        color=None, binwidth=None, binrange=None, bins='auto',
        max_y=None, max_x=None, min_x=None, min_y=None, stats=True):
    """Make a histogram of given parameter."""
    histogram_data = df[col].values

    plt = histplot(
        data=histogram_data,
        bins=bins,
        binwidth=binwidth,
        binrange=binrange,
        color=color)

    # Change axes names
    plt._fig.xaxis.axis_label = xaxis
    plt._fig.yaxis.axis_label = yaxis
    # Change title and add mean as subtitle and horiz. line, if provided
    if isinstance(rounding, int):
        meanv = df[col].mean().round(rounding)
        medianv = df[col].median().round(rounding)
    else:
        meanv = df[col].mean()
        medianv = df[col].median()

    if stats:
        plt._fig.add_layout(
            Title(
                text=f"Mean: {round(meanv, 2)}. Median: {round(medianv, 2)}",
                text_font_size="0.8em"),
            'above'
        )
    plt._fig.add_layout(
        Title(text=title, text_font_size="1.5em"),
        'above'
    )

    # Customize X-axis
    if max_x is not None:
        plt._fig.x_range.end = max_x
    if min_x is not None:
        plt._fig.x_range.start = min_x
    if max_y is not None:
        plt._fig.y_range.end = max_y
    if min_y is not None:
        plt._fig.y_range.start = min_y

    # Tooltips
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [(yaxis, "@top")]
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
        hue=hue,
        marker=False
    )
    # Change axes names
    plt._fig.xaxis.axis_label = xaxis
    plt._fig.yaxis.axis_label = yaxis
    # Set y limit
    if max_y is not None:
        plt._fig.y_range.end = max_y
    # Change title and add mean as subtitle and horiz. line, if provided
    if add_mean is not None:
        plt._fig.add_layout(
            Title(
                text=f"Mean: {round(add_mean, 2)}.",
                text_font_size="0.8em"),
            'above'
        )
    plt._fig.add_layout(
        Title(text=title, text_font_size="1.5em"),
        'above'
    )
    # Tooltips
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [(yaxis, "@y")]
    return plt


def plot_spectra(spectra, sample, cmap=None):
    """Plot the mutation spectra."""
    # Plot change spectrum
    df = spectra[['Change', sample]]
    # Count by change type
    tots = df.groupby('Change').sum().reset_index()
    tots.columns = ['Change', sample]

    # Define appropriate palette
    palette = None
    if cmap:
        palette = [cmap.get(i) for i in sorted(df.get("Change").unique())]

    # Create bar plot
    plot = barplot(
        data=tots,
        x='Change',
        y=sample,
        palette=palette
    )
    # Tooltips
    hover = plot._fig.select(dict(type=HoverTool))
    hover.tooltips = [("Number of sites", "@top")]
    return plot


def plot_profile(df, sample, cmap=None):
    """Plot the mutation profile."""
    # Ensure sorting
    df = df.sort_values(["Change", "Flanks"])

    # Define appropriate palette
    palette = None
    if cmap:
        palette = [cmap.get(i) for i in sorted(df.get("Change").unique())]

    # Create base barplot
    plt = barplot(
        data=df.assign(order=range(df.shape[0])),
        x="Flanks",
        y=sample,
        hue="Change",
        dodge=True,
        palette=palette,
        nested_x=True
    )

    # Sort out axis
    plt._fig.xaxis.major_label_orientation = 1.2
    # Tooltips
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [('Number of sites', "@y")]

    return plt
