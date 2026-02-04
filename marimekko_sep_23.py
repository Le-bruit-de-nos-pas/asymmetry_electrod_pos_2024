import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from plotly.subplots import make_subplots
import math
     

def stacked_bar_width_plot(df, value_cols, width_col, colors=None, **subplot):
    """A stacked column plot with variable bar width.
       :param df
       :param list value_cols: columns of `df`, already normalized (sum=1 for every line).
       :param str width_col: column of `df`, unbounded, used (i) as label (ii) to compute width.
       :param dict subplot: optional figure/row/col
    """
    categories = df.index.to_list()
    width = df[width_col]
    htmpl = width_col + ': %{customdata}, %{y}'
    x = np.cumsum([0] + list(width[:-1]))
    colors = [ "#ffcc66", "#0074b3", "#903149"] # or [None, ] * len(value_cols)


    figure = subplot.pop('figure', go.Figure())
    for colname, color in zip(value_cols, colors):
        figure.add_trace(go.Bar(name=colname, x=x, y=df[colname], text=round(100*df[colname],1), width=width, hovertemplate=htmpl, customdata=width, offset=0, marker_color=color), **subplot)
        figure.update_layout(autosize=False, width=700, height=400)
    return figure.update_xaxes(
        tickvals=x + np.array(width) / 2,
        ticktext=categories,
        range=[0, np.sum(width)], **subplot
    ) \
        .update_yaxes(tickformat=',.0%', range=[0, 1], row=subplot.get('row'), col=subplot.get('col')) \
        .update_layout(barmode='stack', bargap=0)
     


def mekko_plot(df, unit_name=None, colors=None, **subplot):
    """A mekko plot is a normalized stacked column plot with variable bar width.

       :param DataFrame df: already indexed by category, with only numeric columns:
                   there will be one stacked-bar per line (X), labeled after the index, with one bar per column.
       :param str unit_name: used to populate hover.
       :param list colors: color of each column. None for default palette (blue, red, ...).
       The rest (title..) is easily added afterwards.
    """
    # Normalize then defer to stacked_bar_width_plot plot.
    value_cols = df.columns
    w = pd.DataFrame({unit_name: df.sum(axis='columns')})
    w[value_cols] = df.div(w[unit_name], axis='index')

    return stacked_bar_width_plot(w, value_cols, width_col=unit_name, colors=colors, **subplot)

     

# OFFOFF -> ON-ON
df = pd.DataFrame(dict(Worst_Side_OFFOFF=['Equal','Left','Right'],
                       Equal=[7, 7, 18 ],
                       Left=[24 , 230 , 35 ],
                       Right=[9 , 10 , 201]
                       )).set_index('Worst_Side_OFFOFF')

#display(df)

mekko_plot(df, unit_name='some_label')
     

# OFFOFF -> ON-OFF
df = pd.DataFrame(dict(Worst_Side_OFFOFF=['Equal','Left','Right'],
                       Equal=[4, 11, 16 ],
                       Left=[16 , 208 , 67 ],
                       Right=[13 , 30 , 178]
                       )).set_index('Worst_Side_OFFOFF')

#display(df)

mekko_plot(df, unit_name='some_label')
     

# OFFOFF -> OFF-ON
df = pd.DataFrame(dict(Worst_Side_OFFOFF=['Equal','Left','Right'],
                       Equal=[3, 11, 17 ],
                       Left=[23 , 201 , 62 ],
                       Right=[8 , 21 , 191]
                       )).set_index('Worst_Side_OFFOFF')

#display(df)

mekko_plot(df, unit_name='some_label')
     
