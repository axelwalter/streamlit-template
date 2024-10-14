import streamlit as st

import pandas as pd
import numpy as np

from .common.common import *

import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff


from sklearn.preprocessing import StandardScaler


INFO = "💡 Both workflows provide input for the [advanced metabolomics statistis tool](https://metabolomics-statistics.streamlit.app/). Visit [the GitHub repository](https://github.com/axelwalter/streamlit-metabolomics-statistics)."


@st.cache_data
def scale_df(df):
    scaled = pd.DataFrame(
        StandardScaler().fit_transform(df)).set_index(df.index)
    scaled.columns = df.columns
    return scaled


@st.cache_resource
def dendrogram(df):
    fig = ff.create_dendrogram(df, labels=list(df.index))
    fig.update_xaxes(side="top")
    return fig


@st.cache_resource
def heatmap(df):
    fig = px.imshow(
        df,
        y=list(df.index),
        x=list(df.columns),
        # text_auto=True,
        aspect="auto",
        color_continuous_scale="RdBu_r",
        range_color=[-3, 3],
    )

    fig.update_layout(
        autosize=False, width=700, height=1200, xaxis_title="", yaxis_title=""
    )

    # fig.update_yaxes(visible=False)
    fig.update_xaxes(tickangle=35)
    return fig


@st.cache_resource
def fold_change_plot(change):
    fig = px.bar(change)
    fig.update_layout(
        showlegend=False,
        title="log 2 fold change",
        yaxis=dict(title="log 2 fold change"),
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


@st.cache_resource
def mean_intensity_plot(samples, features, mean, std):
    fig = go.Figure()
    for feature in features:
        if not std.empty:
            fig.add_trace(
                go.Bar(
                    x=samples,
                    y=mean.loc[feature],
                    name=feature,
                    error_y=dict(
                        type="data", array=std.loc[feature], visible=True
                    ),
                )
            )
        else:
            fig.add_trace(
                go.Bar(x=samples, y=mean.loc[feature], name=feature))
    fig.update_layout(
        title="Metabolite Intensities",
        yaxis=dict(title="mean intensity"),
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


@st.cache_data
def get_mean_change_std(df, a, b):
    mean = pd.concat(
        [
            df[[col for col in df.columns if a in col]].mean(axis=1),
            df[[col for col in df.columns if b in col]].mean(axis=1),
        ],
        axis=1,
    ).rename(columns={0: a, 1: b})
    change = pd.DataFrame(np.log2(mean[a] / mean[b])).rename(
        columns={0: "log2 fold change"}
    )
    std = pd.concat(
        [
            df[[col for col in df.columns if a in col]].std(axis=1),
            df[[col for col in df.columns if b in col]].std(axis=1),
        ],
        axis=1,
    ).rename(columns={0: a, 1: b})
    return mean, change, std
