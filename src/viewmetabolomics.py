import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from .common.common import show_fig

"""
Functions to generate individual plots.
"""

COLOR_SCALE = [
    (0.00, "rgba(233, 233, 233, 1.0)"),
    (0.01, "rgba(243, 236, 166, 1.0)"),
    (0.1, "rgba(255, 168, 0, 1.0)"),
    (0.2, "rgba(191, 0, 191, 1.0)"),
    (0.4, "rgba(68, 0, 206, 1.0)"),
    (1.0, "rgba(33, 0, 101, 1.0)"),
]


@st.cache_resource
def plot_peak_map_2D(df, cutoff):
    fig = go.Figure()

    ints = np.concatenate([df.loc[index, "intarray"] for index in df.index])
    int_filter = ints > cutoff  # show only ints over threshold
    ints = ints[int_filter]
    mzs = np.concatenate([df.loc[index, "mzarray"] for index in df.index])[int_filter]
    rts = np.concatenate(
        [
            np.full(len(df.loc[index, "mzarray"]), df.loc[index, "RT"])
            for index in df.index
        ]
    )[int_filter]

    sort = np.argsort(ints)
    ints = ints[sort]
    mzs = mzs[sort]
    rts = rts[sort]

    fig.add_trace(
        go.Scattergl(
            name="peaks",
            x=rts,
            y=mzs,
            mode="markers",
            marker_color=ints,
            marker_symbol="square",
        )
    )
    fig.update_layout(
        title="2D peak map",
        xaxis_title="retention time",
        yaxis_title="m/z",
        plot_bgcolor="rgb(255,255,255)",
        height=800,
        width=1000,
    )
    fig.layout.template = "plotly_white"

    fig.update_traces(
        marker_colorscale=COLOR_SCALE,
        hovertext=ints.round(),
        selector=dict(type="scattergl"),
    )

    return fig


def plot_feature_map(df):
    fig = go.Figure()
    # fig.add_traces([go.Scattergl(x=[df.loc[i, "RTstart"], df.loc[i, "RTend"]], y=[df.loc[i, "mz"], df.loc[i, "mz"]], mode="lines", line=dict(color="rgba(41,55,155, 0.5)")) for i in df.index])
    fig.add_trace(
        go.Scattergl(
            name="feature",
            x=df["RT"],
            y=df["mz"],
            mode="markers",
            marker_color=df["intensity"],
            marker_symbol="square",
            marker_size=12,
            customdata=np.stack(
                (
                    df["mz"].round(5),
                    df["intensity"].astype(int),
                    df["RTstart"].astype(int),
                    df["RTend"].astype(int),
                    df["RTend"].astype(int) - df["RTstart"].astype(int),
                    df["fwhm"].round(),
                    df["charge"],
                    df["quality"],
                    df["adduct"],
                ),
                axis=-1,
            ),
            hovertemplate="<b>mz: %{customdata[0]}<br>intensity: %{customdata[1]}<br>RTstart: %{customdata[2]}<br>RTend: %{customdata[3]}<br>RTrange: %{customdata[4]}<br>FWHM: %{customdata[5]}<br>charge: %{customdata[6]}<br>quality: %{customdata[7]}<br>adduct: %{customdata[8]}<br>",
        )
    )

    fig.update_layout(
        title="feature map",
        xaxis_title="retention time",
        yaxis_title="m/z",
        plot_bgcolor="rgb(255,255,255)",
        height=800,
        width=1000,
    )
    fig.layout.template = "plotly_white"
    fig.update_traces(
        showlegend=False,
        marker_colorscale=COLOR_SCALE,
        hovertext="int: "
        + df["intensity"].astype(int).astype(str)
        + " q: "
        + df["quality"].astype(str),
        selector=dict(type="scattergl"),
    )
    return fig


def plot_feature_chromatogram(df):
    fig = px.line(
        x=df.loc[df.index[0], "chrom_rts"],
        y=df.loc[df.index[0], "chrom_intensities"],
    )
    fig.update_layout(
        title=f"monoisotopic peak chromatogram",
        xaxis_title="retention time",
        yaxis_title="intensity (counts per second)",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


def plot_feature_map_alignment(df):
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            name="before aligment",
            x=df["original_rt"],
            y=df["mz"],
            mode="markers",
            marker_symbol="circle",
            marker_size=12,
        )
    )
    fig.add_trace(
        go.Scattergl(
            name="after aligment",
            x=df["RT"],
            y=df["mz"],
            mode="markers",
            marker_symbol="circle",
            marker_size=8,
        )
    )
    fig.update_layout(
        title=f"feature map alignment",
        xaxis_title="retention time",
        yaxis_title="m/z",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


def plot_consensus_map(df, sample):
    fig = go.Figure()

    # define custom data for hovering
    meta_values = [
        df["metabolite"],
        df["mz"].round(5),
        df["RT"].round(),
        df[sample],
        df["charge"],
        df["quality"],
    ]

    hovertemplate = """
<b>name: %{customdata[0]}<br>
mz: %{customdata[1]}<br>
RT: %{customdata[2]}<br>
intensity: %{customdata[3]}<br>
charge: %{customdata[4]}<br>
quality: %{customdata[5]}<br>
"""

    if "adduct" in df.columns:
        meta_values.append(df["adduct"])
        hovertemplate += "adduct: %{customdata[5]}<br>"

    for s in [col for col in df.columns if col.endswith("mzML")]:
        meta_values.append(df[s])
        hovertemplate += (
            s[:-5] + ": %{customdata[" + str(len(meta_values) - 1) + "]}<br>"
        )

    customdata = np.stack(meta_values, axis=-1)

    fig.add_trace(
        go.Scattergl(
            name="feature",
            x=df["RT"],
            y=df["mz"],
            mode="markers",
            marker_color=df[sample],
            marker_symbol="square",
            marker_size=12,
            customdata=customdata,
            hovertemplate=hovertemplate,
        )
    )
    fig.update_layout(
        title=f"consensus map with intensities from {sample}",
        xaxis_title="retention time",
        yaxis_title="m/z",
        showlegend=False,
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.update_traces(
        showlegend=False,
        marker_colorscale=COLOR_SCALE,
        selector=dict(type="scattergl"),
    )

    return fig


def plot_consensus_feature_chromatograms(feature_maps: dict, to_plot: dict, title: str):

    fig = go.Figure()

    for sample, id in to_plot.items():
        df = feature_maps[sample][feature_maps[sample]["feature_id"] == id]

        fig.add_trace(
            go.Scattergl(
                name=sample,
                x=df.loc[df.index[0], "chrom_rts"],
                y=df.loc[df.index[0], "chrom_intensities"],
                mode="lines",
            )
        )
    fig.update_layout(
        showlegend=True,
        title=title,
        xaxis_title="time",
        yaxis_title="intensity (counts per second)",
        legend_title="sample",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


"""
Contains all visualization blocks as modules with streamlit components.
"""


def display_feature_data(feature_maps, spectra, feature_detection_method="FFM"):
    name = st.selectbox("choose feature map", [
                        key for key in feature_maps.keys()], key=feature_detection_method)
    if name:
        df = feature_maps[name]
        c1, c2, c3 = st.columns(3)
        c1.metric("number of features", df.shape[0])
        c2.metric("lowest intensity", int(min(df["intensity"])))
        c3.metric("highest intensity", int(max(df["intensity"])))

        percentile = st.slider("show top intensity features (%)", 10, 100, 10, 10, key=feature_detection_method+"percentile")
        index = int((percentile/10) * (df["intensity"].size/10) - 1)
        int_cutoff = df["intensity"].sort_values(ascending=False).iloc[index]

        df = df[df["intensity"] > int_cutoff]
        df.sort_values(by=["intensity"], inplace=True)
        df.fillna(0, inplace=True)
        # 2D feature map
        fig = plot_feature_map(df)
        show_fig(fig, f"feature-map-{name}")

        selected_feature = st.selectbox(
            "detailed feature view", df["metabolite"])
        feature = df[df["metabolite"] == selected_feature]
        # get start and end RT/mz values to filter peak map
        bb = feature[["RTstart", "RTend", "MZstart", "MZend"]]

        # Define the values for mzmin and mzmax
        mzmin = bb.loc[bb.index[0], "MZstart"]
        mzmax = bb.loc[bb.index[0], "MZend"]

        # # filter df by retention time...
        df_peaks = spectra[name][
            (spectra[name]["RT"] > bb.loc[bb.index[0], "RTstart"])
            & (spectra[name]["RT"] < (bb.loc[bb.index[0], "RTend"]))
        ]

        # Create a function to filter the entries based on the condition
        def filter_entries(row):
            mz_array = row['mzarray']
            int_array = row['intarray']
            filtered_indices = np.where((mz_array > mzmin) & (mz_array < mzmax))
            filtered_mz = mz_array[filtered_indices]
            filtered_int = int_array[filtered_indices]
            return filtered_mz, filtered_int

        # Apply the filter_entries function to the DataFrame and create new columns with the filtered results
        df_peaks[['filtered_mzarray', 'filtered_intarray']] = df_peaks.apply(filter_entries, axis=1, result_type='expand')

        # If you only want the entries that have some values in both 'mzarray' and 'intarray'
        # you can drop rows with empty arrays in the filtered columns
        df_peaks = df_peaks.dropna(subset=['filtered_mzarray', 'filtered_intarray'])

        # Drop the original 'mzarray' and 'intarray' columns if needed
        df_peaks = df_peaks.drop(columns=['mzarray', 'intarray'])

        # Rename to original name
        df_peaks = df_peaks.rename(columns={"filtered_mzarray": "mzarray", "filtered_intarray": "intarray"})


        c1, c2, c3, c4 = st.columns(4)
        c1.metric("m/z", feature["mz"].round(5))
        c2.metric("intensity", feature["intensity"].astype(int))
        c3.metric("charge", feature["charge"])
        c4.metric("adduct", feature.loc[feature.index[0], "adduct"])
        c1.metric("retention time", feature["RT"].round(2))
        c2.metric(
            "retention time range",
            np.round(feature["RTend"] - feature["RTstart"], 2),
        )
        c3.metric("FWHM", feature["fwhm"].round(2))
        if "snr" in feature.columns:
            c4.metric("S/N", feature["snr"].round(2))

        # chromatogram monoisotopic mass
        fig = plot_feature_chromatogram(
            df[df["metabolite"] == selected_feature])
        show_fig(fig, f"chromatogram-{selected_feature}")

        # 2D chromatogram peak map of feature area
        fig = plot_peak_map_2D(df_peaks, 0)
        show_fig(fig, f"feature-area-{selected_feature}")


def display_map_alignement(feature_maps: dict):
    name = st.selectbox(
        "choose feature map",
        [
            key
            for key, value in feature_maps.items()
            if not value["original_rt"].isnull().any()
        ],
    )
    if name:
        df = feature_maps[name]
        c1, c2 = st.columns(2)
        c1.metric(
            "reference map (contains most features)",
            [
                key
                for key, value in feature_maps.items()
                if value["original_rt"].isnull().any()
            ][0],
        )
        c2.metric("mean RT deviation", np.mean(
            df["RT"] - df["original_rt"]).round(2))
        fig = plot_feature_map_alignment(df)
        show_fig(fig, "map-alignment")


def display_consensus_map(df: pd.DataFrame, feature_maps: dict):
    # calculate means of intensities
    df = df.assign(
        means=df[[col for col in df.columns if col.endswith(".mzML")]].mean(
            axis=1)
    )


    # select a sample or all for mean intensities
    sample = st.selectbox(
        "choose sample to highlight",
        ["Show mean intensities"]
        + [col for col in df.columns if col.endswith(".mzML")],
    )

    # select "means" column if selected
    if sample == "Show mean intensities":
        sample = "means"

    percentile = st.slider("show top intensity features (%)", 10, 100, 10, 10, key="consensus-map-percentile")
    top_n = int((percentile/10) * (df[sample].size/10) - 1)

    # sort df by sample value to plot high intensities on top
    df.sort_values(by=[sample], inplace=True)

    # from the sorted dataframe, get only the last n rows from the top n feature intensities
    df = df.tail(top_n)

    fig = plot_consensus_map(df, sample)
    show_fig(fig, "consensus-map")

    # chose a feature for comparison chromatogram plot
    feature = st.selectbox(
        "choose metabolite for intensity comparison",
        [m for m in df["metabolite"]],
    )

    # determine which features in which samples to plot as chromatogram
    to_plot = {
        k: v
        for k, v in df[df["metabolite"] == feature]["sample_to_id"].iloc[0].items()
        if v
    }

    fig = plot_consensus_feature_chromatograms(feature_maps, to_plot, feature)
    show_fig(fig, f"consensus-feature-chromatogram-{feature}")
