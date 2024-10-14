import streamlit as st
import pandas as pd
from pathlib import Path
from src.common import *
from src.stats import *

params = page_setup(page="workflow")

st.title("Statistics")
st.info(INFO)

result_files = {
    "Extracted Ion Chromatograms": Path(
        st.session_state.workspace,
        "extracted-ion-chromatograms",
        "summary.tsv",
    ),
    "Metabolomics": Path(
        st.session_state.workspace, "umetaflow-results", "FeatureMatrix.tsv"
    ),
    "Metabolomics MS1 annotated features": Path(
        st.session_state.workspace,
        "umetaflow-results",
        "MS1-annotations",
        "MS1-annotation.tsv",
    ),
    "Metabolomics MS2 annotated features": Path(
        st.session_state.workspace,
        "umetaflow",
        "MS2-annotations",
        "MS2-annotation.tsv",
    )}

options = []
for name, path in result_files.items():
    if path.exists():
        options.append(name)

cols = st.columns(2)
experiment = cols[0].selectbox("choose data set", options)

if experiment:
    df = pd.read_csv(result_files[experiment],
                     sep="\t").set_index("metabolite")

    samples = set(
        [col[:-5].split("#")[0]
            for col in df.columns if col.endswith(".mzML")]
    )

    tabs = st.tabs(
        ["📊 Sample comparison", "🔥 Clustering & Heatmap", "📁 Feature Matrix"])
    with tabs[0]:
        help = "💡 Samples with equal names before a `#` symbol (replicates) will be grouped together."
        c1, c2 = st.columns(2)
        a = c1.selectbox("samples A", samples, help=help)
        samples.remove(a)
        b = c2.selectbox("samples B", samples, help=help)

        if a != b:
            mean, change, std = get_mean_change_std(df, a, b)

            with st.expander("intermediate results"):
                st.markdown("**Mean intensities**")
                show_table(mean, "mean-intensities")
                st.markdown("**Standard deviations**")
                show_table(std, "standard-deviations")
                st.markdown("**Fold Changes**")
                show_table(change, "fold-changes")

            # show plots
            samples = mean.columns.tolist()
            features = mean.index.tolist()
            fig = mean_intensity_plot(samples, features, mean, std)
            show_fig(fig, "metabolite-intensities")

            fig = fold_change_plot(change)
            show_fig(fig, "fold-changes")

    with tabs[1]:
        df = scale_df(
            df[[col for col in df.columns if col.endswith(".mzML")]]
        )
        fig = dendrogram(df.T)
        show_fig(fig, "dendrogram")

        fig = heatmap(df)
        show_fig(fig, "heatmap")
    with tabs[2]:
        show_table(df, "feature-matrix")
