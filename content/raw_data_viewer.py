from pathlib import Path

import streamlit as st

from src.common.common import page_setup
from src import view


params = page_setup()

st.title("View raw MS data")

# File selection can not be in fragment since it influences the subsequent sections
cols = st.columns(3)
selected_file = cols[0].selectbox(
    "choose file",
    [f.name for f in Path(st.session_state.workspace, "mzML-files").iterdir()],
    key="view_selected_file",
)
if selected_file:
    view.get_df(Path(st.session_state.workspace, "mzML-files", selected_file))

    tabs = st.tabs(
        ["📈 Peak map (MS1)", "📈 Spectra (MS1 + MS2)", "📈 Chromatograms (MS1)"]
    )
    with tabs[0]:
        view.view_peak_map()
    with tabs[1]:
        view.view_spectrum()
    with tabs[2]:
        view.view_bpc_tic()
