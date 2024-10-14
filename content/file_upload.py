from pathlib import Path

import streamlit as st

from src.common.common import *
from src.fileupload import *

params = page_setup()

st.title("mzML Files")
df_path = Path(st.session_state.workspace, "mzML-files.tsv")
mzML_dir = Path(st.session_state.workspace, "mzML-files")

tabs = ["File Upload", "Example Data"]
if st.session_state.location == "local":
    tabs.append("Files from local folder")
elif st.session_state.location == "online":
    tabs.append("Download mzML files")

tabs = st.tabs(tabs)

with tabs[0]:
    with st.form("mzML-upload", clear_on_submit=True):
        files = st.file_uploader(
            "mzML files", accept_multiple_files=(st.session_state.location == "local"))
        _, c2, _ = st.columns(3)
        if c2.form_submit_button("Add files to workspace", use_container_width=True, type="primary"):
            if files:
                save_uploaded_mzML(files)
                update_mzML_df(df_path, mzML_dir).to_csv(df_path, sep="\t", index=False)
                st.rerun()
            else:
                st.error("Nothing to add, please upload file.")

# Example mzML files
with tabs[1]:
    st.markdown("Example data set of bacterial cytosolic fractions. Bacillus subtilis cultures were treated with the antibiotic fosfomycin, which inhibits a step in the biosynthesis of petidoglycan (bacterial cell wall). The major accumulation product is UDP-GlcNAc.")
    _, c2, _ = st.columns(3)
    if c2.button("Load Example Data", type="primary", use_container_width=True):
        load_example_mzML_files()
        update_mzML_df(df_path, mzML_dir).to_csv(df_path, sep="\t", index=False)
        st.rerun()

# Local file upload option: via directory path
if st.session_state.location == "local":
    with tabs[2]:
        # with st.form("local-file-upload"):
        local_mzML_dir = st.text_input(
            "path to folder with mzML files")
        # raw string for file paths
        local_mzML_dir = r"{}".format(local_mzML_dir)
        _, c2, _ = st.columns(3)
        if c2.button("Copy files to workspace", type="primary", use_container_width=True, disabled=(local_mzML_dir == "")):
            copy_local_mzML_files_from_directory(local_mzML_dir)
            update_mzML_df(df_path, mzML_dir).to_csv(df_path, sep="\t", index=False)
            st.rerun()
elif st.session_state.location == "online":
    with tabs[2]:
        c1, c2 = st.columns(2)
        if c1.button("Get all mzML files in workspace as zip file.", type="primary", use_container_width=True):
            zip_buffer = zip_files(Path(st.session_state.workspace, "mzML-files"))
            c2.download_button(
                label="⬇️ Download Now",
                data=zip_buffer,
                file_name=f"mzML_files-{Path(st.session_state.workspace).stem}.zip",
                mime="application/zip",
                type="primary",
                use_container_width=True
            )


df = update_mzML_df(df_path, mzML_dir)

if any(Path(mzML_dir).iterdir()):
    with st.form("select-mzML-files"):
        c1, c2 = st.columns(2)
        # Display all mzML files currently in workspace
        c1.markdown("##### select mzML files for analysis:")
        select_submit = c2.form_submit_button("💾 Save file selection", type="primary", use_container_width=True)
        edited = st.data_editor(df, hide_index=True, use_container_width=True, disabled=["file name"], key="mzML-files-df")
    if select_submit and (st.session_state["mzML-files-df"]["edited_rows"] or st.session_state["mzML-files-df"]["deleted_rows"] or st.session_state["mzML-files-df"]["added_rows"]):
        edited.to_csv(df_path, sep="\t", index=False)

    # Remove files
    with st.form("remove-mzML-files"):
        st.markdown("🗑️ Remove mzML files")
        to_remove = st.multiselect("select mzML files",
                                options=[f.stem for f in sorted(mzML_dir.iterdir())])
        c1, c2 = st.columns(2)
        if c2.form_submit_button("Remove **selected**", use_container_width=True):
            remove_selected_mzML_files(to_remove)
            update_mzML_df(df_path, mzML_dir).to_csv(df_path, sep="\t", index=False)
            st.rerun()

        if c1.form_submit_button("⚠️ Remove **all**", use_container_width=True):
            remove_all_mzML_files()
            update_mzML_df(df_path, mzML_dir).to_csv(df_path, sep="\t", index=False)
            st.rerun()
