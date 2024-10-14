import streamlit as st
from pathlib import Path

if __name__ == '__main__':
    pages = {
        "UmetaFlow App" : [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
            st.Page(Path("content", "file_upload.py"), title="File Upload", icon="📂"),
            st.Page(Path("content", "raw_data_viewer.py"), title="View MS data", icon="👀"),
        ],
        "Untargeted" : [
            st.Page(Path("content", "umetaflow_pyopenms.py"), title="UmetaFlow pyOpenMS", icon="🐍"),
            st.Page(Path("content", "umetaflow_topp.py"), title="UmetaFlow TOPP", icon="🚀"),
        ],
        "Targeted": [
            st.Page(Path("content", "mz_calculator.py"), title="m/z Calculator", icon="📟"),
            st.Page(Path("content", "extracted_ion_chromatograms.py"), title="Extracted Ion Chromatograms", icon="🔍"),
        ],
        "Downstream": [
            st.Page(Path("content", "statistics.py"), title="Statistics", icon="📈"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()