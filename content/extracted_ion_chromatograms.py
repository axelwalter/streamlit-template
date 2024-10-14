import streamlit as st
from src.common import *
from src.eic import *
from src.masscalculator import get_mass, check_formula

from pathlib import Path
import pandas as pd
import numpy as np

params = page_setup(page="workflow")

with st.sidebar:
    st.markdown(HELP)

# c1, c2 = st.columns([0.7, 0.3])
st.title("Extracted Ion Chromatograms")
results_only = st.toggle("view results only")

results_dir = Path(st.session_state.workspace, "extracted-ion-chromatograms")

if not results_only:
    input_table_path = Path(st.session_state.workspace, "EIC-input-table.csv")
    if input_table_path.exists():
        # Read dataframe from path (defined in src.eic)
        df = pd.read_csv(input_table_path)
    else:
        # Load a default example df
        df = pd.DataFrame(
            {"name": [""], "mz": [np.nan], "RT": [np.nan], "peak width": [np.nan]}
        )

    use_mz_calculator_table = st.toggle(
        "Use metabolite table from m/z calculator",
        params["eic_use_mz_table"],
        key="eic_use_mz_table",
    )
    with st.form("eic_form", border=True):
        if not use_mz_calculator_table:
            st.markdown(
                "metabolite table",
                help="""
**Input table for EIC extraction.**

Add metabolites whith their names and m/z values.

Optionally define a retention time value. If no peak width is specified, the default peak width defined in parameter section will be used.
RT and peak width can be entered in seconds or minutes, depending on the time unit setting in parameter section.

New metabolites can be entered by entering sum formula and adduct information as well.

💡 Combine intensities of metabolite variants (e.g. different adducts) using the `#` symbol in the metabolite name. E.g. `GlcNAc#[M+H]+` and `GlcNAc#[M+Na]+`. Make sure to check the **combine variants** box in the result section.

To download the modified table, click on the **Download** button which appears in the top right corner hovering over the table.

To paste a data table from Excel simply select all the cells in Excel, select the top left cell in the metabolite table (turns red) and paste with **Ctrl-V**.
""",
            )
            edited = st.data_editor(df, use_container_width=True, num_rows="dynamic")
            c1, c2, c3 = st.columns(3)
            formula = c1.text_input(
                "sum formula",
                "",
                help="Enter a neutral sum formula for a new compound in the table.",
            ).strip()
            adduct = c2.selectbox(
                "adduct",
                [
                    "[M+H]+",
                    "[M+Na]+",
                    "[M+2H]2+",
                    "[M-H2O+H]+",
                    "[M-H]-",
                    "[M-2H]2-",
                    "[M-H2O-H]-",
                ],
                help="Specify the adduct.",
            )
            v_space(1, c3)
            add_compound_button = c3.form_submit_button(
                "Add Metabolite",
                use_container_width=True,
                help="Calculate m/z from sum formula and adduct and add metabolite to table.",
            )
            if add_compound_button:
                if check_formula(formula):
                    mz = get_mass(formula, adduct)
                    if mz:
                        compound_name = f"{formula}#{adduct}"
                        new_row = pd.DataFrame(
                            {
                                "name": [compound_name],
                                "mz": [mz],
                                "RT": [np.nan],
                                "peak width": [np.nan],
                            }
                        )
                        edited = pd.concat([edited, new_row], ignore_index=True).to_csv(
                            input_table_path, index=False
                        )
                        st.rerun()
                    else:
                        st.warning(
                            "Can not calculate mz of this formula/adduct combination."
                        )
                else:
                    st.warning("Invalid formula.")
        c1, c2, c3 = st.columns(3)
        c1.radio(
            "time unit",
            ["seconds", "minutes"],
            index=["seconds", "minutes"].index(params["eic_time_unit"]),
            key="eic_time_unit",
            help="Retention time unit.",
        )
        c2.number_input(
            "default peak width (seconds)",
            1,
            600,
            params["eic_peak_width"],
            5,
            key="eic_peak_width",
            help="Default value for peak width. Used when a retention time is given without peak width. Adding a peak width in the table will override this setting.",
        )
        c3.number_input(
            "**noise threshold**",
            0,
            1000000,
            params["eic_baseline"],
            100,
            key="eic_baseline",
            help="Peaks below the treshold intensity will not be extracted.",
        )
        # Mass tolerance settings
        c1, c2, c3 = st.columns(3)
        c1.radio(
            "mass tolerance unit",
            ["ppm", "Da"],
            index=["ppm", "Da"].index(params["eic_mz_unit"]),
            key="eic_mz_unit",
        )
        c2.number_input(
            "**mass tolerance ppm**",
            1,
            100,
            params["eic_tolerance_ppm"],
            step=5,
            key="eic_tolerance_ppm",
        )
        c3.number_input(
            "mass tolerance Da",
            0.01,
            10.0,
            params["eic_tolerance_da"],
            0.05,
            key="eic_tolerance_da",
        )
        c1, c2 = st.columns(2)
        save_button = c1.form_submit_button(
            "💾 Save parameters",
            use_container_width=True,
            help="Save selected paramters to your workspace.",
        )
        submitted = c2.form_submit_button(
            "Extract chromatograms", type="primary", use_container_width=True
        )

    v_space(1)

    if save_button:
        if not use_mz_calculator_table:
            edited.to_csv(input_table_path, index=False)
        save_params(params)

    if submitted:
        if not use_mz_calculator_table:
            edited.to_csv(input_table_path, index=False)
        save_params(params)
        df_path = Path(st.session_state.workspace, "mzML-files.tsv")

        if not df_path.exists():
            mzML_files = None
        else:
            df = pd.read_csv(df_path, sep="\t")

            # Filter the DataFrame for files where "use in workflow" is True
            selected_files = df[df["use in workflows"] == True]["file name"].tolist()

            # Construct full file paths
            mzML_files = [
                Path(st.session_state.workspace, "mzML-files", file_name)
                for file_name in selected_files
            ]
        if not mzML_files:
            st.warning("Upload/select some mzML files first!")
        else:
            if use_mz_calculator_table:
                data = pd.read_csv(
                    Path(st.session_state.workspace, "mass-calculator.csv")
                )[["name", "mz", "RT", "peak width"]]
            else:
                data = edited

            if not data.empty:
                extract_chromatograms(
                    results_dir,
                    mzML_files,
                    data,
                    st.session_state["eic_mz_unit"],
                    st.session_state["eic_tolerance_ppm"],
                    st.session_state["eic_tolerance_da"],
                    st.session_state["eic_time_unit"],
                    st.session_state["eic_peak_width"],
                    st.session_state["eic_baseline"],
                )
            else:
                st.error("No input m/z values provided.")

path = Path(results_dir, "summary.tsv")
if path.exists():
    st.checkbox(
        "combine metabolite variants",
        params["eic_combine"],
        help="Combines different variants (e.g. adducts or neutral losses) of a metabolite. Put a `#` with the name first and variant second (e.g. `glucose#[M+H]+` and `glucose#[M+Na]+`)",
        key="eic_combine",
    )
    tabs = st.tabs(
        [
            "📊 Summary",
            "📈 Samples",
            "📈 Metabolites",
            "📁 Chromatogram data",
            "📁 Meta data",
        ]
    )
    with open(Path(results_dir, "run-params.txt"), "r") as f:
        baseline = int(f.readline())
        time_unit = f.readline()

    with tabs[0]:
        if st.session_state["eic_combine"]:
            file_name = "summary-combined.tsv"
        else:
            file_name = "summary.tsv"

        df_auc = pd.read_csv(Path(results_dir, file_name), sep="\t").set_index(
            "metabolite"
        )
        # display the feature matrix and it's bar plot
        fig = get_auc_fig(df_auc)
        show_fig(fig, "xic-summary")
        show_table(df_auc, "feature-matrix-xic")

    with tabs[1]:
        file = st.selectbox(f"mzML file", df_auc.columns)
        c1, c2 = st.columns(2)
        show_bpc = c1.toggle("BPC", True, help="Show base peak chromatogram.")
        show_baseline = c2.toggle(
            "AUC baseline",
            True,
            help="Show baseline used for AUC calculation (from noise threshold parameter).",
        )

        df = pd.read_feather(Path(results_dir, file[:-4] + "ftr"))
        if show_baseline:
            df["AUC baseline"] = [baseline] * df.shape[0]
        if not show_bpc:
            df.drop(columns=["BPC"], inplace=True)
        fig = get_sample_plot(df, file, time_unit)
        show_fig(fig, file)

    with tabs[2]:
        # overlayed EICs for each sample
        if st.session_state["eic_combine"]:
            # read in complete auc file as df_auc to display EICs of different variants
            df_auc = pd.read_csv(Path(results_dir, "summary.tsv"), sep="\t").set_index(
                "metabolite"
            )
        metabolite = st.selectbox("select metabolite", df_auc.index)
        fig = get_metabolite_fig(df_auc, metabolite, time_unit)
        show_fig(fig, f"eic-{metabolite}")

    with tabs[3]:
        with open(os.path.join(results_dir, "chromatograms.zip"), "rb") as fp:
            btn = st.download_button(
                label="Download chromatogram data",
                data=fp,
                file_name=f"chromatogram-data.zip",
                mime="application/zip",
            )

    with tabs[4]:
        md = pd.DataFrame(
            {
                "filename": df_auc.columns,
                "ATTRIBUTE_Sample_Type": ["Sample"] * df_auc.shape[1],
            }
        ).set_index("filename")
        data = st.data_editor(md.T, num_rows="dynamic", use_container_width=True)
        st.download_button(
            "Download Table",
            data.T.to_csv(sep="\t").encode("utf-8"),
            "xic-meta-data.tsv",
        )
