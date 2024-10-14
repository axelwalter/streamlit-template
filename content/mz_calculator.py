import streamlit as st
import pandas as pd
from pathlib import Path
from src.common.common import *
from src.masscalculator import (
    create_compound,
    build_compound,
    save_df,
    HELP,
    validate_dataframe,
)

params = page_setup(page="workflow")

with st.sidebar:
    st.markdown(HELP)

st.title("m/z calculator")
results_only = st.toggle("view results only")

input_table_path = Path(st.session_state.workspace, "mass-calculator.csv")

if not input_table_path.exists():
    pd.DataFrame(
        columns=["name", "sum formula", "adduct", "mz", "RT", "peak width", "comment"]
    ).to_csv(input_table_path, index=False)

tabs = st.tabs(["➕ New", "📟 Combine metabolites", "⬆️ Upload existing table"])

if not results_only:
    with tabs[0]:
        with st.form("new-metabolite-form"):
            c1, c2 = st.columns(2)
            formula = c1.text_input(
                "**sum formula**", "", help="Enter sum formula for new metabolite."
            )
            name = c1.text_input(
                "metabolite name (optional)",
                "",
                help="Will be created automatically if omitted.",
            )
            neutral_loss = c1.text_input(
                "neutral losses (optional)",
                "",
                help="Sum formula of neutral losses (e.g. H2O).",
            )
            add_adduct_info = c1.checkbox(
                "add adduct info to name",
                True,
                help="Always add adduct information to metabolite name separated by `#` if a custom name was chosen.",
            )
            charge = c2.number_input(
                "**charge**",
                -50,
                50,
                1,
                help="Enter charge. Negative numbers for negative ion mode, positive numbers for positive ion mode.",
            )
            c2.markdown(
                "addtional adducts",
                help="Specify adducts except for protons (H) up the number of charges in total, the remaing will be filled with protons (positive mode). In negative mode as the absolute charge number of protons will be removed regardless of specified additional adducts.",
            )
            adducts = c2.data_editor(
                pd.DataFrame({"adduct": ["Na", "K", "HCOOH"], "number": [0, 0, 0]}),
                hide_index=True,
                use_container_width=True,
            )
            add_both_adducts = c2.checkbox(
                "add two entries: protons only **and** with additional adducts",
                False,
                help="Add to entries to table. One contains only addition or loss of protons, the other considers the additional adduct table. Useful to include e.g. always the sodium adduct as well: `metabolite#[M+H]+` and `metabolite#[M+Na]+`.",
            )
            create_compound_button = st.form_submit_button(
                "Add new metabolite",
                use_container_width=True,
                help="Calculate m/z from sum formula and adduct and add metabolite to table.",
            )
        if create_compound_button:
            if add_both_adducts and adducts["number"].sum() > 0:
                save_df(
                    pd.concat(
                        [
                            create_compound(
                                formula,
                                charge,
                                pd.DataFrame({"adduct": [], "number": []}),
                                "",
                                name,
                                True,
                            ),
                            create_compound(
                                formula, charge, adducts, neutral_loss, name, True
                            ),
                        ]
                    ),
                    input_table_path,
                )
            else:
                save_df(
                    create_compound(
                        formula, charge, adducts, neutral_loss, name, add_adduct_info
                    ),
                    input_table_path,
                )

    with tabs[1]:
        with st.form("build-metabolite-form"):
            st.markdown(
                "**metabolites to combine**",
                help="Select from metabolites which are already in the table to combine them into larger molecules from the given numbers.",
            )
            column_types = {"metabolite": "str", "number": "int"}
            builder = st.data_editor(
                pd.DataFrame(columns=column_types.keys()).astype(column_types),
                hide_index=True,
                use_container_width=True,
                num_rows="dynamic",
                column_config={
                    "metabolite": st.column_config.SelectboxColumn(
                        "metabolite",
                        width="large",
                        options=[
                            x
                            for x in pd.read_csv(input_table_path)[
                                ["name", "sum formula"]
                            ]
                            .dropna()["name"]
                            .tolist()
                        ],
                        required=True,
                    ),
                    "number": st.column_config.NumberColumn(
                        "number",
                        width="small",
                        min_value=-100,
                        max_value=100,
                        step=1,
                        required=True,
                        default=1,
                    ),
                },
            )
            c1, c2 = st.columns(2)
            charge = c1.number_input(
                "**charge**",
                -50,
                50,
                1,
                help="Enter charge. Negative numbers for negative ion mode, positive numbers for positive ion mode.",
            )
            name = c1.text_input(
                "metabolite name (optional)",
                "",
                help="Will be created automatically if omitted.",
            )
            elimination = c1.text_input(
                "elimination product (optional)",
                "H2O",
                help="Remove elemination product when combining two metabolites.",
            )
            add_adduct_info = c1.checkbox(
                "add adduct info to name",
                True,
                help="Always add adduct information to metabolite name separated by `#` if a custom name was chosen.",
            )
            c2.markdown(
                "additional adducts",
                help="Specify adducts except for protons (H) up the number of charges in total, the remaing will be filled with protons (positive mode). In negative mode as the absolute charge number of protons will be removed regardless of specified additional adducts.",
            )
            adducts = c2.data_editor(
                pd.DataFrame({"adduct": ["Na", "K", "HCOOH"], "number": [0, 0, 0]}),
                hide_index=True,
                use_container_width=True,
            )
            add_both_adducts = c2.checkbox(
                "add two entries: protons only **and** with additional adducts",
                False,
                help="Add to entries to table. One contains only addition or loss of protons, the other considers the additional adduct table. Useful to include e.g. always the sodium adduct as well: `metabolite#[M+H]+` and `metabolite#[M+Na]+`.",
            )
            build_compound_button = st.form_submit_button(
                "Calculate metabolite",
                use_container_width=True,
                help="Calculate m/z from sum formula and adduct and add metabolite to table.",
            )

        if build_compound_button:
            save_df(
                build_compound(
                    builder,
                    charge,
                    adducts,
                    name,
                    pd.read_csv(input_table_path),
                    elimination,
                    add_adduct_info,
                    add_both_adducts,
                ),
                input_table_path,
            )

    with tabs[2]:
        # Create file uploader
        uploaded_file = st.file_uploader("Upload a CSV file", type="csv")

        if uploaded_file is not None:
            # Read the uploaded file into a DataFrame
            try:
                df = pd.read_csv(uploaded_file)

                # Validate the DataFrame
                if validate_dataframe(df):
                    if st.button(
                        "Replace current table with uploaded table (will delete current data).",
                        type="primary",
                        use_container_width=True,
                    ):
                        # Save the table (customize the path as needed)
                        df.to_csv(input_table_path, index=False)
                        st.session_state["mz_calc_success"] = [
                            f"Uploaded table **{uploaded_file.name}**"
                        ]
                        st.rerun()
                else:
                    st.error(
                        "The uploaded file does not match the required format or data types."
                    )

            except Exception as e:
                st.error(f"An error occurred: {e}")

if "mz_calc_success" in st.session_state:
    for message in st.session_state["mz_calc_success"]:
        st.success(message)
if "mz_calc_error" in st.session_state:
    st.error(st.session_state["mz_calc_error"])

edited = st.data_editor(
    pd.read_csv(
        input_table_path,
        dtype={
            "name": str,
            "sum formula": str,
            "adduct": str,
            "mz": float,
            "RT": float,
            "peak width": float,
            "comment": str,
        },
    ),
    use_container_width=True,
    hide_index=True,
    key="mass-table",
    disabled=("sum formula", "adduct", "mz"),
    num_rows="dynamic",
)

if (
    st.session_state["mass-table"]["edited_rows"]
    or st.session_state["mass-table"]["deleted_rows"]
    or st.session_state["mass-table"]["added_rows"]
):
    if edited["name"].duplicated().any():
        st.error("Metabolite names need to be unique.")
    else:
        edited.to_csv(input_table_path, index=False)
        st.rerun()

save_params(params)
