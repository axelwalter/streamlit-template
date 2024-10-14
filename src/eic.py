import streamlit as st

from .common.common import *

from pathlib import Path

from pyopenms import *

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

HELP = """
This page allows you to extract chromatograms from mzML files and perform various operations on the extracted data. Below is a guide to help you navigate and use the features of this page effectively.

##### 1. Define EIC input metabolite table
- In the metabolite table, you can define metabolite m/z and RT values for EIC extraction.
- Hovering on the edges of the table will allow you to add more rows, or select columns which you can delete with the delete key on your keyboard.
- To download the modified table, click on the "Download" button which appears in the top right corner hovering over the table.
- When executing the workflow the current table will be saved to your workspace, no need to always download and re-upload it!
- 💡 Combine intensities of metabolite variants (e.g. different adducts) using the `#` symbol in the metabolite name. E.g. `GlcNAc#[M+H]+` and `GlcNAc#[M+Na]+`. Make sure to check the **combine variants** box in the result section.

##### 2. Calculate mass and add to table
- You can calculate the mass of a compound using its sum formula and adduct and add it to the table.
- Click on the "Add Metabolite" button to calculate the mass and add the compound to the table.
- Edit the compound name in the table (optional).

##### 3. Chromatogram extraction parameters
- In the "Settings" section, you can customize the parameters for chromatogram extraction.
- Specify the time unit for retention time, default peak width, mass tolerance unit, mass tolerance (in ppm and Da), noise threshold, and whether to combine variants of metabolites.
- Make sure to set these parameters before proceeding with chromatogram extraction.

##### 4. Extract chromatograms
- After setting the parameters, click on the "Extract chromatograms" button to extract the chromatograms from the uploaded mzML files.
- Note: You need to upload/select mzML files before extracting chromatograms.

##### 5. Results
- In the "Results" section, you can explore the extracted chromatogram results.
- The section contains several tabs:
  - "📊 Summary": Displays a summary table of the extracted chromatogram data.
  - "📈 Samples": Allows you to compare and visualize chromatograms for different samples.
  - "📈 Metabolites": Displays overlayed EICs (Extracted Ion Chromatograms) for each metabolite.
  - "📁 Chromatogram data": Provides a download button to download the chromatogram data as a zip file.
  - "📁 Meta data": Displays and allows you to edit the meta data associated with the extracted chromatograms.

Feel free to explore the different features and options on this page to extract and analyze your chromatogram data efficiently.
"""

def extract_chromatograms(results_dir, mzML_files, df_input, mz_unit, mz_ppm, mz_da, time_unit, default_peak_width, baseline):
    with st.status("Extracting chromatograms...") as status:
        # Check for unique index
        if not df_input["name"].is_unique:
            st.error("Metabolite names need to be unique.")
            status.update(label="Error!", state="error", expanded=True)
            return
        
        if any(df_input["name"].isna()):
            st.error("Enter a name for each metabolite.")
            status.update(label="Error!", state="error", expanded=True)
            return

        # Drop all rows without mz value
        df_input = df_input[df_input['mz'].notna()]

        # Get RT times in seconds for extraction
        if time_unit == "minutes":
            df_input["RT"] = df_input["RT"]*60
            df_input["peak width"] = df_input["peak width"]*60

        # Delete and re-create results directory
        reset_directory(Path(results_dir))

        # To make a zip file with tables in tsv format later, create a directory
        tsv_dir = Path(results_dir, "tsv-tables")
        reset_directory(tsv_dir)

        # Create an empty df for AUCs with filenames as columns and mass names as indexes
        df_auc = pd.DataFrame(
            columns=[Path(file).name for file in mzML_files], index=df_input["name"]
        )

        # Iterate over the files and extract chromatograms in a single dataframe per file
        for file in mzML_files:
            st.write(f"Extracting chromatograms from {Path(file).name} ...")
            # Load mzML file into exp
            exp = MSExperiment()
            MzMLFile().load(str(file), exp)

            # Create an empty dataframe to collect chromatogram data in
            df = pd.DataFrame()

            # get BPC and time always for each file
            time = []
            ints = []
            for spec in exp:
                if spec.getMSLevel() == 2:
                    continue
                _, intensities = spec.get_peaks()
                rt = spec.getRT()
                if time_unit == "minutes":
                    rt = rt / 60
                time.append(rt)
                i = int(max(intensities))
                ints.append(i)
            df["time"] = time
            df["BPC"] = ints

            # get EICs
            for mass, name, rt, peak_width in zip(df_input["mz"],
                                                  df_input["name"],
                                                  df_input["RT"],
                                                  df_input["peak width"]):
                if name:
                    metabolite_name = name
                else:
                    metabolite_name = str(mass)

                # If RT information is given, determine RT borders
                if rt:
                    # Custom peak width per metabolite
                    if not np.isnan(peak_width):
                        rt_min = rt-(peak_width/2)
                        rt_max = rt+(peak_width/2)
                    else:
                        rt_min = rt - default_peak_width/2
                        rt_max = rt + default_peak_width/2

                # List with intensity values (len = number of MS1 spectra)
                ints = []
                for spec in exp:
                    # Skip MS2 spectra
                    if spec.getMSLevel() == 2:
                        continue
                    # Skip Spectra which are not within RT range
                    if not np.isnan(rt):
                        if rt_min > spec.getRT() or rt_max < spec.getRT():
                            ints.append(0)
                            continue
                    # Find the highest peak in spec within mz window (Da)
                    if mz_unit == "Da":
                        index_highest_peak_within_window = (
                            spec.findHighestInWindow(
                                mass, mz_da, mz_da
                            )
                        )
                    else:
                        # Find the highest peak in spec within mz window (ppm)
                        ppm_window = float((mz_ppm / 1000000) * mass)
                        index_highest_peak_within_window = (
                            spec.findHighestInWindow(
                                mass, ppm_window, ppm_window)
                        )
                    # Get intensity at highest peak index if hit was found
                    if index_highest_peak_within_window > -1:
                        i = int(
                            spec[
                                index_highest_peak_within_window
                            ].getIntensity()
                        )
                        # Check if higher then baseline and
                        if i > baseline:
                            ints.append(i)
                        else:
                            ints.append(0)
                    # If nothing was found, XIC intensity for this spec is zero
                    else:
                        ints.append(0)

                # Add intensity values to df
                df[metabolite_name] = ints
                # also insert the AUC in the auc dataframe
                df_auc.loc[metabolite_name, Path(file).name] = np.trapz(ints, df["time"])

            # Save to feather dataframe for quick access
            df.to_feather(Path(results_dir, Path(file).stem + ".ftr"))
            # Save as tsv for download option
            df.to_csv(
                Path(tsv_dir, Path(file).stem + ".tsv"), sep="\t", index=False
            )

        # once all files are processed, zip the tsv files and delete their directory
        shutil.make_archive(os.path.join(
            results_dir, "chromatograms"), "zip", tsv_dir)
        shutil.rmtree(tsv_dir)

        # save summary of AUC values to tsv file
        df_auc.index.name = "metabolite"
        df_auc = df_auc.reindex(sorted(df_auc.columns), axis=1)
        df_auc = df_auc[~(df_auc == 0).all(axis=1)]
        if df_auc.empty:
            st.error("No metabolites detected from given input.")
            status.update(label="Error", state="error", expanded=True)
            return
        df_auc.to_csv(Path(results_dir, "summary.tsv"), sep="\t")

        # sum intensities of variants of the same metabolite
        combined = {}
        metabolite_names = list(
            set([c.split("#")[0] for c in df_auc.index]))
        for filename in df_auc.columns:
            aucs = []
            for a in metabolite_names:
                auc = 0
                for b in [
                    b
                    for b in df_auc.index
                    if ((a + "#" in b and b.startswith(a)) or a == b)
                ]:
                    auc += df_auc.loc[b, filename]
                aucs.append(auc)
            combined[filename] = aucs
        df_auc_combined = pd.DataFrame(combined)
        df_auc_combined.set_index(pd.Index(metabolite_names), inplace=True)
        df_auc_combined.index.name = "metabolite"
        df_auc_combined = df_auc_combined.reindex(sorted(df_auc_combined.columns), axis=1)
        df_auc_combined.to_csv(Path(results_dir, "summary-combined.tsv"), sep="\t")

        # Save AUC to text file
        with open(Path(results_dir, "run-params.txt"), "w") as f:
            f.write(f"{baseline}\n{time_unit}")
    # # Re-run to prevent tab jumping back to first tab upon first widget change (streamlit bug)
    st.rerun()


@st.cache_resource
def get_auc_fig(df_auc):
    for col in df_auc.columns:
        df_auc = df_auc.rename(columns={col: col[:-5]})
    fig = px.bar(df_auc.T, barmode="group")
    # fig.update_traces(width=0.2)
    fig.update_layout(
        xaxis_title="",
        yaxis_title="area under curve",
        legend_title="metabolite",
        plot_bgcolor="rgb(255,255,255)"
    )
    fig.layout.template = "plotly_white"
    return fig


@st.cache_resource
def get_metabolite_fig(df_auc, metabolite, time_unit):
    fig = go.Figure()
    for sample in df_auc.columns:
        df = pd.read_feather(Path(st.session_state.workspace,
                             "extracted-ion-chromatograms", sample[:-4] + "ftr"))
        fig.add_trace(
            go.Scattergl(
                name=sample[:-5],
                x=df["time"],
                y=df[[metabolite][0]],
            )
        )
    fig.update_layout(
        title=metabolite,
        xaxis_title=f"time ({time_unit})",
        yaxis_title="intensity (counts per second)",
        legend_title="sample",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig


def get_sample_plot(df, sample, time_unit):
    fig = px.line(df, x="time", y=df.columns)
    fig.update_layout(
        title=sample[:-5],
        xaxis_title=f"time ({time_unit})",
        yaxis_title="intensity (counts per second)",
        legend_title="metabolite",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    return fig
