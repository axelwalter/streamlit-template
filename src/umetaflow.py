import streamlit as st
import shutil
from pathlib import Path
from .common.common import reset_directory
from .core import *
from .sirius import Sirius
from .gnps import GNPSExport
from .dataframes import DataFrames
from .blankremoval import BlankRemoval


class UmetaFlow:
    """
    Implements the steps for UmetaFlow in a reusable way.
    """

    def __init__(self, params, mzML_files, results_dir):
        self.params = params
        self.mzML_files = mzML_files
        self.results_dir = results_dir
        self.interim = Path(results_dir, "interim")
        reset_directory(self.interim)
        self.mzML_dir = Path(self.interim, "mzML_original")
        self.featureXML_dir = Path(self.interim, "FFM")
        self.sirius_ms_dir = ""
        self.consensusXML = Path(self.interim, "FeatureMatrix.consensusXML")
        self.consensus_tsv = Path(self.results_dir, "FeatureMatrix.tsv")
        self.mgf_file = str(Path(self.results_dir, "GNPS", "MS2.mgf"))

    def fetch_raw_data(self):
        reset_directory(self.mzML_dir)
        for file in self.mzML_files:
            shutil.copy(file, self.mzML_dir)

    def feature_detection(self):
        if self.params["ffm_single_traces"]:
            ffm_single_traces = "true"
        else:
            ffm_single_traces = "false"
        FeatureFinderMetabo().run(
            self.mzML_dir,
            os.path.join(self.interim, "FFM"),
            {
                "noise_threshold_int": float(self.params["ffm_noise"]),
                "chrom_fwhm": self.params["ffm_peak_width"],
                "chrom_peak_snr": self.params["ffm_snr"],
                "report_chromatograms": "true",
                "mass_error_ppm": self.params["ffm_mass_error"],
                "remove_single_traces": ffm_single_traces,
                "report_convex_hulls": "true",
                "min_fwhm": self.params["ffm_min_fwhm"],
                "max_fwhm": self.params["ffm_max_fwhm"],
            },
        )
        self.featureXML_dir = Path(self.interim, "FFM")

    def remove_blanks(self):
        blanks = [file + ".mzML" for file in self.params["blank_mzML_files"]]
        BlankRemoval.run(
            self.featureXML_dir,
            blanks,
            self.params["blank_cutoff"],
        )
        # remove blank mzML files
        for file in blanks:
            Path(self.mzML_dir, file).unlink()

    def adduct_detection(self):
        if self.params["ad_ion_mode"] == "negative":
            negative_mode = "true"
        else:
            negative_mode = "false"
        adducts = [line.encode()
                   for line in self.params["ad_adducts"].split("\n")]
        MetaboliteAdductDecharger().run(
            self.featureXML_dir,
            Path(self.interim, "FFM_decharged"),
            {
                "potential_adducts": adducts,
                "charge_min": self.params["ad_charge_min"],
                "charge_max": self.params["ad_charge_max"],
                "charge_span_max": 1,
                "max_neutrals": len(adducts) - 1,
                "negative_mode": negative_mode,
                "retention_max_diff": float(self.params["ad_rt_max_diff"]),
                "retention_max_diff_local": float(self.params["ad_rt_max_diff"]),
            },
        )
        shutil.rmtree(self.featureXML_dir)
        Path(self.interim, "FFM_decharged").rename(Path(self.featureXML_dir))

    def align_feature_maps(self):
        MapAligner().run(
            self.featureXML_dir,
            os.path.join(self.interim, "FFM_aligned"),
            os.path.join(self.interim, "Trafo"),
            {
                "max_num_peaks_considered": -1,
                "superimposer:mz_pair_max_distance": 0.05,
                "pairfinder:distance_MZ:max_difference": self.params["ma_mz_max"],
                "pairfinder:distance_MZ:unit": self.params["ma_mz_unit"],
                "pairfinder:distance_RT:max_difference": float(self.params["ma_rt_max"]),
            },
        )
        self.featureXML_dir = Path(self.interim, "FFM_aligned")

    def feature_maps_to_df(self, requant=False):
        df_path = Path(self.featureXML_dir).name + "_df"
        reset_directory(Path(self.interim, df_path))
        for file in Path(self.featureXML_dir).iterdir():
            DataFrames.featureXML_to_ftr(
                file, Path(self.interim, df_path), requant=requant
            )

    def align_peak_maps(self):
        MapAligner().run(
            self.mzML_dir,
            Path(self.interim, "mzML_aligned"),
            Path(self.interim, "Trafo"),
        )
        self.mzML_dir = Path(self.interim, "mzML_aligned")

    def peak_maps_to_df(self):
        df_path = Path(self.mzML_dir).name + "_df"
        reset_directory(Path(self.interim, df_path))
        for file in Path(self.mzML_dir).iterdir():
            DataFrames.mzML_to_ftr(file, Path(self.interim, df_path))

    def map_MS2(self):
        MapID().run(
            self.mzML_dir,
            self.featureXML_dir,
            Path(self.interim, "FeatureMaps_ID_mapped"),
        )
        self.featureXML_dir = Path(self.interim, "FeatureMaps_ID_mapped")

    def link_feature_maps(self):
        FeatureLinker().run(
            self.featureXML_dir,
            str(self.consensusXML),
            {
                "link:mz_tol": self.params["fl_mz_tol"],
                "link:rt_tol": float(self.params["fl_rt_tol"]),
                "mz_unit": self.params["fl_mz_unit"],
            },
        )

    def consensus_df(self):
        DataFrames().create_consensus_table(
            str(self.consensusXML),
            self.consensus_tsv,
            self.sirius_ms_dir,
        )

    def sirius(self):
        Sirius().run(
            self.mzML_dir,
            self.featureXML_dir,
            Path(self.results_dir, "SIRIUS"),
            "",
            True,
            {"-preprocessing:feature_only": "true"},
        )
        self.sirius_ms_dir = Path(self.results_dir, "SIRIUS", "sirius_files")

    def requantify_alternative(self):
        Requantifier().run(
            str(self.consensusXML),
            self.consensus_tsv,
            self.mzML_dir,
            self.featureXML_dir,
            self.params["ffm_mass_error"],
        )

    def requantify(self):
        FeatureMapHelper.FFMID_library_from_consensus_df(
            self.consensus_tsv, Path(self.interim, "FFMID.tsv")
        )
        FeatureFinderMetaboIdent().run(
            self.mzML_dir,
            str(Path(self.interim, "FFMID")),
            Path(self.interim, "FFMID.tsv"),
            {
                "detect:peak_width": self.params["ffm_peak_width"] * 2,
                "extract:mz_window": self.params["ffm_mass_error"],
                "extract:n_isotopes": 2,
            },
        )
        self.consensusXML = Path(
            self.interim, "FeatureMatrixRequantified.consensusXML")
        self.consensus_tsv = Path(
            self.results_dir, "FeatureMatrixRequantified.tsv")
        self.featureXML_dir = Path(self.interim, "FFMID")

    def export_metadata(self):
        GNPSExport().export_metadata_table_only(
            str(self.consensusXML),
            str(Path(self.results_dir, "MetaData.tsv")),
        )

    def gnps(self):
        # create interim directory for MS2 ids later (they are based on the GNPS mgf file)
        reset_directory(Path(self.interim, "mztab_ms2"))
        GNPSExport().run(
            str(self.consensusXML),
            self.mzML_dir,
            Path(self.results_dir, "GNPS"),
        )
        output_mztab = str(Path(self.interim, "mztab_ms2", "GNPS.mzTab"))
        database = "example-data/ms2-libraries/GNPS-LIBRARY.mgf"
        SpectralMatcher().run(database, self.mgf_file, output_mztab)
        DataFrames().annotate_ms2(
            self.mgf_file,
            output_mztab,
            self.consensus_tsv,
            "GNPS library match",
        )

    def annotate_MS1(self):
        if self.params["ms1_annotation_file"]:
            DataFrames().annotate_ms1(
                self.consensus_tsv,
                self.params["ms1_annotation_file"],
                self.params["annotation_mz_window_ppm"],
                self.params["annoation_rt_window_sec"],
            )
            DataFrames().save_MS_ids(
                self.consensus_tsv,
                Path(self.results_dir, "MS1-annotations"),
                "MS1 annotation",
            )

    def annotate_MS2(self):
        if self.params["ms2_annotation_file"]:
            output_mztab = str(Path(self.interim, "mztab_ms2", "MS2.mzTab"))
            SpectralMatcher().run(
                self.params["ms2_annotation_file"], self.mgf_file, output_mztab
            )
            DataFrames().annotate_ms2(
                self.mgf_file,
                output_mztab,
                self.consensus_tsv,
                "MS2 annotation",
                overwrite_name=True,
            )
            DataFrames().save_MS_ids(
                self.consensus_tsv,
                os.path.join(self.results_dir, "MS2-annotations"),
                "MS2 annotation",
            )

    def make_zip_archives(self):
        if self.sirius_ms_dir:
            path = Path(self.sirius_ms_dir)
            if path.exists():
                shutil.make_archive(
                    Path(self.results_dir, "ExportSIRIUS"), "zip", path)
        path = Path(self.results_dir, "GNPS")
        if path.exists():
            shutil.make_archive(
                Path(self.results_dir, "ExportGNPS"),
                "zip",
                path,
            )

    def additional_data_for_consensus_df(self):
        DataFrames.consensus_df_additional_annotations(
            self.consensus_tsv,
            Path(self.results_dir, "FeatureMatrix.ftr"),
            self.consensusXML,
        )


def run_umetaflow(params, mzML_files, results_dir):
    with st.status("Running UmetaFlow", expanded=True) as status:
        umetaflow = UmetaFlow(params, mzML_files, results_dir)

        st.write("Fetching input mzML files...")
        umetaflow.fetch_raw_data()

        st.write("Detecting features...")
        umetaflow.feature_detection()

        # return if no features have been detected
        if not any(Path(umetaflow.interim, "FFM").iterdir()):
            st.error("No features detected in all files")
            return

        if params["remove_blanks"] and len(params["blank_mzML_files"]) > 0:
            st.write("Removing blank features...")
            umetaflow.remove_blanks()
            if not any(umetaflow.featureXML_dir.iterdir()):
                st.warning(
                    "No samples left after blank removal! Blank samples will not be further processed."
                )
                shutil.rmtree(results_dir)
                return

        if params["use_ad"] and not params["use_requant"]:
            st.write("Determining adducts...")
            umetaflow.adduct_detection()

        # annotate only when necessary
        if (
            params["use_sirius_manual"] or params["annotate_ms2"] or params["use_gnps"]
        ) and not params["use_requant"]:
            umetaflow.feature_maps_to_df()
            st.write("Mapping MS2 data to features...")
            umetaflow.map_MS2()


        if params["use_ma"] and len(list(umetaflow.featureXML_dir.iterdir())) > 1:
            st.write("Aligning feature maps...")
            umetaflow.align_feature_maps()
            umetaflow.feature_maps_to_df()

            st.write("Aligning mzML files...")
            umetaflow.align_peak_maps()
            umetaflow.peak_maps_to_df()
        else:
            umetaflow.feature_maps_to_df()
            umetaflow.peak_maps_to_df()

        # export only sirius ms files to use in the GUI tool
        if params["use_sirius_manual"] and not params["use_requant"]:
            st.write("Exporting files for Sirius...")
            umetaflow.sirius()

        st.write("Linking features...")
        umetaflow.link_feature_maps()
        umetaflow.consensus_df()

        if params["use_requant"]:
            # for requant run FFMID and Feature Linking and optionally Adduct Decharging and Mapping MS2 data
            st.write("Re-quantification...")
            umetaflow.requantify()

            # return if no features have been detected
            if not any(Path(umetaflow.interim, "FFMID").iterdir()):
                st.error("No features detected in all files **after requantification**.")
                shutil.rmtree(results_dir)
                return

            st.write("Exporting re-quantified feature maps for visualization...")
            umetaflow.feature_maps_to_df(requant=True)

            # annotate only when necessary
            if params["use_sirius_manual"] or params["annotate_ms2"] or params["use_gnps"]:
                st.write("Mapping MS2 data to features...")
                umetaflow.map_MS2()

            if params["use_ad"]:
                st.write("Determining adducts...")
                umetaflow.adduct_detection()

            if params["use_sirius_manual"]:
                st.write("Exporting files for Sirius...")
                umetaflow.sirius()

            st.write("Linking re-quantified features...")
            umetaflow.link_feature_maps()
            umetaflow.consensus_df()

        # export metadata
        umetaflow.export_metadata()

        if params["use_gnps"] or params["annotate_ms2"]:
            st.write("Exporting files for GNPS...")
            umetaflow.gnps()

        if params["annotate_ms1"]:
            st.write("Annotating features on MS1 level by m/z and RT")
            umetaflow.annotate_MS1()

        if (
            params["annotate_ms2"] and params["ad_ion_mode"] == "positive"
        ):  # fails with negative mode right now due to wrong charge annotation
            st.write("Annotating features on MS2 level by fragmentation patterns...")
            umetaflow.annotate_MS2()

        umetaflow.additional_data_for_consensus_df()

        # Export files
        umetaflow.make_zip_archives()

        status.update(label="UmetaFlow run complete!", state="complete", expanded=False)
        st.rerun()


METABO = {
    "main": """
This workflow includes the core UmetaFlow pipeline which results in a table of metabolic features.

- The most important parameter are marked as **bold** text. Adjust them according to your instrument.
- All the steps with checkboxes are optional.
""",
    "blank_removal": "Useful to filter out features which are present in blank sample/s or e.g. for differential feature detection to remove features which are present in control, but not in treatment samples.",
    "blank_samples": "The selected samples will be used to calculate avarage feature blank intensities and will not be further processed.",
    "blank_cutoff": "Features that have an intensity ratio below (avagera blank) to (average samples) will be removed. Set low for strict blank removal.",
    "potential_adducts": """
Specify adducts and neutral additions/losses.\n
Format (each in a new line): adducts:charge:probability.\n
The summed up probability for all charged entries needs to be 1.0.\n

Good starting options for positive mode with sodium adduct and water loss:\n
H:+:0.9\n
Na:+:0.1\n
H-2O-1:0:0.4

Good starting options for negative mode with water loss and formic acid addition:\n
H-1:-:1\n
H-2O-1:0:0.5\n
CH2O2:0:0.5
    """,
    "ad_ion_mode": "Carefully adjust settings for each mode. Especially potential adducts and negative min/max charges for negative mode.",
    "ad_charge_min": "e.g. for negative mode -3, for positive mode 1",
    "ad_charge_max": "e.g. for negative mode -1, for positive mode 3",
    "ad_rt_diff": "Groups features with slightly different RT.",
    "re-quant": "Go back into the raw data to re-quantify consensus features that have missing values.",
    "sirius": "Export files for formula and structure predictions. Run Sirius with these pre-processed .ms files, can be found in results -> SIRIUS -> sirius_files.",
    "gnps": "Run GNPS Feature Based Molecular Networking and Ion Identity Molecular Networking with these files, can be found in results -> GNPS.",
    "annotate_gnps": "UmetaFlow contains the complete GNPS library in mgf file format. Check to annotate.",
    "annotate_ms1": "Annotate features on MS1 level with known m/z and retention times values.",
    "annotate_ms1_rt": "Checks around peak apex, e.g. window of 60 s will check left and right 30 s.",
    "annotate_ms2": "Annotate features on MS2 level based on their fragmentation patterns. The library has to be in mgf file format.",
}

HELP = """
This page lets you analyze mass spectrometry data for untargeted metabolomics.

Below are the main steps and options available in the tool:

#### 1. Pre-Processing

This step focuses on feature detection and data filtering.

**Feature Detection**

The most important parameters should be set according to your data:

*Mass Error ppm*: Specify the allowed mass error in parts per million (ppm).

*Noise Threshold*: Set the noise threshold to filter out low-intensity signals.

**Blank Removal**

You can optionally perform blank removal to filter out features that are present in blank samples. Set the blank/sample intensity ratio cutoff for filtering.

**Map Alignment**

If needed, you can perform map alignment to correct retention time differences between samples.

**Adduct Detection**

Optionally, you can perform adduct detection to identify potential ionization adducts and charge states for features.

#### 2. Re-Quantification
In this step, you have the option to re-quantify consensus features with missing values in the raw data.

#### 3. Export Files for SIRIUS and GNPS
You can export files for SIRIUS, which is a tool for formula and structure predictions. Additionally, you can export files for GNPS, a feature-based molecular networking and ion identity molecular networking tool. Optionally, you can annotate features using the GNPS library.

#### 4. Annotation via In-House Library
You can perform MS1 annotation by m/z and retention time. Select a library in TSV format for annotation, and specify the retention time window and m/z window for matching. Optionally, you can perform MS2 annotation via fragmentation patterns. Select an MGF format library for MS2 annotation.

Please use the checkboxes and input fields to customize the parameters according to your specific dataset and analysis needs. Remember to check the information and help provided for each option to understand its impact on your data analysis.
"""
