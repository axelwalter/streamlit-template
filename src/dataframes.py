from pyopenms import *
import pandas as pd
import numpy as np
import os
from pathlib import Path
from .common.common import reset_directory
from pyteomics import mztab, mgf


class DataFrames:
    def create_consensus_table(self, consensusXML_file, table_file, sirius_ms_dir=""):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        df = consensus_map.get_df().drop(["sequence"], axis=1)
        for cf in consensus_map:
            if cf.metaValueExists("best ion"):
                df.insert(
                    4, "adduct", [cf.getMetaValue(
                        "best ion") for cf in consensus_map]
                )
                break
        for cf in consensus_map:
            if cf.metaValueExists("label"):
                df["name"] = [cf.getMetaValue("label") for cf in consensus_map]
                break
        if "adduct" in df.columns:
            df.insert(
                0,
                "metabolite",
                [
                    f"{round(mz, 4)}@{round(rt, 2)}@{adduct}"
                    for mz, rt, adduct in zip(
                        df["mz"].tolist(), df["RT"].tolist(), df["adduct"].tolist()
                    )
                ],
            )
        else:
            df.insert(
                0,
                "metabolite",
                [
                    f"{round(mz, 4)}@{round(rt, 2)}"
                    for mz, rt in zip(df["mz"].tolist(), df["RT"].tolist())
                ],
            )
        not_sample = [
            c
            for c in df.columns
            if c not in ["mz", "RT", "charge", "adduct", "name", "quality"]
        ]
        df[not_sample] = df[not_sample].applymap(
            lambda x: int(round(x, 0)) if isinstance(x, (int, float)) else x
        )
        # annotate original feature Ids which are in the Sirius .ms files
        if sirius_ms_dir:
            print(sirius_ms_dir)
            ms_files = [Path(sirius_ms_dir, file)
                        for file in os.listdir(sirius_ms_dir)]
            map = {
                Path(value.filename).stem: key
                for key, value in consensus_map.getColumnHeaders().items()
            }
            for file in ms_files:
                if file.exists():
                    key = map[file.stem]
                    id_list = []
                    content = file.read_text()
                    for cf in consensus_map:
                        # get a map with map index and feature id for each consensus feature -> get the feature id key exists
                        f_map = {
                            fh.getMapIndex(): fh.getUniqueId()
                            for fh in cf.getFeatureList()
                        }
                        if key in f_map.keys():
                            f_id = str(f_map[key])
                        else:
                            f_id = ""
                        if f_id and f_id in content:
                            id_list.append(f_id)
                        else:
                            id_list.append("")
                    df[file.stem + "_SiriusID"] = id_list
        df_mzML = df[[col for col in df.columns if col.endswith("mzML")]]
        df_mzML = df_mzML.reindex(sorted(df_mzML.columns), axis=1)
        df = pd.concat(
            [df[[col for col in df.columns if not col.endswith(
                "mzML")]], df_mzML],
            axis=1,
        )

        if "second-feature-map.mzML" in df.columns:
            df = df.drop(columns=["second-feature-map.mzML"])
        if str(table_file).endswith("tsv"):
            df.to_csv(table_file, sep="\t")
        elif str(table_file).endswith("ftr"):
            df.reset_index().to_feather(table_file)
        return df

    def FFMID_chroms_to_df(self, featureXML_file, table_file, time_unit="seconds"):
        time_factor = 1
        if time_unit == "minutes":
            time_factor = 60
        fm = FeatureMap()
        FeatureXMLFile().load(featureXML_file, fm)
        chroms = {}
        for f in fm:
            for i, sub in enumerate(f.getSubordinates()):
                name = f.getMetaValue("label") + "_" + str(i + 1)
                chroms[name + "_int"] = [
                    int(y[1]) for y in sub.getConvexHulls()[0].getHullPoints()
                ]
                chroms[name + "_RT"] = [
                    x[0] / time_factor for x in sub.getConvexHulls()[0].getHullPoints()
                ]
        df = pd.DataFrame({key: pd.Series(value)
                          for key, value in chroms.items()})
        if table_file.endswith("tsv"):
            df.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

    def FFMID_auc_to_df(self, featureXML_file, table_file):
        fm = FeatureMap()
        FeatureXMLFile().load(featureXML_file, fm)
        aucs = {}
        for f in fm:
            aucs[f.getMetaValue("label")] = [int(f.getIntensity())]
        df = pd.DataFrame({key: pd.Series(value)
                          for key, value in aucs.items()})
        if table_file.endswith("tsv"):
            df.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

    def FFMID_auc_combined_to_df(self, df_auc_file, table_file):
        if df_auc_file.endswith("tsv"):
            df = pd.read_csv(df_auc_file, sep="\t")
        elif df_auc_file.endswith("ftr"):
            df = pd.read_feather(df_auc_file).drop(columns=["index"])
        aucs_condensed = {}
        for a in set([c.split("#")[0] for c in df.columns]):
            aucs_condensed[a] = 0
            for b in [
                b for b in df.columns if ((a + "#" in b and b.startswith(a)) or a == b)
            ]:
                aucs_condensed[a] += df[b][0]
        df_combined = pd.DataFrame(
            {key: pd.Series(value) for key, value in aucs_condensed.items()}
        )
        if table_file.endswith("tsv"):
            df_combined.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df_combined.reset_index().to_feather(table_file)

    def get_auc_summary(self, df_files, table_file):
        # get a list of auc dataframe file paths (df_files), combine them into a summary (consensus) df
        dfs = []
        indeces = []
        empty = []
        for file in df_files:
            if file.endswith("tsv"):
                df = pd.read_csv(file, sep="\t")
            elif file.endswith("ftr"):
                df = pd.read_feather(file)
                if "index" in df.columns:
                    df.index = df["index"]
                    df = df.drop(columns=["index"])
            sample_name = os.path.basename(file)[:-4].split("AUC")[0]
            if df.empty:
                empty.append(sample_name)
            else:
                dfs.append(df)
                indeces.append(sample_name)
        df = pd.concat(dfs)
        df = df.set_index(pd.Series(indeces))
        df = df.transpose()
        df = df.fillna(0)
        for sample in empty:
            df[sample] = 0
        df = df.applymap(
            lambda x: int(round(x, 0)) if isinstance(x, (int, float)) else x
        )
        df.sort_index(axis=1, inplace=True)
        df.insert(0, "metabolite", df.index)
        df = df.set_index(pd.Series(range(1, len(df) + 1)))
        if table_file.endswith("tsv"):
            df.to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

    def annotate_ms1(self, df_file, library_file, mz_window, rt_window):
        df = pd.read_csv(df_file, sep="\t")
        library = pd.read_csv(library_file, sep="\t")
        df.insert(2, "MS1 annotation", "")

        df["mz"] = df["mz"].astype(float)

        for _, std in library.iterrows():
            delta_Da = abs(mz_window * std["mz"] / 1000000)
            mz_lower = std["mz"] - delta_Da
            mz_upper = std["mz"] + delta_Da
            rt_lower = std["RT"] - rt_window / 2
            rt_upper = std["RT"] + rt_window / 2
            match = df.query(
                "mz > @mz_lower and mz < @mz_upper and RT > @rt_lower and RT < @rt_upper"
            )
            if not match.empty:
                for _, row in match.iterrows():
                    if len(df.loc[df["id"] == row["id"], "MS1 annotation"]) > 1:
                        df.loc[df["id"] == row["id"], "MS1 annotation"] += (
                            ";" + std["name"]
                        )
                    else:
                        df.loc[df["id"] == row["id"],
                               "MS1 annotation"] += std["name"]

        # replace generic metabolite name with actual MS1 annotation
        metabolites = []
        for x, y in zip(df["metabolite"], df["MS1 annotation"]):
            if y and y not in metabolites:
                metabolites.append(y)
            elif y and y in metabolites:
                metabolites.append(y + f"_{metabolites.count(y)}")
            else:
                metabolites.append(x)
        df["metabolite"] = metabolites
        df.to_csv(df_file, sep="\t", index=False)

    def save_MS_ids(self, df_file, ms_dir, column_name):
        df = pd.read_csv(df_file, sep="\t")
        if column_name not in df.columns:
            return
        df = df[df[column_name].notna()]
        filename = column_name.replace(" ", "-") + ".tsv"
        if not df.empty:
            reset_directory(ms_dir)
            df.to_csv(os.path.join(ms_dir, filename), sep="\t", index=False)

    def annotate_ms2(
        self,
        mgf_file,
        output_mztab,
        feature_matrix_df_file,
        match_column_name,
        overwrite_name=False,
    ):
        # clean up the mzTab to a dataframe:
        matches = mztab.MzTab(
            output_mztab, encoding="UTF8", table_format="df"
        ).small_molecule_table
        if matches.empty:
            return
        matches = matches.query(
            "opt_ppm_error <= 10 and opt_ppm_error >= -10 and opt_match_score >= 60"
        )
        matches["opt_spec_native_id"] = matches["opt_spec_native_id"].str.replace(
            r"index=", ""
        )

        # load spectra parameters from mgf_file
        exp = mgf.MGF(
            source=mgf_file,
            use_header=True,
            convert_arrays=2,
            read_charges=True,
            read_ions=False,
            dtype=None,
            encoding=None,
        )
        parameters = []
        for spec in exp:
            parameters.append(spec["params"])
        spectra = pd.DataFrame(parameters)
        spectra["feature_id"] = spectra["feature_id"].str.replace(r"e_", "")

        # annotate matches with feature id, based on scan number
        matches.insert(
            0,
            "feature_id",
            [
                spectra[spectra["scans"].astype(int) == int(scan) + 1][
                    "feature_id"
                ].tolist()[0]
                for scan in matches["opt_spec_native_id"]
            ],
        )

        # annotate features with hits based on feature id
        features = pd.read_csv(feature_matrix_df_file, sep="\t")
        ms2_ids = []
        for f_id in features.id:
            hit = matches[matches["feature_id"] == str(f_id)]
            if hit.shape[0] > 0:
                ms2_ids.append(";".join(hit["description"]))
            else:
                ms2_ids.append("")
        features.insert(2, match_column_name, ms2_ids)

        if overwrite_name:
            # replace generic or MS1 annotated metabolite name with actual MS2 annotation
            metabolites = []
            for x, y in zip(features["metabolite"], features[match_column_name]):
                if y and y not in metabolites:
                    metabolites.append(y)
                elif y and y in metabolites:
                    metabolites.append(y + f"_{metabolites.count(y)}")
                else:
                    metabolites.append(x)
            features["metabolite"] = metabolites
        features.to_csv(feature_matrix_df_file, sep="\t", index=False)

    def mzML_to_ftr(mzML_file_path, ftr_dir):
        exp = MSExperiment()
        MzMLFile().load(str(mzML_file_path), exp)
        df = exp.get_df()
        df.insert(0, "mslevel", [spec.getMSLevel() for spec in exp])
        df.insert(
            0,
            "precursormz",
            [
                spec.getPrecursors()[0].getMZ() if spec.getPrecursors() else 0
                for spec in exp
            ],
        )
        df.to_feather(Path(ftr_dir, mzML_file_path.stem + ".ftr"))

    def featureXML_to_ftr(featureXML_file_path, ftr_dir, requant=False):
        fm = FeatureMap()
        FeatureXMLFile().load(str(featureXML_file_path), fm)
        df = fm.get_df(export_peptide_identifications=False).reset_index()
        df["adduct"] = [f.getMetaValue("dc_charge_adducts") for f in fm]
        df["original_rt"] = [f.getMetaValue("original_RT") for f in fm]

        if requant:
            df["chrom_rts"] = [
                np.array(
                    [
                        x[0]
                        for x in f.getSubordinates()[0]
                        .getConvexHulls()[0]
                        .getHullPoints()
                    ]
                ).astype(np.float64)
                for f in fm
            ]
            df["chrom_intensities"] = [
                np.array(
                    [
                        x[1]
                        for x in f.getSubordinates()[0]
                        .getConvexHulls()[0]
                        .getHullPoints()
                    ]
                ).astype(np.float64)
                for f in fm
            ]
            df["fwhm"] = [f.getMetaValue("model_FWHM") for f in fm]
            df["snr"] = [f.getMetaValue("sn_ratio") for f in fm]

        else:
            df["chrom_rts"] = [
                np.array(f.getMetaValue("chrom_rts").split(
                    ",")).astype(np.float64)
                for f in fm
            ]
            df["chrom_intensities"] = [
                np.array(f.getMetaValue("chrom_intensities").split(",")).astype(
                    np.float64
                )
                for f in fm
            ]
            df["fwhm"] = [f.getMetaValue("FWHM") for f in fm]

        df["metabolite"] = (
            df["mz"].round(5).astype(str) + "@" + df["RT"].round(2).astype(str)
        )
        df.to_feather(Path(ftr_dir, featureXML_file_path.stem + ".ftr"))

    def consensus_df_additional_annotations(
        tsv_file_path, ftr_file_path, consensusXML_file_path
    ):
        df = pd.read_csv(tsv_file_path, sep="\t")

        # load ConsensusMap to extract additional info that should go into the ftr file only
        cm = ConsensusMap()
        ConsensusXMLFile().load(str(consensusXML_file_path), cm)

        # to map feature map file names to feature ids create a map
        map = {
            key: Path(value.filename).name
            for key, value in cm.getColumnHeaders().items()
        }

        # for each cf annotate a dict with key=sample, value=feature unique ID
        df["sample_to_id"] = [
            {
                map[f.getMapIndex()][:-5]: str(f.getUniqueId())
                for f in cf.getFeatureList()
            }
            for cf in cm
        ]

        df.to_feather(ftr_file_path)
