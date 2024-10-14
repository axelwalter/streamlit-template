import os
import csv
from pathlib import Path
import pandas as pd
import numpy as np
from pyopenms import *
from .common.common import reset_directory


class FeatureMapHelper:
    def load_feature_maps(path: str):
        feature_maps = []
        for f in Path(path).glob("*.featureXML"):
            feature_map = FeatureMap()
            FeatureXMLFile().load(str(f), feature_map)
            feature_maps.append(feature_map)
        return feature_maps

    def split_consensus_map(
        self,
        consensusXML_file,
        consensusXML_complete_file="",
        consensusXML_missing_file="",
    ):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        headers = consensus_map.getColumnHeaders()
        complete = ConsensusMap(consensus_map)
        complete.clear(False)
        missing = ConsensusMap(consensus_map)
        missing.clear(False)
        for cf in consensus_map:
            if len(cf.getFeatureList()) < len(headers):
                missing.push_back(cf)
            else:
                complete.push_back(cf)
        if consensusXML_complete_file.endswith("consensusXML"):
            ConsensusXMLFile().store(consensusXML_complete_file, complete)
        if consensusXML_missing_file.endswith("consensusXML"):
            ConsensusXMLFile().store(consensusXML_missing_file, missing)

    # takes a (filtered) ConsensusMap (usually the complete) and reconstructs the FeatureMaps
    def consensus_to_feature_maps(
        self, consensusXML_file, original_featureXML_dir, reconstructed_featureXML_dir
    ):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        to_keep_ids = [
            item
            for sublist in [
                [feature.getUniqueId() for feature in cf.getFeatureList()]
                for cf in consensus_map
            ]
            for item in sublist
        ]
        reset_directory(reconstructed_featureXML_dir)
        for file in os.listdir(original_featureXML_dir):
            feature_map = FeatureMap()
            FeatureXMLFile().load(
                os.path.join(original_featureXML_dir, file), feature_map
            )
            fm_filterd = FeatureMap(feature_map)
            fm_filterd.clear(False)
            for feature in feature_map:
                if feature.getUniqueId() in to_keep_ids:
                    fm_filterd.push_back(feature)
            FeatureXMLFile().store(
                os.path.join(
                    reconstructed_featureXML_dir,
                    os.path.basename(
                        fm_filterd.getMetaValue("spectra_data")[0].decode()
                    )[:-4]
                    + "featureXML",
                ),
                fm_filterd,
            )

    def merge_feature_maps(
        self, featureXML_merged_dir, featureXML_dir_a, featureXML_dir_b
    ):
        reset_directory(featureXML_merged_dir)
        for file_ffm in os.listdir(featureXML_dir_a):
            for file_ffmid in os.listdir(featureXML_dir_b):
                if file_ffm == file_ffmid:
                    fm_ffm = FeatureMap()
                    FeatureXMLFile().load(
                        os.path.join(featureXML_dir_a, file_ffm), fm_ffm
                    )
                    fm_ffmid = FeatureMap()
                    FeatureXMLFile().load(
                        os.path.join(featureXML_dir_b, file_ffm), fm_ffmid
                    )
                    for f in fm_ffmid:
                        fm_ffm.push_back(f)
                    FeatureXMLFile().store(
                        os.path.join(featureXML_merged_dir, file_ffm), fm_ffm
                    )

    # create a library for each FeatureXML file with features that are missing in this particular sample
    # in the given ConsensusXMLFile
    def FFMID_libraries_for_missing_features(self, consensusXML_file, library_dir):
        reset_directory(library_dir)
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        headers = {
            i: header.filename
            for i, header in enumerate(consensus_map.getColumnHeaders().values())
        }
        # write header in all library files here
        all_lib_files = [
            os.path.join(library_dir, headers[i][:-5] + ".tsv") for i in headers.keys()
        ]
        for file in all_lib_files:
            with open(file, "a", newline="") as tsvfile:
                writer = csv.writer(tsvfile, delimiter="\t")
                writer.writerow(
                    [
                        "CompoundName",
                        "SumFormula",
                        "Mass",
                        "Charge",
                        "RetentionTime",
                        "RetentionTimeRange",
                        "IsoDistribution",
                    ]
                )
        for i, cf in enumerate(consensus_map):
            lib_files = [
                file
                for i, file in enumerate(all_lib_files)
                if i not in [f.getMapIndex() for f in cf.getFeatureList()]
            ]
            if lib_files:
                row = [
                    "f_" + str(i),
                    "",
                    cf.getMZ() * cf.getCharge() - (cf.getCharge() * 1.007825),
                    cf.getCharge(),
                    cf.getRT(),
                    0,
                    0,
                ]
                for file in lib_files:
                    with open(file, "a", newline="") as tsvfile:
                        writer = csv.writer(tsvfile, delimiter="\t")
                        writer.writerow(row)

    def FFMID_library_from_consensus_df(consensus_df_path: str, library_path: str):
        df = pd.read_csv(consensus_df_path, sep="\t")
        lib = pd.DataFrame(
            {
                "CompoundName": df["metabolite"],
                "SumFormula": np.full(df.shape[0], ""),
                "Mass": df["mz"] * df["charge"] - (df["charge"] * 1.007825),
                "Charge": df["charge"],
                "RetentionTime": df["RT"],
                "RetentionTimeRange": np.full(df.shape[0], 0),
                "IsoDistribution": np.full(df.shape[0], 0),
            },
            index=range(df.shape[0]),
        )
        lib.to_csv(library_path, sep="\t", index=False)


class FeatureFinderMetabo:
    def run(self, mzML, featureXML, params={}):
        print(params)
        if type(mzML) is list:
            mzML_files = mzML
        elif os.path.isdir(mzML):
            mzML_files = [os.path.join(mzML, file)
                          for file in os.listdir(mzML)]
        else:
            mzML_files = [mzML]
        if not featureXML.endswith(".featureXML"):
            reset_directory(featureXML)

        for mzML_file in mzML_files:
            exp = MSExperiment()
            MzMLFile().load(mzML_file, exp)
            exp.sortSpectra(True)

            mtd = MassTraceDetection()
            mtd_par = mtd.getDefaults()
            for key, value in params.items():
                if key.encode() in mtd_par.keys():
                    mtd_par.setValue(key, value)
            mtd.setParameters(mtd_par)
            mass_traces = []
            # input MSExperiment, empty list for detected mass traces, max_size (if not 0, sets the maximum number of mass traces)
            mtd.run(exp, mass_traces, 0)

            epd = ElutionPeakDetection()
            epd_par = epd.getDefaults()
            for key, value in params.items():
                if key.encode() in epd_par.keys():
                    epd_par.setValue(key, value)
            epd.setParameters(epd_par)
            elution_peaks = []
            # list with all detected mass traces, list of mass traces that represent an elution peak
            epd.detectPeaks(mass_traces, elution_peaks)
            ffm = FeatureFindingMetabo()
            ffm_par = ffm.getDefaults()
            for key, value in params.items():
                if key.encode() in ffm_par.keys():
                    ffm_par.setValue(key, value)
            ffm.setParameters(ffm_par)
            feature_map = FeatureMap()
            chromatograms = []
            # elution peaks, empty FeatureMap, empty list for feature chromatograms
            ffm.run(elution_peaks, feature_map, chromatograms)

            feature_map.setPrimaryMSRunPath([mzML_file.encode()])

            feature_map_chroms = FeatureMap(feature_map)
            feature_map_chroms.clear(False)

            if feature_map.size() == len(chromatograms):
                for f in feature_map:
                    # get the matching monoisotopic chromatogram
                    match = [
                        c[0]
                        for c in chromatograms
                        if int(c[0].getName().split("_")[0]) == f.getUniqueId()
                    ][0]
                    rts, ints = match.get_peaks()
                    f.setMetaValue("chrom_rts", ",".join(rts.astype(str)))
                    f.setMetaValue("chrom_intensities",
                                   ",".join(ints.astype(str)))
                    feature_map_chroms.push_back(f)

            feature_map_chroms.setUniqueIds()

            if feature_map.size() > 0:

                if os.path.isdir(featureXML):
                    FeatureXMLFile().store(
                        os.path.join(
                            featureXML, os.path.basename(
                                mzML_file)[:-4] + "featureXML"
                        ),
                        feature_map_chroms,
                    )
                else:
                    FeatureXMLFile().store(featureXML, feature_map_chroms)


class MapAligner:
    def run(self, input_files, aligned_dir, trafo_dir, params={}):
        aligner = MapAlignmentAlgorithmPoseClustering()
        aligner_par = aligner.getDefaults()
        for key, value in params.items():
            if key.encode() in aligner_par.keys():
                aligner_par.setValue(key, value)
        aligner.setParameters(aligner_par)
        if type(input_files) is list:
            inputs = input_files
        else:
            inputs = os.listdir(input_files)
        if inputs and inputs[0].endswith("featureXML"):
            reset_directory(aligned_dir)
            reset_directory(trafo_dir)
            feature_maps = FeatureMapHelper.load_feature_maps(input_files)
            # store TransformationDescriptions for MapAlignmentTransformer of MSExperiments during Requantification
            transformations = {}
            ref_index = feature_maps.index(
                sorted(feature_maps, key=lambda x: x.size())[-1]
            )
            aligner.setReference(feature_maps[ref_index])
            print(
                "Map Alignment reference map: ",
                feature_maps[ref_index].getMetaValue("spectra_data")[
                    0].decode(),
            )

            for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1:]:
                trafo = TransformationDescription()
                try:
                    # store information on aligmentment in TransformationDescription, RTs in FeatureMap not modified at this point
                    aligner.align(feature_map, trafo)
                    transformer = MapAlignmentTransformer()
                    # FeatureMap, TransformationDescription, bool: keep original RTs as meta value
                    transformer.transformRetentionTimes(
                        feature_map, trafo, True)
                    transformations[
                        feature_map.getMetaValue("spectra_data")[0].decode()
                    ] = trafo
                    TransformationXMLFile().store(
                        os.path.join(
                            trafo_dir,
                            os.path.basename(
                                feature_map.getMetaValue("spectra_data")[
                                    0].decode()
                            )[:-4]
                            + "trafoXML",
                        ),
                        trafo,
                    )
                except RuntimeError:  # WARNING: your map likely has a scaling around inf but your parameters only allow for a maximal scaling of 2
                    pass

            for feature_map in feature_maps:
                print(feature_map.size())
                FeatureXMLFile().store(
                    os.path.join(
                        aligned_dir,
                        os.path.basename(
                            feature_map.getMetaValue("spectra_data")[
                                0].decode()
                        )[:-4]
                        + "featureXML",
                    ),
                    feature_map,
                )

        elif inputs and inputs[0].endswith("mzML") and type(input_files) is not list:
            reset_directory(aligned_dir)
            for file in inputs:
                exp = MSExperiment()
                MzMLFile().load(os.path.join(input_files, file), exp)
                exp.sortSpectra(True)
                if file[:-4] + "trafoXML" not in os.listdir(trafo_dir):
                    MzMLFile().store(os.path.join(aligned_dir, file), exp)
                    continue
                transformer = MapAlignmentTransformer()
                trafo_description = TransformationDescription()
                TransformationXMLFile().load(
                    os.path.join(trafo_dir, file[:-4] + "trafoXML"),
                    trafo_description,
                    False,
                )
                transformer.transformRetentionTimes(
                    exp, trafo_description, True)
                MzMLFile().store(os.path.join(aligned_dir, file), exp)

        elif inputs and inputs[0].endswith("mzML") and type(input_files) is list:
            reset_directory(aligned_dir)
            for file in inputs:
                exp = MSExperiment()
                MzMLFile().load(file, exp)
                exp.sortSpectra(True)
                if file[:-4] + "trafoXML" not in os.listdir(trafo_dir):
                    MzMLFile().store(os.path.join(aligned_dir, Path(file).name), exp)
                    continue
                transformer = MapAlignmentTransformer()
                trafo_description = TransformationDescription()
                TransformationXMLFile().load(
                    os.path.join(trafo_dir, file[:-4] + "trafoXML"),
                    trafo_description,
                    False,
                )
                transformer.transformRetentionTimes(
                    exp, trafo_description, True)
                MzMLFile().store(os.path.join(aligned_dir, Path(file).name), exp)


class FeatureLinker:
    def run(self, featureXML_dir, consensusXML_file, params={}):
        feature_maps = FeatureMapHelper.load_feature_maps(featureXML_dir)
        feature_grouper = FeatureGroupingAlgorithmKD()

        feature_grouper_params = feature_grouper.getDefaults()
        for key, value in params.items():
            if key.encode() in feature_grouper_params.keys():
                feature_grouper_params.setValue(key, value)
        feature_grouper.setParameters(feature_grouper_params)

        consensus_map = ConsensusMap()
        file_descriptions = consensus_map.getColumnHeaders()

        if len(feature_maps) == 1:
            fm = FeatureMap()
            fm.setMetaValue("spectra_data", [b"second-feature-map.mzML"])
            feature_maps.append(fm)
        for i, feature_map in enumerate(feature_maps):
            file_description = file_descriptions.get(i, ColumnHeader())
            file_description.filename = os.path.basename(
                feature_map.getMetaValue("spectra_data")[0].decode()
            )
            file_description.size = feature_map.size()
            file_descriptions[i] = file_description

        feature_grouper.group(feature_maps, consensus_map)
        consensus_map.setColumnHeaders(file_descriptions)

        consensus_map.setUniqueIds()
        ConsensusXMLFile().store(consensusXML_file, consensus_map)
        print(f"ConsensusMap size: {consensus_map.size()}")


class FeatureFinderMetaboIdent:
    def load_library(self, input_file, library_file=""):
        # # input file can be a consensusXML or tsv file
        # if input_file.endswith("consensusXML"):
        #     consensus_map = ConsensusMap()
        #     ConsensusXMLFile().load(input_file, consensus_map)
        #     # Import the consensus tsv table and keep only the columns: RT, mz and charge
        #     library = consensus_map.get_df()[["RT", "mz", "charge"]]
        #     # convert the mz and RT columns to floats and charge to integer for calculations
        #     library["charge"] = pd.to_numeric(library["charge"], downcast="integer")
        #     library["mz"] = pd.to_numeric(library["mz"], downcast="float")
        #     library["RT"] = pd.to_numeric(library["RT"], downcast="float")
        #     library = library.rename(
        #         columns={"RT": "RetentionTime", "charge": "Charge"}
        #     )
        #     # Add a columns named "Mass" and calculate the neutral masses from the charge and mz:
        #     library["Mass"] = 0.0
        #     for i in library.index:
        #         library.at[i, "Mass"] = (
        #             library.loc[i, "mz"] * library.loc[i, "Charge"]
        #         ) - (library.loc[i, "Charge"] * 1.007825)
        #     # drop the mz column
        #     library = library.drop(columns="mz")
        #     library["Charge"] = [[c] for c in library["Charge"]]
        #     library["RetentionTime"] = [[rt] for rt in library["RetentionTime"]]
        #     # add the rest of the columns required for the MetaboliteIdentificationTable and fill with zeros or blanks, except the "Compound Name"
        #     # which, since they are all unknown, can be filled with f_#
        #     library["CompoundName"] = [i for i in range(0, len(library))]
        #     library["CompoundName"] = "f_" + library["CompoundName"].astype(str)
        #     library["SumFormula"] = ""
        #     library["RetentionTimeRange"] = [[0.0] for _ in range(len(library.index))]
        #     library["IsoDistribution"] = [[0.0] for _ in range(len(library.index))]
        #     library = library[
        #         [
        #             "CompoundName",
        #             "SumFormula",
        #             "Mass",
        #             "Charge",
        #             "RetentionTime",
        #             "RetentionTimeRange",
        #             "IsoDistribution",
        #         ]
        #     ]
        # if library_file:
        #     library.to_csv(library_file, sep="\t")
        # metabo_table = []
        # for _, row in library.iterrows():
        #     metabo_table.append(
        #         FeatureFinderMetaboIdentCompound(
        #             row["CompoundName"],
        #             row["SumFormula"],
        #             row["Mass"],
        #             row["Charge"],
        #             row["RetentionTime"],
        #             row["RetentionTimeRange"],
        #             row["IsoDistribution"],
        #         )
        #     )
        # return metabo_table

        # elif input_file.endswith("tsv"):
        metabo_table = []
        with open(input_file, "r") as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            next(tsv_reader)  # skip header
            for row in tsv_reader:
                metabo_table.append(
                    FeatureFinderMetaboIdentCompound(
                        row[0],  # name
                        row[1],  # sum formula
                        float(row[2]),  # mass
                        [int(charge)
                         for charge in row[3].split(",")],  # charges
                        [float(rt) for rt in row[4].split(",")],  # RTs
                        [
                            float(rt_range) for rt_range in row[5].split(",")
                        ],  # RT ranges
                        # isotope distributions
                        [float(iso_distrib)
                         for iso_distrib in row[6].split(",")],
                    )
                )
        return metabo_table

    def run(self, mzML, featureXML, library, params={}):
        if os.path.isdir(mzML):
            mzML_files = [os.path.join(mzML, file)
                          for file in os.listdir(mzML)]
        else:
            mzML_files = [mzML]
        if not featureXML.endswith(".featureXML"):  # -> it is a directory
            reset_directory(featureXML)

        lib_is_dir = False
        if Path(library).is_dir():
            lib_is_dir = True
        elif Path(library).is_file():
            metabo_table = self.load_library(library)
        for mzML_file in mzML_files:
            if Path(library).is_dir():
                metabo_table = self.load_library(
                    Path(library, Path(mzML_file).stem + ".tsv")
                )
            exp = MSExperiment()
            MzMLFile().load(mzML_file, exp)
            ffmid = FeatureFinderAlgorithmMetaboIdent()
            ffmid.setMSData(exp)
            feature_map = FeatureMap()
            ffmid_params = ffmid.getDefaults()
            for key, value in params.items():
                if key.encode() in ffmid_params.keys():
                    ffmid_params.setValue(key, value)
            ffmid.setParameters(ffmid_params)
            if lib_is_dir:
                metabo_table = self.load_library(
                    os.path.join(library, os.path.basename(
                        mzML_file)[:-4] + "tsv")
                )
            # run the FeatureFinderMetaboIdent with the metabo_table and aligned mzML file path
            ffmid.run(metabo_table, feature_map, mzML_file)
            print(feature_map.size())
            feature_map.setUnassignedPeptideIdentifications([])
            feature_map.setProteinIdentifications([])
            # set number of mass traces (for SIRIUS)
            fm_include_mass_traces = FeatureMap(feature_map)
            fm_include_mass_traces.clear(False)

            if feature_map.size() > 0:
                for feature in feature_map:
                    feature.setMetaValue(
                        "num_of_masstraces", ffmid_params[b"extract:n_isotopes"]
                    )
                    fm_include_mass_traces.push_back(feature)
                feature_map = fm_include_mass_traces
                if os.path.isdir(featureXML):
                    FeatureXMLFile().store(
                        os.path.join(
                            featureXML, os.path.basename(
                                mzML_file)[:-4] + "featureXML"
                        ),
                        feature_map,
                    )
                else:
                    FeatureXMLFile().store(featureXML, feature_map)


class MetaboliteAdductDecharger:
    def run(self, fm_dir, fm_decharged_dir, params={}):
        reset_directory(fm_decharged_dir)
        for file in os.listdir(fm_dir):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(fm_dir, file), feature_map)
            mfd = MetaboliteFeatureDeconvolution()
            mdf_par = mfd.getDefaults()
            for key, value in params.items():
                if key.encode() in mdf_par.keys():
                    mdf_par.setValue(key, value)

            mfd.setParameters(mdf_par)

            feature_map_decharged = FeatureMap()
            mfd.compute(
                feature_map, feature_map_decharged, ConsensusMap(), ConsensusMap()
            )
            FeatureXMLFile().store(
                os.path.join(fm_decharged_dir, file), feature_map_decharged
            )


class MapID:
    def run(self, mzML_dir, fm_dir, fm_mapped_dir):
        reset_directory(fm_mapped_dir)
        use_centroid_rt = False
        use_centroid_mz = True
        mapper = IDMapper()
        for mzML_file in os.listdir(mzML_dir):
            exp = MSExperiment()
            MzMLFile().load(os.path.join(mzML_dir, mzML_file), exp)
            for feature_file in os.listdir(fm_dir):
                fm = FeatureMap()
                FeatureXMLFile().load(os.path.join(fm_dir, feature_file), fm)
                if feature_file[:-10] == mzML_file[:-4]:
                    peptide_ids = []
                    protein_ids = []
                    mapper.annotate(
                        fm,
                        peptide_ids,
                        protein_ids,
                        use_centroid_rt,
                        use_centroid_mz,
                        exp,
                    )
                    FeatureXMLFile().store(
                        os.path.join(fm_mapped_dir, feature_file), fm
                    )


class PrecursorCorrector:
    def to_highest_intensity(self, mzML_dir, mzML_corrected_dir):
        mzML_files = os.listdir(mzML_dir)
        reset_directory(mzML_corrected_dir)
        for filename in mzML_files:
            exp = MSExperiment()
            MzMLFile().load(os.path.join(mzML_dir, filename), exp)
            exp.sortSpectra(True)
            delta_mzs = []
            mzs = []
            rts = []
            PrecursorCorrection.correctToHighestIntensityMS1Peak(
                exp, 100.0, True, delta_mzs, mzs, rts
            )
            MzMLFile().store(os.path.join(mzML_corrected_dir, filename), exp)

    def to_nearest_feature(self, mzML_dir, mzML_corrected_dir, featureXML_dir):
        reset_directory(mzML_corrected_dir)
        mzML_files = os.listdir(mzML_dir)
        feature_files = os.listdir(featureXML_dir)
        for mzml_file in mzML_files:
            exp = MSExperiment()
            MzMLFile().load(os.path.join(mzML_dir, mzml_file), exp)
            exp.sortSpectra(True)
            correct = PrecursorCorrection()
            for feature_file in feature_files:
                feature_map_MFD = FeatureMap()
                FeatureXMLFile().load(
                    os.path.join(featureXML_dir, feature_file), feature_map_MFD
                )
                if (
                    os.path.basename(mzml_file)[:-5]
                    == os.path.basename(feature_file)[:-11]
                ):
                    correct.correctToNearestFeature(
                        feature_map_MFD,
                        exp,
                        0.0,
                        100.0,
                        True,
                        False,
                        False,
                        False,
                        3,
                        0,
                    )
                    corrected_file = os.path.join(
                        mzML_corrected_dir, mzml_file)
                    MzMLFile().store(corrected_file, exp)


class SpectralMatcher:
    def run(self, database, mgf_file, output_mztab):
        # Load the MGF file to an MSExperiment format:
        exp = MSExperiment()
        MascotGenericFile().load(mgf_file, exp)

        # Perform spectral matching with a library in MGF format that is located under "resources":
        speclib = MSExperiment()
        MascotGenericFile().load(database, speclib)
        mztab = MzTab()
        out_merged = ""
        MSMS_match = MetaboliteSpectralMatching()
        MSMS_match_par = MSMS_match.getDefaults()
        MSMS_match_par.setValue("merge_spectra", "false")
        MSMS_match.setParameters(MSMS_match_par)
        MSMS_match.run(exp, speclib, mztab, String(out_merged))
        MzTabFile().store(output_mztab, mztab)


class Requantifier:
    def run(
        self,
        consensusXML_file,
        feature_matrix_df_file,
        mzML_dir,
        feautureXML_dir,
        mz_window_ppm,
    ):
        # get consensus df from consensusXML file
        cm = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, cm)

        # load feature_matrix
        df_cm = pd.read_csv(feature_matrix_df_file, sep="\t").set_index("id")

        # to map feature map file names to feature ids create a map
        map = {
            key: Path(value.filename).stem
            for key, value in cm.getColumnHeaders().items()
        }

        # create a database to store all information for requantification with consensus feature ids
        db = {}

        # get total number of files to check if re-quantification is required for a cf
        n_files_total = len(cm.getColumnHeaders().items())

        for cf in cm:
            f_list = cf.getFeatureList()
            if len(f_list) == n_files_total:
                continue  # skip cf that has all values
            # calculate mass delta in Da for upper and lower limits
            mz = cf.getMZ()
            delta_Da = mz_window_ppm * mz / 1000000
            db[cf.getUniqueId()] = {
                "mz_lower": mz - delta_Da,  # lower mz from cf mz
                "mz_upper": mz + delta_Da,  # upper mz from cf mz
                "file_to_id": {
                    map[f.getMapIndex()]: f.getUniqueId() for f in f_list
                },  # map files to feature ids to extract rt values later
                "rt_start": 0,  # here, RT start values will be summed up in the next step and later diveded by number of files
                "rt_end": 0,  # same for RT end values
            }

        # now add rt start end and m/z start end positions to each cf in db -> get from feature map dataframes
        for file in Path(feautureXML_dir).iterdir():
            fm = FeatureMap()
            FeatureXMLFile().load(str(file), fm)
            df = fm.get_df()
            for cf in db.keys():
                map = db[cf]["file_to_id"]
                if file.stem not in map.keys():
                    continue  # if the given file is not part of the consensus feature, skip this step
                f_id = map[file.stem]
                if file.stem in map.keys():
                    db[cf]["rt_start"] += df.loc[str(map[file.stem]),
                                                 "RTstart"]
                    db[cf]["rt_end"] += df.loc[str(map[file.stem]), "RTend"]

        # now that we have all the RT points summed up, devide by number of files in the consensus feature and remove "file_to_id" entry
        for cf in db.keys():
            n = len(db[cf]["file_to_id"])
            db[cf]["rt_start"] /= n
            db[cf]["rt_end"] /= n
            del db[cf]["file_to_id"]

        # with the complete database, iterate over the mzML files to do the actual requantification
        for file in Path(mzML_dir).iterdir():
            # get the mzML file name to link back to the consensus df
            name = str(file.name)
            exp = MSExperiment()
            MzMLFile().load(str(file), exp)
            df = exp.get_df()
            for cf in db.keys():
                rt_start = db[cf]["rt_start"]
                rt_end = db[cf]["rt_end"]
                df_filtered = df.query("RT > @rt_start and RT < @rt_end")
                tic = 0
                # for each matching spectrum extract the intensities between mz boundaries and add them to the TIC
                for _, row in df_filtered.iterrows():
                    tic += sum(
                        row["intarray"][
                            (
                                (row["mzarray"] > db[cf]["mz_lower"])
                                & (row["mzarray"] < db[cf]["mz_upper"])
                            )
                        ]
                    )
                # finally replace the entry in the consensus dataframe at index cf and column name with the re-quantified TIC values
                df_cm.loc[cf, name] = int(tic)

        df_cm.to_csv(feature_matrix_df_file, sep="\t")
