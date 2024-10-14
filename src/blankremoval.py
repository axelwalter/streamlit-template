from pyopenms import *
from pathlib import Path


class BlankRemoval:
    def run(featureXML_dir, blank_files, cutoff=0.3):
        # first we need to link the files (group)
        # 2nd export the dataframe with the individual feature IDs
        feature_maps = []
        for file in Path(featureXML_dir).iterdir():
            feature_map = FeatureMap()
            FeatureXMLFile().load(str(file), feature_map)
            feature_maps.append(feature_map)

        if len(feature_maps) == 1:
            fm = FeatureMap()
            fm.setMetaValue("spectra_data", [b"second-feature-map.mzML"])
            feature_maps.append(fm)

        feature_grouper = FeatureGroupingAlgorithmKD()

        consensus_map = ConsensusMap()
        file_descriptions = consensus_map.getColumnHeaders()

        for i, feature_map in enumerate(feature_maps):
            file_description = file_descriptions.get(i, ColumnHeader())
            file_description.filename = Path(
                feature_map.getMetaValue("spectra_data")[0].decode()
            ).name
            file_description.size = feature_map.size()
            file_descriptions[i] = file_description

        feature_grouper.group(feature_maps, consensus_map)
        consensus_map.setUniqueIds()
        consensus_map.setColumnHeaders(file_descriptions)

        # get intensities as a DataFrame
        df = consensus_map.get_df()
        for cf in consensus_map:
            if cf.metaValueExists("best ion"):
                df["adduct"] = [cf.getMetaValue(
                    "best ion") for cf in consensus_map]
                break
        df["feature_ids"] = [
            [handle.getUniqueId() for handle in cf.getFeatureList()]
            for cf in consensus_map
        ]
        df = df.reset_index()
        df = df.drop(columns=["sequence"])
        df = df.rename(columns={"RT": "RT(s)", "mz": "m/z"})

        # blank filtering (define the blanks== user input)
        df = df.drop(columns=["id", "charge", "quality", "m/z", "RT(s)"])

        samples = df[[col for col in df.columns if col not in blank_files]]
        blanks = df[[col for col in df.columns if col in blank_files]]

        # split samples and blanks
        # Getting mean for every feature in blank and Samples
        avg_blank = blanks.mean(
            axis=1, skipna=False
        )  # set skipna = False do not exclude NA/null values when computing the result.
        avg_samples = samples.drop(columns=["feature_ids"]).mean(axis=1, skipna=False)
        # Getting the ratio of blank vs samples
        ratio_blank_samples = (avg_blank + 1) / (avg_samples + 1)

        # Create an array with boolean values: True (is a real feature, ratio<cutoff) / False (is a blank, background, noise feature, ratio>cutoff)
        is_real_feature = ratio_blank_samples < cutoff
        # Calculating the number of background features and features present (sum(bg_bin) equals number of features to be removed)
        real_features = samples[is_real_feature.values]
        keep_ids = [
            item for sublist in real_features["feature_ids"] for item in sublist
        ]

        for file in Path(featureXML_dir).iterdir():
            if file.stem + ".mzML" in blank_files:
                file.unlink()
            else:
                feature_map = FeatureMap()
                FeatureXMLFile().load(str(file), feature_map)
                fmap_clean = FeatureMap(feature_map)
                fmap_clean.clear(False)
                for f in feature_map:
                    if f.getUniqueId() in keep_ids:
                        fmap_clean.push_back(f)

                FeatureXMLFile().store(str(Path(featureXML_dir, file.name)), fmap_clean)
