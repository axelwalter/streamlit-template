from pyopenms import *
from .common.common import reset_directory
from pathlib import Path
import pandas as pd
import os


class GNPSExport:
    def run(self, consensusXML_file, aligned_mzML_dir, gnps_dir):
        reset_directory(gnps_dir)
        all_mzML_files = [
            str(f)
            for f in Path(aligned_mzML_dir).iterdir()
            if f.is_file() and str(f).endswith("mzML")
        ]
        mzML_files = []
        for mzML in all_mzML_files:
            exp = MSExperiment()
            MzMLFile().load(mzML, exp)
            for spec in exp:
                if spec.getMSLevel() == 2:
                    mzML_files.append(mzML)

        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        filtered_map = ConsensusMap(consensus_map)
        filtered_map.clear(False)
        for feature in consensus_map:
            if feature.getPeptideIdentifications():
                filtered_map.push_back(feature)

        consensusXML_file = "filtered.consensusXML"
        ConsensusXMLFile().store(consensusXML_file, filtered_map)

        # for FFBM
        GNPSMGFFile().store(
            String(consensusXML_file),
            [file.encode() for file in mzML_files],
            String(os.path.join(gnps_dir, "MS2.mgf")),
        )
        GNPSQuantificationFile().store(
            consensus_map, os.path.join(
                gnps_dir, "FeatureQantificationTable.txt")
        )
        GNPSMetaValueFile().store(
            consensus_map, os.path.join(gnps_dir, "MetaValueTable.tsv")
        )

        # for IIMN
        IonIdentityMolecularNetworking().annotateConsensusMap(consensus_map)
        IonIdentityMolecularNetworking().writeSupplementaryPairTable(
            consensus_map, os.path.join(gnps_dir, "SupplementaryPairTable.csv")
        )

        # finally remove filtered.consensusXML file
        os.remove(consensusXML_file)

    def export_metadata_table_only(self, consensusXML_file, metadata_file):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        GNPSMetaValueFile().store(consensus_map, metadata_file)
        df = pd.read_csv(metadata_file, sep="\t", index_col=[0])
        df["ATTRIBUTE_Sample_Type"] = ["Sample"] * len(df["filename"])
        df.drop("ATTRIBUTE_MAPID", inplace=True, axis=1)
        df.to_csv(metadata_file, sep="\t", index=False)
