import os
import collections
import pandas as pd
import rdkit
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig

class Pharmacophore:
    def __init__(self):
        self.sdf = None

    def read_sdf(self, sdf_file: str = None, verbose: bool = True):
        supplier = Chem.SDMolSupplier(sdf_file)
        # extract ROMols into a list
        mol_list = []
        if verbose:
            for mol in tqdm(supplier, desc="Reading Molecules"):
                mol_list.append(mol)
        else:
            for mol in supplier:
                mol_list.append(mol)

        self.sdf = mol_list

        return mol_list


def feature_factory():
    feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    return feature_factory.GetFeatureFamilies()


def molecule_features(mols: list = None, mol_name: list = None):
    molecule_feature_frequencies = []
    for mol in mols:
        features = [feature.GetFamily() for feature in feature_factory.GetFeaturesForMol(mol)]
        feature_frequency = collections.Counter(features)
        molecule_feature_frequencies.append(feature_frequency)

    feature_frequencies_df = (
        pd.DataFrame(
            molecule_feature_frequencies,
            index=[f"Mol{i}" for i, _ in enumerate(mols, 1)]
        )
        .fillna(0)
        .astype(int)
    )

    # reformat table and rename columns to molecule list
    feature_frequencies_df = feature_frequencies_df.transpose()
    feature_frequencies_df = feature_frequencies_df.set_axis(mol_name, axis=1)

    return feature_frequencies_df


def output_features(savepath: str = None, mol: rdkit.Chem.rdchem.Mol = None, sphere_size: float = 0.5):
    """

    :param savepath: str = None
        Must be a file in .pml format.
    :param mol: rdkit.Chem.rdchem.Mol = None
        Input RDKit Mol. 3D conformational form must be generated. Or input must be RDKit Mol read from an sdf file.
    :param sphere_size: float = 0.5
        Set size of spheres.
    :return:
    """
    with open(savepath, "w") as f:
        # todo create class structure?
        feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        features = feature_factory.GetFeaturesForMol(mol)
        print(f"Number of features: {len(features)}")
        for i, feat in enumerate(features):
            type = feat.GetFamily()
            if type == "Acceptor":
                pos = feat.GetPos()  # Feature position
                f.write(
                    f"pseudoatom Acceptor_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=red\n"
                )
            if type == "Donor":
                pos = feat.GetPos()  # Feature position
                f.write(
                    f"pseudoatom Donor_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=marine\n"
                )
            if type == "Hydrophobe":
                pos = feat.GetPos()
                f.write(f"pseudoatom Hydrophobe_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=green\n"
                        )
            if type == 'Aromatic':
                pos = feat.GetPos()
                f.write(f"pseudoatom Aromatic_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=pink\n"
                        )
            if type == 'LumpedHydrophobe':
                pos = feat.GetPos()
                f.write(f"pseudoatom LumpedHydrophobe_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=green\n"
                        )

        f.write("show spheres, Acceptor_*\n")
        f.write("show spheres, Donor_*\n")
        f.write("show spheres, Hydrophobe_*\n")
        f.write("show spheres, Aromatic_*\n")
        f.write("show spheres, LumpedHydrophobe_*\n")
        f.write(f"set sphere_scale, {sphere_size}\n")  # Adjust sphere size in PyMOL
    print(f"Feature visualization script written to {savepath}.")
