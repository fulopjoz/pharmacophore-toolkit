import os
import collections
import pandas as pd
import rdkit
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig
from pharmacophore.constants import feature_factory

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


    def default_features(self):
        """
        A tuple containing default features from RDKit
        :return:
        """
        # feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        phrase = f"Default features from RDKit: \n{feature_factory.GetFeatureFamilies()}"
        return phrase


    def to_df(self, mols: list = None, mol_name: list = None, features: str = 'rdkit'):
        """
        From a list of ROMols and molecule names, create a dataframe of features. Defaults to features from RDKit.
        :param mols:
        :param mol_name:

        :return:
        """
        global feat_factory
        molecule_feature_frequencies = []

        # use default features
        if features == 'rdkit':
            feat_factory = feature_factory
        for mol in mols:
            features = [feature.GetFamily() for feature in feat_factory.GetFeaturesForMol(mol)]
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

    def calc_pharm(self, mol: Chem.Mol = None, features: str = 'rdkit'):
        """
        Generate a list of pharmacophore features and position from a molecule.
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param features: str
            Designate the type of pharmacophore features to calculate from. Defaults to 'rdkit'.
        :return:
        """
        global pharmacophore

        # build feature factory
        if features == 'rdkit':
            feat_factory = feature_factory

            # get list of features for molecule
            features = feat_factory.GetFeaturesForMol(mol)

            # store pharmacophore features
            pharmacophore = []
            for feature in features:
                fam = feature.GetFamily()
                pos = feature.GetPos()
                atom_indices = feature.GetAtomIds()
                p = [fam, atom_indices, pos[0], pos[1], pos[2]]
                pharmacophore.append(p)

        return pharmacophore


    def output_features(self, savepath: str = None, mol: rdkit.Chem.rdchem.Mol = None, sphere_size: float = 0.5):
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
