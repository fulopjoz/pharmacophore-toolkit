import collections
import pandas as pd
import numpy as np
from collections import defaultdict
from tqdm import tqdm
from rdkit import Chem
from pharmacophore.constants import feature_factory, FEATURES


class Pharmacophore:
    def __init__(self):
        self.sdf = None
        self.features = None

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

    def feature_types(self, features="default"):
        """
        A tuple containing default features from RDKit
        :return:
        """
        global phrase
        if features == "default":
            phrase = f"Default features: \n('Donor', 'Acceptor', 'Aromatic', 'Hydrophobe')"
        elif features == "rdkit":
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
            # get default features, count frequency, and append into list
            for mol in mols:
                feats = [feature.GetFamily() for feature in feat_factory.GetFeaturesForMol(mol)]
                feature_frequency = collections.Counter(feats)
                molecule_feature_frequencies.append(feature_frequency)

        # use custom features
        elif features == 'pharmacophore':
            feature_list = []
            # get features from constant file, count pharmacophore type, append to list
            for mol in mols:
                feats = self._calc_pharmacophore(mol)
                for feat in feats:
                    feature_list.append(feat[0])
                feature_frequency = collections.Counter(feature_list)
                molecule_feature_frequencies.append(feature_frequency)
                feature_list = []  # reset list to avoid double counting

        # # for troubleshooting
        # print(molecule_feature_frequencies)

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

    def calc_pharm(self, mol: Chem.Mol = None, features: str = 'default'):
        """
        Generate a list of pharmacophore features and position from a molecule.
        :param mol: Chem.Mol
            A molecule in ROMol format.
        :param features: str
            Designate the type of pharmacophore features to calculate from. Defaults to 'rdkit'.
        :return:
        """
        global pharmacophore

        # calculate pharmacophore by rdkit or pharmacophore dict
        if features == 'rdkit':
            pharmacophore = self._calc_rdkit(mol)
        elif features == 'default':
            pharmacophore = self._calc_pharmacophore(mol)

        # set features to self
        self.features = pharmacophore

        return pharmacophore

    def add_feats(self, mol: Chem.Mol = None, substruct: str = None, type: str = None):
        """
        Add specific features if not present using default or RDKit rules.
        :param mol: Chem.Mol
            Main molecule in ROMol format.
        :param substruct: str
            A SMARTS string for specific substructure to match
        :param type: str
            Set the type of interaction for substructure. Only 'Donor', 'Acceptor', 'Aromatic', and 'Hydrophobe' is
            allowed!
        :return:
        """
        # convert substruct smarts string into ROMol
        query = [Chem.MolFromSmarts(substruct)]

        # find matches
        matches = find_matches(mol, query)

        # check type
        if type == 'Donor' or type == 'Acceptor' or type == 'Aromatic' or type == 'Hydrophobe':
            pass
        else:
            raise ValueError(f"Type {type} not supported! Only 'Donor', 'Acceptor', 'Aromatic', and 'Hydrophobe' is "
                             f"allowed!")

        # obtain specific data from each match
        for match in matches:
            result = [type, match[0], match[1][0], match[1][1], match[1][2]]
            # add custom feats to self.features
            self.features.append(result)

        return self.features

    def output_features(self, features: list = None, savepath: str = None, sphere_size: float = 0.5,
                        color: dict = None):
        """
        Output features as a .pml format for visualization in PyMol.
        :param features: list
            A list containing features, corresponding atom number, and 3D position. Preferably genreated using
            calc_pharm.
        :param savepath: str = None
            Must be a file in .pml format.
        :param sphere_size: float = 0.5
            Set size of spheres.
        :param color: dict
            Set color of pharmacophores
        :return:
        """
        with open(savepath, "w") as f:
            # feat_factory = feature_factory
            # features = feat_factory.GetFeaturesForMol(mol)
            print(f"Number of features: {len(features)}")

            # to give sequential numbering for each group:
            feature_counts = {"Acceptor": 0, "Donor": 0, "Hydrophobe": 0, "Aromatic": 0, "LumpedHydrophobe": 0}

            # get features
            for feat in features:
                feature = feat[0]  # extract feature type
                feature_counts[feature] += 1  # give feature count
                count = feature_counts[feature]  # get current count
                # get feature position
                pos_x = feat[2]
                pos_y = feat[3]
                pos_z = feat[4]

                # Write to file with numbering
                if color is None:
                    color_type = {
                        "Acceptor": "red",
                        "Donor": "marine",
                        "Hydrophobe": "green",
                        "Aromatic": "pink",
                        "LumpedHydrophobe": "green"
                    }[feature]
                else:
                    color_type = color

                f.write(
                    f"pseudoatom {feature}_{count}, pos=[{pos_x}, {pos_y}, {pos_z}], color={color_type}\n"
                )

            f.write("show spheres, Acceptor_*\n")
            f.write("show spheres, Donor_*\n")
            f.write("show spheres, Hydrophobe_*\n")
            f.write("show spheres, Aromatic_*\n")
            f.write("show spheres, LumpedHydrophobe_*\n")
            f.write(f"set sphere_scale, {sphere_size}\n")  # Adjust sphere size in PyMOL
        print(f"Feature visualization script written to {savepath}.")

    def _calc_pharmacophore(self, mol: Chem.Mol = None):
        """
        Calculate pharmacophore features from a molecule using dict from constants
        :param mol: Chem.Mol
            Input molecule in ROMol format.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        global pharmacophore
        # read in custom feature dict
        constant_feats = FEATURES
        matches = {}

        # Identify matches to feature dict
        for key, value in constant_feats.items():
            try:
                query = [Chem.MolFromSmarts(smarts) for smarts in value]
                matches[key] = find_matches(mol, query, verbose=False)
            except:
                pass
        # remove duplicate SMARTS matches
        cleaned_matches = {}
        for key, value in matches.items():
            unique_lists = []
            for list in value:
                if list not in unique_lists:
                    unique_lists.append(list)
            cleaned_matches[key] = unique_lists
        # create list containing pharmacophore and centroid coordinates
        pharmacophore = []
        for key, value in cleaned_matches.items():
            for match in value:
                # extarct feature type, atom match, and XYZ position
                p = [key, match[0], match[1][0], match[1][1], match[1][2]]
                pharmacophore.append(p)

        # remove duplicates for aromatic/hydrophobe set
        atom_indices = set()
        final_pharmacophore = []

        # extract aromatic/hydrophobe matches by atom index adn remove hydrophobe
        for entry in pharmacophore:
            label, atom_indexes, *rest = entry
            index = frozenset(atom_indexes)  # Convert atom indexes to a frozenset so order does not matter

            # Check if the atom indexes are already in the set
            if index not in atom_indices:
                final_pharmacophore.append(entry)
                atom_indices.add(index)
            # keep Aromatic
            elif label == 'Aromatic':
                final_pharmacophore = [e for e in final_pharmacophore if not (e[0] == 'Hydrophobe' and frozenset(e[1]) == index)]
                final_pharmacophore.append(entry)
                atom_indices.add(index)

        return final_pharmacophore

    def _calc_rdkit(self, mol: Chem.Mol = None):
        """
        Calculate pharmacophore features from a molecule using default RDKit methods.
        :param mol: Chem.Mol
            input molecule in ROMol format.
        :return:
            List of pharmacophore type and centroid coordinates.
        """
        global pharmacophore

        # load default RDKit feature factory file
        feat_factory = feature_factory

        # get list of features for molecule
        features = feat_factory.GetFeaturesForMol(mol)

        # store pharmacophore features as list
        pharmacophore = []
        for feature in features:
            fam = feature.GetFamily()
            pos = feature.GetPos()
            atom_indices = feature.GetAtomIds()
            pharmacophore_item = [fam, atom_indices, pos[0], pos[1], pos[2]]
            pharmacophore.append(pharmacophore_item)

        return pharmacophore


def find_matches(mol: Chem.Mol = None, patterns: list[Chem.Mol] = None, verbose=True):
    matches = []
    for pattern in patterns:
        matched = mol.GetSubstructMatches(pattern)
        # get centroid and coordinates for each match
        for m in matched:
            try:
                centroid = _compute_match_centroid(mol, m)
                matches.append([m, centroid])
            except:
                pass
        # output statement if no matches
        if verbose is True:
            if len(matches) == 0:
                output_message = Chem.MolToSmarts(pattern)
                print(f"NO Matches to {output_message}!")
                return matches

    return matches


def _compute_match_centroid(mol, matched_pattern):
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in matched_pattern]
    center = np.mean(positions, axis=0)
    # convert result to float
    center = center.tolist()
    return tuple(center)
