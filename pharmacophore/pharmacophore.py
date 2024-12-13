import os
import rdkit
from rdkit.Chem import AllChem, RDConfig


def feature_factory():
    feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    return feature_factory.GetFeatureFamilies()

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
        feature_factory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        features = feature_factory.GetFeaturesForMol(mol)
        print(f"Number of features: {len(features)}")
        for i, feat in enumerate(features):
            type = feat.GetFamily()
            if type == "Acceptor":
                pos = feat.GetPos()  # Feature position
                f.write(
                    f"pseudoatom feature_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=red\n"
                )
            if type == "Donor":
                pos = feat.GetPos()  # Feature position

                color_hex = (0, 0.9, 0)  # Convert to PyMOL-compatible color format
                f.write(
                    f"pseudoatom feature_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=blue\n"
                )
            if type == "Hydrophobe":
                pos = feat.GetPos()
                f.write(f"pseudoatom feature_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=green\n"
                        )
            if type == 'Aromatic':
                pos = feat.GetPos()
                f.write(f"pseudoatom feature_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=pink\n"
                        )
            if type == 'LumpedHydrophobe':
                pos = feat.GetPos()
                f.write(f"pseudoatom feature_{i}, pos=[{pos.x}, {pos.y}, {pos.z}], color=green\n"
                        )
        f.write("show spheres, feature_*\n")
        f.write(f"set sphere_scale, {sphere_size}\n")  # Adjust sphere size in PyMOL
    print(f"Feature visualization script written to {savepath}.")
