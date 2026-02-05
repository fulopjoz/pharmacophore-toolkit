"""
Hold constant values
"""
import os
import sys
import matplotlib.colors as mcolors
from rdkit.Chem import AllChem, RDConfig


def _find_base_features_file():
    """
    Find the BaseFeatures.fdef file in the RDKit installation.
    Handles both properly configured and misconfigured RDKit installations.
    """
    # Try the official RDConfig path first
    default_path = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    if os.path.exists(default_path):
        return default_path

    # If RDDataDir is relative or incorrect, search for it
    import rdkit
    rdkit_path = os.path.dirname(rdkit.__file__)

    # Common locations relative to rdkit package
    search_paths = [
        # Conda installations
        os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(rdkit_path))),
                     'share', 'RDKit', 'Data', 'BaseFeatures.fdef'),
        # System installations
        os.path.join(rdkit_path, 'Data', 'BaseFeatures.fdef'),
        # Alternative conda structure
        os.path.join(sys.prefix, 'share', 'RDKit', 'Data', 'BaseFeatures.fdef'),
    ]

    for path in search_paths:
        if os.path.exists(path):
            return path

    raise FileNotFoundError(
        f"Could not find BaseFeatures.fdef file. "
        f"RDConfig.RDDataDir points to: {RDConfig.RDDataDir} "
        f"Searched in: {search_paths}"
    )


feature_factory = AllChem.BuildFeatureFactory(_find_base_features_file())

# Feature importance weights based on binding energy contributions
# Source: Langer Presentation 2025, Lipinski's Rule of 5
# H-bonds (~5 kcal/mol) > π-interactions (~3 kcal/mol) > vdW (~1 kcal/mol)
FEATURE_WEIGHTS = {
    'Donor': 2.0,      # Strong H-bond donors (N-H, O-H)
    'Acceptor': 2.0,   # Strong H-bond acceptors (C=O, N, O)
    'Aromatic': 1.5,   # π-π stacking, cation-π interactions
    'Hydrophobe': 1.0, # Weak van der Waals contacts (baseline)
}

# Spatial tolerance scaling per feature type for consensus clustering.
# Multiplied by base tolerance: effective_tolerance = tolerance * SPATIAL_SCALE[type]
# Source: PharmaGist (Schneidman-Duhovny 2008) — hydrophobic features are
# positionally less precise than H-bond features and should cluster more loosely.
SPATIAL_SCALE = {
    'Donor': 0.8,            # Tighter — H-bonds are positionally precise
    'Acceptor': 0.8,         # Tighter — H-bonds are positionally precise
    'Aromatic': 1.0,         # Baseline
    'Hydrophobe': 1.5,       # Looser — position less critical for vdW
    'LumpedHydrophobe': 1.5, # Same as Hydrophobe
    'PosIonizable': 0.8,     # Tighter — charge interactions are precise
}

# Type compatibility matrix for feature matching across pharmacophore models.
# Values: 0.0 = same type, 0.5 = compatible, 1.0 = incompatible.
# Used by Hungarian matching (hungarian_matching.py) and OT scoring (ot_scoring.py).
FEATURE_COMPATIBILITY = {
    ('Donor', 'Acceptor'): 0.5,         # Both H-bond
    ('Donor', 'PosIonizable'): 0.5,     # Both carry H
    ('Acceptor', 'PosIonizable'): 0.7,  # Weak compatibility
    ('Aromatic', 'Hydrophobe'): 0.5,    # Aromatic rings are hydrophobic
    ('Aromatic', 'LumpedHydrophobe'): 0.5,
    ('Hydrophobe', 'LumpedHydrophobe'): 0.0,  # Same concept
}


def get_type_distance(type_a: str, type_b: str) -> float:
    """Get the distance between two feature types from the compatibility matrix.

    Args:
        type_a: First feature type.
        type_b: Second feature type.

    Returns:
        Distance: 0.0 (same), 0.5 (compatible), 1.0 (incompatible).
    """
    if type_a == type_b:
        return 0.0
    key = (type_a, type_b)
    rev_key = (type_b, type_a)
    if key in FEATURE_COMPATIBILITY:
        return FEATURE_COMPATIBILITY[key]
    if rev_key in FEATURE_COMPATIBILITY:
        return FEATURE_COMPATIBILITY[rev_key]
    return 1.0

FEATURES = {
    "Donor": ["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",  # nitrogen in rings or amines
              "[#8!H0&!$([OH][C,S,P]=O)]",  # oxygen atom bonded to hydrogen
              "[#16!H0]"  # sulfur atom
              ],
    "Acceptor": ["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
                 # amines or nitrogen in rings non aromatic rings
                 "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"  # hydroxyl or ether
                 ],
    "Aromatic": ["a1aaaaa1",  # six membered aromatic ring
                 "a1aaaa1",  # five membered aromatic ring
                 "[#6]1[#6]=[#6][#6]=[#7]1",  # fused pyrole ring
                 "[#6]:1:[#6]:[#6]:[#6]:[#6]:[#7]:1"  # fused pyridine ring
                 ],
    "Hydrophobe": ["a1aaaaa1",  # six member aromatic ring
                   "a1aaaa1",  # five member aromatic ring
                   "*~1~*~*~*~*~*~1",  # fix member ring, any atom, any bond
                   "*~1~*~*~*~*~1",  # five member ring, any atom, any bond
                   # methyl, methylene, methine, halogens
                   "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
                   # matches to methyl, methylene, terminal methine
                   "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
                   "[CH3]",  # terminal methyl group
                   "[CH2]~*~!@[*1]",  # terminal methylene group not in a ring
                   ]
}

FEATURE_COLORS = {
    "Donor": (0.2549019607843137, 0.4117647058823529, 0.8823529411764706),  # royalblue
    "Acceptor": (1.0, 0.27058823529411763, 0.0),  # orangered
    "Aromatic": (0.8549019607843137, 0.6470588235294118, 0.12549019607843137),  # goldenrod
    "Hydrophobe": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),  # seagreen
    "LumpedHydrophobe": (0.1803921568627451, 0.5450980392156862, 0.3411764705882353),  # seagreen
    "PosIonizable": (0.0, 0.7490196078431373, 1.0)  # deepskyblue
}

# for rendering py3dmol
INTERACTIVE_COLORS = {
    "Donor": "royalblue",
    "Acceptor": "orangered",
    "Aromatic": "goldenrod",
    "Hydrophobe": "seagreen",
    "LumpedHydrophobe": "seagreen",
    "PosIonizable": "deepskyblue"
}


def color_convert(color: str = None):
    """
    helper funciton to convert color to rgb.
    :param color: str
        Color to convert to rgb. Can be hex or color name.
    :return:
    """
    try:
        # convert color to rgb
        rgb = mcolors.to_rgb(color)
        return rgb
    except:
        raise ValueError(f"{color} is not a valid color!")


if __name__ == "__main__":
    import doctest

    doctest.testmod()
