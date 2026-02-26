"""Shared alignment + feature extraction for pharmacophore distance methods.

Aligns a query molecule to a reference molecule using RDKit shape alignment,
then extracts pharmacophore features from both in the same coordinate frame.

This solves the root cause of broken Hungarian/OT scoring: comparing features
from molecules in different coordinate frames produces meaningless distances.
"""

import logging
from copy import deepcopy
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem.rdShapeAlign import AlignMol

from .pharmacophore import Pharmacophore

logger = logging.getLogger(__name__)


def align_and_extract_features(
    ref_mol: Chem.Mol,
    query_mol: Chem.Mol,
    pharmacophore: Pharmacophore = None,
    opt_param: float = 0.5,
    max_preiters: int = 10,
    max_postiters: int = 30,
) -> Tuple[List, List, float, float]:
    """Align query to reference, then extract features from both.

    AlignMol modifies the probe conformer in-place, so the caller should
    pass a copy if the original must be preserved.

    Args:
        ref_mol: Reference molecule with 3D conformer.
        query_mol: Query molecule with 3D conformer (will be modified in-place).
        pharmacophore: Pharmacophore instance for feature extraction.
            If None, a default one is created.
        opt_param: AlignMol optimization parameter (0=color, 1=shape).
        max_preiters: Maximum pre-optimization iterations.
        max_postiters: Maximum post-optimization iterations.

    Returns:
        Tuple of (ref_features, query_features, shape_tanimoto, color_tanimoto).

    Raises:
        ValueError: If either molecule has no conformer.
    """
    if ref_mol.GetNumConformers() == 0:
        raise ValueError("Reference molecule has no conformer")
    if query_mol.GetNumConformers() == 0:
        raise ValueError("Query molecule has no conformer")

    if pharmacophore is None:
        pharmacophore = Pharmacophore()

    # Add Hs for better alignment (AlignMol benefits from explicit H)
    ref_h = Chem.AddHs(ref_mol, addCoords=True)
    query_h = Chem.AddHs(query_mol, addCoords=True)

    # AlignMol aligns probe (query) onto ref IN-PLACE
    shape_tani, color_tani = AlignMol(
        ref=ref_h,
        probe=query_h,
        probeConfId=0,
        useColors=True,
        opt_param=opt_param,
        max_preiters=max_preiters,
        max_postiters=max_postiters,
    )

    # Extract features — both now in the reference coordinate frame
    ref_features = pharmacophore.calc_pharm(mol=ref_h)
    query_features = pharmacophore.calc_pharm(mol=query_h)

    return ref_features, query_features, shape_tani, color_tani
