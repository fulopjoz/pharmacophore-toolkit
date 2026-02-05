"""Systematic experiments for 3D pharmacophore shape representation.

This module provides a framework for investigating different molecular
representations for pharmacophore-based virtual screening using RDKit's
shape alignment tools.

The goal is to find representations that achieve AUC > 0.80 for
discriminating actives from decoys using 3D shape matching.
"""

from typing import List, Dict, Tuple, Optional, Callable
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol
from rdkit.Chem import rdShapeHelpers
from rdkit.Geometry import Point3D
from sklearn.metrics import roc_auc_score
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol


# Alternative fragment SMARTS sets for experimentation
FRAGMENT_SETS = {
    'original': {
        'Donor': '[NH3]',
        'Acceptor': '[O]',
        'Aromatic': 'c1ccccc1',
        'Hydrophobe': '[CH4]',
        'LumpedHydrophobe': '[CH4]',
        'PosIonizable': '[NH4+]',
    },
    'extended': {
        # Larger fragments with more chemical context
        'Donor': 'NC',           # Methylamine - has C neighbor
        'Acceptor': 'OC',        # Methanol oxygen
        'Aromatic': 'c1ccccc1C', # Toluene - has aliphatic context
        'Hydrophobe': 'CC(C)C',  # Isobutane - larger hydrophobic
        'LumpedHydrophobe': 'CC(C)C',
        'PosIonizable': '[NH3+]C',  # Methylammonium
    },
    'minimal': {
        # Single atoms only (smallest possible)
        'Donor': '[N]',
        'Acceptor': '[O]',
        'Aromatic': '[C]',
        'Hydrophobe': '[C]',
        'LumpedHydrophobe': '[C]',
        'PosIonizable': '[N+]',
    },
    'functional_groups': {
        # Actual functional groups from drugs
        'Donor': 'NC(=O)',           # Amide NH
        'Acceptor': 'C=O',           # Carbonyl
        'Aromatic': 'c1ccc(O)cc1',   # Phenol
        'Hydrophobe': 'C(C)(C)C',    # tert-butyl
        'LumpedHydrophobe': 'C(C)(C)C',
        'PosIonizable': '[NH3+]',    # Protonated amine
    },
}

# Atom mapping for connected scaffold
FEATURE_TO_ATOM = {
    'Donor': (7, 3),           # N with 3 H
    'Acceptor': (8, 0),        # O
    'Aromatic': (6, 0),        # C (aromatic marker)
    'Hydrophobe': (6, 4),      # C with 4 H
    'LumpedHydrophobe': (6, 4),
    'PosIonizable': (7, 4, 1), # N+ with 4 H
}


class ShapeExperiment:
    """Framework for shape representation experiments.

    Provides methods to:
    - Generate consensus pharmacophore from references
    - Score molecules against pharmacophore representations
    - Calculate AUC and other metrics
    - Run comparative experiments

    Example:
        >>> exp = ShapeExperiment(refs, actives, decoys)
        >>> features = exp.generate_consensus()
        >>> pharm_mol = create_connected_scaffold(features)
        >>> result = exp.run_experiment('connected', pharm_mol)
        >>> print(f"AUC: {result['auc']:.3f}")
    """

    def __init__(
        self,
        ref_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        n_conformers: int = 5
    ):
        """Initialize experiment with molecules.

        Args:
            ref_mols: Reference molecules for consensus generation
            actives: Active compounds (should score high)
            decoys: Decoy compounds (should score low)
            n_conformers: Number of conformers per query molecule
        """
        self.ref_mols = [m for m in ref_mols if m is not None]
        self.actives = self._prepare_molecules(actives, n_conformers)
        self.decoys = self._prepare_molecules(decoys, n_conformers)
        self.n_conformers = n_conformers
        self.results: Dict[str, Dict] = {}

    def _prepare_molecules(
        self,
        mols: List[Chem.Mol],
        n_conformers: int
    ) -> List[Chem.Mol]:
        """Ensure molecules have conformers."""
        prepared = []
        for mol in mols:
            if mol is None:
                continue

            # Check if already has conformers
            if mol.GetNumConformers() >= 1:
                prepared.append(mol)
                continue

            # Generate conformers
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)

            if mol_h.GetNumConformers() > 0:
                prepared.append(mol_h)

        return prepared

    def generate_consensus(
        self,
        tolerance: float = 2.0,
        occurrence: float = 0.3
    ) -> List:
        """Generate consensus features from references.

        Args:
            tolerance: Clustering tolerance in Angstroms
            occurrence: Minimum fraction of references with feature

        Returns:
            List of consensus features [type, indices, x, y, z]
        """
        consensus = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence
        )
        return consensus.generate_consensus(self.ref_mols)

    def score_molecules(
        self,
        pharm_mol: Chem.Mol,
        mols: List[Chem.Mol],
        shape_weight: float = 0.6
    ) -> List[float]:
        """Score molecules against pharmacophore mol.

        Args:
            pharm_mol: Pharmacophore representation as RDKit Mol
            mols: Query molecules to score
            shape_weight: Weight for shape vs color (0-1)

        Returns:
            List of best scores for each molecule
        """
        scores = []
        color_weight = 1.0 - shape_weight

        for mol in mols:
            best_score = 0.0
            for conf_id in range(mol.GetNumConformers()):
                try:
                    shape, color = AlignMol(
                        ref=pharm_mol,
                        probe=mol,
                        probeConfId=conf_id,
                        useColors=True,
                        opt_param=0.5
                    )
                    score = shape_weight * shape + color_weight * color
                    best_score = max(best_score, score)
                except Exception:
                    continue
            scores.append(best_score)

        return scores

    def calculate_auc(
        self,
        active_scores: List[float],
        decoy_scores: List[float]
    ) -> float:
        """Calculate ROC-AUC from scores.

        Args:
            active_scores: Scores for active molecules
            decoy_scores: Scores for decoy molecules

        Returns:
            ROC-AUC score (0.5 = random)
        """
        y_true = [1] * len(active_scores) + [0] * len(decoy_scores)
        y_scores = active_scores + decoy_scores

        try:
            return roc_auc_score(y_true, y_scores)
        except ValueError:
            return 0.5

    def run_experiment(
        self,
        name: str,
        pharm_mol: Chem.Mol,
        shape_weight: float = 0.6
    ) -> Dict:
        """Run single experiment and return metrics.

        Args:
            name: Experiment name for tracking
            pharm_mol: Pharmacophore molecule to test
            shape_weight: Weight for shape vs color

        Returns:
            Dictionary with auc, separation, and other metrics
        """
        if pharm_mol is None:
            return {
                'name': name,
                'auc': 0.5,
                'separation': 0.0,
                'active_mean': 0.0,
                'decoy_mean': 0.0,
                'n_atoms': 0,
                'n_bonds': 0,
                'error': 'pharm_mol is None'
            }

        active_scores = self.score_molecules(pharm_mol, self.actives, shape_weight)
        decoy_scores = self.score_molecules(pharm_mol, self.decoys, shape_weight)

        auc = self.calculate_auc(active_scores, decoy_scores)
        separation = np.mean(active_scores) - np.mean(decoy_scores)

        result = {
            'name': name,
            'auc': auc,
            'separation': separation,
            'active_mean': np.mean(active_scores),
            'active_std': np.std(active_scores),
            'decoy_mean': np.mean(decoy_scores),
            'decoy_std': np.std(decoy_scores),
            'n_atoms': pharm_mol.GetNumAtoms(),
            'n_bonds': pharm_mol.GetNumBonds(),
            'n_fragments': len(Chem.GetMolFrags(pharm_mol)),
        }
        self.results[name] = result
        return result


def create_connected_scaffold(
    features: List,
    linker_type: str = 'direct'
) -> Optional[Chem.Mol]:
    """Create connected molecular scaffold from pharmacophore features.

    Strategy: Build minimum spanning tree connecting all feature positions,
    then add bonds along edges to create connected graph.

    Args:
        features: List of [type, indices, x, y, z]
        linker_type: 'direct' (single bond), 'carbon' (CH2), or 'ether' (O)

    Returns:
        Connected RDKit Mol with features at specified positions
    """
    if len(features) < 2:
        return None

    # Extract positions
    positions = np.array([[f[2], f[3], f[4]] for f in features])

    # Calculate distance matrix and MST
    dist_matrix = squareform(pdist(positions))
    mst = minimum_spanning_tree(dist_matrix)

    # Create editable mol
    mol = Chem.RWMol()

    # Add feature atoms
    feature_indices = []
    for i, feat in enumerate(features):
        ftype = feat[0]
        atom_info = FEATURE_TO_ATOM.get(ftype, (6, 0))

        atom = Chem.Atom(atom_info[0])
        if len(atom_info) > 2:
            atom.SetFormalCharge(atom_info[2])
        idx = mol.AddAtom(atom)
        feature_indices.append(idx)

    # Create conformer with feature positions
    conf = Chem.Conformer(len(features))
    for i, feat in enumerate(features):
        pos = Point3D(feat[2], feat[3], feat[4])
        conf.SetAtomPosition(i, pos)

    mol.AddConformer(conf, assignId=True)

    # Add bonds along MST edges
    mst_array = mst.toarray()
    for i in range(len(features)):
        for j in range(i + 1, len(features)):
            if mst_array[i, j] > 0 or mst_array[j, i] > 0:
                mol.AddBond(feature_indices[i], feature_indices[j],
                           Chem.BondType.SINGLE)

    # Skip sanitization - valences may be unusual
    # but shape alignment doesn't require sanitized mol

    return mol.GetMol()


def create_pharmacophore_mol_variant(
    features: List,
    fragment_set: str = 'original'
) -> Optional[Chem.Mol]:
    """Create pharmacophore mol using different fragment sets.

    Args:
        features: List of [type, indices, x, y, z]
        fragment_set: Name of fragment set to use

    Returns:
        RDKit Mol with fragments at feature positions
    """
    if fragment_set not in FRAGMENT_SETS:
        raise ValueError(f"Unknown fragment set: {fragment_set}")

    fragments = FRAGMENT_SETS[fragment_set]
    mol = Chem.RWMol()
    all_atom_positions = []

    for feat in features:
        ftype = feat[0]
        smiles = fragments.get(ftype, '[C]')

        frag = Chem.MolFromSmiles(smiles)
        if frag is None:
            continue

        # Add hydrogens and embed
        frag = Chem.AddHs(frag)
        AllChem.EmbedMolecule(frag, randomSeed=42)

        if frag.GetNumConformers() == 0:
            # Fallback: single atom
            frag = Chem.MolFromSmiles('[C]')
            frag = Chem.AddHs(frag)
            AllChem.EmbedMolecule(frag, randomSeed=42)

        # Get fragment centroid
        frag_conf = frag.GetConformer()
        centroid = np.zeros(3)
        for atom_idx in range(frag.GetNumAtoms()):
            pos = frag_conf.GetAtomPosition(atom_idx)
            centroid += np.array([pos.x, pos.y, pos.z])
        centroid /= frag.GetNumAtoms()

        # Calculate offset to move centroid to feature position
        target_pos = np.array([feat[2], feat[3], feat[4]])
        offset = target_pos - centroid

        # Add atoms with offset positions
        atom_offset = mol.GetNumAtoms()
        for atom in frag.GetAtoms():
            new_idx = mol.AddAtom(atom)
            orig_pos = frag_conf.GetAtomPosition(atom.GetIdx())
            new_pos = np.array([orig_pos.x, orig_pos.y, orig_pos.z]) + offset
            all_atom_positions.append((new_idx, new_pos))

        # Add bonds
        for bond in frag.GetBonds():
            begin_idx = bond.GetBeginAtomIdx() + atom_offset
            end_idx = bond.GetEndAtomIdx() + atom_offset
            mol.AddBond(begin_idx, end_idx, bond.GetBondType())

    # Create conformer with all positions
    if len(all_atom_positions) > 0:
        conf = Chem.Conformer(mol.GetNumAtoms())
        for idx, pos in all_atom_positions:
            conf.SetAtomPosition(idx, Point3D(pos[0], pos[1], pos[2]))
        mol.AddConformer(conf, assignId=True)

    return mol.GetMol()


def add_pseudo_bonds(
    mol: Chem.Mol,
    max_distance: float = 15.0
) -> Chem.Mol:
    """Add pseudo-bonds between disconnected fragments.

    Uses BondType.SINGLE to connect fragments using MST.
    Creates connected graph for shape alignment while
    preserving fragment identities.

    Args:
        mol: Molecule with disconnected fragments
        max_distance: Maximum distance for pseudo-bonds

    Returns:
        Connected molecule with pseudo-bonds
    """
    frags = Chem.GetMolFrags(mol, asMols=False)

    if len(frags) <= 1:
        return mol  # Already connected

    rwmol = Chem.RWMol(mol)
    conf = mol.GetConformer()

    # Get anchor atoms (first atom of each fragment)
    anchors = [frag[0] for frag in frags]
    positions = []
    for a in anchors:
        pos = conf.GetAtomPosition(a)
        positions.append([pos.x, pos.y, pos.z])

    # Calculate MST
    dist_matrix = squareform(pdist(positions))
    mst = minimum_spanning_tree(dist_matrix).toarray()

    # Add pseudo-bonds along MST
    for i in range(len(anchors)):
        for j in range(i + 1, len(anchors)):
            if mst[i, j] > 0 and mst[i, j] <= max_distance:
                rwmol.AddBond(anchors[i], anchors[j], Chem.BondType.SINGLE)

    return rwmol.GetMol()


def score_gaussian_overlap(
    ref_mol: Chem.Mol,
    probe_mol: Chem.Mol
) -> float:
    """Score using Gaussian density overlap instead of grid-based.

    Args:
        ref_mol: Reference molecule
        probe_mol: Probe molecule to score

    Returns:
        Gaussian overlap score
    """
    try:
        overlap = rdShapeHelpers.ComputeGaussianOverlap(ref_mol, probe_mol)
        return overlap
    except Exception:
        return 0.0


def score_shape_protrude(
    ref_mol: Chem.Mol,
    probe_mol: Chem.Mol
) -> float:
    """Score using shape protrusion distance.

    Args:
        ref_mol: Reference molecule
        probe_mol: Probe molecule to score

    Returns:
        Inverted protrusion score (higher = better match)
    """
    try:
        protrude = rdShapeHelpers.ShapeProtrudeDist(ref_mol, probe_mol)
        # Lower protrusion = better match, so invert
        return 1.0 - protrude
    except Exception:
        return 0.0


def score_against_references(
    query_mol: Chem.Mol,
    ref_mols: List[Chem.Mol],
    shape_weight: float = 0.6
) -> float:
    """Score query against best-matching reference molecule.

    Args:
        query_mol: Query molecule to score
        ref_mols: Reference molecules
        shape_weight: Weight for shape vs color

    Returns:
        Best score across all references and conformers
    """
    best_score = 0.0
    color_weight = 1.0 - shape_weight

    for ref in ref_mols:
        if ref is None:
            continue
        for conf_id in range(query_mol.GetNumConformers()):
            try:
                shape, color = AlignMol(
                    ref=ref,
                    probe=query_mol,
                    probeConfId=conf_id,
                    useColors=True,
                    opt_param=0.5
                )
                score = shape_weight * shape + color_weight * color
                best_score = max(best_score, score)
            except Exception:
                continue

    return best_score


def run_all_experiments(
    refs: List[Chem.Mol],
    actives: List[Chem.Mol],
    decoys: List[Chem.Mol],
    tolerance: float = 2.0,
    occurrence: float = 0.3,
    shape_weight: float = 0.6,
    verbose: bool = True
) -> Dict[str, Dict]:
    """Run all shape representation experiments.

    Args:
        refs: Reference molecules
        actives: Active compounds
        decoys: Decoy compounds
        tolerance: Consensus tolerance
        occurrence: Consensus occurrence threshold
        shape_weight: Weight for shape vs color
        verbose: Print progress

    Returns:
        Dictionary of experiment results
    """
    exp = ShapeExperiment(refs, actives, decoys)
    features = exp.generate_consensus(tolerance=tolerance, occurrence=occurrence)

    if verbose:
        print("=" * 60)
        print("SHAPE REPRESENTATION INVESTIGATION")
        print("=" * 60)
        print(f"References: {len(refs)}")
        print(f"Actives: {len(exp.actives)}")
        print(f"Decoys: {len(exp.decoys)}")
        print(f"Consensus features: {len(features)}")
        print(f"Parameters: tol={tolerance}, occ={occurrence}, sw={shape_weight}")
        print()

    # Define experiments
    experiments = [
        ("A. Original disconnected",
         lambda: PharmacophoreToMol.convert(features, enable_color_features=True)),
        ("B. Connected scaffold",
         lambda: create_connected_scaffold(features)),
        ("C. Extended fragments",
         lambda: create_pharmacophore_mol_variant(features, 'extended')),
        ("D. Minimal atoms",
         lambda: create_pharmacophore_mol_variant(features, 'minimal')),
        ("E. Functional groups",
         lambda: create_pharmacophore_mol_variant(features, 'functional_groups')),
        ("F. Pseudo-bonds",
         lambda: add_pseudo_bonds(
             PharmacophoreToMol.convert(features, enable_color_features=True))),
    ]

    results = []
    for name, mol_func in experiments:
        try:
            pharm_mol = mol_func()
            if pharm_mol is not None:
                result = exp.run_experiment(name, pharm_mol, shape_weight)
                results.append(result)
                if verbose:
                    status = "PASS" if result['auc'] > 0.80 else (
                        "WEAK" if result['auc'] > 0.60 else "FAIL")
                    print(f"{name}")
                    print(f"  AUC: {result['auc']:.3f} [{status}]")
                    print(f"  Separation: {result['separation']:.3f}")
                    print(f"  Atoms: {result['n_atoms']}, Bonds: {result['n_bonds']}, "
                          f"Fragments: {result['n_fragments']}")
                    print()
        except Exception as e:
            if verbose:
                print(f"{name}: FAILED - {e}")
                print()

    # Reference control
    if verbose:
        print("G. Reference alignment (CONTROL)")

    active_scores = [score_against_references(m, refs, shape_weight)
                    for m in exp.actives]
    decoy_scores = [score_against_references(m, refs, shape_weight)
                   for m in exp.decoys]
    ref_auc = exp.calculate_auc(active_scores, decoy_scores)
    ref_sep = np.mean(active_scores) - np.mean(decoy_scores)

    ref_result = {
        'name': 'G. Reference (control)',
        'auc': ref_auc,
        'separation': ref_sep,
        'active_mean': np.mean(active_scores),
        'decoy_mean': np.mean(decoy_scores),
    }
    results.append(ref_result)
    exp.results['G. Reference (control)'] = ref_result

    if verbose:
        print(f"  AUC: {ref_auc:.3f} [CONTROL]")
        print(f"  Separation: {ref_sep:.3f}")
        print()

        # Summary table
        print("=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        print(f"{'Experiment':<30} {'AUC':>8} {'Sep':>8} {'Status':>10}")
        print("-" * 60)
        for r in sorted(results, key=lambda x: -x['auc']):
            if 'control' in r['name'].lower():
                status = "CONTROL"
            elif r['auc'] > 0.80:
                status = "PASS"
            elif r['auc'] > 0.60:
                status = "WEAK"
            else:
                status = "FAIL"
            print(f"{r['name']:<30} {r['auc']:>8.3f} {r['separation']:>8.3f} {status:>10}")

    return exp.results
