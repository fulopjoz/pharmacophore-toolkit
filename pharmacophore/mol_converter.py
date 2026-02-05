"""Convert pharmacophore models to RDKit molecules.

This module provides conversion from pharmacophore point clouds to
RDKit Mol objects for shape-based alignment and comparison using
RDKit's shape tools (ShapeProtrudeDist, AlignMol, etc.).

The conversion creates pseudo-molecules with minimal molecular fragments
positioned at pharmacophore feature locations. When enable_color_features=True,
proper molecular fragments are created to support pharmacophore SMARTS matching:
- Donor: NH3 (ammonia) or NH2 fragment with N-H bonds
- Acceptor: Lone oxygen or nitrogen atom
- Aromatic: Benzene ring fragment
- Hydrophobe: CH4 (methane) or CH3 fragment

This enables both shape-based alignment AND pharmacophore color scoring.
"""

from typing import List, Dict, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdGeometry
import numpy as np


class PharmacophoreToMol:
    """Convert pharmacophore features to RDKit Mol object.
    
    This class creates a pseudo-molecule representation of a pharmacophore
    model by placing molecular fragments at feature positions. The fragments
    are chosen to represent the pharmacophore feature types:
    
    - Donor → Ammonia (NH3) or amine fragment
    - Acceptor → Lone oxygen atom
    - Aromatic → Benzene ring fragment
    - Hydrophobe → Methane (CH4) or methyl fragment
    - LumpedHydrophobe → Methane (CH4)
    - PosIonizable → Protonated amine (NH4+)
    
    When enable_color_features=True, the fragments include proper bonding
    to support pharmacophore SMARTS pattern matching for color-based alignment.
    When False, only atoms are created for pure shape-based comparison.
    
    Example:
        >>> from pharmacophore.mol_converter import PharmacophoreToMol
        >>> 
        >>> # Consensus features from PharmacophoreConsensus
        >>> features = [
        ...     ['Donor', (), 1.0, 2.0, 3.0],
        ...     ['Acceptor', (), 4.0, 5.0, 6.0]
        ... ]
        >>> 
        >>> # Convert to RDKit Mol with color features enabled
        >>> mol = PharmacophoreToMol.convert(
        ...     features, 
        ...     name='Consensus',
        ...     enable_color_features=True
        ... )
        >>> print(f"Atoms: {mol.GetNumAtoms()}")  # Will include H atoms
    """
    
    # Mapping from pharmacophore feature types to element symbols (shape-only mode)
    FEATURE_TO_ELEMENT = {
        'Donor': 'N',
        'Acceptor': 'O',
        'Aromatic': 'C',
        'Hydrophobe': 'C',
        'LumpedHydrophobe': 'C',
        'PosIonizable': 'N'
    }
    
    # Mapping from pharmacophore feature types to SMILES fragments (color mode)
    # These fragments are chosen to be recognized by rdShapeAlign's color features
    #
    # Fragment SMILES version history:
    # v1 (<=0.0.4): Acceptor='[O]', Hydrophobe='[CH4]'
    #   - Did not match pharmacophore SMARTS → color scores returned 0.0
    # v2 (>=0.0.5): Acceptor='C=O', Hydrophobe='C1CCCC1'
    #   - Matches SMARTS patterns for color-based alignment scoring
    #   - Breaking change: atom counts per fragment differ
    FRAGMENT_SMILES_VERSION = 2

    FEATURE_TO_SMILES = {
        'Donor': '[NH3]',           # Ammonia - matches [#7!H0] donor SMARTS
        'Acceptor': 'C=O',          # Formaldehyde - matches acceptor SMARTS (v1 used [O])
        'Aromatic': 'c1ccccc1',     # Benzene - matches a1aaaaa1 aromatic SMARTS
        'Hydrophobe': 'C1CCCC1',    # Cyclopentane - matches hydrophobe (v1 used [CH4])
        'LumpedHydrophobe': 'C1CCCC1',
        'PosIonizable': '[NH4+]'    # Ammonium - matches positive ionizable
    }
    
    @classmethod
    def _create_fragment_from_smiles(
        cls,
        smiles: str,
        position: Tuple[float, float, float],
        feature_type: str
    ) -> Chem.Mol:
        """Create a molecular fragment from SMILES at given position.
        
        Args:
            smiles: SMILES string for the fragment.
            position: (x, y, z) coordinates for the fragment center.
            feature_type: Type of pharmacophore feature.
            
        Returns:
            RDKit Mol fragment with 3D coordinates centered at position.
        """
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for {feature_type}: {smiles}")
        
        # Add hydrogens if not already present
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates (small fragment)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Get conformer and center it at the desired position
        conf = mol.GetConformer()
        
        # Calculate current centroid
        centroid = np.array([0.0, 0.0, 0.0])
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            centroid += np.array([pos.x, pos.y, pos.z])
        centroid /= mol.GetNumAtoms()
        
        # Translate fragment to desired position
        target = np.array(position)
        translation = target - centroid
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            new_pos = np.array([pos.x, pos.y, pos.z]) + translation
            conf.SetAtomPosition(i, rdGeometry.Point3D(*new_pos))
        
        return mol
    
    @classmethod
    def _combine_fragments(
        cls,
        fragments: List[Chem.Mol],
        name: str = 'Pharmacophore'
    ) -> Chem.Mol:
        """Combine multiple molecular fragments into a single molecule.
        
        Args:
            fragments: List of RDKit Mol fragments.
            name: Name for the combined molecule.
            
        Returns:
            Combined RDKit Mol with all fragments as disconnected components.
        """
        if not fragments:
            raise ValueError("fragments list cannot be empty")
        
        # Start with first fragment
        combined = Chem.RWMol(fragments[0])
        
        # Add all other fragments
        for frag in fragments[1:]:
            # Add atoms and bonds from fragment
            atom_offset = combined.GetNumAtoms()
            
            for atom in frag.GetAtoms():
                new_idx = combined.AddAtom(atom)
            
            for bond in frag.GetBonds():
                begin_idx = bond.GetBeginAtomIdx() + atom_offset
                end_idx = bond.GetEndAtomIdx() + atom_offset
                bond_type = bond.GetBondType()
                combined.AddBond(begin_idx, end_idx, bond_type)
        
        # Get conformers from all fragments and combine
        combined_mol = combined.GetMol()

        # CRITICAL: Remove existing conformers from RWMol construction
        # RWMol(fragments[0]) copies the first fragment's conformer, but it only
        # has positions for that fragment's atoms. The new atoms added later have
        # default (0,0,0) positions. We must remove this incomplete conformer and
        # create a fresh one with all atom positions correctly set.
        combined_mol.RemoveAllConformers()

        combined_conf = Chem.Conformer(combined_mol.GetNumAtoms())

        atom_idx = 0
        for frag in fragments:
            frag_conf = frag.GetConformer()
            for i in range(frag.GetNumAtoms()):
                pos = frag_conf.GetAtomPosition(i)
                # CRITICAL: Create explicit Point3D to avoid reference issues
                # Passing pos directly from GetAtomPosition() can fail due to
                # RDKit internal reference handling
                combined_conf.SetAtomPosition(
                    atom_idx, rdGeometry.Point3D(pos.x, pos.y, pos.z)
                )
                atom_idx += 1

        combined_mol.AddConformer(combined_conf, assignId=True)
        
        # CRITICAL: Initialize ring info and other molecular properties
        # This is required for SMARTS matching and pharmacophore feature detection
        try:
            Chem.SanitizeMol(combined_mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        except Exception as e:
            # If full sanitization fails, try minimal sanitization
            try:
                Chem.FastFindRings(combined_mol)
            except:
                pass  # Even if this fails, molecule might still be usable for shape
        
        combined_mol.SetProp('_Name', name)
        
        return combined_mol
    
    @classmethod
    def convert(
        cls,
        features: List[List],
        name: str = 'Pharmacophore',
        enable_color_features: bool = True
    ) -> Chem.Mol:
        """Convert pharmacophore features to RDKit Mol.
        
        Creates a pseudo-molecule with molecular fragments positioned at
        pharmacophore feature coordinates. Each feature becomes a small
        molecular fragment chosen to represent the feature type.
        
        When enable_color_features=True, creates proper molecular fragments
        (NH3, benzene, CH4, etc.) with bonds to support pharmacophore SMARTS
        matching for color-based alignment.
        
        When enable_color_features=False, creates only single atoms without
        bonds for pure shape-based alignment.
        
        Args:
            features: List of pharmacophore features, where each feature is:
                [type, atom_indices, x, y, z]
            name: Name to assign to the molecule (default: 'Pharmacophore').
            enable_color_features: If True, create molecular fragments with
                bonds for pharmacophore color matching (default: True).
                If False, create only atoms for shape-only alignment.
        
        Returns:
            RDKit Mol object with fragments at feature positions.
        
        Raises:
            ValueError: If features is empty or contains invalid data.
        
        Example:
            >>> features = [['Donor', (), 1.0, 2.0, 3.0]]
            >>> mol = PharmacophoreToMol.convert(features, enable_color_features=True)
            >>> # Creates NH3 fragment at position (1.0, 2.0, 3.0)
        """
        if not features:
            raise ValueError("features list cannot be empty")
        
        if enable_color_features:
            # Create molecular fragments with proper bonding
            fragments = []
            
            for feature in features:
                if len(feature) < 5:
                    raise ValueError(
                        f"Invalid feature format. Expected at least 5 elements "
                        f"[type, atom_indices, x, y, z], got {len(feature)}"
                    )
                
                feat_type = feature[0]
                position = (float(feature[2]), float(feature[3]), float(feature[4]))
                
                # Get SMILES for this feature type
                smiles = cls.FEATURE_TO_SMILES.get(feat_type)
                
                if smiles is None:
                    raise ValueError(
                        f"Unknown feature type: {feat_type}. "
                        f"Valid types: {list(cls.FEATURE_TO_SMILES.keys())}"
                    )
                
                # Create fragment from SMILES
                fragment = cls._create_fragment_from_smiles(smiles, position, feat_type)
                fragments.append(fragment)
            
            # Combine all fragments into one molecule
            mol = cls._combine_fragments(fragments, name)
            
        else:
            # Create simple atoms-only molecule (shape-only mode)
            mol = cls._convert_atoms_only(features, name)
        
        return mol
    
    @classmethod
    def _convert_atoms_only(
        cls,
        features: List[List],
        name: str = 'Pharmacophore'
    ) -> Chem.Mol:
        """Convert features to atoms-only molecule (no bonds).
        
        This is the original implementation for shape-only alignment.
        
        Args:
            features: List of pharmacophore features.
            name: Name for the molecule.
            
        Returns:
            RDKit Mol with only atoms (no bonds).
        """
        # Create editable molecule
        mol = Chem.RWMol()
        
        # Add atoms for each feature
        for feature in features:
            if len(feature) < 5:
                raise ValueError(
                    f"Invalid feature format. Expected at least 5 elements "
                    f"[type, atom_indices, x, y, z], got {len(feature)}"
                )
            
            feat_type = feature[0]
            x, y, z = float(feature[2]), float(feature[3]), float(feature[4])
            
            # Get element symbol for this feature type
            element = cls.FEATURE_TO_ELEMENT.get(feat_type)
            
            if element is None:
                raise ValueError(
                    f"Unknown feature type: {feat_type}. "
                    f"Valid types: {list(cls.FEATURE_TO_ELEMENT.keys())}"
                )
            
            # Create atom
            atom = Chem.Atom(element)
            mol.AddAtom(atom)
        
        # Convert to non-editable Mol
        mol = mol.GetMol()
        
        # Create conformer with 3D coordinates
        conformer = Chem.Conformer(len(features))
        
        for idx, feature in enumerate(features):
            x, y, z = float(feature[2]), float(feature[3]), float(feature[4])
            conformer.SetAtomPosition(idx, rdGeometry.Point3D(x, y, z))
        
        # Add conformer to molecule
        mol.AddConformer(conformer, assignId=True)
        
        # Set molecule name
        mol.SetProp('_Name', name)
        
        return mol
    
    @classmethod
    def convert_with_metadata(
        cls,
        features: List[List],
        name: str = 'Pharmacophore',
        add_feature_labels: bool = True
    ) -> Chem.Mol:
        """Convert features to Mol with additional metadata.
        
        Similar to convert() but adds feature type information as atom
        properties for visualization and analysis.
        
        Args:
            features: List of pharmacophore features.
            name: Name to assign to the molecule.
            add_feature_labels: If True, add feature type as atom property
                (default: True).
        
        Returns:
            RDKit Mol object with atoms and metadata.
        
        Example:
            >>> features = [['Donor', (), 1.0, 2.0, 3.0]]
            >>> mol = PharmacophoreToMol.convert_with_metadata(features)
            >>> atom = mol.GetAtomWithIdx(0)
            >>> atom.GetProp('FeatureType')
            'Donor'
        """
        mol = cls.convert(features, name=name)
        
        if add_feature_labels:
            for idx, feature in enumerate(features):
                feat_type = feature[0]
                atom = mol.GetAtomWithIdx(idx)
                atom.SetProp('FeatureType', feat_type)
        
        return mol
    
    @classmethod
    def validate_for_shape_alignment(cls, mol: Chem.Mol) -> bool:
        """Validate that a molecule is suitable for shape alignment.
        
        Checks that the molecule has:
        1. At least one atom
        2. A 3D conformer
        3. Valid 3D coordinates
        
        Args:
            mol: RDKit Mol object to validate.
        
        Returns:
            True if molecule is valid for shape alignment.
        
        Raises:
            ValueError: If molecule is invalid with description of issue.
        
        Example:
            >>> features = [['Donor', (), 1.0, 2.0, 3.0]]
            >>> mol = PharmacophoreToMol.convert(features)
            >>> PharmacophoreToMol.validate_for_shape_alignment(mol)
            True
        """
        if mol is None:
            raise ValueError("Molecule is None")
        
        if mol.GetNumAtoms() == 0:
            raise ValueError("Molecule has no atoms")
        
        if not mol.GetNumConformers():
            raise ValueError("Molecule has no conformer")
        
        conformer = mol.GetConformer()
        
        if not conformer.Is3D():
            raise ValueError("Conformer is not 3D")
        
        # Check for valid coordinates
        for i in range(mol.GetNumAtoms()):
            pos = conformer.GetAtomPosition(i)
            if not all(abs(coord) < 1e6 for coord in [pos.x, pos.y, pos.z]):
                raise ValueError(
                    f"Atom {i} has invalid coordinates: "
                    f"({pos.x}, {pos.y}, {pos.z})"
                )
        
        return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
