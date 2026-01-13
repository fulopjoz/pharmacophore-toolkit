"""Tests for interactive consensus visualization in pharmacophore/draw.py."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from pharmacophore import Pharmacophore
from pharmacophore.draw import View
from pharmacophore.constants import INTERACTIVE_COLORS


@pytest.fixture
def sample_molecules():
    """Create sample aligned molecules for testing."""
    smiles_list = [
        "NCCc1c[nH]c2ccc(O)cc12",
        "NCCc1ccc(O)c(O)c1",
        "NCC(O)c1ccc(O)c(O)c1"
    ]
    
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    mols = [Chem.AddHs(m) for m in mols]
    
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    for m in mols:
        AllChem.EmbedMolecule(m, ps)
    
    for i in range(1, len(mols)):
        alignment = rdMolAlign.GetO3A(mols[i], mols[0])
        alignment.Align()
    
    return mols


@pytest.fixture
def sample_consensus_models(sample_molecules):
    """Generate consensus models for testing."""
    pharm = Pharmacophore()
    models = pharm.generate_consensus_models(
        mols=sample_molecules,
        tolerance=2.0,
        occurrence_threshold=0.5,
        linkage='average',
        model_set='standard'
    )
    return models


class TestViewConsensusInitialization:
    """Test View class initialization with consensus models."""
    
    def test_view_init_basic(self):
        """Test basic View initialization."""
        v = View()
        assert v.mol is None
        assert v.pharmacophore == 'default'
        assert v.consensus_models is None
    
    def test_view_init_with_consensus_models(self, sample_consensus_models):
        """Test View initialization with consensus models."""
        v = View(consensus_models=sample_consensus_models)
        assert v.consensus_models == sample_consensus_models
        assert 'strict' in v.consensus_models
        assert 'moderate' in v.consensus_models
        assert 'relaxed' in v.consensus_models


class TestViewConsensusParameters:
    """Test parameter validation for view_consensus method."""
    
    def test_view_consensus_requires_consensus_models(
        self,
        sample_molecules
    ):
        """Test that consensus_models must be provided."""
        v = View()
        with pytest.raises(ValueError, match="consensus_models must be"):
            v.view_consensus(mols=sample_molecules)
    
    def test_view_consensus_mol_names_length_mismatch(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test error when mol_names length doesn't match mols."""
        v = View()
        with pytest.raises(ValueError, match="mol_names length"):
            v.view_consensus(
                mols=sample_molecules,
                mol_names=['Mol1', 'Mol2'],
                consensus_models=sample_consensus_models
            )
    
    def test_view_consensus_default_mol_names(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that default mol_names are generated correctly."""
        v = View()
        
        v.mols_consensus = None
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models
            )
        except Exception:
            pass
        
        expected_names = [
            "Molecule 1", "Molecule 2", "Molecule 3"
        ]
        assert v.mol_names_consensus == expected_names
    
    def test_view_consensus_custom_mol_names(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test custom mol_names are stored correctly."""
        v = View()
        custom_names = ['Serotonin', 'Dopamine', 'Norepinephrine']
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                mol_names=custom_names,
                consensus_models=sample_consensus_models
            )
        except Exception:
            pass
        
        assert v.mol_names_consensus == custom_names


class TestViewConsensusWidgetCreation:
    """Test widget creation in view_consensus."""
    
    def test_checkboxes_created_for_molecules(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that checkboxes are created for each molecule."""
        v = View()
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models
            )
        except Exception:
            pass
        
        assert len(v.mol_checkboxes_consensus) == len(sample_molecules)
        
        for checkbox in v.mol_checkboxes_consensus:
            assert checkbox.value is True
    
    def test_model_dropdown_created(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that model dropdown is created with correct options."""
        v = View()
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models,
                default_model='moderate'
            )
        except Exception:
            pass
        
        dropdown = v.model_dropdown_consensus
        assert set(dropdown.options) == {'strict', 'moderate', 'relaxed'}
        assert dropdown.value == 'moderate'
    
    def test_consensus_visibility_checkbox_created(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that consensus visibility checkbox is created."""
        v = View()
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models
            )
        except Exception:
            pass
        
        assert v.show_consensus_checkbox.value is True


class TestViewConsensusColorHandling:
    """Test color handling in view_consensus."""
    
    def test_default_colors_used(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that default colors are used when none provided."""
        v = View()
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models
            )
        except Exception:
            pass
        
        assert v.colors_consensus == INTERACTIVE_COLORS
    
    def test_custom_colors_converted(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that custom colors are properly converted."""
        v = View()
        custom_colors = {
            'Donor': 'red',
            'Acceptor': 'blue',
            'Aromatic': 'green',
            'Hydrophobe': 'yellow'
        }
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models,
                color=custom_colors
            )
        except Exception:
            pass
        
        assert 'Donor' in v.colors_consensus
        assert isinstance(v.colors_consensus['Donor'], (str, tuple))


class TestViewConsensusStateManagement:
    """Test that view_consensus properly manages state."""
    
    def test_instance_variables_stored(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test that all parameters are stored as instance variables."""
        v = View()
        
        try:
            v.view_consensus(
                mols=sample_molecules,
                consensus_models=sample_consensus_models,
                labels=False,
                window=(1000, 800)
            )
        except Exception:
            pass
        
        assert v.mols_consensus == sample_molecules
        assert v.consensus_models_consensus == sample_consensus_models
        assert v.labels_consensus is False
        assert v.window_consensus == (1000, 800)
    
    def test_consensus_models_from_instance(
        self,
        sample_molecules,
        sample_consensus_models
    ):
        """Test using consensus_models from instance variable."""
        v = View(consensus_models=sample_consensus_models)
        
        try:
            v.view_consensus(mols=sample_molecules)
        except Exception:
            pass
        
        assert v.consensus_models_consensus == sample_consensus_models


class TestViewBackwardCompatibility:
    """Test that old view() method still works."""
    
    def test_old_view_method_exists(self):
        """Test that view() method still exists."""
        v = View()
        assert hasattr(v, 'view')
        assert callable(v.view)
    
    def test_old_render_method_exists(self):
        """Test that _render() method still exists."""
        v = View()
        assert hasattr(v, '_render')
        assert callable(v._render)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
