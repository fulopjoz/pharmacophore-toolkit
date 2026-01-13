# Pharmacophore Toolkit - AI Agent Instructions

## Project Overview

**Pharmacophore-Toolkit** is a RDKit-based library for building pharmacophore models from crystal structures, docking poses, or SMILES strings. It generates 3D models (PyMOL `.pml` scripts) and 2D visualizations, with a major focus on **deterministic consensus pharmacophore generation** from multiple aligned molecules.

**Core Architecture**: Feature extraction → Clustering → Consensus model generation → RDKit molecule conversion for shape alignment

## Essential Knowledge

### 1. Dual Pharmacophore Systems (Critical to Understand!)

The codebase has **TWO DISTINCT** pharmacophore generation approaches:

#### Single-Molecule Pharmacophores (`pharmacophore.py`)
- **Purpose**: Extract features from ONE molecule
- **Main class**: `Pharmacophore`
- **Key methods**: `calc_pharm()`, `output_features()`, `draw_pharm()`
- **Features**: Donor, Acceptor, Aromatic, Hydrophobe (defined in `constants.py`)

#### Consensus Pharmacophores (`consensus.py` - NEW)
- **Purpose**: Find common features across MULTIPLE aligned molecules
- **Main class**: `PharmacophoreConsensus`
- **Algorithm**: Agglomerative hierarchical clustering (deterministic, NOT DBSCAN)
- **Parameters**: `tolerance` (spatial, Å) + `occurrence_threshold` (frequency, 0.0-1.0)
- **Output**: Features + RDKit Mol objects for shape-based alignment

**Why This Matters**: When modifying code, understand which system you're working with. They have different inputs, outputs, and use cases.

### 2. Feature Representation Format

All pharmacophore features use this standardized format:
```python
[feature_type, atom_indices, x, y, z]
# Examples:
['Donor', (5, 6), 1.234, 5.678, 9.012]  # Single molecule
['Acceptor', (), 2.345, 6.789, 1.234]    # Consensus (no atom mapping)
```

- **Tuple indices**: Empty `()` for consensus features (centroid positions)
- **Coordinates**: Always 3D Cartesian (Angstroms)
- **Types**: Match keys in `constants.FEATURES` and `constants.FEATURE_COLORS`

### 3. Consensus Model Generation Workflow

```
Aligned 3D molecules → Extract features → Group by type → 
Cluster spatially (AgglomerativeClustering) → Calculate centroids →
Filter by occurrence → Convert to RDKit Mol
```

**Key Implementation Files**:
- `pharmacophore/consensus.py`: Clustering logic
- `pharmacophore/mol_converter.py`: Feature → RDKit Mol conversion
- `pharmacophore/pharmacophore.py`: User-facing API (`generate_consensus_models()`)

### 4. Determinism is Non-Negotiable

The consensus system was specifically redesigned to be **100% deterministic**:
- Uses `AgglomerativeClustering` (NOT DBSCAN which was non-deterministic)
- Same input → Same output every time
- Critical test: `test_determinism()` in `tests/test_consensus.py`

**When modifying clustering code**: Always run determinism tests. Any randomness breaks the core design.

### 5. RDKit Mol Conversion (Shape Alignment)

**Why it exists**: RDKit's shape tools (`ShapeProtrudeDist`, `AlignMol`) require Mol objects, not coordinate lists.

**How it works** (`mol_converter.py`):
- **NEW (Oct 2025)**: Creates molecular fragments (NH3, benzene, CH4) at feature positions
- **Fragment-based approach**: Proper bonds and ring systems for SMARTS matching
- **Dual mode**:
  - `enable_color_features=True` (default): Molecular fragments with bonds → supports color scoring
  - `enable_color_features=False`: Atoms only (legacy) → shape-only alignment
- Element mapping: Donor→NH3, Acceptor→O, Aromatic→benzene, Hydrophobe→CH4
- **CRITICAL**: Sanitization initializes RingInfo required for SMARTS matching
- **Validate with**: `PharmacophoreToMol.validate_for_shape_alignment(mol)`

**Color Feature Support**:
```python
# NEW: Full shape+color alignment
features = [['Donor', (), 1.0, 2.0, 3.0]]
mol = PharmacophoreToMol.convert(features, enable_color_features=True)  
# Creates NH3 fragment → matches SMARTS [#7!H0]

from rdkit.Chem.rdShapeAlign import AlignMol
shape, color = AlignMol(ref=consensus_mol, probe=query_mol, useColors=True)
# Color score now works (0.0-1.0)
```

### 6. Interactive 3D Visualization (NEW - Oct 2025)

**Purpose**: Interactive Jupyter widget for visualizing aligned molecules with consensus pharmacophore overlay

**Implementation** (`pharmacophore/draw.py`):
- **Class**: `View`
- **Main method**: `view_consensus()`
- **Legacy method**: `view()` (backward compatible, single molecule dropdown)

**Key Features**:
- Multi-select molecule display with checkboxes
- Consensus model switching (strict/moderate/relaxed) via dropdown
- Consensus pharmacophore visibility toggle
- Real-time 3D rendering with py3Dmol
- Molecule color differentiation for overlays

**Usage Example**:
```python
from pharmacophore.draw import View

v = View()
v.view_consensus(
    mols=aligned_molecules,
    mol_names=['Serotonin', 'Dopamine', 'Norepinephrine'],
    consensus_models=models,
    default_model='moderate',
    labels=True,
    window=(900, 700)
)
```

**Widget Components**:
- **Molecule checkboxes**: Select which molecules to display
- **Model dropdown**: Switch between consensus model variants
- **Consensus toggle**: Show/hide consensus pharmacophore spheres
- **Rendering function**: `_render_consensus()` updates viewer on widget changes

**Technical Details**:
- Uses `ipywidgets` for interactive controls
- Uses `py3Dmol` for 3D molecular visualization
- Molecules rendered with distinct colors (gray, lightblue, lightgreen, etc.)
- Pharmacophore spheres color-coded by feature type
- Automatic view zooming to fit all visible objects

**Testing**: See `tests/test_interactive_view.py` for comprehensive widget tests

### 7. Shape and Color Tanimoto Scoring

**Combo Tanimoto Calculation**:
```python
shape_tanimoto, color_tanimoto = AlignMol(ref, probe, useColors=True)
combo_tanimoto = shape_tanimoto + color_tanimoto  # Simple sum, range 0-2
```

**Important**: Combo Tanimoto is NOT weighted by `opt_param`. It is always the simple sum of shape and color scores.

## Code Conventions

### Style Standards
- **PEP8 compliant** with Black formatter
- **Google Python Style Guide** for docstrings
- **Type hints required** for all function signatures
- **Comprehensive docstrings**: Args, Returns, Raises, Examples

### Feature Definition Pattern
When adding custom features:
```python
CUSTOM_FEATURES = {
    "FeatureType": ["smarts1", "smarts2", ...]
}
pharm = Pharmacophore(features=CUSTOM_FEATURES)
```

### API Design Pattern
- **Legacy method preserved**: `consensus_pharm()` exists but deprecated
- **New method**: `generate_consensus_models()` returns dict of {model_name: (features, mol)}
- **Always emit deprecation warnings** when keeping old APIs

## Testing Requirements

### Must-Run Tests Before Committing
```bash
pytest tests/test_consensus.py::TestPharmacophoreConsensus::test_determinism
pytest tests/test_consensus.py::TestPharmacophoreIntegration::test_determinism_across_runs
```

### Coverage Expectations
- All consensus code: >90% coverage
- **Parametric validation tests**: Check invalid inputs raise `ValueError`
- **Integration tests**: Test full pipeline (molecules → features → mol conversion)

## Development Workflows

### Adding New Feature Types
1. Update `constants.FEATURES` dict with SMARTS patterns
2. Update `constants.FEATURE_COLORS` for visualization
3. Update `mol_converter.FEATURE_TO_ELEMENT` mapping
4. Add test in `tests/test_consensus.py::test_feature_to_element_mapping`

### Modifying Clustering Algorithm
1. Changes must maintain determinism (verify with `test_determinism`)
2. Update `PharmacophoreConsensus._cluster_features()`
3. Valid linkage types: 'average', 'complete', 'single', 'ward'
4. **Do not** replace with DBSCAN or other non-deterministic methods

### Visualization Outputs
- **PyMOL**: `.pml` scripts via `output_features()` (sphere rendering)
- **2D images**: `draw_pharm()` in `draw.py` (uses CairoSVG)
- **Interactive 3D**: `View().view()` using py3Dmol in Jupyter notebooks
- **Tutorials**: Update Jupyter notebooks in `tutorials/` with working examples

## Common Pitfalls

1. **Forgetting 3D conformers**: Consensus requires molecules with `GetNumConformers() > 0`
2. **Mixing aligned/unaligned molecules**: Consensus assumes pre-aligned molecules (use `rdMolAlign.GetO3A`)
3. **Modifying feature list format**: Breaking `[type, indices, x, y, z]` structure breaks everything
4. **Using DBSCAN**: Was removed for non-determinism; use AgglomerativeClustering
5. **Hardcoded paths**: Demo files use `data/` directory; tutorials use `tutorials/data/`

## Key Files Reference

- `pharmacophore/pharmacophore.py`: Main API entry point
- `pharmacophore/consensus.py`: Deterministic clustering implementation
- `pharmacophore/mol_converter.py`: Feature-to-molecule conversion
- `pharmacophore/constants.py`: Feature definitions and color palettes
- `pharmacophore/draw.py`: 2D/3D visualization
- `tests/test_consensus.py`: Comprehensive test suite (28 tests)
- `tutorials/consensus_pharmacophore_tutorial.ipynb`: Complete usage examples

## Documentation Links

- Main README: Basic usage, installation, examples
- `README_CONSENSUS_UPGRADE.md`: Consensus system architecture
- `CONSENSUS_PHARMACOPHORE_IMPLEMENTATION_PLAN.md`: Technical design document (1341 lines!)
- `IMPLEMENTATION_COMPLETE.md`: Implementation summary and checklist
