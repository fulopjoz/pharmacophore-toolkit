# PheSA: An Open-Source Tool for Pharmacophore-Enhanced Shape Alignment

**Authors**: J Wahl  
**Year**: 2024  
**Source**: Google Scholar search  
**Topic**: Shape+Color Scoring

---

## Summary

PheSA is an open-source pharmacophore-enhanced shape alignment tool that combines
shape-based 3D similarity with pharmacophore feature matching for virtual screening.
This is highly relevant to our toolkit's combo scoring approach.

## Key Findings

- Shape alignment alone (ROCS-like) is insufficient for selectivity
- Adding pharmacophore features ("color") dramatically improves enrichment
- Open-source implementation allows for customization and integration

## Relevance to Pharmacophore-Toolkit

**Direct Application**:
- Our `PharmacophoreToMol.convert()` with `enable_color_features=True` does similar
- Validates our fragment-based approach (NH3, benzene, CH4 for features)
- Confirms combo scoring (shape+color) is industry standard

**Implementation Comparison**:
| Feature | PheSA | Our Toolkit |
|---------|-------|-------------|
| Shape Alignment | Yes | Yes (AlignMol) |
| Color Features | Yes | Yes (fragments) |
| Deterministic | Unknown | Yes ✅ |
| RDKit Integration | Unclear | Yes |

## Methods

- Shape overlay using Gaussian functions
- Feature matching using pharmacophore definitions
- Combo scoring: weighted sum of shape and color Tanimoto

## Notes

- **TO DO**: Download paper and compare exact algorithm
- Check if their color implementation differs from ours
- Look for benchmark results on DUD-E or similar

---

**Status**: 📋 To review  
**Priority**: ⭐ High - directly comparable to our approach
