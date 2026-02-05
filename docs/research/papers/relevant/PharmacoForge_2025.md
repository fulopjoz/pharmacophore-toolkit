# PharmacoForge: Pharmacophore Generation with Diffusion Models

**Authors**: EL Flynn, R Shah, I Dunn, R Aggarwal  
**Year**: 2025  
**Source**: Google Scholar search  
**Topic**: ML-based Pharmacophore Generation, Diffusion Models

---

## Summary

PharmacoForge uses diffusion models for pharmacophore generation, representing
a cutting-edge ML approach to the traditional clustering-based methods. This
could be an alternative to our deterministic agglomerative clustering.

## Key Findings

- Diffusion models can generate pharmacophore hypotheses from ligand sets
- Learns feature placement distributions rather than hard clustering
- Can generate multiple plausible pharmacophore models

## Relevance to Pharmacophore-Toolkit

**Potential Application**:
- Alternative to hierarchical clustering for consensus generation
- Could generate diverse pharmacophore hypotheses for ensemble screening
- May handle sparse/noisy data better than clustering

**Key Difference from Our Approach**:
| Aspect | Our Method | PharmacoForge |
|--------|------------|---------------|
| Approach | Deterministic clustering | Learned generation |
| Output | Single consensus | Multiple hypotheses |
| Determinism | ✅ Yes | ❌ Stochastic |
| Training | Not required | Requires training |
| Speed | Fast (ms) | Slower (inference) |

## Why We Keep Hierarchical Clustering

- **Determinism**: Non-negotiable for reproducibility
- **No Training**: Works out-of-box on any ligand set
- **Interpretable**: Clear spatial tolerance → cluster membership
- **Speed**: Essential for large-scale screening

## Methods (Inferred)

- Diffusion process: Noise → pharmacophore coordinates
- Conditioning: Input ligand features guide generation
- Sampling: Multiple hypotheses via different seeds

## Notes

- **Valuable for**: Research comparison, method benchmarking
- **Not for our core**: Violates determinism requirement
- Could use as validation: do our models match diffusion output?

---

**Status**: 📋 To review  
**Priority**: Medium - interesting but conflicts with our determinism principle
