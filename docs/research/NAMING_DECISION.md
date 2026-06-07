# Method naming decision (2026-06-07)

**Chosen:** **PRISM** — *Pharmacophore-Resolved Importance-weighted Similarity Model* (was the internal codename `s3_3d`). Metaphor: a prism resolves white light into a spectrum, as the method resolves ROCS "color" into its 6 feature types and reweights them by learned actives-vs-decoys importance over consensus templates in 3D. To apply to **both code** (`s3_3d`→`prism`, `s3_3d_fixed`→`prism_fixed`, `s3_3d_esp`→`prism_esp`) **and paper**.

**Reserved backup:** **REFRACT** (*REweighted Feature-Resolved Adaptive Color Tanimoto*) — a zero-collision name in the same optical metaphor, to use if the comp-bio reuse of "PRISM" (protein-interaction matching / NRPS genome-mining) becomes a problem.

**Rejected:** SHARC (collides phonetically/visually with SHAFTS, the benchmarked tool, and "sharts"); CHROMA (taken — Generate Biomedicines' 2023 generative protein model); IRIS / MOSAIC (overloaded acronyms). Note: scite keyword search showed no PRISM/REFRACT-named tool dominates ligand-based shape VS, but it ranks by domain terms, not exact name — so this is a niche-collision check, not a formal prior-art search.
