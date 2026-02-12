"""Greedy forward feature selection for consensus pharmacophores.

Sequential forward selection (Whitney 1971) adapted for pharmacophore
feature subsets. Uses Pharm2D scoring (15 ms/mol) as the fast objective
for feature-subset evaluation, making the search O(N^2) instead of the
O(2^N) combinatorial explosion of exhaustive enumeration.

Feature selection applies to ``scoring_mode='consensus_mol'`` and
``scoring_mode='pharm2d'`` — in ``scoring_mode='reference'`` queries
align directly to reference molecules and consensus features are not
used for scoring.

Example:
    >>> from pharmacophore.greedy_selector import GreedyFeatureSelector
    >>> selector = GreedyFeatureSelector(ref_mols, actives, decoys)
    >>> result = selector.select(tolerance=2.0, occurrence=0.3)
    >>> print(f"Selected {len(result.selected_features)} features, AUC={result.best_auc:.4f}")
"""

from dataclasses import dataclass, field
from typing import List, Optional
import logging
import time
import itertools
import numpy as np
from rdkit import Chem

from .consensus import PharmacophoreConsensus
from .screening_metrics import calculate_all_metrics

logger = logging.getLogger(__name__)


@dataclass
class FeatureSelectionResult:
    """Result of greedy feature selection.

    Attributes:
        selected_features: Final selected pharmacophore features.
        selected_indices: Indices into the original consensus feature list.
        selection_history: List of (step, added_index, auc, bedroc) tuples.
        individual_scores: Per-feature AUC attribution from pair analysis.
        best_auc: AUC of the final selected subset.
        best_bedroc: BEDROC of the final selected subset.
        n_evaluations: Total Pharm2D evaluations performed.
        wall_time_sec: Total wall-clock time.
    """
    selected_features: List = field(default_factory=list)
    selected_indices: List[int] = field(default_factory=list)
    selection_history: List[tuple] = field(default_factory=list)
    individual_scores: dict = field(default_factory=dict)
    best_auc: float = 0.5
    best_bedroc: float = 0.0
    n_evaluations: int = 0
    wall_time_sec: float = 0.0


class GreedyFeatureSelector:
    """Greedy forward selection of pharmacophore features.

    Algorithm:
    1. Generate consensus pharmacophore from references.
    2. Score all C(N,2) feature pairs via Pharm2D.
    3. Attribute each pair's AUC to both features (max-pooling).
    4. Greedy forward: start with the best feature, add the one that
       improves AUC most at each step.
    5. Stop when improvement < convergence_threshold or max_features reached.

    Args:
        reference_mols: Aligned reference molecules.
        actives: Active compounds (Mol objects).
        decoys: Decoy compounds (Mol objects).
        random_state: Seed for reproducibility.
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        random_state: int = 42,
    ):
        from rdkit.Chem import AllChem

        self.reference_mols = reference_mols
        self.actives = actives
        self.decoys = decoys
        self.random_state = random_state

        self.y_true = [1] * len(actives) + [0] * len(decoys)

        # Ensure all molecules have 3D conformers for spatial scoring
        self.all_mols = []
        for mol in list(actives) + list(decoys):
            if mol.GetNumConformers() > 0:
                self.all_mols.append(mol)
            else:
                mol_h = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = random_state
                params.numThreads = 1
                AllChem.EmbedMolecule(mol_h, params=params)
                if mol_h.GetNumConformers() > 0:
                    self.all_mols.append(mol_h)
                else:
                    # Last resort: keep original (will score 0)
                    self.all_mols.append(mol)

    def _score_subset_spatial(
        self, features: List, all_mols: List[Chem.Mol],
    ) -> dict:
        """Score a feature subset using Pharm2D with feature-type masking.

        Computes Pharm2D fingerprints for references and queries, then
        masks out bits that involve feature types NOT present in the
        selected consensus subset. This makes Pharm2D scoring aware of
        which pharmacophore features we've selected.

        The Gobbi_Pharm2D factory encodes 6 feature types. We identify
        which types appear in the selected consensus features and only
        count fingerprint bits that involve those types.

        Args:
            features: Consensus feature subset to evaluate.
            all_mols: All molecules (actives + decoys).

        Returns:
            Dict of screening metrics.
        """
        from rdkit import DataStructs
        from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

        n_consensus = len(features)
        if n_consensus == 0:
            return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                    'ef_5': 0.0, 'ef_10': 0.0}

        factory = Gobbi_Pharm2D.factory

        # Map our feature types to Gobbi family abbreviations
        # AG=Anion, AR=Aromatic, BG=BasicGroup, HA=HBAcceptor,
        # HD=HBDonor, LH=LipophilicHeavy, RR=Ring, X=Other
        TYPE_TO_GOBBI = {
            'Donor': 'HD',
            'Acceptor': 'HA',
            'Aromatic': 'AR',
            'Hydrophobe': 'LH',
            'PosIonizable': 'BG',
            'NegIonizable': 'AG',
        }

        # Identify selected Gobbi families
        selected_gobbi = set()
        for f in features:
            gobbi = TYPE_TO_GOBBI.get(f[0])
            if gobbi:
                selected_gobbi.add(gobbi)

        # Build allowed bits: bits where ALL families are in selected set
        # Bit descriptions: "FAM1 FAM2 |bins|..." or "FAM1 FAM2 FAM3 |bins|..."
        allowed_bits = set()
        for bit_idx in range(factory.GetSigSize()):
            try:
                desc = factory.GetBitDescription(bit_idx)
                families = desc.split('|')[0].strip().split()
                if families and all(f in selected_gobbi for f in families):
                    allowed_bits.add(bit_idx)
            except Exception:
                continue

        if not allowed_bits:
            return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                    'ef_5': 0.0, 'ef_10': 0.0}

        # Generate reference fingerprints (masked)
        ref_on_bits = []
        for ref in self.reference_mols:
            try:
                fp = Generate.Gen2DFingerprint(ref, factory)
                masked = set(fp.GetOnBits()) & allowed_bits
                ref_on_bits.append(masked)
            except Exception:
                continue

        if not ref_on_bits:
            return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                    'ef_5': 0.0, 'ef_10': 0.0}

        # Score each molecule
        scores = []
        for mol in all_mols:
            try:
                fp = Generate.Gen2DFingerprint(mol, factory)
                mol_bits = set(fp.GetOnBits()) & allowed_bits

                # Max Tanimoto to any reference
                best_sim = 0.0
                for rb in ref_on_bits:
                    intersection = len(mol_bits & rb)
                    union = len(mol_bits | rb)
                    sim = intersection / union if union > 0 else 0.0
                    best_sim = max(best_sim, sim)
                scores.append(best_sim)
            except Exception:
                scores.append(0.0)

        try:
            return calculate_all_metrics(self.y_true, scores)
        except Exception:
            return {'roc_auc': 0.5, 'bedroc': 0.0, 'ef_1': 0.0,
                    'ef_5': 0.0, 'ef_10': 0.0}

    def select(
        self,
        tolerance: float = 2.0,
        occurrence: float = 0.3,
        max_features: int = 8,
        convergence_threshold: float = 0.001,
        verbose: bool = True,
    ) -> FeatureSelectionResult:
        """Run greedy forward feature selection.

        Args:
            tolerance: Consensus clustering spatial tolerance (Angstroms).
            occurrence: Minimum feature occurrence fraction.
            max_features: Maximum features to select.
            convergence_threshold: Stop when AUC improvement < this.
            verbose: Log progress.

        Returns:
            FeatureSelectionResult with selected features and history.
        """
        start_time = time.time()
        n_evals = 0

        # Step 1: Generate consensus pharmacophore
        consensus = PharmacophoreConsensus(
            tolerance=tolerance,
            occurrence_threshold=occurrence,
        )
        all_features = consensus.generate_consensus(self.reference_mols)

        if len(all_features) < 2:
            return FeatureSelectionResult(
                selected_features=all_features,
                selected_indices=list(range(len(all_features))),
                best_auc=0.5,
                wall_time_sec=time.time() - start_time,
            )

        n_features = len(all_features)
        if verbose:
            logger.info("Consensus has %d features, starting selection", n_features)

        # Step 2: Score all pairs to get per-feature attribution
        individual_scores = {i: 0.0 for i in range(n_features)}

        if n_features <= 2:
            # With <= 2 features, no pair analysis needed
            for i in range(n_features):
                individual_scores[i] = 0.5
        else:
            for i, j in itertools.combinations(range(n_features), 2):
                subset = [all_features[i], all_features[j]]
                metrics = self._score_subset_spatial(subset, self.all_mols)
                n_evals += 1
                auc = metrics['roc_auc']
                # Max-pool attribution
                individual_scores[i] = max(individual_scores[i], auc)
                individual_scores[j] = max(individual_scores[j], auc)

        # Step 3: Greedy forward selection
        # Start with the feature that has the highest individual score
        ranked = sorted(individual_scores.items(), key=lambda x: -x[1])
        selected_indices = [ranked[0][0]]
        remaining = set(range(n_features)) - set(selected_indices)

        # Evaluate starting feature
        current_subset = [all_features[i] for i in selected_indices]
        current_metrics = self._score_subset_spatial(
            current_subset, self.all_mols
        )
        n_evals += 1
        current_auc = current_metrics['roc_auc']
        current_bedroc = current_metrics.get('bedroc', 0.0)

        history = [(0, selected_indices[0], current_auc, current_bedroc)]

        if verbose:
            logger.info(
                "Step 0: feature %d, AUC=%.4f",
                selected_indices[0], current_auc,
            )

        # Greedy loop
        step = 1
        while remaining and len(selected_indices) < max_features:
            best_new_idx = None
            best_new_auc = current_auc

            for candidate in sorted(remaining):
                trial = selected_indices + [candidate]
                trial_features = [all_features[i] for i in trial]
                metrics = self._score_subset_spatial(
                    trial_features, self.all_mols
                )
                n_evals += 1
                trial_auc = metrics['roc_auc']

                if trial_auc > best_new_auc:
                    best_new_auc = trial_auc
                    best_new_idx = candidate
                    best_new_bedroc = metrics.get('bedroc', 0.0)

            # Check convergence
            improvement = best_new_auc - current_auc
            if best_new_idx is None or improvement < convergence_threshold:
                if verbose:
                    logger.info(
                        "Converged at step %d (improvement=%.6f < %.6f)",
                        step, improvement, convergence_threshold,
                    )
                break

            selected_indices.append(best_new_idx)
            remaining.discard(best_new_idx)
            current_auc = best_new_auc
            current_bedroc = best_new_bedroc
            history.append((step, best_new_idx, current_auc, current_bedroc))

            if verbose:
                logger.info(
                    "Step %d: +feature %d, AUC=%.4f (+%.4f)",
                    step, best_new_idx, current_auc, improvement,
                )
            step += 1

        selected_features = [all_features[i] for i in selected_indices]
        wall_time = time.time() - start_time

        return FeatureSelectionResult(
            selected_features=selected_features,
            selected_indices=selected_indices,
            selection_history=history,
            individual_scores=individual_scores,
            best_auc=current_auc,
            best_bedroc=current_bedroc,
            n_evaluations=n_evals,
            wall_time_sec=wall_time,
        )
