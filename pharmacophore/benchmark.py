"""Benchmarking framework for pharmacophore-based virtual screening.

This module provides tools for systematically comparing different
screening methods, measuring performance, and generating visualizations.

Example:
    >>> from pharmacophore.benchmark import ScreeningBenchmark
    >>>
    >>> benchmark = ScreeningBenchmark(reference_mols, actives, decoys)
    >>> benchmark.run_method('combo', optimizer.score_molecule)
    >>> benchmark.run_method('shape_only', shape_scorer)
    >>> print(benchmark.comparison_table())
"""

import time
import warnings
import numpy as np
import pandas as pd
from typing import List, Dict, Callable, Optional, Tuple, Any
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdShapeAlign import AlignMol

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from .screening_metrics import calculate_all_metrics, enrichment_factor
from .mol_converter import PharmacophoreToMol
from .consensus import PharmacophoreConsensus
from .pharm2d_scoring import Pharm2DScorer


class ScreeningBenchmark:
    """Benchmark framework for comparing virtual screening methods.

    This class provides infrastructure for:
    - Running multiple screening methods with timing
    - Calculating comprehensive metrics for each method
    - Generating comparison tables
    - Creating visualization plots

    Attributes:
        reference_mols: List of reference molecules for consensus.
        actives: List of active compounds.
        decoys: List of decoy compounds.
        results: Dict of method results.

    Example:
        >>> benchmark = ScreeningBenchmark(refs, actives, decoys)
        >>> benchmark.run_combo_scoring()
        >>> benchmark.run_shape_only_scoring()
        >>> table = benchmark.comparison_table()
        >>> print(table)
    """

    def __init__(
        self,
        reference_mols: List[Chem.Mol],
        actives: List[Chem.Mol],
        decoys: List[Chem.Mol],
        n_conformers: int = 10,
        random_state: int = 42,
        max_align_iters: int = 50,
        use_feature_weights: bool = False
    ):
        """Initialize benchmark with molecules.

        Args:
            reference_mols: Aligned reference molecules for consensus.
            actives: Active compounds to score.
            decoys: Decoy compounds to score.
            n_conformers: Conformers per molecule (default: 10).
            random_state: Random seed for reproducibility.
            max_align_iters: Maximum AlignMol iterations (default: 50).
                20=fast, 50=balanced, 300=accurate (Langer HAT).
            use_feature_weights: Weight features by importance (default: False).
                Donor/Acceptor=2.0x, Aromatic=1.5x, Hydrophobe=1.0x.
        """
        self.reference_mols = reference_mols
        self.actives = actives
        self.decoys = decoys
        self.n_conformers = n_conformers
        self.random_state = random_state
        self.max_align_iters = max_align_iters
        self.use_feature_weights = use_feature_weights

        self.results: Dict[str, Dict] = {}
        self._conformer_cache: Dict[str, List[Chem.Mol]] = {}
        self._pharmacophore_mol: Optional[Chem.Mol] = None

        # Precompute labels
        self.y_true = [1] * len(actives) + [0] * len(decoys)
        self.all_mols = actives + decoys

    def _get_pharmacophore_mol(
        self,
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5
    ) -> Chem.Mol:
        """Generate or retrieve cached pharmacophore mol."""
        if self._pharmacophore_mol is None:
            consensus = PharmacophoreConsensus(
                tolerance=tolerance,
                occurrence_threshold=occurrence_threshold
            )
            features = consensus.generate_consensus(self.reference_mols)
            self._pharmacophore_mol = PharmacophoreToMol.convert(
                features,
                name='Consensus',
                enable_color_features=True
            )
        return self._pharmacophore_mol

    def _get_conformers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Get or generate conformers with caching."""
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = None

        if smiles and smiles in self._conformer_cache:
            return self._conformer_cache[smiles]

        conformers = []
        try:
            mol_h = Chem.AddHs(mol)
            conf_ids = AllChem.EmbedMultipleConfs(
                mol_h,
                numConfs=self.n_conformers,
                randomSeed=self.random_state,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )

            if not conf_ids:
                AllChem.EmbedMolecule(mol_h, randomSeed=self.random_state)
                if mol_h.GetNumConformers() > 0:
                    conformers.append(mol_h)
            else:
                for conf_id in conf_ids:
                    conf_mol = Chem.Mol(mol_h)
                    conf_mol.RemoveAllConformers()
                    conf_mol.AddConformer(mol_h.GetConformer(conf_id), assignId=True)
                    conformers.append(conf_mol)
        except Exception:
            if mol.GetNumConformers() > 0:
                conformers.append(mol)

        if smiles:
            self._conformer_cache[smiles] = conformers

        return conformers

    def _score_single_molecule(
        self,
        mol: Chem.Mol,
        pharm_mol: Chem.Mol,
        use_colors: bool = True,
        shape_weight: float = 0.5,
        color_weight: float = 0.5
    ) -> Tuple[float, float, float]:
        """Score a molecule against pharmacophore.

        Returns:
            Tuple of (shape_score, color_score, weighted_combo).
        """
        conformers = self._get_conformers(mol)

        if not conformers:
            return 0.0, 0.0, 0.0

        best_combo = 0.0
        best_shape = 0.0
        best_color = 0.0

        for conf_mol in conformers:
            try:
                shape, color = AlignMol(
                    ref=pharm_mol,
                    probe=conf_mol,
                    useColors=use_colors,
                    opt_param=0.5,
                    max_preiters=self.max_align_iters  # Control pre-optimization iterations
                )
                combo = shape_weight * shape + color_weight * color

                if combo > best_combo:
                    best_combo = combo
                    best_shape = shape
                    best_color = color
            except Exception:
                continue

        return best_shape, best_color, best_combo

    def run_method(
        self,
        name: str,
        scoring_func: Callable[[Chem.Mol], float],
        description: str = ""
    ) -> Dict[str, Any]:
        """Run a custom scoring method with timing.

        Args:
            name: Method name for results.
            scoring_func: Function that takes a molecule and returns a score.
            description: Optional description of the method.

        Returns:
            Dict with timing and metrics.
        """
        scores = []
        start_time = time.time()

        for mol in self.all_mols:
            try:
                score = scoring_func(mol)
            except Exception:
                score = 0.0
            scores.append(score)

        elapsed = time.time() - start_time

        # Calculate metrics
        metrics = calculate_all_metrics(self.y_true, scores)

        result = {
            'name': name,
            'description': description,
            'time_sec': elapsed,
            'time_per_mol_ms': elapsed / len(self.all_mols) * 1000,
            'scores': scores,
            **metrics
        }

        self.results[name] = result
        return result

    def run_combo_scoring(
        self,
        name: str = "Combo (shape+color)",
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5
    ) -> Dict[str, Any]:
        """Run combo scoring (shape + color, current default method)."""
        warnings.warn(
            "This method uses anti-discriminative PharmacophoreToMol scoring "
            "(AUC < 0.5). Use run_reference_alignment_scoring() or "
            "run_pharm2d_scoring() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        pharm_mol = self._get_pharmacophore_mol(tolerance, occurrence_threshold)

        def scorer(mol: Chem.Mol) -> float:
            _, _, combo = self._score_single_molecule(mol, pharm_mol, use_colors=True)
            return combo

        return self.run_method(name, scorer, "Shape + Color Tanimoto (0-2)")

    def run_shape_only_scoring(
        self,
        name: str = "Shape Only",
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5
    ) -> Dict[str, Any]:
        """Run shape-only scoring (no color features)."""
        warnings.warn(
            "This method uses anti-discriminative PharmacophoreToMol scoring "
            "(AUC < 0.5). Use run_reference_alignment_scoring() or "
            "run_pharm2d_scoring() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        pharm_mol = self._get_pharmacophore_mol(tolerance, occurrence_threshold)

        def scorer(mol: Chem.Mol) -> float:
            shape, _, _ = self._score_single_molecule(mol, pharm_mol, use_colors=False)
            return shape

        return self.run_method(name, scorer, "Pure 3D shape similarity (0-1)")

    def run_color_weighted_scoring(
        self,
        name: str = "Color Weighted (0.3/0.7)",
        shape_weight: float = 0.3,
        color_weight: float = 0.7,
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5
    ) -> Dict[str, Any]:
        """Run color-weighted scoring."""
        warnings.warn(
            "This method uses anti-discriminative PharmacophoreToMol scoring "
            "(AUC < 0.5). Use run_reference_alignment_scoring() or "
            "run_pharm2d_scoring() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        pharm_mol = self._get_pharmacophore_mol(tolerance, occurrence_threshold)

        def scorer(mol: Chem.Mol) -> float:
            _, _, combo = self._score_single_molecule(
                mol, pharm_mol, use_colors=True,
                shape_weight=shape_weight, color_weight=color_weight
            )
            return combo

        desc = f"Weighted: {shape_weight}×shape + {color_weight}×color"
        return self.run_method(name, scorer, desc)

    def _score_against_references(
        self,
        mol: Chem.Mol,
        use_colors: bool = True,
        shape_weight: float = 0.5,
        color_weight: float = 0.5
    ) -> float:
        """Score molecule by aligning to best-matching reference ligand.

        This method aligns query molecules to the actual reference ligands
        (connected molecules) rather than the disconnected pharmacophore mol.
        This provides meaningful shape overlap for alignment.

        Args:
            mol: Query molecule to score.
            use_colors: Whether to use color (pharmacophore) features.
            shape_weight: Weight for shape component.
            color_weight: Weight for color component.

        Returns:
            Best combo score across all references and conformers.
        """
        conformers = self._get_conformers(mol)

        if not conformers:
            return 0.0

        best_score = 0.0

        for ref_mol in self.reference_mols:
            for conf_mol in conformers:
                try:
                    shape, color = AlignMol(
                        ref=ref_mol,
                        probe=conf_mol,
                        useColors=use_colors,
                        opt_param=0.5
                    )
                    combo = shape_weight * shape + color_weight * color

                    if combo > best_score:
                        best_score = combo
                except Exception:
                    continue

        return best_score

    def run_reference_alignment_scoring(
        self,
        name: str = "Reference Alignment",
        use_colors: bool = True,
        shape_weight: float = 0.5,
        color_weight: float = 0.5
    ) -> Dict[str, Any]:
        """Run scoring by aligning to reference ligands.

        This is the RECOMMENDED method. It aligns query molecules to the
        actual reference ligands (connected molecules) rather than the
        disconnected pharmacophore mol, providing meaningful shape overlap.

        Args:
            name: Method name for results.
            use_colors: Whether to use color (pharmacophore) features.
            shape_weight: Weight for shape component (0.5 default).
            color_weight: Weight for color component (0.5 default).

        Returns:
            Dict with timing and metrics.
        """
        def scorer(mol: Chem.Mol) -> float:
            return self._score_against_references(
                mol, use_colors, shape_weight, color_weight
            )

        desc = "Align to reference ligands (recommended)"
        return self.run_method(name, scorer, desc)

    def run_pharm2d_scoring(
        self,
        name: str = "Pharm2D Tanimoto",
        method: str = "tanimoto"
    ) -> Dict[str, Any]:
        """Run Pharm2D fingerprint-based scoring.

        This is the BEST performing method for discrimination. Pharm2D
        fingerprints encode pharmacophore feature PAIRS with binned distances,
        providing excellent separation of actives from property-matched decoys.

        Expected AUC: ~0.90 (vs 0.47 for shape-only methods)

        Args:
            name: Method name for results.
            method: Scoring method - 'tanimoto' (default) or 'consensus'.
                - 'tanimoto': Score by max similarity to any reference
                - 'consensus': Score against consensus bits (present in 60%+ refs)

        Returns:
            Dict with timing and metrics.
        """
        # Initialize Pharm2D scorer with reference molecules
        pharm2d_scorer = Pharm2DScorer(self.reference_mols)

        if method == "consensus":
            def scorer(mol: Chem.Mol) -> float:
                return pharm2d_scorer.score_consensus(mol)
            desc = "Pharm2D fingerprint - consensus bits"
        else:
            def scorer(mol: Chem.Mol) -> float:
                return pharm2d_scorer.score(mol)
            desc = "Pharm2D fingerprint - best reference Tanimoto"

        return self.run_method(name, scorer, desc)

    def run_all_standard_methods(
        self,
        tolerance: float = 2.0,
        occurrence_threshold: float = 0.5,
        verbose: bool = True,
        include_reference_alignment: bool = True,
        include_pharm2d: bool = True,
        include_deprecated: bool = False
    ) -> Dict[str, Dict]:
        """Run all standard scoring methods.

        Args:
            tolerance: Consensus tolerance parameter.
            occurrence_threshold: Consensus occurrence parameter.
            verbose: Print progress messages.
            include_reference_alignment: Include reference alignment method
                (recommended, default True).
            include_pharm2d: Include Pharm2D fingerprint method
                (BEST performance, default True).
            include_deprecated: Include deprecated methods (combo, shape,
                color_weighted) that have anti-discriminative AUC < 0.5
                (default False).

        Returns:
            Dict of all results.
        """
        # Build method list
        methods = []

        # Pharm2D (BEST performance) - run first
        if include_pharm2d:
            methods.append(("Pharm2D Tanimoto", "pharm2d"))

        # Reference alignment (recommended) - run second
        if include_reference_alignment:
            methods.append(("Reference Alignment", "ref_align"))

        # Deprecated methods (anti-discriminative, AUC < 0.5)
        if include_deprecated:
            methods.extend([
                ("Combo (shape+color)", "combo"),
                ("Shape Only", "shape"),
                ("Color Weighted (0.3/0.7)", "color_weighted"),
            ])

        for name, method_type in methods:
            if verbose:
                print(f"Running {name}...", end=" ", flush=True)

            if method_type == "pharm2d":
                self.run_pharm2d_scoring(name)
            elif method_type == "ref_align":
                self.run_reference_alignment_scoring(name)
            elif method_type == "combo":
                self.run_combo_scoring(name, tolerance, occurrence_threshold)
            elif method_type == "shape":
                self.run_shape_only_scoring(name, tolerance, occurrence_threshold)
            elif method_type == "color_weighted":
                self.run_color_weighted_scoring(name)

            if verbose:
                result = self.results[name]
                print(f"AUC={result['roc_auc']:.3f}, "
                      f"EF@1%={result['ef_1']:.1f}, "
                      f"time={result['time_per_mol_ms']:.1f}ms/mol")

        return self.results

    def comparison_table(self) -> pd.DataFrame:
        """Generate comparison table of all methods.

        Returns:
            DataFrame with one row per method and columns for metrics.
        """
        if not self.results:
            return pd.DataFrame()

        rows = []
        for name, result in self.results.items():
            rows.append({
                'Method': name,
                'ROC-AUC': result['roc_auc'],
                'BEDROC': result['bedroc'],
                'EF@1%': result['ef_1'],
                'EF@5%': result['ef_5'],
                'EF@10%': result['ef_10'],
                'Youden J': result['youden_j'],
                'Threshold': result['optimal_threshold'],
                'Sensitivity': result['sensitivity'],
                'Specificity': result['specificity'],
                'Time (ms/mol)': result['time_per_mol_ms']
            })

        df = pd.DataFrame(rows)
        df = df.sort_values('ROC-AUC', ascending=False).reset_index(drop=True)
        return df

    def plot_roc_curves(
        self,
        methods: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (8, 6)
    ) -> Optional['Figure']:
        """Plot ROC curves for all or selected methods.

        Args:
            methods: List of method names to plot (None = all).
            figsize: Figure size in inches.

        Returns:
            Matplotlib Figure or None if matplotlib unavailable.
        """
        if not HAS_MATPLOTLIB:
            print("Matplotlib not available for plotting")
            return None

        from sklearn.metrics import roc_curve

        if methods is None:
            methods = list(self.results.keys())

        fig, ax = plt.subplots(figsize=figsize)

        colors = plt.cm.tab10(np.linspace(0, 1, len(methods)))

        for i, name in enumerate(methods):
            if name not in self.results:
                continue

            result = self.results[name]
            scores = result['scores']

            fpr, tpr, _ = roc_curve(self.y_true, scores)
            auc = result['roc_auc']

            ax.plot(fpr, tpr, color=colors[i],
                    label=f"{name} (AUC={auc:.3f})", linewidth=2)

        ax.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')
        ax.set_xlabel('False Positive Rate', fontsize=12)
        ax.set_ylabel('True Positive Rate', fontsize=12)
        ax.set_title('ROC Curves Comparison', fontsize=14)
        ax.legend(loc='lower right')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def plot_enrichment_curves(
        self,
        methods: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (8, 6)
    ) -> Optional['Figure']:
        """Plot enrichment curves (cumulative actives vs % screened).

        Args:
            methods: List of method names to plot (None = all).
            figsize: Figure size in inches.

        Returns:
            Matplotlib Figure or None if matplotlib unavailable.
        """
        if not HAS_MATPLOTLIB:
            return None

        if methods is None:
            methods = list(self.results.keys())

        fig, ax = plt.subplots(figsize=figsize)

        n_total = len(self.y_true)
        n_actives = sum(self.y_true)
        percentages = np.linspace(0, 100, 101)

        colors = plt.cm.tab10(np.linspace(0, 1, len(methods)))

        for i, name in enumerate(methods):
            if name not in self.results:
                continue

            result = self.results[name]
            scores = result['scores']

            # Sort by score descending
            sorted_indices = np.argsort(scores)[::-1]
            sorted_labels = np.array(self.y_true)[sorted_indices]

            # Cumulative actives
            cumulative_actives = np.cumsum(sorted_labels)
            cumulative_pct = cumulative_actives / n_actives * 100

            # Sample at percentage points
            x_vals = np.arange(1, n_total + 1) / n_total * 100
            ax.plot(x_vals, cumulative_pct, color=colors[i],
                    label=name, linewidth=2)

        # Random baseline
        ax.plot([0, 100], [0, 100], 'k--', linewidth=1, label='Random')

        ax.set_xlabel('% of Library Screened', fontsize=12)
        ax.set_ylabel('% of Actives Found', fontsize=12)
        ax.set_title('Enrichment Curves', fontsize=14)
        ax.legend(loc='lower right')
        ax.set_xlim([0, 100])
        ax.set_ylim([0, 100])
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def plot_score_distributions(
        self,
        methods: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (12, 5)
    ) -> Optional['Figure']:
        """Plot score distributions for actives vs decoys.

        Args:
            methods: List of method names to plot (None = all).
            figsize: Figure size in inches.

        Returns:
            Matplotlib Figure or None if matplotlib unavailable.
        """
        if not HAS_MATPLOTLIB:
            return None

        if methods is None:
            methods = list(self.results.keys())

        n_methods = len(methods)
        if n_methods == 0:
            return None

        fig, axes = plt.subplots(1, n_methods, figsize=figsize)
        if n_methods == 1:
            axes = [axes]

        n_actives = sum(self.y_true)

        for i, name in enumerate(methods):
            if name not in self.results:
                continue

            result = self.results[name]
            scores = np.array(result['scores'])

            active_scores = scores[:n_actives]
            decoy_scores = scores[n_actives:]

            ax = axes[i]
            parts = ax.violinplot([active_scores, decoy_scores],
                                   positions=[1, 2], showmeans=True)

            # Color the violins
            parts['bodies'][0].set_facecolor('green')
            parts['bodies'][0].set_alpha(0.7)
            parts['bodies'][1].set_facecolor('red')
            parts['bodies'][1].set_alpha(0.7)

            # Add threshold line
            thresh = result['optimal_threshold']
            ax.axhline(y=thresh, color='blue', linestyle='--',
                      linewidth=1.5, label=f'Thresh={thresh:.2f}')

            ax.set_xticks([1, 2])
            ax.set_xticklabels(['Actives', 'Decoys'])
            ax.set_title(f'{name}\nAUC={result["roc_auc"]:.3f}', fontsize=10)
            ax.set_ylabel('Score')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig

    def plot_timing_comparison(
        self,
        methods: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (10, 5)
    ) -> Optional['Figure']:
        """Plot timing comparison bar chart.

        Args:
            methods: List of method names to plot (None = all).
            figsize: Figure size in inches.

        Returns:
            Matplotlib Figure or None if matplotlib unavailable.
        """
        if not HAS_MATPLOTLIB:
            return None

        if methods is None:
            methods = list(self.results.keys())

        fig, ax = plt.subplots(figsize=figsize)

        names = []
        times = []
        aucs = []

        for name in methods:
            if name not in self.results:
                continue
            result = self.results[name]
            names.append(name)
            times.append(result['time_per_mol_ms'])
            aucs.append(result['roc_auc'])

        x = np.arange(len(names))
        width = 0.35

        bars = ax.bar(x, times, width, color='steelblue')

        # Add AUC labels on bars
        for i, (bar, auc) in enumerate(zip(bars, aucs)):
            height = bar.get_height()
            ax.annotate(f'AUC={auc:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3), textcoords="offset points",
                       ha='center', va='bottom', fontsize=9)

        ax.set_ylabel('Time (ms/molecule)', fontsize=12)
        ax.set_title('Timing Comparison', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        return fig

    def plot_all(
        self,
        figsize: Tuple[int, int] = (14, 10)
    ) -> Optional['Figure']:
        """Create 4-panel figure with all visualizations.

        Args:
            figsize: Figure size in inches.

        Returns:
            Matplotlib Figure or None if matplotlib unavailable.
        """
        if not HAS_MATPLOTLIB:
            return None

        fig = plt.figure(figsize=figsize)

        # ROC curves
        ax1 = fig.add_subplot(2, 2, 1)
        self._plot_roc_on_ax(ax1)

        # Enrichment curves
        ax2 = fig.add_subplot(2, 2, 2)
        self._plot_enrichment_on_ax(ax2)

        # Score distributions (simplified)
        ax3 = fig.add_subplot(2, 2, 3)
        self._plot_distributions_on_ax(ax3)

        # Timing
        ax4 = fig.add_subplot(2, 2, 4)
        self._plot_timing_on_ax(ax4)

        plt.tight_layout()
        return fig

    def _plot_roc_on_ax(self, ax):
        """Plot ROC curves on given axes."""
        from sklearn.metrics import roc_curve

        colors = plt.cm.tab10(np.linspace(0, 1, len(self.results)))

        for i, (name, result) in enumerate(self.results.items()):
            fpr, tpr, _ = roc_curve(self.y_true, result['scores'])
            ax.plot(fpr, tpr, color=colors[i],
                    label=f"{name} ({result['roc_auc']:.3f})", linewidth=2)

        ax.plot([0, 1], [0, 1], 'k--', linewidth=1)
        ax.set_xlabel('FPR')
        ax.set_ylabel('TPR')
        ax.set_title('ROC Curves')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    def _plot_enrichment_on_ax(self, ax):
        """Plot enrichment curves on given axes."""
        n_total = len(self.y_true)
        n_actives = sum(self.y_true)
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.results)))

        for i, (name, result) in enumerate(self.results.items()):
            sorted_indices = np.argsort(result['scores'])[::-1]
            sorted_labels = np.array(self.y_true)[sorted_indices]
            cumulative = np.cumsum(sorted_labels) / n_actives * 100
            x_vals = np.arange(1, n_total + 1) / n_total * 100
            ax.plot(x_vals, cumulative, color=colors[i], label=name, linewidth=2)

        ax.plot([0, 100], [0, 100], 'k--', linewidth=1)
        ax.set_xlabel('% Screened')
        ax.set_ylabel('% Actives Found')
        ax.set_title('Enrichment Curves')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    def _plot_distributions_on_ax(self, ax):
        """Plot score distributions on given axes."""
        n_actives = sum(self.y_true)
        data = []
        labels = []

        for name, result in self.results.items():
            scores = np.array(result['scores'])
            active_mean = np.mean(scores[:n_actives])
            decoy_mean = np.mean(scores[n_actives:])
            data.append([active_mean, decoy_mean])
            labels.append(name)

        x = np.arange(len(labels))
        width = 0.35

        active_means = [d[0] for d in data]
        decoy_means = [d[1] for d in data]

        ax.bar(x - width/2, active_means, width, label='Actives', color='green', alpha=0.7)
        ax.bar(x + width/2, decoy_means, width, label='Decoys', color='red', alpha=0.7)

        ax.set_ylabel('Mean Score')
        ax.set_title('Score Distributions')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_timing_on_ax(self, ax):
        """Plot timing comparison on given axes."""
        names = list(self.results.keys())
        times = [self.results[n]['time_per_mol_ms'] for n in names]

        ax.barh(names, times, color='steelblue')
        ax.set_xlabel('Time (ms/mol)')
        ax.set_title('Timing Comparison')
        ax.grid(True, alpha=0.3, axis='x')

    def clear_cache(self):
        """Clear conformer cache to free memory."""
        self._conformer_cache.clear()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
