"""ECFP6 fingerprint classification baseline for CCR2 active/decoy discrimination.

Establishes a machine-learning baseline using standard molecular fingerprints
to contextualize the pharmacophore-based virtual screening approaches.

Setup:
  - Positives: 5 reference ligands + 74 known actives = 79
  - Negatives: 500 decoys = 500
  - Features: ECFP6 (Morgan radius=3, 2048 bits)
  - Models: Random Forest, XGBoost (default hyperparameters)
  - Evaluation:
      1. Stratified 5-fold CV (random split — kept for comparison)
      2. Scaffold-split 5-fold CV (GroupKFold on Murcko scaffolds)
      3. Scaffold holdout 80/20 (single deterministic split)

Design choices:
  - class_weight='balanced' / scale_pos_weight for imbalance (1:6.3 ratio)
  - Stratified folds preserve class ratio in every split
  - Scaffold splits ensure no Murcko scaffold appears in both train and test
  - No hyperparameter tuning — this is a floor estimate
  - Predicted probabilities used for ranking metrics (BEDROC, EF)

Usage:
    python experiments/ecfp6_baseline.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier, rdFingerprintGenerator
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import (
    roc_auc_score, accuracy_score, precision_score,
    recall_score, f1_score, matthews_corrcoef,
)
import xgboost as xgb

from pharmacophore.screening_metrics import bedroc, enrichment_factor

SEED = 42
N_BITS = 2048
RADIUS = 3  # ECFP6 = Morgan radius 3
N_FOLDS = 5


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dataset():
    """Load CCR2 references, actives, and decoys.

    Returns:
        all_smiles, labels, n_refs, n_actives, n_decoys, precomputed_murcko
        where precomputed_murcko is a list of scaffold SMILES (None for
        refs/decoys whose scaffolds are computed on the fly).
    """
    base = Path(__file__).parent.parent / 'tutorials' / 'data'

    # References from SDF → SMILES
    ref_smiles = []
    supplier = SDMolSupplier(str(base / 'CCR2_reference_ligands.sdf'), removeHs=False)
    for mol in supplier:
        if mol:
            ref_smiles.append(Chem.MolToSmiles(mol))

    # Actives (with pre-computed Murcko column)
    actives_df = pd.read_csv(base / 'actives_ccr2_N75.csv')
    active_smiles = actives_df['SMILES'].tolist()
    active_murcko = actives_df['SMILES_Murcko'].tolist()

    # Decoys
    decoys_df = pd.read_csv(base / 'decoys_ccr2_N500.csv')
    decoy_smiles = decoys_df['Smiles'].tolist()

    # Combine: refs + actives = positive, decoys = negative
    positive_smiles = ref_smiles + active_smiles
    all_smiles = positive_smiles + decoy_smiles
    labels = [1] * len(positive_smiles) + [0] * len(decoy_smiles)

    # Pre-computed Murcko: None for refs and decoys (computed on the fly)
    precomputed_murcko = [None] * len(ref_smiles) + active_murcko + [None] * len(decoy_smiles)

    return all_smiles, np.array(labels), len(ref_smiles), len(active_smiles), len(decoy_smiles), precomputed_murcko


# ---------------------------------------------------------------------------
# Scaffold utilities
# ---------------------------------------------------------------------------

def compute_scaffolds(smiles_list, precomputed_murcko=None):
    """Compute Murcko scaffolds for all molecules.

    Uses pre-computed values from the actives CSV where available.
    Falls back to RDKit MurckoScaffoldSmiles for refs/decoys.
    Uncomputable scaffolds get unique 'singleton_N' placeholders.
    """
    scaffolds = []
    singleton_count = 0

    for i, smi in enumerate(smiles_list):
        # Use pre-computed scaffold if available and not NaN
        if precomputed_murcko is not None and precomputed_murcko[i] is not None:
            val = precomputed_murcko[i]
            if isinstance(val, str) and val.strip():
                scaffolds.append(val)
                continue

        # Compute scaffold from SMILES
        try:
            scaf = MurckoScaffoldSmiles(smiles=smi, includeChirality=False)
            if scaf:
                scaffolds.append(scaf)
            else:
                scaffolds.append(f'singleton_{singleton_count}')
                singleton_count += 1
        except Exception:
            scaffolds.append(f'singleton_{singleton_count}')
            singleton_count += 1

    return scaffolds


def scaffold_group_ids(scaffolds):
    """Map scaffold SMILES to integer group IDs for GroupKFold."""
    scaffold_to_id = {}
    next_id = 0
    groups = []

    for scaf in scaffolds:
        if scaf not in scaffold_to_id:
            scaffold_to_id[scaf] = next_id
            next_id += 1
        groups.append(scaffold_to_id[scaf])

    return np.array(groups), scaffold_to_id


# ---------------------------------------------------------------------------
# Fingerprint generation
# ---------------------------------------------------------------------------

def smiles_to_ecfp6(smiles_list, radius=RADIUS, n_bits=N_BITS):
    """Generate ECFP6 fingerprints as numpy bit vectors.

    Returns:
        X: numpy array of shape (n_valid, n_bits)
        valid_mask: boolean array indicating which SMILES were valid
    """
    fpgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=n_bits,
    )

    fps = []
    valid_mask = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fp = fpgen.GetFingerprintAsNumPy(mol)
            fps.append(fp)
            valid_mask.append(True)
        else:
            valid_mask.append(False)

    X = np.array(fps, dtype=np.float32)
    return X, np.array(valid_mask)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

def get_models(n_pos, n_neg):
    """Return dict of classifiers with proper imbalance handling."""
    imbalance_ratio = n_neg / n_pos

    models = {
        'RandomForest': RandomForestClassifier(
            n_estimators=500,
            class_weight='balanced',
            random_state=SEED,
            n_jobs=-1,
        ),
        'XGBoost': xgb.XGBClassifier(
            n_estimators=500,
            scale_pos_weight=imbalance_ratio,
            eval_metric='logloss',
            random_state=SEED,
            n_jobs=-1,
            verbosity=0,
        ),
    }
    return models


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_fold(y_true, y_prob):
    """Compute all metrics for a single fold."""
    y_pred = (y_prob >= 0.5).astype(int)

    metrics = {
        'roc_auc': roc_auc_score(y_true, y_prob),
        'accuracy': accuracy_score(y_true, y_pred),
        'precision': precision_score(y_true, y_pred, zero_division=0),
        'recall': recall_score(y_true, y_pred, zero_division=0),
        'f1': f1_score(y_true, y_pred, zero_division=0),
        'mcc': matthews_corrcoef(y_true, y_pred),
        'bedroc': bedroc(y_true, y_prob),
        'ef_1': enrichment_factor(y_true, y_prob, percentage=0.01),
        'ef_5': enrichment_factor(y_true, y_prob, percentage=0.05),
    }
    return metrics


def run_cv(X, y, models, n_folds=N_FOLDS):
    """Run stratified k-fold CV for all models."""
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=SEED)

    all_results = {}

    for name, clf in models.items():
        fold_metrics = []

        for fold_i, (train_idx, test_idx) in enumerate(skf.split(X, y)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            clf_clone = type(clf)(**clf.get_params())
            clf_clone.fit(X_train, y_train)

            y_prob = clf_clone.predict_proba(X_test)[:, 1]
            metrics = evaluate_fold(y_test, y_prob)
            fold_metrics.append(metrics)

            print(f"  {name} fold {fold_i+1}/{n_folds}: "
                  f"AUC={metrics['roc_auc']:.4f}  "
                  f"BEDROC={metrics['bedroc']:.4f}  "
                  f"EF@1%={metrics['ef_1']:.2f}")

        # Aggregate
        agg = {}
        for key in fold_metrics[0]:
            vals = [m[key] for m in fold_metrics]
            agg[f'{key}_mean'] = np.mean(vals)
            agg[f'{key}_std'] = np.std(vals)
        agg['fold_details'] = fold_metrics

        all_results[name] = agg
        print()

    return all_results


def run_scaffold_cv(X, y, groups, models, n_folds=N_FOLDS):
    """Run scaffold-split k-fold CV using GroupKFold.

    GroupKFold guarantees no scaffold appears in both train and test.
    Folds with zero positives in test are skipped and reported.
    """
    gkf = GroupKFold(n_splits=n_folds)

    all_results = {}

    for name, clf in models.items():
        fold_metrics = []

        for fold_i, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            n_pos_test = int(y_test.sum())
            n_neg_test = len(y_test) - n_pos_test

            if n_pos_test == 0:
                print(f"  {name} fold {fold_i+1}/{n_folds}: "
                      f"SKIPPED (0 positives in test, {n_neg_test} negatives)")
                continue

            clf_clone = type(clf)(**clf.get_params())
            clf_clone.fit(X_train, y_train)

            y_prob = clf_clone.predict_proba(X_test)[:, 1]
            metrics = evaluate_fold(y_test, y_prob)
            fold_metrics.append(metrics)

            print(f"  {name} fold {fold_i+1}/{n_folds}: "
                  f"AUC={metrics['roc_auc']:.4f}  "
                  f"BEDROC={metrics['bedroc']:.4f}  "
                  f"EF@1%={metrics['ef_1']:.2f}  "
                  f"(test: {n_pos_test}+ / {n_neg_test}-)")

        if not fold_metrics:
            all_results[name] = {
                f'{k}_mean': 0.0 for k in ['roc_auc', 'accuracy', 'precision',
                                             'recall', 'f1', 'mcc', 'bedroc',
                                             'ef_1', 'ef_5']
            }
            all_results[name].update({
                f'{k}_std': 0.0 for k in ['roc_auc', 'accuracy', 'precision',
                                            'recall', 'f1', 'mcc', 'bedroc',
                                            'ef_1', 'ef_5']
            })
            all_results[name]['fold_details'] = []
            all_results[name]['n_valid_folds'] = 0
        else:
            agg = {}
            for key in fold_metrics[0]:
                vals = [m[key] for m in fold_metrics]
                agg[f'{key}_mean'] = np.mean(vals)
                agg[f'{key}_std'] = np.std(vals)
            agg['fold_details'] = fold_metrics
            agg['n_valid_folds'] = len(fold_metrics)
            all_results[name] = agg

        print()

    return all_results


def run_scaffold_holdout(X, y, groups, models, seed=SEED):
    """Single 80/20 scaffold-group holdout split.

    Greedily assigns scaffold groups to the test set until reaching ~20%.
    Deterministic via np.random.RandomState(seed).
    """
    rng = np.random.RandomState(seed)

    # Count samples per group
    unique_groups = np.unique(groups)
    group_sizes = {g: int(np.sum(groups == g)) for g in unique_groups}

    # Shuffle groups deterministically
    group_order = rng.permutation(unique_groups)

    target_test_size = int(0.2 * len(y))
    test_groups = set()
    current_test_size = 0

    for g in group_order:
        if current_test_size >= target_test_size:
            break
        test_groups.add(g)
        current_test_size += group_sizes[g]

    test_mask = np.array([g in test_groups for g in groups])
    train_mask = ~test_mask

    X_train, X_test = X[train_mask], X[test_mask]
    y_train, y_test = y[train_mask], y[test_mask]

    n_pos_test = int(y_test.sum())
    n_neg_test = len(y_test) - n_pos_test
    n_pos_train = int(y_train.sum())
    n_neg_train = len(y_train) - n_pos_train

    print(f"  Train: {len(y_train)} ({n_pos_train}+ / {n_neg_train}-)")
    print(f"  Test:  {len(y_test)} ({n_pos_test}+ / {n_neg_test}-)")
    print(f"  Test scaffold groups: {len(test_groups)}/{len(unique_groups)}")

    all_results = {}

    if n_pos_test == 0:
        print("  WARNING: 0 positives in test set — cannot compute metrics")
        for name in models:
            all_results[name] = {k: 0.0 for k in [
                'roc_auc', 'accuracy', 'precision', 'recall',
                'f1', 'mcc', 'bedroc', 'ef_1', 'ef_5']}
        return all_results

    for name, clf in models.items():
        clf_clone = type(clf)(**clf.get_params())
        clf_clone.fit(X_train, y_train)

        y_prob = clf_clone.predict_proba(X_test)[:, 1]
        metrics = evaluate_fold(y_test, y_prob)
        all_results[name] = metrics

        print(f"  {name}: AUC={metrics['roc_auc']:.4f}  "
              f"BEDROC={metrics['bedroc']:.4f}  "
              f"EF@1%={metrics['ef_1']:.2f}  "
              f"MCC={metrics['mcc']:.3f}")

    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"\n{'='*70}")
    print(f"ECFP6 CLASSIFICATION BASELINE — CCR2 Active vs Decoy")
    print(f"{'='*70}")

    # Load data
    all_smiles, labels, n_refs, n_actives, n_decoys, precomputed_murcko = load_dataset()
    n_pos = n_refs + n_actives
    n_neg = n_decoys

    print(f"  Positives: {n_refs} refs + {n_actives} actives = {n_pos}")
    print(f"  Negatives: {n_neg} decoys")
    print(f"  Imbalance ratio: 1:{n_neg/n_pos:.1f}")
    print(f"  Fingerprint: ECFP6 (Morgan r={RADIUS}, {N_BITS} bits)")
    print(f"  Seed: {SEED}")
    print(f"{'='*70}\n")

    # Generate fingerprints
    print("Generating ECFP6 fingerprints...")
    X, valid_mask = smiles_to_ecfp6(all_smiles)
    y = labels[valid_mask]
    print(f"  Valid molecules: {valid_mask.sum()}/{len(all_smiles)}")
    print(f"  Feature matrix: {X.shape}")
    print(f"  Bit density: {X.mean():.4f} ({X.sum(axis=1).mean():.0f} bits ON per molecule)\n")

    # Compute scaffolds (filter to valid molecules)
    valid_smiles = [s for s, v in zip(all_smiles, valid_mask) if v]
    valid_murcko = [m for m, v in zip(precomputed_murcko, valid_mask) if v]
    scaffolds = compute_scaffolds(valid_smiles, valid_murcko)
    groups, scaffold_to_id = scaffold_group_ids(scaffolds)
    n_unique = len(scaffold_to_id)
    print(f"Scaffold analysis: {n_unique} unique Murcko scaffolds\n")

    # Models
    models = get_models(n_pos, n_neg)

    # ===================================================================
    # EVALUATION 1: Random Stratified 5-Fold CV
    # ===================================================================
    print(f"{'='*70}")
    print(f"EVALUATION 1: Random Stratified {N_FOLDS}-Fold CV")
    print(f"{'='*70}\n")

    random_results = run_cv(X, y, models)

    # ===================================================================
    # EVALUATION 2: Scaffold-Split 5-Fold CV
    # ===================================================================
    print(f"{'='*70}")
    print(f"EVALUATION 2: Scaffold-Split {N_FOLDS}-Fold CV (GroupKFold)")
    print(f"{'='*70}\n")

    scaffold_results = run_scaffold_cv(X, y, groups, models)

    # ===================================================================
    # EVALUATION 3: Scaffold Holdout 80/20
    # ===================================================================
    print(f"{'='*70}")
    print(f"EVALUATION 3: Scaffold Holdout 80/20")
    print(f"{'='*70}\n")

    holdout_results = run_scaffold_holdout(X, y, groups, models)

    # ===================================================================
    # COMBINED RESULTS TABLE
    # ===================================================================
    print(f"\n{'='*70}")
    print(f"COMBINED RESULTS")
    print(f"{'='*70}\n")

    header = (f"{'Model':<15} {'Split':<15} {'ROC-AUC':>14} {'BEDROC':>14} {'MCC':>12}")
    print(header)
    print('-' * len(header))

    for name in models:
        # Random CV
        rr = random_results[name]
        print(f"{name:<15} {'Random CV':<15} "
              f"{rr['roc_auc_mean']:>5.4f}\u00b1{rr['roc_auc_std']:.3f} "
              f"{rr['bedroc_mean']:>5.4f}\u00b1{rr['bedroc_std']:.3f} "
              f"{rr['mcc_mean']:>5.3f}\u00b1{rr['mcc_std']:.2f}")

        # Scaffold CV
        sr = scaffold_results[name]
        n_valid = sr.get('n_valid_folds', N_FOLDS)
        print(f"{'':<15} {'Scaffold CV':<15} "
              f"{sr['roc_auc_mean']:>5.4f}\u00b1{sr['roc_auc_std']:.3f} "
              f"{sr['bedroc_mean']:>5.4f}\u00b1{sr['bedroc_std']:.3f} "
              f"{sr['mcc_mean']:>5.3f}\u00b1{sr['mcc_std']:.02f}"
              f"  ({n_valid}/{N_FOLDS} folds)")

        # Scaffold holdout
        hr = holdout_results[name]
        print(f"{'':<15} {'Scaffold Hold':<15} "
              f"{'':>2}{hr['roc_auc']:>5.4f}{'':>8} "
              f"{'':>2}{hr['bedroc']:>5.4f}{'':>8} "
              f"{'':>2}{hr['mcc']:>5.3f}")

        print()

    # ===================================================================
    # PERFORMANCE GAP ANALYSIS
    # ===================================================================
    print(f"{'='*70}")
    print(f"PERFORMANCE GAP: Random vs Scaffold")
    print(f"{'='*70}\n")

    for name in models:
        rand_auc = random_results[name]['roc_auc_mean']
        scaf_auc = scaffold_results[name]['roc_auc_mean']
        gap = rand_auc - scaf_auc

        label = "LEAK DETECTED" if gap > 0.05 else "MINOR" if gap > 0.01 else "NEGLIGIBLE"
        print(f"  {name}: Random AUC={rand_auc:.4f}, "
              f"Scaffold AUC={scaf_auc:.4f}, Gap={gap:+.4f} ({label})")

    print()

    # Context: compare with pharmacophore approaches
    print(f"{'='*70}")
    print(f"CONTEXT: Pharmacophore Approach Results (from benchmark)")
    print(f"{'='*70}\n")

    results_dir = Path(__file__).parent / 'results'
    benchmark_files = sorted(results_dir.glob('comparison_*.json'))
    if benchmark_files:
        with open(benchmark_files[-1]) as f:
            bench = json.load(f)

        print(f"{'Approach':<30} {'ROC-AUC':>8} {'BEDROC':>8}")
        print('-' * 50)
        for key, r in bench.items():
            print(f"{r['approach']:<30} {r['best_auc']:>8.4f} {r['best_bedroc']:>8.4f}")
        print('-' * 50)

    # Add ECFP6 rows for comparison
    for name in models:
        rr = random_results[name]
        sr = scaffold_results[name]
        print(f"{'ECFP6 + ' + name:<30} "
              f"{rr['roc_auc_mean']:>8.4f} {rr['bedroc_mean']:>8.4f}  (random CV)")
        print(f"{'  scaffold CV':<30} "
              f"{sr['roc_auc_mean']:>8.4f} {sr['bedroc_mean']:>8.4f}  (scaffold CV)")

    print()

    # Save results
    results_dir.mkdir(exist_ok=True)

    save_data = {
        'config': {
            'fingerprint': 'ECFP6',
            'radius': RADIUS,
            'n_bits': N_BITS,
            'n_folds': N_FOLDS,
            'seed': SEED,
            'n_positives': int(n_pos),
            'n_negatives': int(n_neg),
            'n_unique_scaffolds': n_unique,
        },
        'random_cv': {},
        'scaffold_cv': {},
        'scaffold_holdout': {},
    }

    for name in models:
        rr = random_results[name]
        save_data['random_cv'][name] = {
            k: float(v) for k, v in rr.items() if k != 'fold_details'
        }

        sr = scaffold_results[name]
        sr_save = {k: float(v) for k, v in sr.items()
                   if k not in ('fold_details', 'n_valid_folds')}
        sr_save['n_valid_folds'] = sr.get('n_valid_folds', N_FOLDS)
        save_data['scaffold_cv'][name] = sr_save

        hr = holdout_results[name]
        save_data['scaffold_holdout'][name] = {
            k: float(v) for k, v in hr.items()
        }

    json_file = results_dir / f'ecfp6_baseline_{timestamp}.json'
    with open(json_file, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"Results saved: {json_file}")


if __name__ == '__main__':
    main()
