#!/usr/bin/env env python
"""
Phase 2: Central Composite Design (CCD) for Response Surface Modeling

Based on Phase 1 findings:
- Tolerance: narrow to [0.8, 1.5] Å (center at 1.15)
- Threshold: narrow to [0.3, 0.6] (center at 0.45)
- Linkage: fix at 'complete' (best from Phase 1)

CCD Design:
- 2 factors (tolerance, threshold)
- 5 levels per factor (factorial + axial + center points)
- ~13 experiments + 5 center point replicates = 18 total runs
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, List, Tuple, Any
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from pyDOE3 import ccdesign
from scipy.optimize import minimize
import warnings

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from experiments.parameter_sweep import evaluate_model
from experiments.experiment_logger import ExperimentLogger


def create_ccd_design(center_tolerance=1.15, center_threshold=0.45, 
                      range_tolerance=0.35, range_threshold=0.15):
    """
    Create Central Composite Design for 2 factors.
    
    Args:
        center_tolerance: Center point for tolerance (Å)
        center_threshold: Center point for threshold
        range_tolerance: Half-range for tolerance (axial distance)
        range_threshold: Half-range for threshold (axial distance)
    
    Returns:
        DataFrame with experiment design matrix
    """
    # Generate CCD design (2 factors, circumscribed design)
    # Returns coded values in [-1, 1] range
    # center=(nf, nc) where nf=factorial points, nc=center point replicates
    design = ccdesign(2, center=(1, 5), alpha='orthogonal', face='ccc')
    
    # Convert to DataFrame
    df_design = pd.DataFrame(design, columns=['tolerance_coded', 'threshold_coded'])
    
    # Decode to actual parameter values
    df_design['tolerance'] = center_tolerance + df_design['tolerance_coded'] * range_tolerance
    df_design['threshold'] = center_threshold + df_design['threshold_coded'] * range_threshold
    
    # Clip to valid ranges
    df_design['tolerance'] = df_design['tolerance'].clip(0.8, 1.5)
    df_design['threshold'] = df_design['threshold'].clip(0.3, 0.6)
    
    # Add run order (randomized for blocking effects)
    df_design['run_order'] = np.random.permutation(len(df_design)) + 1
    df_design = df_design.sort_values('run_order').reset_index(drop=True)
    
    # Add experiment metadata
    df_design['linkage'] = 'complete'
    df_design['experiment_id'] = [f"CCD_{i+1:03d}" for i in range(len(df_design))]
    
    return df_design


def evaluate_model(
    tolerance: float,
    threshold: float,
    linkage: str = 'complete',
    n_refs: int = 5,
    n_actives: int = 74,
    n_decoys: int = 499,
    verbose: bool = False
) -> Dict[str, Any]:
    """
    Evaluate consensus model with given parameters on CCR2 dataset.
    
    Args:
        tolerance: Spatial clustering tolerance (Å)
        threshold: Occurrence threshold (0.0-1.0)
        linkage: Linkage method for clustering
        n_refs: Number of reference molecules to use
        n_actives: Number of actives to use
        n_decoys: Number of decoys to use
        verbose: Print detailed output
    
    Returns:
        Dict with performance metrics
    """
    from experiments.parameter_sweep import (
        load_ccr2_dataset, score_molecule, evaluate_model as eval_fn
    )
    
    # Load dataset
    ref_mols, actives, decoys = load_ccr2_dataset()
    
    # Subset if requested
    ref_mols = ref_mols[:n_refs]
    actives = actives[:n_actives]
    decoys = decoys[:n_decoys]
    
    # Evaluate
    return eval_fn(ref_mols, actives, decoys, tolerance, threshold, linkage)


def run_ccd_experiments(design_matrix, dataset='CCR2', n_refs=5, 
                        n_actives=74, n_decoys=499):
    """
    Execute CCD experiments and log results.
    
    Args:
        design_matrix: DataFrame from create_ccd_design()
        dataset: Dataset name
        n_refs: Number of reference molecules
        n_actives: Number of active compounds
        n_decoys: Number of decoy compounds
    
    Returns:
        DataFrame with results
    """
    # Initialize logger
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = ExperimentLogger(
        experiment_name=f"ccd_tolerance_threshold_{timestamp}",
        output_dir=Path(__file__).parent.parent / "docs/research/experiments/results"
    )
    
    results = []
    
    print(f"\n{'='*80}")
    print(f"PHASE 2: Central Composite Design (CCD) Experiments")
    print(f"{'='*80}")
    print(f"Total runs: {len(design_matrix)}")
    print(f"Dataset: {dataset} ({n_refs} refs, {n_actives} actives, {n_decoys} decoys)")
    print(f"Fixed parameter: linkage='complete'")
    print(f"{'='*80}\n")
    
    for idx, row in design_matrix.iterrows():
        tolerance = row['tolerance']
        threshold = row['threshold']
        linkage = row['linkage']
        exp_id = row['experiment_id']
        run_num = idx + 1
        
        print(f"\n[{run_num}/{len(design_matrix)}] {exp_id}")
        print(f"  Tolerance: {tolerance:.3f} Å")
        print(f"  Threshold: {threshold:.3f}")
        print(f"  Linkage: {linkage}")
        
        try:
            # Run evaluation
            metrics = evaluate_model(
                tolerance=tolerance,
                threshold=threshold,
                linkage=linkage,
                n_refs=n_refs,
                n_actives=n_actives,
                n_decoys=n_decoys,
                verbose=False
            )
            
            # Add parameters to metrics
            metrics['tolerance'] = tolerance
            metrics['threshold'] = threshold
            metrics['linkage'] = linkage
            metrics['experiment_id'] = exp_id
            metrics['run_order'] = row['run_order']
            metrics['tolerance_coded'] = row['tolerance_coded']
            metrics['threshold_coded'] = row['threshold_coded']
            
            # Log to experiment logger
            logger.log_run(
                run_id=run_num,
                parameters={'tolerance': tolerance, 'threshold': threshold, 'linkage': linkage},
                metrics=metrics,
                metadata={'experiment_id': exp_id, 'phase': 'Phase2_CCD'}
            )
            
            results.append(metrics)
            
            print(f"  ✓ ROC-AUC: {metrics['roc_auc']:.4f}")
            print(f"  ✓ EF@1%: {metrics.get('ef_1', 0.0):.2f}")
            print(f"  ✓ Features: {metrics['n_features']}")
            
        except Exception as e:
            print(f"  ✗ ERROR: {str(e)}")
            warnings.warn(f"Experiment {exp_id} failed: {e}")
            continue
    
    # Save results
    logger.save_csv()
    logger.save_json()
    
    df_results = pd.DataFrame(results)
    
    # Generate summary report
    summary = f"""# Phase 2: CCD Response Surface Results

## Experiment Overview
- **Date**: {datetime.now().strftime("%Y-%m-%d %H:%M")}
- **Design**: Central Composite Design (2 factors)
- **Total Runs**: {len(df_results)}/{len(design_matrix)}
- **Dataset**: {dataset}

## Parameter Ranges
- **Tolerance**: {df_results['tolerance'].min():.3f} - {df_results['tolerance'].max():.3f} Å
- **Threshold**: {df_results['threshold'].min():.3f} - {df_results['threshold'].max():.3f}
- **Linkage**: complete (fixed)

## Performance Summary
- **Best ROC-AUC**: {df_results['roc_auc'].max():.4f}
- **Median ROC-AUC**: {df_results['roc_auc'].median():.4f}
- **Worst ROC-AUC**: {df_results['roc_auc'].min():.4f}

## Top 5 Configurations
"""
    top5 = df_results.nlargest(5, 'roc_auc')[['experiment_id', 'tolerance', 'threshold', 
                                                'roc_auc', 'n_features']]
    if 'ef_1' in df_results.columns:
        top5 = df_results.nlargest(5, 'roc_auc')[['experiment_id', 'tolerance', 'threshold',
                                                    'roc_auc', 'ef_1', 'n_features']]
    summary += top5.to_markdown(index=False)
    
    logger.generate_summary_report()
    
    # Save custom summary
    summary_path = logger.results_dir / f"PHASE2_CCD_SUMMARY_{datetime.now().strftime('%Y%m%d')}.md"
    with open(summary_path, 'w') as f:
        f.write(summary)
    
    print(f"\n{'='*80}")
    print(f"✓ CCD experiments complete!")
    print(f"✓ Results saved to: {logger.results_csv}")
    print(f"✓ Summary saved to: {summary_path}")
    print(f"{'='*80}\n")
    
    return df_results


def fit_response_surface(df_results):
    """
    Fit quadratic response surface model with interaction terms.
    
    Model: AUC = β0 + β1*x1 + β2*x2 + β11*x1^2 + β22*x2^2 + β12*x1*x2
    where x1=tolerance, x2=threshold
    
    Args:
        df_results: DataFrame from run_ccd_experiments()
    
    Returns:
        dict with model coefficients and statistics
    """
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score, mean_squared_error
    
    # Prepare design matrix (coded values for better numerical stability)
    X = df_results[['tolerance_coded', 'threshold_coded']].values
    y = df_results['roc_auc'].values
    
    # Create quadratic terms
    X_quad = np.column_stack([
        np.ones(len(X)),           # Intercept
        X[:, 0],                    # x1 (tolerance)
        X[:, 1],                    # x2 (threshold)
        X[:, 0]**2,                 # x1^2
        X[:, 1]**2,                 # x2^2
        X[:, 0] * X[:, 1]           # x1*x2 (interaction)
    ])
    
    # Fit model
    model = LinearRegression(fit_intercept=False)
    model.fit(X_quad, y)
    
    # Predictions
    y_pred = model.predict(X_quad)
    
    # Statistics
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    adj_r2 = 1 - (1 - r2) * (len(y) - 1) / (len(y) - X_quad.shape[1] - 1)
    
    # Coefficient names
    coef_names = ['Intercept', 'Tolerance', 'Threshold', 
                  'Tolerance²', 'Threshold²', 'Tolerance×Threshold']
    
    # Build results
    results = {
        'coefficients': dict(zip(coef_names, model.coef_)),
        'r2': r2,
        'adj_r2': adj_r2,
        'rmse': rmse,
        'predictions': y_pred,
        'residuals': y - y_pred,
        'model': model
    }
    
    print(f"\n{'='*80}")
    print(f"RESPONSE SURFACE MODEL (Quadratic)")
    print(f"{'='*80}")
    print(f"R² = {r2:.4f}")
    print(f"Adjusted R² = {adj_r2:.4f}")
    print(f"RMSE = {rmse:.4f}")
    print(f"\nCoefficients (coded scale):")
    for name, coef in results['coefficients'].items():
        print(f"  {name:20s}: {coef:+.6f}")
    print(f"{'='*80}\n")
    
    return results


def optimize_response_surface(model_results, center_tolerance=1.15, 
                                center_threshold=0.45, 
                                range_tolerance=0.35, 
                                range_threshold=0.15):
    """
    Find optimum on response surface using numerical optimization.
    
    Args:
        model_results: dict from fit_response_surface()
        center_tolerance, center_threshold: CCD center points
        range_tolerance, range_threshold: CCD ranges
    
    Returns:
        dict with optimal parameters and predicted response
    """
    model = model_results['model']
    
    def objective(x_coded):
        """Negative AUC to minimize (optimization finds minimum)"""
        X_quad = np.array([[1, x_coded[0], x_coded[1], 
                            x_coded[0]**2, x_coded[1]**2, 
                            x_coded[0]*x_coded[1]]])
        return -model.predict(X_quad)[0]
    
    # Optimize in coded space [-1.5, 1.5] (slightly beyond factorial points)
    result = minimize(objective, x0=[0, 0], method='L-BFGS-B',
                      bounds=[(-1.5, 1.5), (-1.5, 1.5)])
    
    # Decode to actual parameters
    optimal_coded = result.x
    optimal_tolerance = center_tolerance + optimal_coded[0] * range_tolerance
    optimal_threshold = center_threshold + optimal_coded[1] * range_threshold
    
    # Clip to valid ranges
    optimal_tolerance = np.clip(optimal_tolerance, 0.8, 1.5)
    optimal_threshold = np.clip(optimal_threshold, 0.3, 0.6)
    
    # Predicted response
    predicted_auc = -result.fun
    
    print(f"\n{'='*80}")
    print(f"OPTIMAL PARAMETERS (Response Surface)")
    print(f"{'='*80}")
    print(f"Tolerance:  {optimal_tolerance:.4f} Å")
    print(f"Threshold:  {optimal_threshold:.4f}")
    print(f"Linkage:    complete")
    print(f"\nPredicted ROC-AUC: {predicted_auc:.4f}")
    print(f"{'='*80}\n")
    
    return {
        'optimal_tolerance': optimal_tolerance,
        'optimal_threshold': optimal_threshold,
        'optimal_linkage': 'complete',
        'predicted_auc': predicted_auc,
        'coded_values': optimal_coded
    }


def run_confirmation_experiments(optimal_params, n_replicates=5):
    """
    Run confirmation experiments at predicted optimum.
    
    Args:
        optimal_params: dict from optimize_response_surface()
        n_replicates: Number of replicate runs
    
    Returns:
        DataFrame with confirmation results
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = ExperimentLogger(
        experiment_name=f"ccd_confirmation_{timestamp}",
        output_dir=Path(__file__).parent.parent / "docs/research/experiments/results"
    )
    
    tolerance = optimal_params['optimal_tolerance']
    threshold = optimal_params['optimal_threshold']
    linkage = optimal_params['optimal_linkage']
    
    print(f"\n{'='*80}")
    print(f"CONFIRMATION EXPERIMENTS")
    print(f"{'='*80}")
    print(f"Parameters: tolerance={tolerance:.4f}, threshold={threshold:.4f}, linkage={linkage}")
    print(f"Replicates: {n_replicates}")
    print(f"{'='*80}\n")
    
    results = []
    for rep in range(1, n_replicates + 1):
        print(f"[{rep}/{n_replicates}] Running confirmation replicate...")
        
        metrics = evaluate_model(
            tolerance=tolerance,
            occurrence_threshold=threshold,
            linkage=linkage,
            n_refs=5,
            n_actives=74,
            n_decoys=499,
            verbose=False
        )
        
        metrics['replicate'] = rep
        metrics['tolerance'] = tolerance
        metrics['threshold'] = threshold
        metrics['linkage'] = linkage
        
        logger.log_run(
            run_id=rep,
            parameters={'tolerance': tolerance, 'threshold': threshold, 'linkage': linkage},
            metrics=metrics,
            metadata={'replicate': rep, 'phase': 'Phase2_Confirmation'}
        )
        
        results.append(metrics)
        print(f"  ✓ ROC-AUC: {metrics['roc_auc']:.4f}")
    
    logger.save_csv()
    logger.save_json()
    
    df_confirm = pd.DataFrame(results)
    
    # Statistical summary
    mean_auc = df_confirm['roc_auc'].mean()
    std_auc = df_confirm['roc_auc'].std()
    ci_95 = 1.96 * std_auc / np.sqrt(n_replicates)
    
    print(f"\n{'='*80}")
    print(f"CONFIRMATION RESULTS")
    print(f"{'='*80}")
    print(f"Mean ROC-AUC:     {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"95% CI:           [{mean_auc - ci_95:.4f}, {mean_auc + ci_95:.4f}]")
    print(f"Predicted:        {optimal_params['predicted_auc']:.4f}")
    print(f"Prediction Error: {abs(mean_auc - optimal_params['predicted_auc']):.4f}")
    print(f"{'='*80}\n")
    
    # Add to summary
    summary = f"""# Phase 2: Confirmation Experiments

## Optimal Parameters
- **Tolerance**: {tolerance:.4f} Å
- **Threshold**: {threshold:.4f}
- **Linkage**: {linkage}

## Confirmation Results (n={n_replicates})
- **Mean ROC-AUC**: {mean_auc:.4f} ± {std_auc:.4f}
- **95% CI**: [{mean_auc - ci_95:.4f}, {mean_auc + ci_95:.4f}]
- **Range**: {df_confirm['roc_auc'].min():.4f} - {df_confirm['roc_auc'].max():.4f}

## Model Validation
- **Predicted AUC**: {optimal_params['predicted_auc']:.4f}
- **Actual AUC**: {mean_auc:.4f}
- **Prediction Error**: {abs(mean_auc - optimal_params['predicted_auc']):.4f}
- **Model Status**: {'✓ VALIDATED' if abs(mean_auc - optimal_params['predicted_auc']) < 0.05 else '⚠ NEEDS REVIEW'}
"""
    
    logger.generate_summary_report()
    
    # Save custom confirmation summary
    summary_path = logger.results_dir / f"PHASE2_CONFIRMATION_SUMMARY_{datetime.now().strftime('%Y%m%d')}.md"
    with open(summary_path, 'w') as f:
        f.write(summary)
    
    return df_confirm, {'mean_auc': mean_auc, 'std_auc': std_auc, 'ci_95': ci_95}


def plot_response_surface(df_results, model_results, optimal_params, 
                          output_dir=None):
    """
    Create 3D surface plot and 2D contour plot.
    
    Args:
        df_results: DataFrame with CCD results
        model_results: dict from fit_response_surface()
        optimal_params: dict from optimize_response_surface()
        output_dir: Directory to save plots
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "docs/research/experiments/plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create prediction grid
    tol_range = np.linspace(df_results['tolerance'].min(), 
                            df_results['tolerance'].max(), 50)
    thr_range = np.linspace(df_results['threshold'].min(), 
                            df_results['threshold'].max(), 50)
    
    # Center points for coding
    center_tol = (df_results['tolerance'].max() + df_results['tolerance'].min()) / 2
    center_thr = (df_results['threshold'].max() + df_results['threshold'].min()) / 2
    range_tol = (df_results['tolerance'].max() - df_results['tolerance'].min()) / 2
    range_thr = (df_results['threshold'].max() - df_results['threshold'].min()) / 2
    
    TOL_grid, THR_grid = np.meshgrid(tol_range, thr_range)
    
    # Code grid values
    TOL_coded = (TOL_grid - center_tol) / range_tol
    THR_coded = (THR_grid - center_thr) / range_thr
    
    # Predict on grid
    model = model_results['model']
    AUC_grid = np.zeros_like(TOL_grid)
    for i in range(TOL_grid.shape[0]):
        for j in range(TOL_grid.shape[1]):
            X_quad = np.array([[1, TOL_coded[i,j], THR_coded[i,j], 
                               TOL_coded[i,j]**2, THR_coded[i,j]**2, 
                               TOL_coded[i,j]*THR_coded[i,j]]])
            AUC_grid[i,j] = model.predict(X_quad)[0]
    
    # 3D Surface Plot
    fig1 = plt.figure(figsize=(12, 9))
    ax1 = fig1.add_subplot(111, projection='3d')
    
    surf = ax1.plot_surface(TOL_grid, THR_grid, AUC_grid, 
                            cmap='viridis', alpha=0.7, edgecolor='none')
    
    # Plot experimental points
    ax1.scatter(df_results['tolerance'], df_results['threshold'], 
                df_results['roc_auc'], c='red', s=100, marker='o', 
                label='Experimental', edgecolors='black', linewidths=1.5)
    
    # Plot optimum
    ax1.scatter([optimal_params['optimal_tolerance']], 
                [optimal_params['optimal_threshold']], 
                [optimal_params['predicted_auc']], 
                c='gold', s=300, marker='*', edgecolors='black', 
                linewidths=2, label='Optimum')
    
    ax1.set_xlabel('Tolerance (Å)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Threshold', fontsize=12, fontweight='bold')
    ax1.set_zlabel('ROC-AUC', fontsize=12, fontweight='bold')
    ax1.set_title('Response Surface: ROC-AUC vs Parameters', 
                  fontsize=14, fontweight='bold', pad=20)
    ax1.legend(fontsize=10)
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    
    plt.tight_layout()
    fig1.savefig(output_dir / "ccd_response_surface_3d.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved 3D surface plot: {output_dir / 'ccd_response_surface_3d.png'}")
    
    # 2D Contour Plot
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    
    contour = ax2.contourf(TOL_grid, THR_grid, AUC_grid, levels=20, cmap='viridis')
    contour_lines = ax2.contour(TOL_grid, THR_grid, AUC_grid, levels=10, 
                                 colors='white', alpha=0.4, linewidths=1)
    ax2.clabel(contour_lines, inline=True, fontsize=8)
    
    # Plot experimental points
    ax2.scatter(df_results['tolerance'], df_results['threshold'], 
                c=df_results['roc_auc'], s=200, marker='o', 
                edgecolors='white', linewidths=2, cmap='viridis',
                vmin=AUC_grid.min(), vmax=AUC_grid.max())
    
    # Plot optimum
    ax2.scatter([optimal_params['optimal_tolerance']], 
                [optimal_params['optimal_threshold']], 
                c='gold', s=500, marker='*', edgecolors='black', linewidths=3)
    ax2.annotate(f"Optimum\nAUC={optimal_params['predicted_auc']:.4f}", 
                 xy=(optimal_params['optimal_tolerance'], 
                     optimal_params['optimal_threshold']),
                 xytext=(10, 10), textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.8),
                 fontsize=10, fontweight='bold')
    
    ax2.set_xlabel('Tolerance (Å)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Threshold', fontsize=12, fontweight='bold')
    ax2.set_title('Response Surface Contour Plot', fontsize=14, fontweight='bold')
    
    cbar = fig2.colorbar(contour, ax=ax2)
    cbar.set_label('ROC-AUC', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    fig2.savefig(output_dir / "ccd_response_surface_contour.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved contour plot: {output_dir / 'ccd_response_surface_contour.png'}")
    
    plt.close('all')


def main():
    """Main execution pipeline for Phase 2."""
    
    print("\n" + "="*80)
    print("PHASE 2: RESPONSE SURFACE METHODOLOGY")
    print("Central Composite Design (CCD) for Parameter Optimization")
    print("="*80 + "\n")
    
    # Step 1: Create CCD design
    print("Step 1: Creating CCD design matrix...")
    design_matrix = create_ccd_design(
        center_tolerance=1.15,      # Midpoint of [0.8, 1.5]
        center_threshold=0.45,      # Midpoint of [0.3, 0.6]
        range_tolerance=0.35,       # Half-range
        range_threshold=0.15        # Half-range
    )
    print(f"✓ Created {len(design_matrix)} design points")
    print("\nDesign Matrix Preview:")
    print(design_matrix[['experiment_id', 'tolerance', 'threshold', 'linkage']].head(10))
    
    # Step 2: Run CCD experiments
    print("\n" + "="*80)
    print("Step 2: Running CCD experiments...")
    print("="*80)
    df_results = run_ccd_experiments(design_matrix)
    
    # Step 3: Fit response surface model
    print("\n" + "="*80)
    print("Step 3: Fitting quadratic response surface model...")
    print("="*80)
    model_results = fit_response_surface(df_results)
    
    # Step 4: Optimize response surface
    print("\n" + "="*80)
    print("Step 4: Optimizing response surface...")
    print("="*80)
    optimal_params = optimize_response_surface(model_results)
    
    # Step 5: Create visualizations
    print("\n" + "="*80)
    print("Step 5: Creating response surface visualizations...")
    print("="*80)
    plot_response_surface(df_results, model_results, optimal_params)
    
    # Step 6: Run confirmation experiments
    print("\n" + "="*80)
    print("Step 6: Running confirmation experiments...")
    print("="*80)
    df_confirm, confirm_stats = run_confirmation_experiments(optimal_params, n_replicates=5)
    
    # Final summary
    print("\n" + "="*80)
    print("PHASE 2 COMPLETE: RESPONSE SURFACE OPTIMIZATION")
    print("="*80)
    print(f"\nOptimal Parameters:")
    print(f"  Tolerance:  {optimal_params['optimal_tolerance']:.4f} Å")
    print(f"  Threshold:  {optimal_params['optimal_threshold']:.4f}")
    print(f"  Linkage:    {optimal_params['optimal_linkage']}")
    print(f"\nConfirmed Performance:")
    print(f"  Mean ROC-AUC: {confirm_stats['mean_auc']:.4f} ± {confirm_stats['std_auc']:.4f}")
    print(f"  95% CI: [{confirm_stats['mean_auc'] - confirm_stats['ci_95']:.4f}, "
          f"{confirm_stats['mean_auc'] + confirm_stats['ci_95']:.4f}]")
    print(f"\nModel Quality:")
    print(f"  R² = {model_results['r2']:.4f}")
    print(f"  Adjusted R² = {model_results['adj_r2']:.4f}")
    print(f"  RMSE = {model_results['rmse']:.4f}")
    print("="*80 + "\n")
    
    return df_results, model_results, optimal_params, df_confirm


if __name__ == "__main__":
    df_results, model_results, optimal_params, df_confirm = main()
