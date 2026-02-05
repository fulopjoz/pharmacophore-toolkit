"""Experiment logging utilities for DOE studies.

Provides automatic logging, result saving, and report generation
for consensus pharmacophore optimization experiments.
"""

import json
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional
import pandas as pd


class ExperimentLogger:
    """Logger for DOE experiments with automatic file management."""
    
    def __init__(
        self,
        experiment_name: str,
        output_dir: Path,
        description: str = ""
    ):
        """Initialize experiment logger.
        
        Args:
            experiment_name: Name of experiment (e.g., 'tolerance_sweep')
            output_dir: Directory to save results
            description: Optional experiment description
        """
        self.experiment_name = experiment_name
        self.output_dir = Path(output_dir)
        self.description = description
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create output directories
        self.results_dir = self.output_dir / "results"
        self.plots_dir = self.output_dir / "plots"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Result storage
        self.results: List[Dict[str, Any]] = []
        self.metadata: Dict[str, Any] = {
            'experiment_name': experiment_name,
            'description': description,
            'start_time': datetime.now().isoformat(),
            'timestamp': self.timestamp
        }
        
        print(f"📊 Experiment Logger initialized: {experiment_name}")
        print(f"   Output: {self.results_dir}")
    
    def log_run(
        self,
        run_id: int,
        parameters: Dict[str, Any],
        metrics: Dict[str, float],
        metadata: Optional[Dict[str, Any]] = None
    ):
        """Log a single experimental run.
        
        Args:
            run_id: Unique run identifier
            parameters: Parameter settings (tolerance, threshold, etc.)
            metrics: Performance metrics (roc_auc, ef_1, etc.)
            metadata: Optional additional metadata
        """
        result = {
            'run_id': run_id,
            'timestamp': datetime.now().isoformat(),
            **parameters,
            **metrics
        }
        
        if metadata:
            result['metadata'] = metadata
        
        self.results.append(result)
        
        # Print summary
        print(f"  Run {run_id}: ", end="")
        for key, val in parameters.items():
            print(f"{key}={val} ", end="")
        print(f"→ ROC-AUC={metrics.get('roc_auc', 0):.4f}")
    
    def get_dataframe(self) -> pd.DataFrame:
        """Get results as pandas DataFrame."""
        if not self.results:
            return pd.DataFrame()
        
        # Flatten results (exclude nested metadata)
        flat_results = []
        for r in self.results:
            flat = {k: v for k, v in r.items() if k != 'metadata'}
            flat_results.append(flat)
        
        return pd.DataFrame(flat_results)
    
    def save_csv(self, filename: Optional[str] = None) -> Path:
        """Save results to CSV file.
        
        Args:
            filename: Optional custom filename
        
        Returns:
            Path to saved CSV file
        """
        if filename is None:
            filename = f"{self.experiment_name}_{self.timestamp}.csv"
        
        csv_path = self.results_dir / filename
        df = self.get_dataframe()
        df.to_csv(csv_path, index=False)
        
        print(f"💾 Results saved: {csv_path}")
        return csv_path
    
    def save_json(self, filename: Optional[str] = None) -> Path:
        """Save results to JSON file.
        
        Args:
            filename: Optional custom filename
        
        Returns:
            Path to saved JSON file
        """
        if filename is None:
            filename = f"{self.experiment_name}_{self.timestamp}.json"
        
        json_path = self.results_dir / filename
        
        data = {
            'metadata': self.metadata,
            'results': self.results
        }
        
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"💾 Results saved: {json_path}")
        return json_path
    
    def generate_summary_report(self, filename: Optional[str] = None) -> Path:
        """Generate markdown summary report.
        
        Args:
            filename: Optional custom filename
        
        Returns:
            Path to saved markdown file
        """
        if filename is None:
            filename = f"{self.experiment_name}_{self.timestamp}_summary.md"
        
        md_path = self.results_dir / filename
        df = self.get_dataframe()
        
        if df.empty:
            print("⚠️  No results to summarize")
            return md_path
        
        with open(md_path, 'w') as f:
            # Header
            f.write(f"# {self.experiment_name.replace('_', ' ').title()}\n\n")
            f.write(f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
            f.write(f"**Description**: {self.description}\n")
            f.write(f"**Total Runs**: {len(df)}\n\n")
            
            # Best result
            if 'roc_auc' in df.columns:
                best_idx = df['roc_auc'].idxmax()
                best = df.loc[best_idx]
                
                f.write("## 🏆 Best Configuration\n\n")
                f.write(f"- **ROC-AUC**: {best['roc_auc']:.4f}\n")
                
                if 'ef_1' in df.columns:
                    f.write(f"- **EF@1%**: {best['ef_1']:.2f}\n")
                if 'bedroc' in df.columns:
                    f.write(f"- **BEDROC**: {best['bedroc']:.4f}\n")
                if 'n_features' in df.columns:
                    f.write(f"- **Features**: {int(best['n_features'])}\n")
                
                # Parameters
                f.write("\n**Parameters:**\n")
                param_cols = [c for c in df.columns if c in ['tolerance', 'threshold', 'linkage', 'clustering_method']]
                for col in param_cols:
                    if col in best:
                        f.write(f"- **{col.title()}**: {best[col]}\n")
                
                f.write("\n")
            
            # Summary statistics
            f.write("## 📊 Summary Statistics\n\n")
            
            metric_cols = ['roc_auc', 'ef_1', 'ef_5', 'bedroc', 'time_per_mol_ms']
            available_metrics = [c for c in metric_cols if c in df.columns]
            
            if available_metrics:
                summary = df[available_metrics].describe()
                f.write("```\n")
                f.write(summary.to_string())
                f.write("\n```\n\n")
            
            # Full results table
            f.write("## 📋 Full Results\n\n")
            
            # Select key columns for display
            display_cols = ['run_id']
            display_cols.extend([c for c in ['tolerance', 'threshold', 'linkage', 'clustering_method'] if c in df.columns])
            display_cols.extend([c for c in ['roc_auc', 'ef_1', 'bedroc', 'time_per_mol_ms'] if c in df.columns])
            
            display_df = df[display_cols]
            f.write(display_df.to_markdown(index=False, floatfmt=".4f"))
            f.write("\n")
        
        print(f"📝 Summary report: {md_path}")
        return md_path
    
    def finalize(self):
        """Finalize experiment: save all results and generate report."""
        self.metadata['end_time'] = datetime.now().isoformat()
        self.metadata['total_runs'] = len(self.results)
        
        # Save in multiple formats
        self.save_csv()
        self.save_json()
        self.generate_summary_report()
        
        print(f"\n✅ Experiment '{self.experiment_name}' complete!")
        print(f"   Total runs: {len(self.results)}")
        print(f"   Output: {self.results_dir}")


def load_experiment_results(json_path: Path) -> Dict[str, Any]:
    """Load experiment results from JSON file.
    
    Args:
        json_path: Path to JSON results file
    
    Returns:
        Dict with metadata and results
    """
    with open(json_path) as f:
        data = json.load(f)
    return data


def merge_experiment_results(
    experiment_paths: List[Path],
    output_path: Optional[Path] = None
) -> pd.DataFrame:
    """Merge results from multiple experiments.
    
    Args:
        experiment_paths: List of paths to CSV or JSON files
        output_path: Optional path to save merged results
    
    Returns:
        DataFrame with merged results
    """
    dfs = []
    
    for path in experiment_paths:
        if path.suffix == '.csv':
            df = pd.read_csv(path)
        elif path.suffix == '.json':
            data = load_experiment_results(path)
            df = pd.DataFrame(data['results'])
        else:
            continue
        
        df['source_file'] = path.name
        dfs.append(df)
    
    merged = pd.concat(dfs, ignore_index=True)
    
    if output_path:
        merged.to_csv(output_path, index=False)
        print(f"💾 Merged results saved: {output_path}")
    
    return merged


if __name__ == "__main__":
    # Example usage
    logger = ExperimentLogger(
        experiment_name="test_experiment",
        output_dir=Path("docs/research/experiments"),
        description="Test of experiment logger"
    )
    
    # Simulate some runs
    for i in range(5):
        logger.log_run(
            run_id=i+1,
            parameters={'tolerance': 1.0 + i*0.2, 'threshold': 0.5},
            metrics={'roc_auc': 0.7 + i*0.02, 'ef_1': 3.0 + i*0.5}
        )
    
    logger.finalize()
