from .pharmacophore import *
from .draw import *
from .consensus import PharmacophoreConsensus
from .mol_converter import PharmacophoreToMol
from .optimization import ConsensusOptimizer, optimize_consensus
from .screening_metrics import (
    youden_index, bedroc, rie, confusion_metrics,
    enrichment_factor, screening_report, calculate_all_metrics
)
from .benchmark import ScreeningBenchmark
from .pharm2d_scoring import Pharm2DScorer, score_molecules_pharm2d
from .pharm3d_scoring import Pharm3DScorer, score_molecules_pharm3d
from .rdshape_optimizer import (
    ReferenceEnsembleScorer,
    ConsensusPharmacophoreScorer,
    score_molecules,
    quick_benchmark,
)

# Caching infrastructure
from .caching import (
    PharmacophoreCache,
    ConsensusCacheKey,
    hash_features,
    hash_reference_mols,
    create_cache_key,
)

# DOE optimization (optional - requires scikit-optimize)
try:
    from .doe_optimization import PharmacophoreOptimizer, optimize_pharmacophore
except ImportError:
    pass  # scikit-optimize not installed

# Auto-optimization (optional - requires scikit-optimize)
try:
    from .auto_optimizer import (
        AutoPharmacophoreOptimizer,
        auto_optimize_pharmacophore,
    )
except ImportError:
    pass  # scikit-optimize not installed

# Combo optimizer (optional - requires scikit-optimize)
try:
    from .combo_optimizer import ComboPharmacophoreOptimizer
except ImportError:
    pass  # scikit-optimize not installed

# Unified evaluation framework
from .evaluation import (
    EvaluationConfig,
    EvaluationResult,
    UnifiedEvaluator,
)

# Multi-objective optimization (optional - requires Optuna)
try:
    from .optuna_optimizer import OptunaPharmacophoreOptimizer
except ImportError:
    pass  # Optuna not installed

# HypoGen 3-phase optimization
from .hypogen_optimizer import HypoGenOptimizer

# New modules (Phase 2-5 enhancements)
from .hungarian_matching import match_features, pharmacophore_distance
from .ot_scoring import wasserstein_pharmacophore_distance, wasserstein_similarity
from .ensemble_consensus import EnsembleConsensus
from .evaluation import compute_sdbw

# Strategy selector (optional - requires scikit-optimize)
try:
    from .auto_strategy import StrategySelector
except ImportError:
    pass

__version__ = "0.0.5"
