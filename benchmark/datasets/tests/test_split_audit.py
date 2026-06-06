"""
TDD tests for benchmark/datasets/split.py and audit.py.

Run with:
    ~/miniconda3/bin/python -m pytest benchmark/datasets/tests/test_split_audit.py -v
or:
    ~/miniconda3/bin/python benchmark/datasets/tests/test_split_audit.py
"""

import sys
import os
import unittest

# Ensure the benchmark/datasets directory is on the path regardless of cwd.
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATASETS_DIR = os.path.dirname(_HERE)
if _DATASETS_DIR not in sys.path:
    sys.path.insert(0, _DATASETS_DIR)

from split import scaffold_split, random_split
from audit import decoy_bias_audit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_generic_scaffold(smiles: str) -> str:
    """Return the generic Bemis-Murcko scaffold SMILES for a molecule."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(generic)


# ---------------------------------------------------------------------------
# Carefully chosen tiny datasets with KNOWN scaffolds
# ---------------------------------------------------------------------------

# Two distinct Bemis-Murcko generic scaffolds:
#   SCAFFOLD_A: benzene ring (phenyl / monosubstituted)
#   SCAFFOLD_B: naphthalene ring

SCAFFOLD_A_SMILES = [
    "c1ccccc1",           # benzene itself
    "Cc1ccccc1",          # toluene
    "CCc1ccccc1",         # ethylbenzene
    "Oc1ccccc1",          # phenol
]

SCAFFOLD_B_SMILES = [
    "c1ccc2ccccc2c1",     # naphthalene itself
    "Cc1ccc2ccccc2c1",    # 2-methylnaphthalene
    "CCc1ccc2ccccc2c1",   # 2-ethylnaphthalene
    "Oc1ccc2ccccc2c1",    # 2-naphthol
]

# For scaffold split: 4 actives (2 from each scaffold), 4 decoys (2 from each)
ACTIVES = SCAFFOLD_A_SMILES[:2] + SCAFFOLD_B_SMILES[:2]
DECOYS  = SCAFFOLD_A_SMILES[2:] + SCAFFOLD_B_SMILES[2:]

# For bias audit: deliberately biased decoys = same molecules as actives
BIASED_DECOYS = list(ACTIVES)

# Unrelated decoys: aliphatic chains — no ring scaffolds → very different from ring actives
UNBIASED_DECOYS = [
    "CCCCCCCC",   # octane
    "CCCCCCCCO",  # 1-octanol
    "CCCCCCN",    # hexylamine
    "CCCCC(=O)O", # pentanoic acid
]


# ---------------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------------

class TestScaffoldSplit(unittest.TestCase):
    """(a) scaffold_split must produce ZERO shared generic scaffolds between
    the train and test partitions."""

    def test_zero_shared_scaffolds_actives(self):
        result = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        train_act = result["train_actives"]
        test_act  = result["test_actives"]
        self.assertGreater(len(train_act), 0, "train_actives must be non-empty")
        self.assertGreater(len(test_act),  0, "test_actives must be non-empty")
        train_scaffolds = {_get_generic_scaffold(s) for s in train_act}
        test_scaffolds  = {_get_generic_scaffold(s) for s in test_act}
        shared = train_scaffolds & test_scaffolds
        self.assertEqual(shared, set(),
                         f"Scaffold leak: {shared} appear in both train and test actives")

    def test_zero_shared_scaffolds_decoys(self):
        result = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        train_dec = result["train_decoys"]
        test_dec  = result["test_decoys"]
        self.assertGreater(len(train_dec), 0, "train_decoys must be non-empty")
        self.assertGreater(len(test_dec),  0, "test_decoys must be non-empty")
        train_scaffolds = {_get_generic_scaffold(s) for s in train_dec}
        test_scaffolds  = {_get_generic_scaffold(s) for s in test_dec}
        shared = train_scaffolds & test_scaffolds
        self.assertEqual(shared, set(),
                         f"Scaffold leak: {shared} appear in both train and test decoys")

    def test_all_molecules_accounted_for(self):
        result = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        all_act = set(result["train_actives"]) | set(result["test_actives"])
        all_dec = set(result["train_decoys"])  | set(result["test_decoys"])
        self.assertEqual(all_act, set(ACTIVES), "All actives must appear in one split")
        self.assertEqual(all_dec, set(DECOYS),  "All decoys must appear in one split")

    def test_deterministic(self):
        r1 = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        r2 = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        self.assertEqual(r1["train_actives"], r2["train_actives"])
        self.assertEqual(r1["test_actives"],  r2["test_actives"])

    def test_returns_expected_keys(self):
        result = scaffold_split(ACTIVES, DECOYS)
        for key in ("train_actives", "test_actives", "train_decoys", "test_decoys"):
            self.assertIn(key, result, f"Missing key: {key}")


class TestRandomSplit(unittest.TestCase):
    """(b) random_split CAN share scaffolds (no scaffold constraint)."""

    def test_all_molecules_accounted_for(self):
        result = random_split(ACTIVES, DECOYS, test_frac=0.5, seed=0)
        all_act = set(result["train_actives"]) | set(result["test_actives"])
        all_dec = set(result["train_decoys"])  | set(result["test_decoys"])
        self.assertEqual(all_act, set(ACTIVES))
        self.assertEqual(all_dec, set(DECOYS))

    def test_random_split_may_share_scaffolds(self):
        """With enough repetitions at least one seed shares a scaffold — demonstrating
        random_split has no scaffold-disjoint guarantee."""
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        # Use a set where all molecules share a single scaffold type so random
        # split will almost certainly place some in train and some in test.
        all_smi = SCAFFOLD_A_SMILES + SCAFFOLD_B_SMILES  # 8 molecules, 2 scaffolds
        found_shared = False
        for seed in range(50):
            result = random_split(all_smi, [], test_frac=0.5, seed=seed)
            train_scaffolds = {_get_generic_scaffold(s) for s in result["train_actives"] if s}
            test_scaffolds  = {_get_generic_scaffold(s) for s in result["test_actives"]  if s}
            if train_scaffolds & test_scaffolds:
                found_shared = True
                break
        self.assertTrue(found_shared,
                        "random_split should share scaffolds across splits in at least one seed")

    def test_deterministic(self):
        r1 = random_split(ACTIVES, DECOYS, test_frac=0.25, seed=7)
        r2 = random_split(ACTIVES, DECOYS, test_frac=0.25, seed=7)
        self.assertEqual(r1["train_actives"], r2["train_actives"])


class TestDecoyBiasAudit(unittest.TestCase):
    """(c) audit flags biased sets and passes unbiased ones."""

    def test_biased_set_flagged(self):
        """Decoys that are identical to actives must be 'biased (analog leakage)'."""
        result = decoy_bias_audit(ACTIVES, BIASED_DECOYS)
        self.assertIn("verdict", result)
        self.assertIn("biased", result["verdict"].lower(),
                      f"Expected biased verdict, got: {result['verdict']!r}")

    def test_unbiased_set_passes(self):
        """Aliphatic-chain decoys vs ring actives must report 'unbiased'."""
        result = decoy_bias_audit(ACTIVES, UNBIASED_DECOYS)
        self.assertIn("verdict", result)
        self.assertEqual(result["verdict"], "unbiased",
                         f"Expected unbiased verdict, got: {result['verdict']!r}")

    def test_audit_returns_required_keys(self):
        result = decoy_bias_audit(ACTIVES, DECOYS)
        for key in ("mean_max_tc", "median_max_tc", "p95_max_tc",
                    "fraction_above_0_35", "verdict",
                    "actives_mean_mw", "decoys_mean_mw",
                    "actives_mean_logp", "decoys_mean_logp",
                    "actives_mean_hbd", "decoys_mean_hbd",
                    "actives_mean_hba", "decoys_mean_hba"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_audit_statistics_are_numeric(self):
        result = decoy_bias_audit(ACTIVES, DECOYS)
        for key in ("mean_max_tc", "median_max_tc", "p95_max_tc", "fraction_above_0_35"):
            self.assertIsInstance(result[key], float,
                                  f"Key {key} should be float, got {type(result[key])}")

    def test_biased_high_fraction_above_threshold(self):
        """When decoys == actives, fraction_above_0_35 must be 1.0."""
        result = decoy_bias_audit(ACTIVES, BIASED_DECOYS)
        self.assertAlmostEqual(result["fraction_above_0_35"], 1.0, places=5)


# ---------------------------------------------------------------------------
# Entry point for direct execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite  = unittest.TestSuite()
    for cls in (TestScaffoldSplit, TestRandomSplit, TestDecoyBiasAudit):
        suite.addTests(loader.loadTestsFromTestCase(cls))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
