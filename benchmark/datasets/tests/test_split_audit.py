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

from split import scaffold_split, scaffold_split_indices, random_split
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

    # --- regression tests for the 2026-06-06 code review ------------------- #

    def test_acyclic_molecules_are_separate_groups(self):
        """#2: distinct acyclic molecules must NOT collapse into one '' group.

        The old _generic_scaffold returned '' for every ring-free molecule, so
        all acyclic compounds moved together as one indivisible block.
        """
        from split import _generic_scaffold
        acyclics = ["CCCCN", "CCOCCO", "CC(C)CC(=O)O", "CCCCCCCC"]
        keys = {_generic_scaffold(s) for s in acyclics}
        self.assertEqual(
            len(keys), len(acyclics),
            f"Each distinct acyclic molecule must get its own scaffold key; got {keys}")

    def test_cross_list_scaffold_consistency(self):
        """#7: a generic scaffold shared by actives AND decoys must land on the
        SAME split side for both lists (no scaffold crosses train/test at all).
        """
        actives = ["c1ccccc1", "c1ccc2ccccc2c1"]      # benzene, naphthalene
        decoys  = ["Cc1ccccc1", "Cc1ccc2ccccc2c1"]    # benzene, naphthalene
        r = scaffold_split(actives, decoys, test_frac=0.5, seed=1)
        act_benzene_side = "test" if "c1ccccc1"  in r["test_actives"] else "train"
        dec_benzene_side = "test" if "Cc1ccccc1" in r["test_decoys"]  else "train"
        self.assertEqual(
            act_benzene_side, dec_benzene_side,
            "A scaffold shared by actives and decoys must be on the same split side")

    def test_dominant_scaffold_does_not_empty_partition(self):
        """#1: a dominant scaffold must not be able to empty a partition for
        ANY seed (the old random-order greedy could dump everything into test).
        """
        benzene = ["c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
                   "CCCc1ccccc1", "CCCCc1ccccc1", "CCCCCc1ccccc1"]
        others  = ["c1ccc2ccccc2c1", "c1ccncc1", "c1ccoc1", "c1ccsc1"]
        actives = ["c1ccccc1", "c1ccc2ccccc2c1"]      # shares benzene + naphthalene
        decoys  = benzene + others                    # benzene dominates (6 of 10)
        for seed in range(8):
            r = scaffold_split(actives, decoys, test_frac=0.5, seed=seed)
            self.assertGreater(len(r["train_decoys"]), 0, f"seed {seed}: train_decoys empty")
            self.assertGreater(len(r["test_decoys"]),  0, f"seed {seed}: test_decoys empty")
            self.assertGreater(len(r["train_actives"]), 0, f"seed {seed}: train_actives empty")
            self.assertGreater(len(r["test_actives"]),  0, f"seed {seed}: test_actives empty")

    def test_split_indices_consistent_with_smiles_api(self):
        """scaffold_split_indices must induce exactly the same partition as the
        SMILES-returning scaffold_split (one tested core, two views)."""
        smiles = list(ACTIVES) + list(DECOYS)
        labels = [1] * len(ACTIVES) + [0] * len(DECOYS)
        train_idx, test_idx = scaffold_split_indices(smiles, labels, test_frac=0.5, seed=42)
        # No index appears in both partitions; together they cover everything.
        self.assertEqual(set(train_idx) & set(test_idx), set())
        self.assertEqual(set(train_idx) | set(test_idx), set(range(len(smiles))))
        # Same molecule partition as the SMILES API.
        ref = scaffold_split(ACTIVES, DECOYS, test_frac=0.5, seed=42)
        test_act_by_idx = {smiles[i] for i in test_idx if labels[i] == 1}
        self.assertEqual(test_act_by_idx, set(ref["test_actives"]))

    def test_split_indices_no_scaffold_crosses_partition(self):
        """The index split must keep every scaffold on a single side."""
        smiles = list(ACTIVES) + list(DECOYS)
        labels = [1] * len(ACTIVES) + [0] * len(DECOYS)
        train_idx, test_idx = scaffold_split_indices(smiles, labels, test_frac=0.5, seed=7)
        train_sc = {_get_generic_scaffold(smiles[i]) for i in train_idx}
        test_sc  = {_get_generic_scaffold(smiles[i]) for i in test_idx}
        self.assertEqual(train_sc & test_sc, set())

    def test_stratify_balances_each_class(self):
        """Stratified split must put ~test_frac of BOTH actives and decoys in test
        (the pooled split can starve one class -- e.g. 3/74 test actives)."""
        actives = ["C" * n + "N" for n in range(2, 14)]   # 12 distinct acyclic actives
        decoys  = ["C" * n + "O" for n in range(2, 14)]   # 12 distinct acyclic decoys
        smiles = actives + decoys
        labels = [1] * len(actives) + [0] * len(decoys)
        train_idx, test_idx = scaffold_split_indices(smiles, labels, test_frac=0.25,
                                                     seed=3, stratify=True)
        n_te_act = sum(labels[i] == 1 for i in test_idx)
        n_te_dec = sum(labels[i] == 0 for i in test_idx)
        self.assertEqual(n_te_act, 3, "stratify should give ~25% of actives in test")
        self.assertEqual(n_te_dec, 3, "stratify should give ~25% of decoys in test")

    def test_split_deterministic_across_hashseed(self):
        """The split must depend only on `seed`, not on PYTHONHASHSEED — set
        iteration order over scaffold keys must not leak into the partition."""
        import subprocess
        code = (
            "import sys; sys.path.insert(0, %r)\n"
            "from split import scaffold_split_indices\n"
            "act=['C'*n+'N' for n in range(2,32)]; dec=['C'*n+'O' for n in range(2,32)]\n"
            "sm=act+dec; lab=[1]*30+[0]*30\n"
            "tr,te=scaffold_split_indices(sm,lab,test_frac=0.25,seed=42,stratify=True)\n"
            "print(','.join(map(str,sorted(te))))\n"
        ) % _DATASETS_DIR
        outs = set()
        for hs in ("0", "1", "2", "7"):
            env = dict(os.environ, PYTHONHASHSEED=hs)
            r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
            outs.add(r.stdout.strip())
        self.assertEqual(len(outs), 1, f"split varies with PYTHONHASHSEED: {outs}")

    def test_raises_when_a_class_cannot_be_split(self):
        """A dataset where all actives share ONE scaffold cannot be scaffold-split
        with actives in both train and test -> must fail loudly, not silently
        return a degenerate (empty-partition) split.
        """
        actives = ["c1ccccc1", "Cc1ccccc1"]                  # both benzene
        decoys  = ["c1ccc2ccccc2c1", "Cc1ccc2ccccc2c1"]      # both naphthalene
        with self.assertRaises(ValueError):
            scaffold_split(actives, decoys, test_frac=0.5, seed=0)


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

    def test_empty_actives_not_silently_unbiased(self):
        """#3: with no actives the audit cannot judge bias; it must NOT report
        'unbiased' (which would pass a dataset it never actually evaluated).
        """
        result = decoy_bias_audit([], list(DECOYS))
        self.assertNotEqual(
            result["verdict"], "unbiased",
            "Empty actives must not yield an 'unbiased' verdict")
        self.assertEqual(result["verdict"], "insufficient_data")

    def test_all_unparseable_actives_not_silently_unbiased(self):
        """All-unparseable actives behave like empty actives -> insufficient_data."""
        result = decoy_bias_audit(["not_a_smiles", "@@@@"], list(DECOYS))
        self.assertEqual(result["verdict"], "insufficient_data")


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
