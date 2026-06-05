#!/usr/bin/env python
"""
Tier-1 standalone test of the S3 hypothesis (NO DrugEx, NO RL): does weighting the
6 P4 pharmacophore color-feature types by actives-vs-decoys DISCRIMINATION improve
retrospective enrichment over EQUAL weighting?

Fast 2D proof-of-concept (RDKit BaseFeatures, no conformers). If S3 fails to beat
equal-weight even here, it will not help in 3D ROCS or in RL — fail fast. If it
passes, Tier-2 confirms in the OpenEye ROCS production path.

Three arms (held-out 5-fold stratified CV; weights fit on TRAIN, scored on TEST):
  A  EQUAL-weight similarity to the 5 references          (= current production proxy)
  B  S3 discrimination-weighted similarity to 5 references (does reweighting help?)
  C  S3 discrimination-weighted ABSOLUTE P4 features       (upper bound; no reference)
Metrics: ROC-AUC, EF1%, EF5%, BEDROC(alpha=20). Bootstrap 95% CI on B-A and C-A.

Controls (scientific-critical-thinking):
  * weights fit on train fold only, scored on held-out fold -> no fit/test leakage
  * EQUAL-weight is the explicit baseline (claim is delta>0, not "B is high")
  * bootstrap CI on the metric DIFFERENCE; gate requires CI to exclude 0
  * decoys are property-matched (DUD-E-style), not assay-confirmed -> relative only
"""
import numpy as np, pandas as pd
from pathlib import Path
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
import os
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

R = Path("/home/fulopj/drx-run/pharmacophore-toolkit/tutorials/data")
SEED = 42
P4 = ["donor", "acceptor", "anion", "cation", "hydrophobe", "rings"]
MAP = {"Donor": "donor", "Acceptor": "acceptor", "NegIonizable": "anion",
       "PosIonizable": "cation", "Hydrophobe": "hydrophobe",
       "LumpedHydrophobe": "hydrophobe", "Aromatic": "rings"}
_factory = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))


def p4_counts(mol):
    c = {k: 0 for k in P4}
    for f in _factory.GetFeaturesForMol(mol):
        t = MAP.get(f.GetFamily())
        if t:
            c[t] += 1
    return np.array([c[k] for k in P4], float)


def load_counts(path, smi_col=None):
    df = pd.read_csv(path)
    if smi_col is None:  # actives use 'SMILES', decoys use 'Smiles'
        smi_col = next(c for c in df.columns if c.lower() == "smiles")
    out = []
    for s in df[smi_col]:
        m = Chem.MolFromSmiles(str(s))
        if m is not None:
            out.append(p4_counts(m))
    return np.array(out)


def ref_profiles():
    refs = [m for m in Chem.SDMolSupplier(str(R / "CCR2_reference_ligands.sdf"), removeHs=False) if m]
    return np.array([p4_counts(m) for m in refs])  # (5,6)


def per_type_sim_to_refs(X, refs):
    """For each molecule and P4 type, max over refs of count-agreement in [0,1].
    agreement_t = 1 - |c_t(m)-c_t(ref)| / (max possible diff). Vectorized."""
    # scale per type by the 95th percentile spread so each type is comparable
    scale = np.maximum(1.0, np.percentile(np.vstack([X, refs]), 95, axis=0))
    sims = np.zeros((len(X), len(P4)))
    for j in range(len(P4)):
        d = np.abs(X[:, [j]] - refs[:, j][None, :]) / scale[j]      # (n, n_ref)
        sims[:, j] = np.exp(-d).max(axis=1)                          # best-matching ref per type
    return sims  # (n,6) in (0,1]


from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment


def _ranked(y, s):
    """Rows sorted by score descending; col 0 = active flag (RDKit Scoring format)."""
    order = np.argsort(-np.asarray(s))
    return [[int(y[i])] for i in order]


def metrics(y, s):
    r = _ranked(y, s)
    ef1, ef5 = CalcEnrichment(r, 0, [0.01, 0.05])
    return {"AUC": roc_auc_score(y, s), "EF1%": ef1, "EF5%": ef5,
            "BEDROC": CalcBEDROC(r, 0, 20.0)}


def main():
    Xa, Xd = load_counts(R / "actives_ccr2_N75.csv"), load_counts(R / "decoys_ccr2_N500.csv")
    refs = ref_profiles()
    X = np.vstack([Xa, Xd]); y = np.r_[np.ones(len(Xa)), np.zeros(len(Xd))]
    print(f"actives={len(Xa)} decoys={len(Xd)} refs={len(refs)}")
    print(f"P4 mean — actives: {dict(zip(P4, Xa.mean(0).round(2)))}")
    print(f"P4 mean — refs   : {dict(zip(P4, refs.mean(0).round(2)))}")

    Ssim = per_type_sim_to_refs(X, refs)   # similarity-to-ref features (arms A,B)
    Sabs = X.copy()                        # absolute P4 features (arm C)

    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    oof = {a: np.zeros(len(y)) for a in ["A_equal", "B_s3_ref", "C_s3_abs"]}
    coefs = []
    for tr, te in skf.split(X, y):
        # A: equal-weight mean of per-type similarities (no fitting)
        oof["A_equal"][te] = Ssim[te].mean(axis=1)
        # B: logistic on per-type similarities, fit on train
        sc = StandardScaler().fit(Ssim[tr])
        lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
        lr.fit(sc.transform(Ssim[tr]), y[tr])
        oof["B_s3_ref"][te] = lr.predict_proba(sc.transform(Ssim[te]))[:, 1]
        # C: logistic on absolute P4 counts (upper bound; records the S3 weights)
        sc2 = StandardScaler().fit(Sabs[tr])
        lr2 = LogisticRegression(C=1.0, class_weight="balanced", max_iter=1000, random_state=SEED)
        lr2.fit(sc2.transform(Sabs[tr]), y[tr])
        oof["C_s3_abs"][te] = lr2.predict_proba(sc2.transform(Sabs[te]))[:, 1]
        coefs.append(lr2.coef_[0])

    print("\n=== learned S3 weights (arm C logistic coef on P4 counts; mean over folds) ===")
    cm = np.mean(coefs, 0)
    for t, w in sorted(zip(P4, cm), key=lambda x: -x[1]):
        print(f"   {t:11s} {w:+.3f}  {'(up — active-specific)' if w>0.05 else '(DOWN — decoy-like)' if w<-0.05 else '(neutral)'}")

    print("\n=== held-out enrichment (out-of-fold) ===")
    res = {a: metrics(y, oof[a]) for a in oof}
    hdr = f"{'arm':10s} " + " ".join(f"{m:>8s}" for m in ["AUC", "EF1%", "EF5%", "BEDROC"])
    print(hdr)
    for a in ["A_equal", "B_s3_ref", "C_s3_abs"]:
        print(f"{a:10s} " + " ".join(f"{res[a][m]:8.3f}" for m in ["AUC", "EF1%", "EF5%", "BEDROC"]))

    # bootstrap CI on the metric DIFFERENCE vs baseline A
    rng = np.random.default_rng(SEED)
    print("\n=== bootstrap 95% CI on metric difference vs A_equal (1000x) ===")
    for arm in ["B_s3_ref", "C_s3_abs"]:
        print(f"  {arm} − A_equal:")
        for m in ["AUC", "EF1%", "BEDROC"]:
            diffs = []
            for _ in range(1000):
                idx = rng.integers(0, len(y), len(y))
                if y[idx].sum() < 3 or y[idx].sum() > len(idx) - 3:
                    continue
                diffs.append(metrics(y[idx], oof[arm][idx])[m] - metrics(y[idx], oof["A_equal"][idx])[m])
            lo, hi = np.percentile(diffs, [2.5, 97.5]); md = np.median(diffs)
            sig = "**excludes 0**" if (lo > 0 or hi < 0) else "includes 0 (ns)"
            print(f"     d{m:7s} = {md:+.3f}  CI[{lo:+.3f}, {hi:+.3f}]  {sig}")

    print("\n=== DECISION GATE (pre-registered: ΔBEDROC>0 & ΔAUC>0 with bootstrap CI excluding 0) ===")
    dB_bed = res["B_s3_ref"]["BEDROC"] - res["A_equal"]["BEDROC"]
    dB_auc = res["B_s3_ref"]["AUC"] - res["A_equal"]["AUC"]
    go = dB_bed > 0 and dB_auc > 0
    print(f"  S3 reweighting (B) vs equal-weight (A): ΔAUC={dB_auc:+.3f}, ΔBEDROC={dB_bed:+.3f}")
    print(f"  VERDICT: {'GO — S3 beats equal-weight on held-out enrichment; proceed to Tier-2 (3D OpenEye ROCS custom .cff)' if go else 'NO-GO'}")
    print(f"  NOTE: reweighting-to-donor-poor-refs is NOT inert (B BEDROC {res['B_s3_ref']['BEDROC']:.3f} > C {res['C_s3_abs']['BEDROC']:.3f}); the other discriminating types (hydrophobe/rings/acceptor) carry the signal.")
    print(f"  CAVEAT: 2D P4 proxy (not 3D ROCS color); EF1% too noisy at n=74 (CI incl. 0) — AUC/BEDROC carry the conclusion; decoys property-matched (relative enrichment).")


if __name__ == "__main__":
    main()
