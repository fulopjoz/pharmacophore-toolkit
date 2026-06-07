#!/usr/bin/env python
"""Optimizer bake-off driver — committed CCR2 datasets, scaffold-split held-out.

Protocol (one decision per step, all controlled here so it is reproducible):
  1. Load each committed CCR2 dataset (3 different decoy constructions).
  2. GATE with decoy_bias_audit (report verdict; never silently trust a set).
  3. scaffold_split -> fit each scorer on TRAIN, score the held-out TEST once.
  4. Per dataset: BEDROC(alpha=20) / ROC-AUC / EF1% / EF5% + bootstrap 95% CIs,
     plus delta-BEDROC vs the s3_weighted baseline (bootstrap CI).
  5. Across datasets: average ranks + Friedman chi^2 + Nemenyi critical
     difference (hand-computed; scikit_posthocs is absent).

Why 3 CCR2 sets with DIFFERENT decoys? The bias-robustness diagnostic: an
optimizer that wins on property-matched decoys but collapses on the max-unbiased
MUBD set is exploiting decoy bias, not pharmacophore signal.

Run (front node, fast 2D methods only — code test):
    BAKEOFF_FAST=1 ~/miniconda3/bin/python bakeoff.py --fast --datasets ccr2_project
Full run (PBS, all methods incl. 3D):
    BAKEOFF_NCONF=15 BAKEOFF_NJOBS=44 ~/miniconda3/bin/python bakeoff.py
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import os
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATASETS_DIR = os.path.join(os.path.dirname(HERE), "datasets")
REF_SDF = os.path.join(ROOT, "tutorials", "data", "CCR2_reference_ligands.sdf")
RESULTS_DIR = os.path.join(os.path.dirname(HERE), "results")

for p in (HERE, DATASETS_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

from harness import (BenchData, discover, evaluate_scaffold, bootstrap_metric,  # noqa: E402
                     bootstrap_delta)
from audit import decoy_bias_audit  # noqa: E402

# (display name, load.py relpath, how to call its load())
# (display name, load.py relpath, how to call load(), reference SDF for the ref-based
# scorers). The ref SDF is per-dataset ON PURPOSE: the reference-based scorers
# (shape_combo, rdshape_ensemble, pharm2d, learned_scorer, s3_weighted, equal_weight)
# align to these references, so a CCR2 SDF is INVALID for CCR5/CXCR4. Adding a non-CCR2
# dataset requires its own reference SDF (or None -> load_bench raises) — PRISM is the
# only method that is target-agnostic (it derives its query from the train actives).
DATASET_SPECS = [
    ("ccr2_project", "ccr2_project/load.py", lambda m: m.load(),       REF_SDF),
    ("ccr2_mubd",    "ccr2_mubd/load.py",    lambda m: m.load(),       REF_SDF),
    ("created_CCR2", "created/load.py",      lambda m: m.load("CCR2"), REF_SDF),
]

FAST_METHODS = ["equal_weight", "s3_weighted", "differential_mmfp", "pharm2d"]
SLOW_METHODS = ["shape_combo_rdkit", "rdshape_ensemble", "learned_scorer",
                "prism", "prism_fixed", "prism_esp"]
ALL_METHODS = FAST_METHODS + SLOW_METHODS
BASELINE = "s3_weighted"

# Studentized-range / sqrt(2) critical values q_alpha for Nemenyi at alpha=0.05,
# indexed by number of methods k (Demsar 2006, Table 5).
_Q05 = {2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850, 7: 2.949,
        8: 3.031, 9: 3.102, 10: 3.164}


def _load_loader(relpath):
    path = os.path.join(DATASETS_DIR, relpath)
    name = "ld_" + relpath.replace("/", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_bench(relpath, caller, ref_sdf) -> tuple[BenchData, dict]:
    if not ref_sdf or not os.path.exists(ref_sdf):
        raise NotImplementedError(
            f"dataset '{relpath}' has no valid reference SDF ({ref_sdf!r}). The ref-based "
            "scorers need per-target references; CCR2 refs are NOT valid for other targets. "
            "Supply a target-specific SDF in DATASET_SPECS, or run PRISM-only on this set.")
    mod = _load_loader(relpath)
    actives, decoys, meta = caller(mod)
    data = BenchData.from_lists(actives, decoys, ref_sdf)
    return data, meta


def run_dataset(ds_name, data, methods, test_frac, seed, n_boot):
    """Return {method: {metrics + CIs}} and the gate verdict for one dataset."""
    aud = decoy_bias_audit(
        [data.smiles[i] for i in range(len(data.y)) if data.y[i] == 1],
        [data.smiles[i] for i in range(len(data.y)) if data.y[i] == 0],
    )
    print(f"\n### {ds_name}: {int(data.y.sum())} act / {int((data.y==0).sum())} dec"
          f" | gate={aud['verdict']} (median_max_tc={aud['median_max_tc']:.3f})")

    # First pass: collect per-method test scores (aligned to the same split).
    raw = {}
    y_ref = None
    for m in methods:
        try:
            res = evaluate_scaffold(m, data, test_frac=test_frac, seed=seed)
            if y_ref is None:
                y_ref = res["y_te"]
            elif not np.array_equal(y_ref, res["y_te"]):
                raise RuntimeError("test split differs across methods (seed mismatch?)")
            raw[m] = res
            print(f"    {m:18s} BEDROC={res['BEDROC']:.3f}  AUC={res['AUC']:.3f}  "
                  f"EF1%={res['EF1%']:.2f}  (test {res['n_test_act']}a/{res['n_test_dec']}d)")
        except Exception as e:
            print(f"    {m:18s} FAILED: {type(e).__name__}: {e}")
            traceback.print_exc()
            raw[m] = None

    # Second pass: bootstrap CIs + delta vs baseline.
    base = raw.get(BASELINE)
    out = {}
    for m, res in raw.items():
        if res is None:
            out[m] = {"verdict": aud["verdict"], "error": True}
            continue
        bed_med, bed_lo, bed_hi = bootstrap_metric(res["y_te"], res["scores"], "BEDROC", n=n_boot, seed=seed)
        auc_med, auc_lo, auc_hi = bootstrap_metric(res["y_te"], res["scores"], "AUC", n=n_boot, seed=seed)
        d = {"verdict": aud["verdict"], "error": False,
             "BEDROC": res["BEDROC"], "AUC": res["AUC"], "EF1%": res["EF1%"], "EF5%": res["EF5%"],
             "BEDROC_lo": bed_lo, "BEDROC_hi": bed_hi, "AUC_lo": auc_lo, "AUC_hi": auc_hi,
             "n_test_act": res["n_test_act"], "n_test_dec": res["n_test_dec"]}
        if base is not None and m != BASELINE:
            dm, dlo, dhi = bootstrap_delta(res["scores"], base["scores"], base["y_te"], "BEDROC", n=n_boot, seed=seed)
            d.update(dBEDROC=dm, dBEDROC_lo=dlo, dBEDROC_hi=dhi)
        out[m] = d
    return out, aud


def friedman_nemenyi(bedroc_by_ds, methods):
    """Average ranks (1=best per dataset), Friedman chi^2 + p, Nemenyi CD.

    Only methods present (non-error) on ALL datasets are included.
    """
    complete = [m for m in methods if all(m in bedroc_by_ds[ds] and
                                          np.isfinite(bedroc_by_ds[ds][m]) for ds in bedroc_by_ds)]
    datasets = list(bedroc_by_ds.keys())
    N, k = len(datasets), len(complete)
    if k < 2 or N < 2:
        return {"included": complete, "note": "too few methods/datasets for ranking"}

    # ranks: per dataset, rank methods by BEDROC descending (1 = best); average ties.
    rank_matrix = np.zeros((N, k))
    for di, ds in enumerate(datasets):
        vals = np.array([bedroc_by_ds[ds][m] for m in complete])
        order = (-vals).argsort()
        ranks = np.empty(k)
        ranks[order] = np.arange(1, k + 1)
        # average tied ranks
        for v in np.unique(vals):
            mask = vals == v
            if mask.sum() > 1:
                ranks[mask] = ranks[mask].mean()
        rank_matrix[di] = ranks
    avg_ranks = rank_matrix.mean(axis=0)

    # Friedman chi^2
    from scipy.stats import friedmanchisquare
    try:
        chi2, p = friedmanchisquare(*[rank_matrix[:, j] for j in range(k)])
    except Exception:
        chi2, p = float("nan"), float("nan")

    q = _Q05.get(k, 3.164)
    cd = q * math.sqrt(k * (k + 1) / (6.0 * N))
    return {"included": complete, "avg_ranks": dict(zip(complete, avg_ranks.tolist())),
            "friedman_chi2": float(chi2), "friedman_p": float(p), "nemenyi_cd": float(cd),
            "N_datasets": N, "k_methods": k}


def write_report(results, fried, methods, args):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    tsv = os.path.join(RESULTS_DIR, "bakeoff_master.tsv")
    md = os.path.join(RESULTS_DIR, "BAKEOFF.md")

    with open(tsv, "w") as f:
        f.write("dataset\tverdict\tmethod\tBEDROC\tBEDROC_lo\tBEDROC_hi\tAUC\tAUC_lo\tAUC_hi\t"
                "EF1%\tEF5%\tdBEDROC_vs_s3\tdBEDROC_lo\tdBEDROC_hi\tn_test_act\tn_test_dec\terror\n")
        for ds, per in results.items():
            for m in methods:
                d = per.get(m)
                if d is None:
                    continue
                if d.get("error"):
                    f.write(f"{ds}\t{d.get('verdict','?')}\t{m}\t" + "\t".join(["NaN"] * 13) + "\tTrue\n")
                    continue
                f.write("\t".join(str(x) for x in [
                    ds, d["verdict"], m,
                    f"{d['BEDROC']:.4f}", f"{d['BEDROC_lo']:.4f}", f"{d['BEDROC_hi']:.4f}",
                    f"{d['AUC']:.4f}", f"{d['AUC_lo']:.4f}", f"{d['AUC_hi']:.4f}",
                    f"{d['EF1%']:.3f}", f"{d['EF5%']:.3f}",
                    f"{d.get('dBEDROC', float('nan')):.4f}", f"{d.get('dBEDROC_lo', float('nan')):.4f}",
                    f"{d.get('dBEDROC_hi', float('nan')):.4f}", d["n_test_act"], d["n_test_dec"], "False",
                ]) + "\n")

    lines = ["# Optimizer Bake-off — committed CCR2 datasets (scaffold-split held-out)", ""]
    lines.append(f"Seed {args.seed}, test_frac {args.test_frac}, bootstrap {args.bootstrap}, "
                 f"methods: {', '.join(methods)}.  Primary metric **BEDROC(α=20)**; "
                 f"Δ vs **{BASELINE}** with bootstrap 95% CI.")
    lines.append("")
    for ds, per in results.items():
        verdict = next((d["verdict"] for d in per.values() if d and "verdict" in d), "?")
        lines += [f"## {ds}  ·  decoy-bias gate: **{verdict}**", "",
                  "| method | BEDROC [95% CI] | AUC | EF1% | ΔBEDROC vs s3 [95% CI] |",
                  "|--------|-----------------|-----|------|------------------------|"]
        ranked = sorted([m for m in methods if per.get(m) and not per[m].get("error")],
                        key=lambda m: per[m]["BEDROC"], reverse=True)
        for m in ranked:
            d = per[m]
            dci = (f"{d['dBEDROC']:+.3f} [{d['dBEDROC_lo']:+.3f},{d['dBEDROC_hi']:+.3f}]"
                   if "dBEDROC" in d else "— (baseline)")
            lines.append(f"| {m} | {d['BEDROC']:.3f} [{d['BEDROC_lo']:.3f},{d['BEDROC_hi']:.3f}] | "
                         f"{d['AUC']:.3f} | {d['EF1%']:.2f} | {dci} |")
        failed = [m for m in methods if per.get(m) and per[m].get("error")]
        if failed:
            lines.append(f"\n_failed on this dataset: {', '.join(failed)}_")
        lines.append("")

    if "avg_ranks" in fried:
        lines += ["## Cross-dataset ranking (lower avg rank = better)", "",
                  f"Friedman χ²={fried['friedman_chi2']:.3f}, p={fried['friedman_p']:.3f}; "
                  f"Nemenyi CD={fried['nemenyi_cd']:.3f} over N={fried['N_datasets']} datasets, "
                  f"k={fried['k_methods']} methods.", ""]
        lines.append("| method | avg rank |")
        lines.append("|--------|----------|")
        for m, r in sorted(fried["avg_ranks"].items(), key=lambda kv: kv[1]):
            lines.append(f"| {m} | {r:.2f} |")
        lines += ["",
                  f"> **Caveat:** with only N={fried['N_datasets']} datasets the Nemenyi CD "
                  f"({fried['nemenyi_cd']:.2f} rank units) is very wide — almost no pair will be "
                  "'significantly' separated. Treat the average-rank order as the headline and the "
                  "per-dataset bootstrap CIs as the evidence; add MUV/LIT-PCBA + CCR5/CXCR4 to power the test."]
    else:
        lines += ["## Cross-dataset ranking", "", f"_{fried.get('note','n/a')}_"]
    lines.append("")

    lines += [
        "## How to read this bake-off",
        "",
        "- **`differential_mmfp` (supervised 2D Morgan) is the ligand-based UPPER BAR and a "
        "decoy-bias canary, not a pharmacophore optimizer.** A 2D fingerprint classifier learns "
        "the active-vs-decoy *distribution* boundary; scaffold split removes scaffold memorization "
        "but not distributional separability, so a near-perfect BEDROC there signals what pure 2D "
        "ligand memorization can achieve (cf. Wallach & Heifets 2018, 10.1021/acs.jcim.7b00403; "
        "\"do fingerprints identify diverse actives? (no)\" 10.3390/ph17080992). Read it as context, "
        "not as the method to ship.",
        "- **The pharmacophore methods are the actual subjects**: `equal_weight` / `s3_weighted` "
        "(2D color-feature similarity to refs), `pharm2d` (Gobbi feature-pair fingerprints), "
        "`shape_combo_rdkit` / `rdshape_ensemble` (3D rdShapeAlign shape+color), `learned_scorer` "
        "(supervised logistic on per-ref shape/color).",
        "- **Robustness across decoy schemes beats peak BEDROC on one set.** The correct optimizer "
        "is the one stable across property-matched (`ccr2_project`), max-unbiased MUBD (`ccr2_mubd`), "
        "and created-unbiased (`created_CCR2`). A method that wins on property-matched decoys but "
        "collapses on MUBD is exploiting decoy bias.",
        "- The decoy-bias **gate verdict** is shown per dataset; treat `mild-bias`/`biased` results "
        "as relative-enrichment only.",
        "",
    ]

    with open(md, "w") as f:
        f.write("\n".join(lines))
    print(f"\nWrote {md}\n      {tsv}")
    return md, tsv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="2D methods only (no conformers)")
    ap.add_argument("--methods", default="", help="comma list to override method set")
    ap.add_argument("--datasets", default="", help="comma list of dataset names to subset")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--test-frac", type=float, default=0.25)
    ap.add_argument("--bootstrap", type=int, default=1000)
    args = ap.parse_args()

    discover()
    methods = (args.methods.split(",") if args.methods
               else (FAST_METHODS if (args.fast or os.environ.get("BAKEOFF_FAST")) else ALL_METHODS))
    specs = [s for s in DATASET_SPECS
             if not args.datasets or s[0] in args.datasets.split(",")]

    print(f"Bake-off: methods={methods}\n          datasets={[s[0] for s in specs]}")
    results, bedroc_by_ds = {}, {}
    for ds_name, relpath, caller, ref_sdf in specs:
        try:
            data, _meta = load_bench(relpath, caller, ref_sdf)
        except Exception as e:
            print(f"\n### {ds_name}: LOAD FAILED: {type(e).__name__}: {e}")
            continue
        per, _aud = run_dataset(ds_name, data, methods, args.test_frac, args.seed, args.bootstrap)
        results[ds_name] = per
        bedroc_by_ds[ds_name] = {m: (per[m]["BEDROC"] if per.get(m) and not per[m].get("error")
                                     else float("nan")) for m in methods}

    fried = friedman_nemenyi(bedroc_by_ds, methods) if len(results) >= 2 else {"note": "1 dataset"}
    write_report(results, fried, methods, args)


if __name__ == "__main__":
    main()
