#!/usr/bin/env python
"""Aggregate the multi-seed bake-off: read every per-seed bakeoff_master.tsv under
results/seeds/seed_*/, and report BEDROC mean +/- across-seed 95% CI per (dataset,
method), plus the PAIRED comparison deltas (same split per seed -> paired) that the
single-split run cannot give honestly:
  - prism - prism_fixed   (does the LEARNED weighting beat fixed weights?)
  - prism_esp - prism      (does electrostatics add signal, or overfit?)

Across-seed variance captures the split-choice noise that a single-split bootstrap CI
misses. (It does NOT power the cross-dataset Friedman test — seeds of one dataset are
correlated resamples, not independent targets.)

Usage: python aggregate_seeds.py [results/seeds]
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(HERE), "results", "seeds")
OUT_MD = os.path.join(os.path.dirname(SEEDS_DIR), "BAKEOFF_MULTISEED.md")

COMPARISONS = [("prism", "prism_fixed"), ("prism_esp", "prism")]


def _load():
    rows = []
    for tsv in sorted(glob.glob(os.path.join(SEEDS_DIR, "seed_*", "bakeoff_master.tsv"))):
        seed = os.path.basename(os.path.dirname(tsv)).replace("seed_", "")
        df = pd.read_csv(tsv, sep="\t")
        df = df[df["error"].astype(str) != "True"].copy()
        df["seed"] = seed
        df["BEDROC"] = pd.to_numeric(df["BEDROC"], errors="coerce")
        rows.append(df[["seed", "dataset", "method", "BEDROC", "verdict"]])
    if not rows:
        sys.exit(f"no per-seed TSVs under {SEEDS_DIR}")
    return pd.concat(rows, ignore_index=True)


def _ci95(vals):
    """Mean and 95% CI of the mean (t-free normal approx) across seeds."""
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan"), 0
    m, sd = float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else 0.0
    sem = sd / max(1, len(v)) ** 0.5
    return m, m - 1.96 * sem, m + 1.96 * sem, len(v)


def main():
    df = _load()
    datasets = list(dict.fromkeys(df["dataset"]))
    methods = list(dict.fromkeys(df["method"]))
    n_seeds = df["seed"].nunique()

    lines = ["# Multi-seed bake-off — BEDROC across split seeds", "",
             f"{n_seeds} scaffold-split seeds; BEDROC(alpha=20) mean +/- across-seed 95% CI "
             "of the mean. Across-seed variance = split-choice noise (the single-split "
             "bootstrap CI misses it). Does NOT power the Friedman test (correlated resamples).", ""]

    for ds in datasets:
        sub = df[df["dataset"] == ds]
        verdict = sub["verdict"].iloc[0] if len(sub) else "?"
        lines += [f"## {ds}  ·  gate: **{verdict}**  ·  seeds={sub['seed'].nunique()}", "",
                  "| method | BEDROC mean [95% CI] | std | n_seeds |",
                  "|--------|----------------------|-----|---------|"]
        stats = {}
        for m in methods:
            vals = sub[sub["method"] == m]["BEDROC"].values
            mean, lo, hi, n = _ci95(vals)
            stats[m] = (mean, vals)
            if n:
                sd = float(np.std(vals, ddof=1)) if n > 1 else 0.0
                lines.append(f"| {m} | {mean:.3f} [{lo:.3f},{hi:.3f}] | {sd:.3f} | {n} |")
        lines.append("")
        # paired comparison deltas (per seed -> paired across seeds)
        lines.append("**Paired deltas across seeds:**")
        piv = sub.pivot_table(index="seed", columns="method", values="BEDROC")
        for a, b in COMPARISONS:
            if a in piv and b in piv:
                d = (piv[a] - piv[b]).dropna().values
                if len(d):
                    md, lo, hi, n = _ci95(d)
                    frac = float(np.mean(np.asarray(d) > 0))
                    lines.append(f"- `{a} - {b}`: {md:+.3f} [{lo:+.3f},{hi:+.3f}] "
                                 f"(>0 in {frac*100:.0f}% of {n} seeds)")
        lines.append("")

    with open(OUT_MD, "w") as f:
        f.write("\n".join(lines))
    print(f"wrote {OUT_MD}")
    print("\n".join(lines[:40]))


if __name__ == "__main__":
    main()
