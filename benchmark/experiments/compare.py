"""
Master comparison driver: score every registered method on the SAME held-out CCR2
split and rank them against our S3 solution.

Output: results/comparison_master.tsv + comparison_report.md, with each method's
held-out AUC/BEDROC/EF and ΔBEDROC-vs-S3 + bootstrap CI. A method "beats" S3 only
if ΔBEDROC > 0 with a CI excluding 0 (the Tier-1 decision gate).

Run: python compare.py
"""
from __future__ import annotations

import os

import harness

BASELINE = "s3_weighted"
DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def run_comparison(methods=None, data=None):
    """Evaluate each method (held-out CV) and attach ΔBEDROC/ΔAUC vs S3 with bootstrap CIs."""
    harness.discover()
    methods = methods or sorted(harness.REGISTRY)
    data = data or harness.BenchData.load_default()

    oof = {m: harness.evaluate_oof(m, data) for m in methods}
    base = oof.get(BASELINE)
    rows = []
    for m in methods:
        met = harness.metrics(data.y, oof[m])
        row = {"method": m, **{k: round(v, 4) for k, v in met.items()}}
        if base is not None and m != BASELINE:
            bedroc_lo = 0.0
            for metric in ("BEDROC", "AUC"):
                md, lo, hi = harness.bootstrap_delta(oof[m], base, data.y, metric, n=1000)
                row[f"d{metric}_vs_s3"] = round(md, 4)
                row[f"d{metric}_ci"] = f"[{lo:+.3f},{hi:+.3f}]"
                if metric == "BEDROC":
                    bedroc_lo = lo
            # beats S3 iff the bootstrap CI of the BEDROC gain strictly excludes 0
            row["beats_s3"] = bool(bedroc_lo > 0)
        rows.append(row)
    rows.sort(key=lambda r: -r["BEDROC"])
    return rows


def main(methods=None, outdir=None):
    outdir = outdir or DEFAULT_OUT
    os.makedirs(outdir, exist_ok=True)
    rows = run_comparison(methods=methods)

    cols = ["method", "AUC", "BEDROC", "EF1%", "EF5%",
            "dBEDROC_vs_s3", "dBEDROC_ci", "dAUC_vs_s3", "dAUC_ci", "beats_s3"]
    with open(os.path.join(outdir, "comparison_master.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    with open(os.path.join(outdir, "comparison_report.md"), "w") as f:
        f.write("# CCR2 Enrichment — Methods Comparison (held-out 5-fold CV)\n\n")
        f.write("Metric: BEDROC(α=20) primary; AUC/EF secondary. Δ vs our **s3_weighted** "
                "baseline with bootstrap 95% CI. A method beats S3 only if ΔBEDROC>0 & CI excludes 0.\n\n")
        f.write("| method | AUC | BEDROC | EF1% | ΔBEDROC vs S3 | CI | beats S3? |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['method']} | {r['AUC']} | {r['BEDROC']} | {r['EF1%']} | "
                    f"{r.get('dBEDROC_vs_s3','—')} | {r.get('dBEDROC_ci','—')} | "
                    f"{'✅' if r.get('beats_s3') else ('baseline' if r['method']==BASELINE else 'no')} |\n")
    return rows


if __name__ == "__main__":
    for r in main():
        print(r)
