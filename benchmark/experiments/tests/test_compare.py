"""Test the master comparison driver."""
import os

import compare
import harness


def test_run_comparison_returns_row_per_method():
    rows = compare.run_comparison(methods=["equal_weight", "s3_weighted", "differential_mmfp"])
    names = {r["method"] for r in rows}
    assert {"equal_weight", "s3_weighted", "differential_mmfp"} <= names
    for r in rows:
        assert {"AUC", "BEDROC", "EF1%"} <= set(r)
        # every non-baseline row carries a delta-vs-S3 with a bootstrap CI
        if r["method"] != "s3_weighted":
            assert "dBEDROC_vs_s3" in r and "dBEDROC_ci" in r


def test_writes_master_table_and_report(tmp_path):
    out = tmp_path / "out"
    compare.main(methods=["equal_weight", "s3_weighted"], outdir=str(out))
    assert (out / "comparison_master.tsv").exists()
    assert (out / "comparison_report.md").exists()
