"""End-to-end test for the spatial transcriptomics demo.

Runs on the real public Visium section with reduced permutations, so the spatial neighbour
graph, Moran's I and neighbourhood enrichment are actually exercised.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("sp") / "results"
    proc = subprocess.run(
        [PY, str(REPO / "scripts" / "analyze.py"), "--outdir", str(out),
         "--n-hvg", "800", "--n-pcs", "20", "--n-perms", "10"],
        capture_output=True, text=True, check=True,
    )
    return out, proc.stdout


def test_domains_found(run):
    out, _ = run
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_spots"] > 1000
    assert summary["n_clusters"] >= 4


def test_moran_sorted(run):
    out, _ = run
    moran = pd.read_csv(out / "tables" / "moranI.csv", index_col=0)
    assert "I" in moran.columns
    assert len(moran) > 100
    assert moran["I"].is_monotonic_decreasing, "Moran's I table must be sorted for the report"


def test_neighbourhood_matrix(run):
    out, _ = run
    z = pd.read_csv(out / "tables" / "neighbourhood_zscore.csv", index_col=0)
    summary = json.loads((out / "summary.json").read_text())
    assert z.shape == (summary["n_clusters"], summary["n_clusters"])


def test_report_and_figures(run):
    out, _ = run
    html = (out / "report.html").read_text()
    assert "Spatial transcriptomics report" in html
    assert "Moran" in html
    for fig in ["01_spot_qc.png", "02_spatial_clusters.png",
                "04_spatially_variable_genes.png", "05_neighbourhood_enrichment.png"]:
        assert (out / "figures" / fig).exists(), fig
