#!/usr/bin/env python3
"""Spatial transcriptomics analysis of a real 10x Visium dataset.

Dataset: 10x Genomics Visium, mouse brain sagittal anterior (public). ~2,700 spots on a
tissue section, each with a barcoded capture area and an (x, y) position — so every cluster
can be put back on the section instead of only on a UMAP.

Pipeline: QC -> normalisation -> clustering -> spatial plots -> Moran's I (spatially variable
genes) -> neighbourhood enrichment -> HTML report.

Usage:
    python scripts/analyze.py --outdir results [--sample V1_Mouse_Brain_Sagittal_Anterior]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from jinja2 import Template

sc.settings.verbosity = 1
plt.rcParams.update({"figure.dpi": 130, "font.size": 9})


def log(msg: str) -> None:
    print(f"[spatial] {msg}", flush=True)


def load(sample: str):
    log(f"loading Visium sample {sample} (public 10x dataset, cached after first run)")
    adata = sc.datasets.visium_sge(sample_id=sample)
    adata.var_names_make_unique()
    # total_counts / n_genes_by_counts are added by calculate_qc_metrics in qc()
    log(f"loaded {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    return adata


def qc(adata, outdir: Path, min_counts: int, min_genes: int):
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    axes[0].hist(adata.obs["total_counts"], bins=60, color="#4C72B0")
    axes[0].axvline(min_counts, color="r", ls="--", lw=1)
    axes[0].set_xlabel("UMI counts per spot")
    axes[1].hist(adata.obs["n_genes_by_counts"], bins=60, color="#55A868")
    axes[1].axvline(min_genes, color="r", ls="--", lw=1)
    axes[1].set_xlabel("genes per spot")
    fig.suptitle(f"Spot QC before filtering ({adata.n_obs:,} spots)")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "01_spot_qc.png")
    plt.close(fig)

    before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    log(f"QC: {before:,} -> {adata.n_obs:,} spots, {adata.n_vars:,} genes")
    return adata, before


def cluster(adata, resolution: float, n_pcs: int, n_hvg: int):
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="cluster", flavor="igraph", n_iterations=2)
    log(f"clustered into {adata.obs['cluster'].nunique()} domains at resolution {resolution}")
    return adata


def spatial_plots(adata, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    sc.pl.spatial(adata, color="cluster", ax=axes[0], show=False, frameon=False,
                  title="Spatial domains", legend_loc="right margin", size=1.35)
    sc.pl.spatial(adata, color="total_counts", ax=axes[1], show=False, frameon=False,
                  title="UMI counts per spot", cmap="viridis", size=1.35)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "02_spatial_clusters.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4.6))
    sc.pl.umap(adata, color="cluster", ax=ax, show=False, frameon=False, legend_loc="on data",
               legend_fontsize=7, title="Clusters (UMAP, no spatial info)")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "03_umap_clusters.png")
    plt.close(fig)


def spatially_variable_genes(adata, outdir: Path, n_top: int = 4):
    sq.gr.spatial_neighbors(adata)
    sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[adata.var["highly_variable"]],
                           n_perms=100, n_jobs=2)
    moran = adata.uns["moranI"].copy()
    moran.index.name = "gene"
    moran.sort_values("I", ascending=False).to_csv(outdir / "tables" / "moranI.csv")

    top = moran.sort_values("I", ascending=False).head(n_top).index.tolist()
    fig, axes = plt.subplots(1, n_top, figsize=(3.1 * n_top, 3.4))
    for ax, gene in zip(np.atleast_1d(axes), top):
        sc.pl.spatial(adata, color=gene, ax=ax, show=False, frameon=False, size=1.2,
                      title=f"{gene} (I={moran.loc[gene, 'I']:.2f})", cmap="magma")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "04_spatially_variable_genes.png")
    plt.close(fig)
    log(f"Moran's I: top genes {top}")
    return moran, top


def neighbourhood(adata, outdir: Path):
    sq.gr.nhood_enrichment(adata, cluster_key="cluster")
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    sq.pl.nhood_enrichment(adata, cluster_key="cluster", ax=ax, title="Neighbourhood enrichment")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "05_neighbourhood_enrichment.png")
    plt.close(fig)
    z = pd.DataFrame(adata.uns["cluster_nhood_enrichment"]["zscore"],
                     index=sorted(adata.obs["cluster"].unique()),
                     columns=sorted(adata.obs["cluster"].unique()))
    z.to_csv(outdir / "tables" / "neighbourhood_zscore.csv")
    return z


TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Spatial transcriptomics report</title>
<style>
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:980px;margin:40px auto;padding:0 20px;color:#1b2430;line-height:1.6}
h1{border-bottom:3px solid #2ec4b6;padding-bottom:10px}
h2{margin-top:34px;border-bottom:1px solid #e3e8ef;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0}
th,td{border:1px solid #e3e8ef;padding:7px 9px;text-align:left}
th{background:#f6f8fb}
img{max-width:100%;border:1px solid #e3e8ef;border-radius:6px;margin:8px 0}
code{background:#f2f4f7;padding:2px 5px;border-radius:3px;font-size:13px}
.meta{color:#5b6b7c;font-size:14px}
</style></head><body>
<h1>Spatial transcriptomics report — 10x Visium</h1>
<p class="meta">Sample: {{ sample }} (public 10x dataset) &middot; {{ n_spots }} spots after QC &middot;
{{ n_genes }} genes &middot; {{ n_clusters }} spatial domains &middot; scanpy {{ version }} + squidpy {{ sq_version }}
&middot; {{ date }}</p>

<h2>1. Spots and QC</h2>
<table><tr><th>Step</th><th>Spots</th></tr>
<tr><td>Loaded</td><td>{{ n_before }}</td></tr>
<tr><td>After filter (≥{{ min_counts }} UMIs, ≥{{ min_genes }} genes)</td><td>{{ n_spots }}</td></tr>
</table>
<img src="figures/01_spot_qc.png" alt="spot QC">
<p>Spots below the count and gene thresholds are removed before normalisation. On a Visium
section this mostly drops tissue-edge spots and off-tissue captures, not biology.</p>

<h2>2. Spatial domains</h2>
<img src="figures/02_spatial_clusters.png" alt="spatial clusters">
<img src="figures/03_umap_clusters.png" alt="UMAP">
<p>Clustering uses expression only; the left panel puts the same labels back on the section.
Domains that look contiguous on tissue but overlap on the UMAP (or the reverse) are the
interesting cases — that difference is the reason to run both views.</p>

<h2>3. Spatially variable genes (Moran's I)</h2>
<img src="figures/04_spatially_variable_genes.png" alt="SVGs">
<table><tr><th>Gene</th><th>Moran's I</th><th>p-value</th></tr>
{% for r in svgs %}<tr><td>{{ r.gene }}</td><td>{{ '%.3f'|format(r.I) }}</td><td>{{ '%.2e'|format(r.pval_norm) }}</td></tr>{% endfor %}
</table>
<p>Moran's I computed on the {{ n_hvg }} highly variable genes with a spatial neighbour graph
({{ n_perms }} permutations). High I means expression is spatially structured rather than
randomly scattered. Full table in <code>results/tables/moranI.csv</code>.</p>

<h2>4. Neighbourhood enrichment</h2>
<img src="figures/05_neighbourhood_enrichment.png" alt="neighbourhood enrichment">
<p>z-scores of domain-domain adjacency: positive means those domains sit next to each other more
often than expected under label permutation. This is a property of the section, not of the
clustering quality, so it is reported as an observation.</p>

<h2>5. Methods paragraph</h2>
<p>Visium spots with fewer than {{ min_counts }} UMIs or {{ min_genes }} detected genes were
excluded. Counts were library-size normalised and log-transformed, {{ n_hvg }} highly variable
genes were selected, and PCA ({{ n_pcs }} components) was followed by Leiden clustering
(resolution {{ resolution }}) on a k-nearest-neighbour graph. Spatial structure was quantified
with Moran's I on the highly variable genes ({{ n_perms }} permutations), and domain adjacency
was tested with neighbourhood enrichment on the Visium hex grid.</p>

<h2>6. Files</h2>
<ul>
<li><code>results/tables/moranI.csv</code> — Moran's I for all highly variable genes</li>
<li><code>results/tables/neighbourhood_zscore.csv</code> — domain adjacency z-scores</li>
<li><code>results/adata_spatial.h5ad</code> — processed object with spatial coordinates</li>
<li><code>results/figures/</code> — spot QC, spatial domains, SVGs, enrichment</li>
</ul>
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", default="V1_Mouse_Brain_Sagittal_Anterior")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--min-counts", type=int, default=500)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--resolution", type=float, default=0.5)
    ap.add_argument("--n-perms", type=int, default=100)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)

    adata = load(args.sample)
    adata, before = qc(adata, outdir, args.min_counts, args.min_genes)
    adata = cluster(adata, args.resolution, args.n_pcs, args.n_hvg)
    spatial_plots(adata, outdir)
    moran, top = spatially_variable_genes(adata, outdir)
    neighbourhood(adata, outdir)

    svgs = moran.loc[top].reset_index().rename(columns={"index": "gene"}).to_dict("records")
    ctx = dict(
        sample=args.sample, n_spots=int(adata.n_obs), n_genes=int(adata.n_vars),
        n_before=int(before), n_clusters=int(adata.obs["cluster"].nunique()),
        min_counts=args.min_counts, min_genes=args.min_genes, n_hvg=args.n_hvg,
        n_pcs=args.n_pcs, resolution=args.resolution, n_perms=args.n_perms,
        version=sc.__version__, sq_version=sq.__version__,
        date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"), svgs=svgs,
    )
    (outdir / "report.html").write_text(Template(TEMPLATE).render(**ctx))
    adata.write_h5ad(outdir / "adata_spatial.h5ad")
    (outdir / "summary.json").write_text(json.dumps(
        {"sample": args.sample, "n_spots": ctx["n_spots"], "n_genes": ctx["n_genes"],
         "n_clusters": ctx["n_clusters"], "n_before": ctx["n_before"]}, indent=2))
    log(f"done: {ctx['n_clusters']} domains, report -> {outdir/'report.html'}")


if __name__ == "__main__":
    main()
