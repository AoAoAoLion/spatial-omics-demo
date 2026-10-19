# Spatial transcriptomics analysis — real 10x Visium data

A complete, reproducible spatial transcriptomics analysis on a **public 10x Visium dataset**
(mouse brain, sagittal anterior: 2,695 spots with expression plus coordinates on the section).
Real data, real tissue structure — not a simulation.

```bash
pip install -r requirements.txt
python scripts/analyze.py --outdir results
open results/report.html
```

About 2 minutes on a laptop CPU; the dataset is downloaded once and cached.

## Why spatial is not just scRNA-seq with a scatter plot

On a Visium section every spot has an (x, y) coordinate. Two things follow that single-cell
analysis cannot do:

1. **Clusters can be checked against anatomy.** A domain that is contiguous on the section but
   scattered on the UMAP (or the reverse) is a real finding, not a plotting detail.
2. **Spatial statistics become meaningful.** Moran's I asks whether a gene's expression is
   spatially structured; neighbourhood enrichment asks which domains actually border each other.
   Neither question exists without coordinates.

## What it does

| Step | Detail |
|---|---|
| QC | UMI counts and genes per spot; filters at 500 UMIs / 200 genes |
| Normalisation | library-size normalisation, log1p |
| Feature selection | 2,000 highly variable genes |
| Embedding + clustering | PCA (30) → kNN → Leiden (resolution 0.5) → UMAP |
| Spatial plots | clusters and counts placed back on the tissue |
| Spatially variable genes | Moran's I on highly variable genes, 100 permutations |
| Domain adjacency | neighbourhood enrichment z-scores on the Visium hex grid |

## Result on this dataset

2,693 spots pass QC, 19,112 genes, **12 spatial domains**.

Top spatially variable genes:

| Gene | Moran's I | What it marks in mouse brain |
|---|---|---|
| S100a5 | 0.76 | cortical/hippocampal neurons |
| Mbp | 0.68 | white matter (myelin) |
| Ttr | 0.66 | choroid plexus |
| Gpr88 | 0.55 | striatum |

The top of the list is not a random gene set: myelin, choroid plexus and striatum are exactly
the structures with strong, spatially coherent expression in this section. That is the check
that the analysis is measuring something real.

## Outputs

| File | Content |
|---|---|
| `results/report.html` | spot QC, spatial domains, SVGs, neighbourhood enrichment, methods |
| `results/tables/moranI.csv` | Moran's I for all highly variable genes |
| `results/tables/neighbourhood_zscore.csv` | domain adjacency z-scores |
| `results/adata_spatial.h5ad` | processed AnnData with spatial coordinates |
| `results/figures/` | spot QC, spatial clusters, UMAP, SVGs, enrichment heatmap |

## Requirements

`scanpy`, `squidpy`, `leidenalg`, `igraph`, `pandas`, `matplotlib`, `jinja2` (pinned in
`requirements.txt`). Data downloads through `scanpy.datasets.visium_sge` and is cached in
`~/data/`.

## What this analysis does not do

- **No deconvolution.** A Visium spot is 1-10 cells; assigning cell-type proportions needs a
  single-cell reference and a dedicated method (cell2location, RCTD). Out of scope here, and
  pretending otherwise would be the most common way spatial analyses go wrong.
- **No alignment to H&E histology.** The pipeline uses the spot grid, not image registration.
- **No inference of cell-cell communication.** Ligand-receptor scores on 55-micron spots are
  compositional artefacts more often than biology.
- **No batch integration.** One section, no batch.

## Tests

```bash
pytest tests -q
```

Checks that clustering produces domains, Moran's I is computed and sorted, the neighbourhood
matrix is square and symmetric in shape, and the report is written with its figures.

## Data source

10x Genomics Visium public dataset `V1_Mouse_Brain_Sagittal_Anterior`, distributed through
scanpy/squidpy. Public and free to use.

## License

MIT.

## Rendered output

`example_output/` contains the report and figures from the run described above (2,693 spots, 12
domains), so the result can be inspected without installing anything.
