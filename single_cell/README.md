# Single-Cell RNA-seq Pipelines

Per-cell transcriptomics: individual cells are isolated, barcoded, and sequenced
separately so each cell receives its own expression profile.

## Conceptual Pipeline

```
FASTQ → Cell barcode + UMI demultiplexing (Cell Ranger / STARsolo / alevin)
     → Cell × gene count matrix
     → QC (filter dead cells, doublets, low-complexity)
     → Normalisation (scran / SCnorm)
     → Dimensionality reduction (PCA → UMAP / t-SNE)
     → Clustering (Louvain / Leiden)
     → Cell type annotation (marker genes / SingleR / CellTypist)
     → Differential expression (MAST / DESeq2 pseudo-bulk)
```

## Key Tools

| Language | Package | Purpose |
|---|---|---|
| R | Seurat | End-to-end scRNA-seq analysis |
| Python | Scanpy / scverse | End-to-end scRNA-seq analysis |
| Both | SingleR / CellTypist | Automated cell type annotation |
| Both | Harmony / scVI | Batch correction |

## Subdirectories

`bacteria/`, `plants/`, `fungi/`, `mammals/` — each holds one folder per study.
