# transcriptomics

A modality-first collection of RNA-seq pipelines spanning bacteria, plants, fungi, and mammals.

## Repository Map

```
transcriptomics/
├── bulk/           — Reference-genome aligned, population-level expression
├── single_cell/    — Per-cell expression (10x, Smart-seq2, etc.)
├── spatial/        — Spatially resolved expression (Visium, Xenium, etc.)
└── shared/         — Tool-agnostic utilities used across all modalities
```

## Modality Overview

| Modality | Key Tools | Organisms Covered |
|---|---|---|
| Bulk RNA-seq | Bowtie2/HISAT2/STAR, featureCounts, DESeq2/edgeR | Bacteria, Plants, Fungi, Mammals |
| Single-Cell | STARsolo/Cell Ranger, Seurat/Scanpy | Plants, Fungi, Mammals |
| Spatial | Space Ranger, Squidpy, Giotto | Plants, Fungi, Mammals |

## Navigation

Each modality folder contains organism subfolders; each organism subfolder contains
one directory per study/dataset, named descriptively (e.g. `ecoli_ale`, `arabidopsis_drought`).

Every study folder has its own `README.md` with: study context, GEO/SRA accession,
organism + reference genome, and a numbered description of the scripts.

## Conventions

- Scripts are numbered by execution order: `01_`, `02_`, ..., `11_`, etc.
- Shell scripts use `set -euo pipefail` and accept `-h/--help`.
- Python scripts are self-contained; dependencies listed in the module docstring.
- `.gitkeep` files mark placeholder directories awaiting future studies.
