# Spatial Transcriptomics Pipelines

Spatially resolved transcriptomics: gene expression is profiled while preserving
the physical location of cells or spots within a tissue section.

## Conceptual Pipeline

```
Tissue section image + FASTQ
  → Spot/cell detection + barcode mapping (Space Ranger / Banksy)
  → Spot × gene count matrix with (x, y) coordinates
  → QC + normalisation
  → Spatially variable gene detection
  → Cell type deconvolution (RCTD / Stereoscope / cell2location)
  → Spatial domain segmentation
  → Ligand–receptor interaction mapping (CellChat / COMMOT)
```

## Platforms

| Platform | Resolution | Notes |
|---|---|---|
| 10x Visium | ~55 µm spots | Most widely adopted; Space Ranger for preprocessing |
| 10x Xenium | Single-cell | In situ, targeted panel |
| Slide-seq v2 | ~10 µm beads | Near single-cell resolution |
| MERFISH / seqFISH | Single-cell | Imaging-based, no sequencing step |

## Subdirectories

`bacteria/`, `plants/`, `fungi/`, `mammals/` — each holds one folder per study.
