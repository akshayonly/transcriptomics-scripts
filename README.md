# transcriptomics-scripts

A modality-first collection of transcriptomics pipelines spanning bacteria, plants, fungi, and mammals.

---

## Pipeline Status

### Bulk RNA-seq

Reference-genome aligned, population-level expression profiling.
Tools: Bowtie2 / HISAT2 / STAR · featureCounts · DESeq2 / edgeR

| Organism | Study / Dataset | Status | Notes |
|---|---|:---:|---|
| Bacteria | [E. coli ndh ALE](bulk/bacteria/ecoli_ale/) | `complete` | delta-ndh knockout + 4 evolved lineages · GSE242875 |
| Plants | — | `planned` | |
| Fungi | — | `planned` | |
| Mammals | — | `planned` | |

### Single-Cell RNA-seq

Per-cell expression profiling with barcode demultiplexing and clustering.
Tools: STARsolo / Cell Ranger · Seurat / Scanpy · Harmony

| Organism | Study / Dataset | Status | Notes |
|---|---|:---:|---|
| Bacteria | — | `planned` | |
| Plants | — | `planned` | |
| Fungi | — | `planned` | |
| Mammals | — | `planned` | |

### Spatial Transcriptomics

Spatially resolved expression with tissue morphology context.
Tools: Space Ranger · Squidpy / Giotto · cell2location

| Organism | Study / Dataset | Status | Notes |
|---|---|:---:|---|
| Plants | — | `planned` | |
| Fungi | — | `planned` | |
| Mammals | — | `planned` | |

### Shared Utilities

Modality-agnostic helpers reusable across any pipeline.

| Utility | Status | Notes |
|---|:---:|---|
| [GEO / SRA download tools](shared/geo_utils/) | `complete` | Works for any GEO accession |
| [FastQC + MultiQC wrapper](shared/qc/) | `complete` | Organism- and modality-agnostic |

---

## Repository Structure

```
transcriptomics-scripts/
├── bulk/
│   ├── bacteria/ecoli_ale/     complete
│   ├── plants/                 planned
│   ├── fungi/                  planned
│   └── mammals/                planned
├── single_cell/
│   ├── bacteria/               planned
│   ├── plants/                 planned
│   ├── fungi/                  planned
│   └── mammals/                planned
├── spatial/
│   ├── plants/                 planned
│   ├── fungi/                  planned
│   └── mammals/                planned
└── shared/
    ├── geo_utils/              complete
    └── qc/                     complete
```

## Modality Overview

| Modality | Key Tools | Organisms |
|---|---|---|
| Bulk RNA-seq | Bowtie2 / HISAT2 / STAR · featureCounts · DESeq2 / edgeR | Bacteria, Plants, Fungi, Mammals |
| Single-Cell RNA-seq | STARsolo / Cell Ranger · Seurat / Scanpy | Plants, Fungi, Mammals |
| Spatial Transcriptomics | Space Ranger · Squidpy / Giotto | Plants, Fungi, Mammals |

## Navigation

Each modality folder contains organism subfolders. Each organism subfolder holds one directory per study or dataset, named descriptively (e.g. `ecoli_ale`, `arabidopsis_drought`).

Every study folder contains its own `README.md` with study context, GEO/SRA accession, reference genome, and a numbered description of the scripts.

## Conventions

- Scripts are numbered by execution order: `01_`, `02_`, ..., `11_`, etc.
- Shell scripts use `set -euo pipefail` and accept `-h/--help`.
- Python scripts are self-contained; dependencies are listed in the module docstring.
- `.gitkeep` files mark placeholder directories awaiting future studies.
- Update the status tables above when adding a new study folder.
