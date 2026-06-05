# Bulk RNA-seq Pipelines

Population-level transcriptomics: RNA is extracted from a tissue or culture,
sequenced, aligned to a reference genome, and quantified per gene.

## Conceptual Pipeline

```
FASTQ → QC (FastQC/MultiQC) → Alignment (Bowtie2 / HISAT2 / STAR)
     → Quantification (featureCounts / HTSeq / Salmon / kallisto)
     → Normalisation & DE (DESeq2 / edgeR / limma-voom)
```

## Aligner Notes by Organism

| Organism | Recommended Aligner | Reason |
|---|---|---|
| Bacteria | Bowtie2 | Prokaryote — no splicing, compact genome |
| Plants / Fungi / Mammals | HISAT2 or STAR | Eukaryote — splice-junction-aware required |

## Subdirectories

| Folder | Contents |
|---|---|
| `bacteria/` | Prokaryotic bulk-seq (E. coli, etc.) |
| `plants/` | Plant bulk-seq (Arabidopsis, maize, etc.) |
| `fungi/` | Fungal bulk-seq (S. cerevisiae, etc.) |
| `mammals/` | Mammalian bulk-seq (human, mouse, etc.) |
