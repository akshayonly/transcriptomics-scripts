# Shared Utilities

Scripts and helpers that are modality-agnostic and reusable across any pipeline in this repo.

## Contents

| Folder | Purpose |
|---|---|
| `geo_utils/` | NCBI GEO/SRA metadata fetching (`01_geo_to_sra.sh`, `02_download_geo_processed.sh`) |
| `qc/` | FastQC + MultiQC wrappers (`04_fastqc.sh`) |

## Usage

Scripts here can be called directly from any study folder:

```bash
bash ../../shared/geo_utils/01_geo_to_sra.sh -g GSE242875 -r SraRunTable.csv -s SampleInfo.csv
bash ../../shared/qc/04_fastqc.sh -d data/raw_fastq -o qc -t 8
```

Or copy them into the study folder if you need to modify them for a specific organism.
