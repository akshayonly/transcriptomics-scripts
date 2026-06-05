# GEO / SRA Utilities

Scripts for fetching metadata and raw/processed data from NCBI GEO and SRA.

- `01_geo_to_sra.sh` — Convert a GEO accession to SRA Run Table + Sample Info CSV
- `02_download_fastq.sh` — Download SRA runs and extract to FASTQ
- `02_download_geo_processed.sh` — Mirror the GEO supplementary (processed) file directory

These scripts are organism- and modality-agnostic: any GEO-deposited dataset works.
