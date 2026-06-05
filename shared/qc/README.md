# Quality Control Utilities

Generic QC wrappers compatible with any organism or modality.

- `04_fastqc.sh` — Runs FastQC on all `.fastq.gz` files in a directory, then
  aggregates results with MultiQC into a single `multiqc_report.html`.

Works for bulk, single-cell, and spatial FASTQ files equally.
