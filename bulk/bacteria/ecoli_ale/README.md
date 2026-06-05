# E. coli ndh Adaptive Laboratory Evolution — Bulk RNA-seq

## Study Link
https://doi.org/10.1128/spectrum.02225-23

## Study Context

Transcriptomic profiling of *E. coli* K-12 MG1655 Δ*ndh* (NADH dehydrogenase II knockout)
and four independently evolved lineages (A–D) from an adaptive laboratory evolution (ALE) experiment.

- **Organism:** *Escherichia coli* K-12 MG1655
- **Reference genome:** GCF_000005845.2 (ASM584v2) — auto-downloaded by `05_align.sh`
- **GEO accession:** GSE242875
- **Design:** 1 parent strain + 4 evolved lineages × 2 biological replicates = 10 samples
- **Library type:** Paired-end, strand-specific (reverse-stranded, `-s 2`)

## Pipeline Scripts

| Step | Script | Tool(s) | Output |
|---|---|---|---|
| 01 | `01_geo_to_sra.sh` | NCBI E-utilities, curl | SraRunTable.csv, SampleInfo.csv |
| 02a | `02_download_fastq.sh` | prefetch, fasterq-dump | \*.fastq.gz |
| 02b | `02_download_geo_processed.sh` | wget/curl | Supplementary files (alternative) |
| 03 | `03_check_fastq.sh` | zgrep | R1/R2 read-count parity report |
| 04 | `04_fastqc.sh` | FastQC, MultiQC | multiqc_report.html |
| 05 | `05_align.sh` | Bowtie2, samtools | \*.sorted.bam, \*.bai |
| 06 | `06_test_strandedness.sh` | featureCounts, samtools | strand_config.txt |
| 07 | `07_count.sh` | featureCounts | counts_raw.txt |
| 11 | `11_analyze_counts.py` | pandas, sklearn, seaborn | PCA, heatmap, correlation PDFs |

## Quick Start

```bash
# 1. Fetch metadata
bash 01_geo_to_sra.sh -g GSE242875 -r SraRunTable.csv -s SampleInfo.csv

# 2. Download reads
bash 02_download_fastq.sh -i SraRunTable.csv -o data/raw_fastq -t 8

# 3. Validate
bash 03_check_fastq.sh -d data/raw_fastq

# 4. QC
bash 04_fastqc.sh -d data/raw_fastq -o qc -t 8

# 5. Align
bash 05_align.sh -i data/raw_fastq -o alignments -r ref -t 8

# 6. Test strandedness (pick one BAM to test)
bash 06_test_strandedness.sh -b alignments/SRR26027886.sorted.bam -g ref/ecoli_k12.gff -o strandness

# 7. Count
bash 07_count.sh -b alignments -g ref/ecoli_k12.gff -o counts -c strandness/strand_config.txt

# 8. Analyse
python 11_analyze_counts.py
```

## Dependencies

```
sra-tools (prefetch, fasterq-dump)
fastqc, multiqc
bowtie2, samtools
subread (featureCounts)
python: pandas, numpy, matplotlib, seaborn, scipy, scikit-learn
```
