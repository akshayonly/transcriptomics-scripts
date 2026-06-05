#!/usr/bin/env bash
# =============================================================================
# Script: 04_fastqc.sh
# Description: Generates sequence quality reports for raw FASTQ files using 
#              FastQC and aggregates them into a single report via MultiQC.
# =============================================================================

# Strict mode for robust error handling
set -euo pipefail

# Default variables
FASTQ_DIR=""
OUT_DIR=""
THREADS=4

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <FASTQ_DIRECTORY> -o <QC_OUTPUT_DIRECTORY> [OPTIONS]

This script runs FastQC on all compressed FASTQ files in a given directory and 
compiles the results into a single interactive MultiQC report.

Required Arguments:
  -i, --input DIR          Directory containing the raw *.fastq.gz files
  -o, --outdir DIR         Directory to save the FastQC and MultiQC reports

Optional Arguments:
  -t, --threads INT        Number of processing threads to use (Default: 4)
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -i data/raw_fastq -o results/fastqc -t 8
EOF
}

check_dependencies() {
    for cmd in fastqc multiqc; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not installed or not in PATH." >&2
            echo "Try: conda install -c bioconda fastqc multiqc" >&2
            exit 1
        fi
    done
}

# =============================================================================
# Argument Parsing
# =============================================================================

if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) usage; exit 0 ;;
        -i|--input) FASTQ_DIR="$2"; shift ;;
        -o|--outdir) OUT_DIR="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        *) echo "Error: Unknown parameter passed: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$FASTQ_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: Missing required arguments (-i and -o)." >&2
    usage
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: Input directory '$FASTQ_DIR' does not exist." >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies

# Safely check for fastq.gz files to prevent globbing errors
shopt -s nullglob
FASTQ_FILES=("$FASTQ_DIR"/*.fastq.gz)
shopt -u nullglob # Turn off immediately after use

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "Error: No '.fastq.gz' files found in '$FASTQ_DIR'." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "────────────────────────────────────────────────────────────"
echo "Step 1/2: Running FastQC on ${#FASTQ_FILES[@]} files..."
echo "Using $THREADS threads. This may take a while depending on file size."
echo "────────────────────────────────────────────────────────────"

# FastQC natively accepts multiple files at once.
# By passing the array, we avoid bash "argument list too long" errors.
fastqc "${FASTQ_FILES[@]}" \
    --outdir "$OUT_DIR" \
    --threads "$THREADS" \
    --quiet

echo "✓ FastQC completed successfully."
echo ""

echo "────────────────────────────────────────────────────────────"
echo "Step 2/2: Aggregating results with MultiQC..."
echo "────────────────────────────────────────────────────────────"

multiqc "$OUT_DIR" \
    --outdir "$OUT_DIR" \
    --filename multiqc_report \
    --quiet

echo "✓ MultiQC completed successfully."
echo "────────────────────────────────────────────────────────────"
echo ""

# Post-run guidance
REPORT_PATH="${OUT_DIR}/multiqc_report.html"

echo "🎉 Done! Open the following file in your web browser to review results:"
echo "   -->  $REPORT_PATH"
echo ""
echo "🔍 Key RNA-Seq Metrics to Check:"
echo "  - Per base sequence quality : Expect PASS (Phred > 28 across most of the read)"
echo "  - Adapter content           : Expect PASS (If FAIL, trimming is required in next step)"
echo "  - GC content                : Expect ~51% (Standard for E. coli K-12)"
echo "  - Sequence Duplication      : High is NORMAL for RNA-Seq (Highly expressed genes get duplicated)"