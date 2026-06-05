#!/usr/bin/env bash
# =============================================================================
# Script: 02_download_fastq.sh
# Description: Downloads SRA runs and extracts them to FASTQ format.
# =============================================================================

# Strict mode for robust error handling
set -euo pipefail

# Default variables
SRA_FILE=""
OUTDIR=""
THREADS=4

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <SRA_RUN_TABLE> -o <OUTPUT_DIR> [OPTIONS]

This script processes an SRA Run Table CSV, prefetches the raw runs, extracts 
them into FASTQ formats dynamically based on layout, and compresses them.

Required Arguments:
  -i, --input FILE         Path to the SraRunTable.csv file
  -o, --outdir DIR         Directory where output FASTQ files should be saved

Optional Arguments:
  -t, --threads INT        Number of processing threads to use (Default: 4)
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -i oxyr_sra.csv -o data/raw_fastq -t 8
EOF
}

check_dependencies() {
    for cmd in prefetch fasterq-dump; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not installed (sra-tools)." >&2
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
        -i|--input) SRA_FILE="$2"; shift ;;
        -o|--outdir) OUTDIR="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        *) echo "Error: Unknown parameter passed: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$SRA_FILE" || -z "$OUTDIR" ]]; then
    echo "Error: Missing required arguments (-i and -o)." >&2
    usage
    exit 1
fi

if [[ ! -f "$SRA_FILE" ]]; then
    echo "Error: Input file '$SRA_FILE' does not exist." >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies
mkdir -p "$OUTDIR"

# Determine the best compression utility and syntax available
if command -v pigz &>/dev/null; then
    GZIP_CMD="pigz -f -p $THREADS"
    echo "Using multi-threaded compressor: pigz ($THREADS threads)"
else
    GZIP_CMD="gzip -f"
    echo "Using fallback single-threaded compressor: gzip"
fi

# Locate the exact column indices dynamically from the CSV header
HEADER=$(head -n 1 "$SRA_FILE")

get_col_idx() {
    local col_name=$1
    echo "$HEADER" | tr ',' '\n' | grep -n "^${col_name}$" | cut -d':' -f1
}

RUN_IDX=$(get_col_idx "Run")
LAYOUT_IDX=$(get_col_idx "LibraryLayout")

if [[ -z "$RUN_IDX" || -z "$LAYOUT_IDX" ]]; then
    echo "Error: Could not find required columns ('Run' or 'LibraryLayout') in CSV header." >&2
    exit 1
fi

echo "Starting download queue..."

# Process data entries, skipping the header row safely
tail -n +2 "$SRA_FILE" | while IFS=, read -r -a row || [[ -n "${row[0]}" ]]; do
    # Extract data using the discovered column indices (adjusting for 0-indexed bash arrays)
    SRR="${row[$((RUN_IDX-1))]}"
    LAYOUT="${row[$((LAYOUT_IDX-1))]}"
    
    # Strip carriage returns if file has Windows line endings
    SRR=$(echo "$SRR" | tr -d '\r')
    LAYOUT=$(echo "$LAYOUT" | tr -d '\r' | tr '[:lower:]' '[:upper:]')

    echo "────────────────────────────────────────────────────────────"
    echo "Processing: $SRR ($LAYOUT-end)"

    # Establish conditional target checks depending on structural layouts
    if [[ "$LAYOUT" == "PAIRED" ]]; then
        TARGETS=("${OUTDIR}/${SRR}_1.fastq.gz" "${OUTDIR}/${SRR}_2.fastq.gz")
    else
        TARGETS=("${OUTDIR}/${SRR}.fastq.gz")
    fi

    # Check if files already exist
    SKIP=true
    for target in "${TARGETS[@]}"; do
        if [[ ! -f "$target" ]]; then
            SKIP=false
        fi
    done

    if [ "$SKIP" = true ]; then
        echo "  [skip] Outputs already exist for $SRR."
        continue
    fi

    # Step 1: Prefetch raw SRA file
    echo "  -> Prefetching $SRR..."
    prefetch "$SRR" --output-directory "$OUTDIR"

    # Step 2: Convert SRA to FASTQ formats via fasterq-dump
    # --split-3 cleanly handles both single and paired data seamlessly
    echo "  -> Extracting FASTQ reads..."
    if ! fasterq-dump "$OUTDIR/$SRR/$SRR.sra" \
        --outdir "$OUTDIR" \
        --split-3 \
        --threads "$THREADS" \
        --progress; then
        echo "Error: fasterq-dump failed on $SRR." >&2
        exit 1
    fi

    # Step 3: Compress extracted output structures
    echo "  -> Compressing reads..."
    if [[ "$LAYOUT" == "PAIRED" ]]; then
        # Check if single file fell out due to mismatched pairs
        if [[ -f "${OUTDIR}/${SRR}.fastq" ]]; then
            $GZIP_CMD "${OUTDIR}/${SRR}.fastq"
        fi
        if [[ -f "${OUTDIR}/${SRR}_1.fastq" ]]; then
            $GZIP_CMD "${OUTDIR}/${SRR}_1.fastq"
        fi
        if [[ -f "${OUTDIR}/${SRR}_2.fastq" ]]; then
            $GZIP_CMD "${OUTDIR}/${SRR}_2.fastq"
        fi
    else
        if [[ -f "${OUTDIR}/${SRR}.fastq" ]]; then
            $GZIP_CMD "${OUTDIR}/${SRR}.fastq"
        fi
    fi

    # Step 4: Clean up intermediate cache directories
    rm -rf "${OUTDIR:?}/${SRR}"
    echo "  [OK] Finished processing $SRR"

done

echo "────────────────────────────────────────────────────────────"
echo "All downloads completed. Output location: $OUTDIR"
ls -lh "$OUTDIR"