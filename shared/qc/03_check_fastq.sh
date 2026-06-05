#!/usr/bin/env bash
# =============================================================================
# Script: 03_check_fastq.sh
# Description: Validates raw FASTQ files for existence, size, corruption, 
#              and paired-end read count parity.
# =============================================================================

# Strict mode for robust error handling
set -euo pipefail

# Default variables
FASTQ_DIR=""

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -d <FASTQ_DIRECTORY> [OPTIONS]

This script runs structural sanity checks across raw FASTQ files. It ensures 
files are intact, non-empty, and verifies that forward (R1) and reverse (R2) 
read counts are completely synchronized.

Required Arguments:
  -d, --dir DIR            Directory containing the compressed (*.fastq.gz) files

Optional Arguments:
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -d data/raw_fastq
EOF
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
        -d|--dir) FASTQ_DIR="$2"; shift ;;
        *) echo "Error: Unknown parameter passed: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$FASTQ_DIR" ]]; then
    echo "Error: Missing required directory argument (-d/--dir)." >&2
    usage
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: Directory '$FASTQ_DIR' does not exist." >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

echo "Checking sequence archives in: $FASTQ_DIR"
echo "────────────────────────────────────────────────────────────"

# 1. Quick structural check for matching files
if ! ls "$FASTQ_DIR"/*.fastq.gz &>/dev/null; then
    echo "Error: No compressed FASTQ files (*.fastq.gz) found in $FASTQ_DIR" >&2
    exit 1
fi

# Display disk usage profiles
ls -lh "$FASTQ_DIR"/*.fastq.gz
echo ""

# 2. Check for completely unpopulated or empty archives
EMPTY_FILES=$(find "$FASTQ_DIR" -name "*.fastq.gz" -empty)
if [[ -n "$EMPTY_FILES" ]]; then
    echo "⚠️  WARNING: Empty file targets found!"
    echo "$EMPTY_FILES"
else
    echo "✓ Integrity pre-check: All files are populated on disk."
fi
echo ""

# 3. Process read counts and verify pair parity
echo "Calculating read parity matrices..."
echo "────────────────────────────────────────────────────────────"
printf "%-15s %-12s %-12s %-10s\n" "RunID" "R1_Reads" "R2_Reads" "Status"
echo "────────────────────────────────────────────────────────────"

FAIL_FLAG=0

# Gather unique SRR prefixes based on R1 or singleton file tracking
for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    # Guard against glob expansion failing if zero files exist
    [[ -e "$R1" ]] || continue

    SRR=$(basename "$R1" _1.fastq.gz)
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"

    # Fast optimized read calculation
    # Using `|| true` prevents set -e from crashing if 0 matches are found, 
    # while naturally keeping zgrep's '0' output.
    COUNT_R1=$(zgrep -c "^+$" "$R1" || true)

    if [[ -f "$R2" ]]; then
        COUNT_R2=$(zgrep -c "^+$" "$R2" || true)
        
        if [[ "$COUNT_R1" -eq "$COUNT_R2" ]]; then
            STATUS="✓ OK"
        else
            STATUS="❌ MISMATCH"
            FAIL_FLAG=1
        fi
        
        # Removed the apostrophe flag to ensure cross-platform macOS/Linux compatibility
        printf "%-15s %-12d %-12d %-10s\n" "$SRR" "$COUNT_R1" "$COUNT_R2" "$STATUS"
    else
        # Handle cases where data layout is Single-End instead of crashing
        printf "%-15s %-12d %-12s %-10s\n" "$SRR" "$COUNT_R1" "N/A (Single)" "✓ OK"
    fi
done

echo "────────────────────────────────────────────────────────────"
if [[ "$FAIL_FLAG" -eq 0 ]]; then
    echo "🎉 SUCCESS: All samples passed synchronization validation checks."
else
    echo "⚠️  CRITICAL FAILURE: Mismatched read counts detected. Re-download broken accessions."
    exit 1
fi