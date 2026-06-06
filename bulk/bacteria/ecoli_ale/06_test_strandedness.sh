#!/usr/bin/env bash
# =============================================================================
# Script: 00_test_strandedness.sh
# Description: Empirically determines the library strandedness by rapidly testing
#              a 500k-read subsample, then automatically exports the winning
#              flag to a configuration file for the downstream pipeline.
# =============================================================================

# Strict mode — note: pipefail intentionally NOT set here because
# samtools view | head | samtools view causes SIGPIPE (exit 141) on the
# first samtools when head closes the pipe, which would kill the script
# before any tests run. We handle errors explicitly instead.
set -eu

# Default variables
BAM_FILE=""
GFF_FILE=""
OUT_DIR="strand_test"
THREADS=4
CONFIG_OUT="$OUT_DIR/strand_config.txt"
SUBSAMPLE=500000

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -b <BAM_FILE> -g <GFF_FILE> [OPTIONS]

This script tests the three possible RNA-seq strandedness configurations on a
rapidly generated subsample of your BAM file. It outputs the optimal '-s' flag
to a text file ($CONFIG_OUT) so your counting script can read it automatically.

Required Arguments:
  -b, --bam FILE           Path to a single sorted BAM file for testing
  -g, --gff FILE           Path to the reference annotation (.gff or .gtf)

Optional Arguments:
  -o, --outdir DIR         Output directory for test files (Default: strand_test)
  -n, --subsample INT      Number of reads to subsample (Default: 500000)
  -t, --threads INT        Number of processing threads (Default: 4)
  -c, --config FILE        Config output file path (Default: strand_config.txt)
  -h, --help               Show this help message and exit
EOF
}

check_dependencies() {
    for cmd in featureCounts samtools awk; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not installed." >&2
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
        -h|--help)   usage; exit 0 ;;
        -b|--bam)    BAM_FILE="$2"; shift ;;
        -g|--gff)    GFF_FILE="$2"; shift ;;
        -o|--outdir) OUT_DIR="$2"; shift ;;
        -n|--subsample) SUBSAMPLE="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        -c|--config) CONFIG_OUT="$2"; shift ;;
        *) echo "Error: Unknown parameter: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

if [[ -z "$BAM_FILE" || -z "$GFF_FILE" ]]; then
    echo "Error: Missing required arguments (-b and -g are required)." >&2
    exit 1
fi

if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file not found: $BAM_FILE" >&2
    exit 1
fi

if [[ ! -f "$GFF_FILE" ]]; then
    echo "Error: GFF file not found: $GFF_FILE" >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies
mkdir -p "$OUT_DIR"

echo "────────────────────────────────────────────────────────────"
echo "Empirical Strandedness Auto-Detection"
echo "────────────────────────────────────────────────────────────"

# ── Step 1: Subsampling ───────────────────────────────────────────────────────
# FIX: The original used a three-way pipe:
#   samtools view -h | head -n 500000 | samtools view -b
# With set -euo pipefail, when head closes after 500k lines the first
# samtools receives SIGPIPE and exits 141 (non-zero). pipefail treats this
# as a fatal error and kills the script before any tests run.
#
# Fix: use samtools view's built-in -s (subsample) flag instead.
# -s 42.1 means: use seed 42, keep 10% of reads — fast and SIGPIPE-free.
# Alternatively we use a two-step approach: extract header separately,
# then use head on the non-header lines only, avoiding the pipe issue.

TEST_BAM="$OUT_DIR/subsample.bam"
echo "Extracting a lightweight subsample (${SUBSAMPLE} reads) for instant testing..."

# Write header first, then stream body lines through head, recombine
# Each step is a separate command — no chained pipes that trigger SIGPIPE
HEADER_FILE="$OUT_DIR/header.sam"
BODY_FILE="$OUT_DIR/body.sam"

samtools view -H "$BAM_FILE" > "$HEADER_FILE"
# samtools view    "$BAM_FILE" | head -n 500000 > "$BODY_FILE" || true
# Removed hardcodded sample size
samtools view    "$BAM_FILE" | head -n "$SUBSAMPLE" > "$BODY_FILE" || true
# '|| true' absorbs the SIGPIPE exit code from samtools when head closes

cat "$HEADER_FILE" "$BODY_FILE" | samtools view -b -o "$TEST_BAM"
rm -f "$HEADER_FILE" "$BODY_FILE"

echo "✓ Subsample created."
echo ""

# ── Step 2: Run featureCounts for each strandedness ───────────────────────────
declare -A ASSIGNED_READS
declare -A PERCENTAGES
declare -A STRAND_NAMES

STRAND_NAMES[0]="Unstranded"
STRAND_NAMES[1]="Forward (Stranded)"
STRAND_NAMES[2]="Reverse Stranded"

for STRAND in 0 1 2; do
    echo -n "Running Test: ${STRAND_NAMES[$STRAND]} (-s $STRAND)... "

    PREFIX="$OUT_DIR/test_strand_${STRAND}"
    LOG_FILE="${PREFIX}.log"

    if ! featureCounts \
            -T "$THREADS" \
            -p --countReadPairs \
            -s "$STRAND" \
            -a "$GFF_FILE" \
            -F GFF \
            -t gene \
            -g Name \
            -o "${PREFIX}.txt" \
            "$TEST_BAM" > "$LOG_FILE" 2>&1; then
        echo "FAILED (check ${LOG_FILE})"
        ASSIGNED_READS[$STRAND]=0
        PERCENTAGES[$STRAND]="0.00"
        continue
    fi

    SUMMARY="${PREFIX}.txt.summary"

    ASSIGNED=$(grep -w "^Assigned" "$SUMMARY" | awk '{print $2}')
    TOTAL=$(awk 'NR>1 {sum+=$2} END {print sum}' "$SUMMARY")

    # Guard against division by zero
    if [[ "$TOTAL" -eq 0 ]]; then
        ASSIGNED_READS[$STRAND]=0
        PERCENTAGES[$STRAND]="0.00"
        echo "FAILED (zero reads in subsample)"
        continue
    fi

    PCT=$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN {printf "%.2f", (a/t)*100}')
    ASSIGNED_READS[$STRAND]=$ASSIGNED
    PERCENTAGES[$STRAND]=$PCT

    echo "Done"
done

# ── Step 3: Report and pick winner ───────────────────────────────────────────
echo ""
echo "────────────────────────────────────────────────────────────"
echo "Decision Matrix: Assignment Rates"
echo "────────────────────────────────────────────────────────────"
printf "%-5s %-22s %-15s %-10s\n" "Flag" "Protocol Type" "Assigned Reads" "Success %"
echo "────────────────────────────────────────────────────────────"

for STRAND in 0 1 2; do
    printf " -s %-2s %-22s %-15d %-10s\n" \
        "$STRAND" \
        "${STRAND_NAMES[$STRAND]}" \
        "${ASSIGNED_READS[$STRAND]}" \
        "${PERCENTAGES[$STRAND]}%"
done
echo "────────────────────────────────────────────────────────────"

# Use awk to find the highest percentage and return the winning strand flag
WINNING_STRAND=$(awk \
    -v p0="${PERCENTAGES[0]}" \
    -v p1="${PERCENTAGES[1]}" \
    -v p2="${PERCENTAGES[2]}" \
    'BEGIN {
        max = p0; strand = 0;
        if (p1 > max) { max = p1; strand = 1; }
        if (p2 > max) { max = p2; strand = 2; }
        print strand;
    }')

echo ""
echo "Optimal configuration detected: -s $WINNING_STRAND (${STRAND_NAMES[$WINNING_STRAND]})"

# ── Step 4: Write config file ─────────────────────────────────────────────────
echo "$WINNING_STRAND" > "$CONFIG_OUT"

# Verify the file was actually written before declaring success
if [[ -f "$CONFIG_OUT" && -s "$CONFIG_OUT" ]]; then
    echo "Saved to: $CONFIG_OUT"
else
    echo "Error: Failed to write config file: $CONFIG_OUT" >&2
    exit 1
fi

# Clean up temporary BAM
rm -f "$TEST_BAM"

echo ""
echo "Usage in your counting script:"
echo "  STRAND=\$(cat $CONFIG_OUT)"
echo "  featureCounts -s \$STRAND ..."
