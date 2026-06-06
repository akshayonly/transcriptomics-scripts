#!/usr/bin/env bash
# =============================================================================
# Script: 06_test_strandedness.sh
# Description: Empirically determines the library strandedness by testing
#              a subsample across multiple BAM files, then exports the consensus
#              winning flag to a configuration file for the downstream pipeline.
#
# Usage:
#   bash scripts/06_test_strandedness.sh -b data/align -g data/ref/ecoli_k12.gff
#   bash scripts/06_test_strandedness.sh -b data/align -g data/ref/ecoli_k12.gff -n 4 -s 750000
#   bash scripts/06_test_strandedness.sh -b data/align/SRR001.sorted.bam -g data/ref/ecoli_k12.gff
# =============================================================================

# Strict mode — pipefail intentionally NOT set: samtools view | head triggers
# SIGPIPE (exit 141) when head closes the pipe, which would kill the script.
# Errors are handled explicitly instead.
set -eu

# =============================================================================
# Default Variables
# =============================================================================
BAM_INPUT=""          # path to a single BAM file OR a directory of BAMs
GFF_FILE=""
OUT_DIR="strand_test"
THREADS=4
CONFIG_OUT=""         # set after OUT_DIR is resolved
MAX_BAMS=4            # -n: how many BAMs to test (picks first N alphabetically)
SUBSAMPLE=500000      # -s: reads to subsample per BAM

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -b <BAM_FILE|BAM_DIR> -g <GFF_FILE> [OPTIONS]

Empirically determines RNA-seq library strandedness by subsampling one or more
BAM files and testing featureCounts assignment rates for -s 0, 1, and 2.
Reports per-BAM results, flags any inconsistencies, and writes the consensus
winning flag to a config file consumed by 07_count.sh.

Required Arguments:
  -b, --bam PATH           Path to a single sorted BAM file OR a directory
                           containing *.sorted.bam files
  -g, --gff FILE           Path to the reference annotation (.gff or .gtf)

Optional Arguments:
  -n, --num-bams INT       Number of BAM files to test when -b is a directory
                           (Default: 4; set to 0 to test all BAMs)
  -s, --subsample INT      Number of reads to subsample per BAM (Default: 500000)
  -o, --outdir DIR         Output directory for test files (Default: strand_test)
  -t, --threads INT        Number of processing threads (Default: 4)
  -c, --config FILE        Config output file path
                           (Default: <outdir>/strand_config.txt)
  -h, --help               Show this help message and exit

Examples:
  # Test 4 BAMs from a directory with default 500k subsample
  $(basename "$0") -b data/align -g data/ref/ecoli_k12.gff

  # Test 6 BAMs with a larger 750k subsample
  $(basename "$0") -b data/align -g data/ref/ecoli_k12.gff -n 6 -s 750000

  # Single BAM mode (original behaviour)
  $(basename "$0") -b data/align/SRR001.sorted.bam -g data/ref/ecoli_k12.gff

  # Test all BAMs in a directory
  $(basename "$0") -b data/align -g data/ref/ecoli_k12.gff -n 0
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

# Subsample a BAM to N reads, write to a temp BAM
# Usage: subsample_bam <input.bam> <output.bam> <n_reads>
subsample_bam() {
    local INPUT_BAM="$1"
    local OUTPUT_BAM="$2"
    local N_READS="$3"

    local HEADER_FILE="${OUTPUT_BAM%.bam}.header.sam"
    local BODY_FILE="${OUTPUT_BAM%.bam}.body.sam"

    samtools view -H "$INPUT_BAM" > "$HEADER_FILE"
    samtools view    "$INPUT_BAM" | head -n "$N_READS" > "$BODY_FILE" || true
    # '|| true' absorbs SIGPIPE exit code when head closes the pipe

    cat "$HEADER_FILE" "$BODY_FILE" | samtools view -b -o "$OUTPUT_BAM"
    rm -f "$HEADER_FILE" "$BODY_FILE"
}

# Run featureCounts for all three strand modes on a given BAM
# Populates global arrays ASSIGNED_READS and PERCENTAGES
# Usage: run_strand_tests <test.bam> <prefix>
run_strand_tests() {
    local TEST_BAM="$1"
    local PREFIX="$2"

    for STRAND in 0 1 2; do
        local FC_PREFIX="${PREFIX}_strand_${STRAND}"
        local LOG_FILE="${FC_PREFIX}.log"

        echo -n "    Testing ${STRAND_NAMES[$STRAND]} (-s $STRAND)... "

        if ! featureCounts \
                -T "$THREADS" \
                -p --countReadPairs \
                -s "$STRAND" \
                -a "$GFF_FILE" \
                -F GFF \
                -t gene \
                -g Name \
                -o "${FC_PREFIX}.txt" \
                "$TEST_BAM" > "$LOG_FILE" 2>&1; then
            echo "FAILED (check ${LOG_FILE})"
            ASSIGNED_READS[$STRAND]=0
            PERCENTAGES[$STRAND]="0.00"
            continue
        fi

        local SUMMARY="${FC_PREFIX}.txt.summary"
        local ASSIGNED TOTAL PCT

        ASSIGNED=$(grep -w "^Assigned" "$SUMMARY" | awk '{print $2}')
        TOTAL=$(awk 'NR>1 {sum+=$2} END {print sum}' "$SUMMARY")

        if [[ "$TOTAL" -eq 0 ]]; then
            ASSIGNED_READS[$STRAND]=0
            PERCENTAGES[$STRAND]="0.00"
            echo "FAILED (zero reads in subsample)"
            continue
        fi

        PCT=$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN {printf "%.2f", (a/t)*100}')
        ASSIGNED_READS[$STRAND]=$ASSIGNED
        PERCENTAGES[$STRAND]=$PCT
        echo "Done (${PCT}%)"
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
        -h|--help)      usage; exit 0 ;;
        -b|--bam)       BAM_INPUT="$2";   shift ;;
        -g|--gff)       GFF_FILE="$2";    shift ;;
        -o|--outdir)    OUT_DIR="$2";     shift ;;
        -t|--threads)   THREADS="$2";     shift ;;
        -c|--config)    CONFIG_OUT="$2";  shift ;;
        -n|--num-bams)  MAX_BAMS="$2";    shift ;;
        -s|--subsample) SUBSAMPLE="$2";   shift ;;
        *) echo "Error: Unknown parameter: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Set config output default now that OUT_DIR is resolved
if [[ -z "$CONFIG_OUT" ]]; then
    CONFIG_OUT="$OUT_DIR/strand_config.txt"
fi

# Validate required arguments
if [[ -z "$BAM_INPUT" || -z "$GFF_FILE" ]]; then
    echo "Error: Missing required arguments (-b and -g are required)." >&2
    exit 1
fi

if [[ ! -e "$BAM_INPUT" ]]; then
    echo "Error: BAM path not found: $BAM_INPUT" >&2
    exit 1
fi

if [[ ! -f "$GFF_FILE" ]]; then
    echo "Error: GFF file not found: $GFF_FILE" >&2
    exit 1
fi

if [[ ! "$SUBSAMPLE" =~ ^[0-9]+$ ]] || [[ "$SUBSAMPLE" -lt 1 ]]; then
    echo "Error: --subsample must be a positive integer (got '$SUBSAMPLE')." >&2
    exit 1
fi

if [[ ! "$MAX_BAMS" =~ ^[0-9]+$ ]]; then
    echo "Error: --num-bams must be a non-negative integer (got '$MAX_BAMS')." >&2
    exit 1
fi

# =============================================================================
# Resolve BAM list
# =============================================================================

declare -a BAM_LIST

if [[ -f "$BAM_INPUT" ]]; then
    # Single BAM file passed directly
    BAM_LIST=("$BAM_INPUT")
elif [[ -d "$BAM_INPUT" ]]; then
    # Directory: collect all *.sorted.bam files
    shopt -s nullglob
    ALL_BAMS=("$BAM_INPUT"/*.sorted.bam)
    shopt -u nullglob

    if [[ ${#ALL_BAMS[@]} -eq 0 ]]; then
        echo "Error: No *.sorted.bam files found in: $BAM_INPUT" >&2
        exit 1
    fi

    # Sort alphabetically for reproducibility, then slice to MAX_BAMS
    # MAX_BAMS=0 means use all
    IFS=$'\n' SORTED_BAMS=($(sort <<<"${ALL_BAMS[*]}")); unset IFS

    if [[ "$MAX_BAMS" -eq 0 || "$MAX_BAMS" -ge "${#SORTED_BAMS[@]}" ]]; then
        BAM_LIST=("${SORTED_BAMS[@]}")
    else
        BAM_LIST=("${SORTED_BAMS[@]:0:$MAX_BAMS}")
    fi
else
    echo "Error: -b must be a BAM file or directory of BAMs." >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies
mkdir -p "$OUT_DIR"

declare -A STRAND_NAMES
STRAND_NAMES[0]="Unstranded"
STRAND_NAMES[1]="Forward (Stranded)"
STRAND_NAMES[2]="Reverse Stranded"

# Accumulators for consensus vote across BAMs
declare -A VOTE_TOTALS
VOTE_TOTALS[0]=0
VOTE_TOTALS[1]=0
VOTE_TOTALS[2]=0

# Per-BAM results table: BAM_RESULTS[i]="SRR  pct0  pct1  pct2  winner"
declare -a BAM_RESULTS

N_TESTED=${#BAM_LIST[@]}
TOTAL_BAMS_IN_DIR=1
[[ -d "$BAM_INPUT" ]] && TOTAL_BAMS_IN_DIR=$(find "$BAM_INPUT" -name "*.sorted.bam" | wc -l)

echo "════════════════════════════════════════════════════════════"
echo " Empirical Strandedness Auto-Detection"
echo "════════════════════════════════════════════════════════════"
echo " BAM source    : $BAM_INPUT"
echo " BAMs available: $TOTAL_BAMS_IN_DIR"
echo " BAMs to test  : $N_TESTED"
echo " Subsample size: $(printf "%'d" $SUBSAMPLE) reads per BAM"
echo " GFF           : $GFF_FILE"
echo " Threads       : $THREADS"
echo " Output dir    : $OUT_DIR"
echo "════════════════════════════════════════════════════════════"
echo ""

# ── Per-BAM testing loop ──────────────────────────────────────────────────────
for i in "${!BAM_LIST[@]}"; do
    BAM="${BAM_LIST[$i]}"
    SAMPLE=$(basename "$BAM" .sorted.bam)
    BAM_NUM=$((i + 1))

    echo "────────────────────────────────────────────────────────────"
    echo " BAM ${BAM_NUM}/${N_TESTED}: ${SAMPLE}"
    echo "────────────────────────────────────────────────────────────"

    # Subsample
    TEST_BAM="$OUT_DIR/${SAMPLE}_subsample.bam"
    echo "  Subsampling ${SUBSAMPLE} reads..."
    subsample_bam "$BAM" "$TEST_BAM" "$SUBSAMPLE"
    echo "  ✓ Subsample ready"
    echo ""

    # Run all three strand tests
    declare -A ASSIGNED_READS
    declare -A PERCENTAGES

    run_strand_tests "$TEST_BAM" "$OUT_DIR/${SAMPLE}"

    # Determine winner for this BAM
    BAM_WINNER=$(awk \
        -v p0="${PERCENTAGES[0]}" \
        -v p1="${PERCENTAGES[1]}" \
        -v p2="${PERCENTAGES[2]}" \
        'BEGIN {
            max = p0; strand = 0;
            if (p1 > max) { max = p1; strand = 1; }
            if (p2 > max) { max = p2; strand = 2; }
            print strand;
        }')

    # Accumulate vote
    VOTE_TOTALS[$BAM_WINNER]=$(( VOTE_TOTALS[$BAM_WINNER] + 1 ))

    # Store result row for summary table
    BAM_RESULTS+=("${SAMPLE}|${PERCENTAGES[0]}|${PERCENTAGES[1]}|${PERCENTAGES[2]}|${BAM_WINNER}")

    echo ""
    echo "  → Winner for ${SAMPLE}: -s ${BAM_WINNER} (${STRAND_NAMES[$BAM_WINNER]})"

    # Clean up subsample BAM
    rm -f "$TEST_BAM"

    unset ASSIGNED_READS
    unset PERCENTAGES
    echo ""
done

# ── Per-BAM Summary Table ─────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════"
echo " Per-BAM Results"
echo "════════════════════════════════════════════════════════════"
printf "%-20s %12s %12s %12s %10s\n" "Sample" "-s 0 %" "-s 1 %" "-s 2 %" "Winner"
echo "────────────────────────────────────────────────────────────"

for ROW in "${BAM_RESULTS[@]}"; do
    IFS='|' read -r SNAME P0 P1 P2 WIN <<< "$ROW"
    printf "%-20s %11s%% %11s%% %11s%% %10s\n" \
        "$SNAME" "$P0" "$P1" "$P2" "-s ${WIN} (${STRAND_NAMES[$WIN]})"
done

# ── Consistency Check ─────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════"
echo " Consensus Vote"
echo "════════════════════════════════════════════════════════════"
printf "  -s 0 (Unstranded)       : %d vote(s)\n" "${VOTE_TOTALS[0]}"
printf "  -s 1 (Forward Stranded) : %d vote(s)\n" "${VOTE_TOTALS[1]}"
printf "  -s 2 (Reverse Stranded) : %d vote(s)\n" "${VOTE_TOTALS[2]}"
echo ""

# Check for inconsistency (more than one strand flag received votes)
NONZERO_VOTES=0
for S in 0 1 2; do
    [[ "${VOTE_TOTALS[$S]}" -gt 0 ]] && NONZERO_VOTES=$(( NONZERO_VOTES + 1 ))
done

if [[ "$NONZERO_VOTES" -gt 1 ]]; then
    echo "  ⚠  WARNING: BAMs do not agree on strandedness."
    echo "     This may indicate:"
    echo "       - Mixed library preps across datasets (check GSE accessions)"
    echo "       - A corrupted or low-depth BAM skewing one result"
    echo "       - A genuine protocol difference between samples"
    echo "     Review the per-BAM table above before proceeding."
    echo "     The majority-vote winner will be written to config, but"
    echo "     verify manually before running 07_count.sh."
    echo ""
fi

# ── Final Consensus Winner ────────────────────────────────────────────────────
CONSENSUS_STRAND=$(awk \
    -v v0="${VOTE_TOTALS[0]}" \
    -v v1="${VOTE_TOTALS[1]}" \
    -v v2="${VOTE_TOTALS[2]}" \
    'BEGIN {
        max = v0; strand = 0;
        if (v1 > max) { max = v1; strand = 1; }
        if (v2 > max) { max = v2; strand = 2; }
        print strand;
    }')

echo "  Consensus: -s ${CONSENSUS_STRAND} (${STRAND_NAMES[$CONSENSUS_STRAND]})"

# ── Write Config ──────────────────────────────────────────────────────────────
echo "$CONSENSUS_STRAND" > "$CONFIG_OUT"

if [[ -f "$CONFIG_OUT" && -s "$CONFIG_OUT" ]]; then
    echo "  Saved to: $CONFIG_OUT"
else
    echo "Error: Failed to write config file: $CONFIG_OUT" >&2
    exit 1
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Usage in your counting script:"
echo "   bash scripts/07_count.sh -b alignments -g $GFF_FILE \\"
echo "     -o counts -c $CONFIG_OUT"
echo "════════════════════════════════════════════════════════════"
