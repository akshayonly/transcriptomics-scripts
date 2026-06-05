#!/usr/bin/env bash
# =============================================================================
# Script: 06_count.sh
# Description: Counts mapped reads against genomic features using featureCounts,
#              automatically loading the optimal strand setting if available.
#
# Usage:
#   bash scripts/06_count.sh -b alignments -g ref/ecoli_k12.gff -o counts
#   bash scripts/06_count.sh -b alignments -g ref/ecoli_k12.gff -o counts -c strandness/strand_config.txt
#   bash scripts/06_count.sh -b alignments -g ref/ecoli_k12.gff -o counts -s 0
# =============================================================================

set -euo pipefail

# Default variables
ALN_DIR=""
GFF_FILE=""
OUT_DIR=""
THREADS=4
FEATURE="gene"
ATTR="Name"
STRAND=""
CONFIG_FILE=""      # explicit path to strand_config.txt — set via -c flag

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -b <ALN_DIR> -g <GFF_FILE> -o <OUT_DIR> [OPTIONS]

Counts mapped RNA-seq reads per gene using featureCounts.
Automatically loads strand configuration from a config file if present.

Required Arguments:
  -b, --bamdir DIR         Directory containing *.sorted.bam files
  -g, --gff FILE           Path to the reference annotation file (.gff or .gtf)
  -o, --outdir DIR         Directory to save the count matrix

Optional Arguments:
  -c, --config FILE        Path to strand_config.txt produced by 00_test_strandedness.sh
                           (Default: looks in <bamdir>/strand_config.txt, then working dir)
  -s, --strand INT         Strandedness override: 0=unstranded, 1=stranded, 2=reverse
                           (Overrides -c and auto-detection if explicitly provided)
  -t, --threads INT        Number of processing threads (Default: 4)
  -f, --feature STRING     Feature type to count in GFF (Default: gene)
  -a, --attribute STRING   Attribute key to use as gene ID (Default: Name)
  -h, --help               Show this help message and exit

Examples:
  # Auto-detect strand from strandedness output directory
  $(basename "$0") -b alignments -g ref/ecoli_k12.gff -o counts -c strandness/strand_config.txt

  # Explicit strand override
  $(basename "$0") -b alignments -g ref/ecoli_k12.gff -o counts -s 0

  # Custom feature and attribute (e.g. for a different organism/annotation)
  $(basename "$0") -b alignments -g ref/genome.gff -o counts -c strandness/strand_config.txt -f gene -a locus_tag
EOF
}

check_dependencies() {
    if ! command -v featureCounts &> /dev/null; then
        echo "Error: 'featureCounts' (subread package) is not installed." >&2
        echo "Install: conda install -c bioconda subread" >&2
        exit 1
    fi
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
        -b|--bamdir)    ALN_DIR="$2";     shift ;;
        -g|--gff)       GFF_FILE="$2";    shift ;;
        -o|--outdir)    OUT_DIR="$2";     shift ;;
        -c|--config)    CONFIG_FILE="$2"; shift ;;
        -s|--strand)    STRAND="$2";      shift ;;
        -t|--threads)   THREADS="$2";     shift ;;
        -f|--feature)   FEATURE="$2";     shift ;;
        -a|--attribute) ATTR="$2";        shift ;;
        *) echo "Error: Unknown parameter: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# ── Validate required arguments ───────────────────────────────────────────────
if [[ -z "$ALN_DIR" || -z "$GFF_FILE" || -z "$OUT_DIR" ]]; then
    echo "Error: -b, -g, and -o are all required." >&2
    usage
    exit 1
fi

if [[ ! -d "$ALN_DIR" ]]; then
    echo "Error: BAM directory not found: $ALN_DIR" >&2
    exit 1
fi

if [[ ! -f "$GFF_FILE" ]]; then
    echo "Error: Annotation file not found: $GFF_FILE" >&2
    exit 1
fi

# ── Validate strandedness value if explicitly provided ────────────────────────
if [[ -n "$STRAND" && ! "$STRAND" =~ ^[012]$ ]]; then
    echo "Error: -s must be 0, 1, or 2. Got: '$STRAND'" >&2
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies
mkdir -p "$OUT_DIR"

ALN_DIR="${ALN_DIR%/}"

# ── Strand configuration ──────────────────────────────────────────────────────
# Priority order:
#   1. Explicit -s flag (user override — highest priority)
#   2. Explicit -c config file path
#   3. Auto-search: <bamdir>/strand_config.txt
#   4. Auto-search: ./strand_config.txt (working directory — least preferred)
#   5. Default to 0 (unstranded) with a warning

if [[ -n "$STRAND" ]]; then
    # Priority 1: explicit -s override
    echo "User override: using -s $STRAND"

elif [[ -n "$CONFIG_FILE" ]]; then
    # Priority 2: explicit -c path provided
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "Error: Config file not found: $CONFIG_FILE" >&2
        echo "Run scripts/00_test_strandedness.sh first to generate it." >&2
        exit 1
    fi
    STRAND=$(tr -d '[:space:]' < "$CONFIG_FILE")
    echo "Loaded strand config from: $CONFIG_FILE  →  -s $STRAND"

elif [[ -f "$ALN_DIR/strand_config.txt" ]]; then
    # Priority 3: config lives alongside BAM files
    STRAND=$(tr -d '[:space:]' < "$ALN_DIR/strand_config.txt")
    echo "Auto-found strand config: $ALN_DIR/strand_config.txt  →  -s $STRAND"

elif [[ -f "strand_config.txt" ]]; then
    # Priority 4: config in working directory (least preferred — warn user)
    STRAND=$(tr -d '[:space:]' < "strand_config.txt")
    echo "Warning: Found strand_config.txt in working directory." >&2
    echo "  Consider using -c to specify the path explicitly." >&2
    echo "  Loaded: -s $STRAND"

else
    # Priority 5: nothing found — default with clear warning
    STRAND=0
    echo "Warning: No strand config found and -s not provided." >&2
    echo "  Defaulting to -s 0 (Unstranded)." >&2
    echo "  Run scripts/00_test_strandedness.sh to auto-detect, then pass:" >&2
    echo "  -c <outdir>/strand_config.txt" >&2
fi

# ── Validate the value loaded from file is 0, 1, or 2 ────────────────────────
if [[ ! "$STRAND" =~ ^[012]$ ]]; then
    echo "Error: Invalid strand value '$STRAND' — must be 0, 1, or 2." >&2
    echo "Check the contents of your strand_config.txt." >&2
    exit 1
fi

# ── Collect BAM files ─────────────────────────────────────────────────────────
shopt -s nullglob
BAM_FILES=("$ALN_DIR"/*.sorted.bam)
shopt -u nullglob

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "Error: No *.sorted.bam files found in: $ALN_DIR" >&2
    exit 1
fi

# ── Summary before running ────────────────────────────────────────────────────
echo "────────────────────────────────────────────────────────────"
echo "Read Quantification"
echo "────────────────────────────────────────────────────────────"
echo "BAM directory      : $ALN_DIR"
echo "BAM files found    : ${#BAM_FILES[@]}"
echo "Reference GFF      : $GFF_FILE"
echo "Feature type       : $FEATURE"
echo "Gene attribute     : $ATTR"
echo "Strandedness (-s)  : $STRAND"
echo "Output directory   : $OUT_DIR"
echo "────────────────────────────────────────────────────────────"

# ── Run featureCounts ─────────────────────────────────────────────────────────
if ! featureCounts \
        -T  "$THREADS" \
        -p  \
        --countReadPairs \
        -s  "$STRAND" \
        -a  "$GFF_FILE" \
        -F  GFF \
        -t  "$FEATURE" \
        -g  "$ATTR" \
        -o  "$OUT_DIR/counts_raw.txt" \
        "${BAM_FILES[@]}"; then
    echo "Error: featureCounts failed. Check output above." >&2
    exit 1
fi

# ── Output summary ────────────────────────────────────────────────────────────
echo ""
echo "────────────────────────────────────────────────────────────"
echo "Quantification complete"
echo "Count matrix : $OUT_DIR/counts_raw.txt"
echo "────────────────────────────────────────────────────────────"

SUMMARY_FILE="$OUT_DIR/counts_raw.txt.summary"
if [[ -f "$SUMMARY_FILE" ]]; then
    echo "Assignment summary:"
    echo ""
    column -t "$SUMMARY_FILE"
fi