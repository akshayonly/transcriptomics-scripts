#!/usr/bin/env bash
# =============================================================================
# Script: 05_align.sh
# Description: Downloads reference genome/annotation, builds a Bowtie2 index, 
#              and aligns raw FASTQ files, outputting sorted, indexed BAMs.
# =============================================================================

# Strict mode for robust error handling
# pipefail is CRUCIAL here to catch bowtie2 failures before they hit samtools
set -euo pipefail

# Default variables
REF_DIR="ref"
FASTQ_DIR="data/raw_fastq"
ALN_DIR="alignments"
THREADS=4

# Hardcoded reference URLs (can be parameterized if needed, but standard for E.coli)
FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz"

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <FASTQ_DIR> -o <ALN_DIR> -r <REF_DIR> [OPTIONS]

This script handles the core genome alignment workflow. It automatically fetches 
the E. coli K-12 MG1655 reference genome, builds a Bowtie2 index, aligns raw reads,
and streams the output directly into compressed, sorted BAM files.

Required Arguments:
  -i, --input DIR          Directory containing raw *.fastq.gz files
  -o, --outdir DIR         Directory to output alignments (.bam, .bai, .log)
  -r, --refdir DIR         Directory to store reference genome and indexes

Optional Arguments:
  -t, --threads INT        Number of processing threads to use (Default: 4)
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -i data/raw_fastq -o alignments -r ref -t 8
EOF
}

check_dependencies() {
    for cmd in wget gunzip bowtie2 bowtie2-build samtools; do
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
        -h|--help) usage; exit 0 ;;
        -i|--input) FASTQ_DIR="$2"; shift ;;
        -o|--outdir) ALN_DIR="$2"; shift ;;
        -r|--refdir) REF_DIR="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        *) echo "Error: Unknown parameter passed: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$FASTQ_DIR" || -z "$ALN_DIR" || -z "$REF_DIR" ]]; then
    echo "Error: Missing required arguments (-i, -o, or -r)." >&2
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
mkdir -p "$ALN_DIR" "$REF_DIR"

INDEX_BASE="$REF_DIR/ecoli_k12"
FASTA_FILE="${INDEX_BASE}.fasta"
GFF_FILE="${INDEX_BASE}.gff"

# ── Step 1: Download Reference Genome ────────────────────────────────────────
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Downloading E. coli K-12 MG1655 reference genome..."
    wget -q --show-progress -c "$FASTA_URL" -O "${FASTA_FILE}.gz"
    gunzip -f "${FASTA_FILE}.gz"
    echo "Genome Check: $(grep -c '^>' "$FASTA_FILE") sequence(s) found."
else
    echo "✓ Reference genome already exists."
fi

# ── Step 2: Download Annotation (GFF) ────────────────────────────────────────
if [[ ! -f "$GFF_FILE" ]]; then
    echo "Downloading annotation (GFF)..."
    wget -q --show-progress -c "$GFF_URL" -O "${GFF_FILE}.gz"
    gunzip -f "${GFF_FILE}.gz"
    echo "Annotation Check: $(grep -v '^#' "$GFF_FILE" | awk '$3=="gene"' | wc -l) gene records found."
else
    echo "✓ Annotation already exists."
fi

# ── Step 3: Build Bowtie2 Index ──────────────────────────────────────────────
# Check for the primary index file to determine if build is needed
if [[ ! -f "${INDEX_BASE}.1.bt2" ]]; then
    echo "Building Bowtie2 index (this may take a few minutes)..."
    bowtie2-build --threads "$THREADS" --quiet "$FASTA_FILE" "$INDEX_BASE"
    echo "✓ Index built successfully."
else
    echo "✓ Bowtie2 index already exists."
fi

echo "────────────────────────────────────────────────────────────"
echo "Starting Alignment Queue..."
echo "────────────────────────────────────────────────────────────"

# ── Step 4: Align Samples ────────────────────────────────────────────────────
shopt -s nullglob
FASTQ_PAIRS=("$FASTQ_DIR"/*_1.fastq.gz)
shopt -u nullglob

if [[ ${#FASTQ_PAIRS[@]} -eq 0 ]]; then
    echo "Error: No paired-end FASTQ files (*_1.fastq.gz) found in $FASTQ_DIR" >&2
    exit 1
fi

for R1 in "${FASTQ_PAIRS[@]}"; do
    SRR=$(basename "$R1" _1.fastq.gz)
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    OUT_BAM="${ALN_DIR}/${SRR}.sorted.bam"
    LOG_FILE="${ALN_DIR}/${SRR}.log"

    # Ensure the R2 file actually exists before starting
    if [[ ! -f "$R2" ]]; then
         echo "Error: Missing R2 partner for $SRR. Skipping." >&2
         continue
    fi

    # Skip if final BAM and its index already exist
    if [[ -f "$OUT_BAM" && -f "${OUT_BAM}.bai" ]]; then
        echo "  [skip] $SRR already aligned and indexed."
        continue
    fi

    echo "▶ Aligning: $SRR"

    # The Core Pipeline: Bowtie2 -> Samtools Sort
    # By using a pipe, we avoid writing massive intermediate SAM files to disk
    if ! bowtie2 -x "$INDEX_BASE" \
                 -1 "$R1" \
                 -2 "$R2" \
                 --threads "$THREADS" \
                 --no-unal \
                 2> "$LOG_FILE" | \
         samtools sort -@ "$THREADS" -m 2G -o "$OUT_BAM" -; then
         
         echo "❌ Error: Alignment pipeline failed for $SRR. Cleaning up partial files." >&2
         rm -f "$OUT_BAM" "${OUT_BAM}.bai"
         exit 1
    fi

    # Index the newly created BAM file
    samtools index -@ "$THREADS" "$OUT_BAM"

    # Extract and display the alignment rate
    RATE=$(grep 'overall alignment rate' "$LOG_FILE" || echo "Rate not found")
    echo "  ✓ $RATE"
done

# ── Step 5: Summary ──────────────────────────────────────────────────────────
echo "────────────────────────────────────────────────────────────"
echo "Alignment Summary (Expect >90% for pure E. coli isolates):"
printf "%-15s %-20s\n" "Sample (SRR)" "Alignment Rate"
echo "────────────────────────────────────────────────────────────"

for LOG in "$ALN_DIR"/*.log; do
    [[ -e "$LOG" ]] || break
    SRR=$(basename "$LOG" .log)
    RATE=$(grep 'overall alignment rate' "$LOG" | awk '{print $1}')
    printf "%-15s %-20s\n" "$SRR" "${RATE:-N/A}"
done