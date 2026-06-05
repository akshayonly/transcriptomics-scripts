#!/usr/bin/env bash
# =============================================================================
# Script: 02_download_geo_processed.sh
# Description: Downloads all processed/supplementary files for a GEO Series
#              accession from NCBI's canonical FTP tree.
#
#              NCBI builds the FTP "range" directory by replacing the last
#              three digits of the accession number with "nnn", e.g.
#                GSE242875 -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242875/suppl/
#
#              This is far more reliable than the /geo/download/?acc=...&format=file
#              endpoint, which 404s for series that have no bundled tar archive.
# =============================================================================

set -euo pipefail

# ----- Defaults --------------------------------------------------------------
GEO=""
OUT_DIR="."
BUNDLE=false        # --bundle: grab the single supplementary .tar instead of the dir
RETRIES=3

# ----- Logging helpers -------------------------------------------------------
err()  { printf 'Error: %s\n' "$*" >&2; }
info() { printf '%s\n' "$*"; }

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -g <GEO_ACCESSION> [-o <OUTPUT_DIR>] [--bundle]

Downloads the supplementary (processed) files for a GEO Series from NCBI's
FTP tree. By default it mirrors every file in the series' suppl/ directory.

Required Arguments:
  -g, --geo ACCESSION      GEO Series accession (e.g., GSE242875)

Optional Arguments:
  -o, --out-dir DIR        Directory to save downloads (default: current dir)
      --bundle             Fetch the single <ACC>_RAW.tar bundle instead of
                           mirroring the suppl/ directory
  -r, --retries N          Download retry attempts (default: ${RETRIES})
  -h, --help               Show this help message and exit

Examples:
  $(basename "$0") -g GSE242875 -o ./processed_data
  $(basename "$0") -g GSE242875 -o ./processed_data --bundle
EOF
}

check_dependencies() {
    if ! command -v wget >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1; then
        err "Neither 'wget' nor 'curl' is installed or in PATH."
        exit 1
    fi
}

# Build the FTP "range" subdirectory name from an accession.
#   GSE242875 -> GSE242nnn   (last 3 digits -> nnn)
#   GSE5290    -> GSE5nnn
#   GSE939     -> GSEnnn     (<= 3 digits)
geo_range_dir() {
    local acc="$1"
    local prefix="${acc%%[0-9]*}"   # GSE
    local num="${acc#"$prefix"}"    # 242875
    if [[ -z "$num" || ! "$num" =~ ^[0-9]+$ ]]; then
        err "Could not parse a numeric ID from accession '$acc'."
        exit 1
    fi
    if (( ${#num} > 3 )); then
        printf '%s%snnn' "$prefix" "${num:0:${#num}-3}"
    else
        printf '%snnn' "$prefix"
    fi
}

# Mirror an FTP directory using wget (preferred).
download_dir_wget() {
    local url="$1" dest="$2"
    # -e robots=off: NCBI's robots.txt is "Disallow: /", which would otherwise
    # make wget skip every file and download only robots.txt.
    wget -e robots=off --no-verbose --tries="$RETRIES" -c \
         -r -np -nd -R "index.html*,robots.txt" \
         -P "$dest" "$url"
}

# Fallback: scrape the directory listing with curl and fetch each file.
download_dir_curl() {
    local url="$1" dest="$2" index f
    local -a files

    index="$(curl -fsSL --retry "$RETRIES" "$url")" || {
        err "Could not list directory: $url"
        return 1
    }

    # Pull relative filenames out of the HTML index, dropping the parent-dir
    # link, sort-order query links (?C=...), and any absolute paths.
    mapfile -t files < <(
        printf '%s' "$index" \
        | grep -oE 'href="[^"]+"' \
        | sed -E 's/href="([^"]+)"/\1/' \
        | grep -vE '^[?]|/$|^/|^\.\.'
    )

    if (( ${#files[@]} == 0 )); then
        err "No supplementary files found at: $url"
        return 1
    fi

    for f in "${files[@]}"; do
        info "  -> ${f}"
        curl -fSL --retry "$RETRIES" -C - -o "${dest}/${f}" "${url}${f}"
    done
}

# Fetch the single bundled supplementary tar (the old behaviour).
download_bundle() {
    local geo="$1" dest="$2"
    local url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${geo}&format=file"
    local out="${dest}/${geo}_RAW.tar"
    info " URL: $url"
    if command -v curl >/dev/null 2>&1; then
        curl -fSL --retry "$RETRIES" -C - -o "$out" "$url"
    else
        wget --no-verbose --tries="$RETRIES" -c -O "$out" "$url"
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
    case "$1" in
        -h|--help)    usage; exit 0 ;;
        -g|--geo)     GEO="${2:?-g/--geo requires a value}"; shift ;;
        -o|--out-dir) OUT_DIR="${2:?-o/--out-dir requires a value}"; shift ;;
        -r|--retries) RETRIES="${2:?-r/--retries requires a value}"; shift ;;
        --bundle)     BUNDLE=true ;;
        --)           shift; break ;;
        *) err "Unknown parameter: $1"; usage; exit 1 ;;
    esac
    shift
done

# ----- Validate --------------------------------------------------------------
if [[ -z "$GEO" ]]; then
    err "Missing required argument -g/--geo."
    usage
    exit 1
fi

GEO="${GEO^^}"   # normalise to uppercase
if [[ ! "$GEO" =~ ^GSE[0-9]+$ ]]; then
    err "Accession '$GEO' doesn't look like a GSE Series id (expected e.g. GSE242875)."
    exit 1
fi

if [[ ! "$RETRIES" =~ ^[0-9]+$ ]]; then
    err "--retries must be a non-negative integer (got '$RETRIES')."
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies
mkdir -p "$OUT_DIR"

# Stage downloads in a temp dir so a failure never leaves partial files behind.
STAGING="$(mktemp -d "${TMPDIR:-/tmp}/geo_${GEO}.XXXXXX")"
cleanup() { rm -rf "$STAGING"; }
trap cleanup EXIT

echo "============================================================"
echo " Fetching processed data for: $GEO"
echo " Destination: $OUT_DIR"

if [[ "$BUNDLE" == true ]]; then
    echo " Mode: single supplementary bundle (.tar)"
    echo "============================================================"
    download_bundle "$GEO" "$STAGING"
else
    RANGE="$(geo_range_dir "$GEO")"
    SUPPL_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${RANGE}/${GEO}/suppl/"
    echo " Mode: mirror suppl/ directory"
    echo " URL:  $SUPPL_URL"
    echo "============================================================"

    if command -v wget >/dev/null 2>&1; then
        download_dir_wget "$SUPPL_URL" "$STAGING"
    else
        download_dir_curl "$SUPPL_URL" "$STAGING"
    fi
fi

# Anything to move?
if [[ -z "$(ls -A "$STAGING")" ]]; then
    err "No files were downloaded for $GEO. Check that the accession is correct"
    err "and that it has supplementary files at the URL above."
    exit 1
fi

# Atomically move staged files into the output directory.
mv "$STAGING"/* "$OUT_DIR"/

echo ""
echo "Success! Files saved to: $OUT_DIR"
ls -lh "$OUT_DIR"

# Hint for any .tar archives that landed.
if compgen -G "${OUT_DIR}/*.tar" >/dev/null; then
    echo ""
    echo "To extract a tar archive, run: tar -xvf <file>.tar -C \"$OUT_DIR\""
fi
