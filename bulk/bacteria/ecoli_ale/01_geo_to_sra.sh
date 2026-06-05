#!/usr/bin/env bash
# =============================================================================
# Script: 01_geo_to_sra.sh
# Description: Fetches SRA run tables and Sample information for a GEO accession.
# =============================================================================

# Strict mode for robust error handling
set -euo pipefail

# Default variables
BASE_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO=""
RUN_TABLE=""
SAMPLE_INFO=""
API_KEY=""

# =============================================================================
# Functions
# =============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") -g <GEO_ACCESSION> -r <RUN_TABLE_OUT> -s <SAMPLE_INFO_OUT> [OPTIONS]

This script retrieves the SRA Run Table and Sample Information for a given GEO 
accession using NCBI E-utilities.

Required Arguments:
  -g, --geo ACCESSION      GEO Accession ID (e.g., GSE242875)
  -r, --run-table FILE     Output filename for the SRA Run Table CSV
  -s, --sample-info FILE   Output filename for the Sample Info CSV

Optional Arguments:
  -k, --api-key KEY        NCBI API Key to increase rate limits from 3 to 10 requests/sec
  -u, --base-url URL       Base URL for NCBI E-utilities 
                           (Default: https://eutils.ncbi.nlm.nih.gov/entrez/eutils)
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -g GSE242875 -r SraRunTable.csv -s SampleInfo.csv -k YOUR_API_KEY
EOF
}

check_dependencies() {
    for cmd in curl python3; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not installed or not in PATH." >&2
            exit 1
        fi
    done
}

# =============================================================================
# Argument Parsing
# =============================================================================

# If no arguments provided, show help
if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) usage; exit 0 ;;
        -g|--geo) GEO="$2"; shift ;;
        -r|--run-table) RUN_TABLE="$2"; shift ;;
        -s|--sample-info) SAMPLE_INFO="$2"; shift ;;
        -k|--api-key) API_KEY="$2"; shift ;;
        -u|--base-url) BASE_URL="$2"; shift ;;
        *) echo "Error: Unknown parameter passed: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$GEO" || -z "$RUN_TABLE" || -z "$SAMPLE_INFO" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
    exit 1
fi

# =============================================================================
# Main Execution
# =============================================================================

check_dependencies

echo "Looking up GEO accession: $GEO"

# Format the API key parameter if provided
API_PARAM=""
if [[ -n "$API_KEY" ]]; then
    API_PARAM="&api_key=${API_KEY}"
    echo "Using NCBI API Key for increased rate limits."
fi

# Step 1: GEO accession → internal GEO ID
GEO_ID=$(curl -sg "${BASE_URL}/esearch.fcgi?db=gds&term=${GEO}[GSE]&retmode=json${API_PARAM}" \
    | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    ids = data.get('esearchresult', {}).get('idlist', [])
    valid_ids = [i for i in ids if len(i) == 9]
    if not valid_ids:
        print('Error: No valid 9-digit internal GEO ID found.', file=sys.stderr)
        sys.exit(1)
    print(valid_ids[0])
except Exception as e:
    print(f'Error parsing Step 1 JSON: {e}', file=sys.stderr)
    sys.exit(1)
")
echo "GEO internal ID: $GEO_ID"

# Step 2: GEO ID → internal SRA IDs
SRA_IDS=$(curl -sg "${BASE_URL}/elink.fcgi?dbfrom=gds&db=sra&id=${GEO_ID}&retmode=json${API_PARAM}" \
    | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    links = data['linksets'][0]['linksetdbs'][0]['links']
    if not links:
        sys.exit(1)
    print(','.join(links))
except Exception as e:
    print(f'Error parsing Step 2 JSON (No SRA links found?): {e}', file=sys.stderr)
    sys.exit(1)
")
echo "SRA internal IDs found."

# Step 3: SRA IDs → SraRunTable.csv (SRR accessions + run metadata)
curl -fsg "${BASE_URL}/efetch.fcgi?db=sra&id=${SRA_IDS}&rettype=runinfo&retmode=text${API_PARAM}" \
  -o "$RUN_TABLE"
echo "Saved: $RUN_TABLE"

# Step 4: GEO Series → SampleInfo.csv (GSM accession + sample title/genotype)
curl -sg "${BASE_URL}/esummary.fcgi?db=gds&id=${GEO_ID}&retmode=json${API_PARAM}" \
    | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print('LibraryName,LibraryDesc')
    for uid, record in data.get('result', {}).items():
        if uid == 'uids': continue
        for sample in record.get('samples', []):
            print(f\"{sample.get('accession', '')},{sample.get('title', '')}\")
except Exception as e:
    print(f'Error parsing Step 4 JSON: {e}', file=sys.stderr)
    sys.exit(1)
" > "$SAMPLE_INFO"
echo "Saved: $SAMPLE_INFO"

echo ""
echo "=== $RUN_TABLE (first 2 rows) ==="
head -2 "$RUN_TABLE"
echo ""
echo "=== $SAMPLE_INFO (first 2 rows) ==="
head -2 "$SAMPLE_INFO"