#!/usr/bin/env bash
# =============================================================================
# Script: 01_geo_to_sra.sh
# Description: Fetches SRA run tables and Sample information for a GEO accession.
#
# Fixes vs original:
#   - Bug 1 (Step 2): linksetdbs can be empty for series with no direct SRA
#     links at the series level. Original did [0] unconditionally → IndexError
#     → script exited before writing SampleInfo.csv. Fixed with a safe check
#     and a fallback: fetch SRA IDs via individual GSM→SRA elinks.
#   - Bug 2 (Step 4): some large/SuperSeries return samples=[] in esummary.
#     Fixed by falling back to a separate elink gds→gds call to retrieve
#     GSM-level records when the samples list is empty.
# =============================================================================

set -euo pipefail

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

Retrieves the SRA Run Table and Sample Information for a given GEO accession
using NCBI E-utilities.

Required Arguments:
  -g, --geo ACCESSION      GEO Accession ID (e.g., GSE242875)
  -r, --run-table FILE     Output filename for the SRA Run Table CSV
  -s, --sample-info FILE   Output filename for the Sample Info CSV

Optional Arguments:
  -k, --api-key KEY        NCBI API Key (increases rate limit from 3 to 10 req/sec)
  -u, --base-url URL       Base URL for NCBI E-utilities
                           (Default: https://eutils.ncbi.nlm.nih.gov/entrez/eutils)
  -h, --help               Show this help message and exit

Example:
  $(basename "$0") -g GSE242875 -r SraRunTable.csv -s SampleInfo.csv
  $(basename "$0") -g GSE135516 -r SraRunTable.csv -s SampleInfo.csv -k YOUR_API_KEY
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

if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)        usage; exit 0 ;;
        -g|--geo)         GEO="$2"; shift ;;
        -r|--run-table)   RUN_TABLE="$2"; shift ;;
        -s|--sample-info) SAMPLE_INFO="$2"; shift ;;
        -k|--api-key)     API_KEY="$2"; shift ;;
        -u|--base-url)    BASE_URL="$2"; shift ;;
        *) echo "Error: Unknown parameter: $1" >&2; usage; exit 1 ;;
    esac
    shift
done

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

API_PARAM=""
if [[ -n "$API_KEY" ]]; then
    API_PARAM="&api_key=${API_KEY}"
    echo "Using NCBI API Key for increased rate limits."
fi

# ── Step 1: GEO accession → internal GEO ID ───────────────────────────────────
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

# ── Step 2: GEO ID → SRA IDs ──────────────────────────────────────────────────
# FIX (Bug 1): linksetdbs can be empty for series without direct SRA links at
# the series level (e.g. large studies, SuperSeries). Original code did
# linksetdbs[0] unconditionally — IndexError → sys.exit(1) → script died here,
# so SampleInfo.csv was never written.
#
# Fix: check if linksetdbs is non-empty before indexing. If empty, fall back to
# fetching SRA IDs via individual GSM accessions linked to the series.

echo "Fetching SRA IDs..."

ELINK_RESPONSE=$(curl -sg \
    "${BASE_URL}/elink.fcgi?dbfrom=gds&db=sra&id=${GEO_ID}&retmode=json${API_PARAM}")

SRA_IDS=$(echo "$ELINK_RESPONSE" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    linksetdbs = data.get('linksets', [{}])[0].get('linksetdbs', [])

    # Safe check — original code did [0] here without guarding
    if not linksetdbs:
        print('', end='')   # signal empty to the shell
        sys.exit(0)

    links = linksetdbs[0].get('links', [])
    if not links:
        print('', end='')
        sys.exit(0)

    print(','.join(links))
except Exception as e:
    print(f'Error parsing Step 2 JSON: {e}', file=sys.stderr)
    sys.exit(1)
")

# Fallback: if series-level elink returned nothing, fetch SRA IDs via GSM records
if [[ -z "$SRA_IDS" ]]; then
    echo "No direct SRA links at series level — falling back to GSM-level SRA lookup..."

    # Get all GSM IDs linked to this series
    GSM_IDS=$(curl -sg \
        "${BASE_URL}/elink.fcgi?dbfrom=gds&db=gds&id=${GEO_ID}&retmode=json${API_PARAM}" \
        | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    linksetdbs = data.get('linksets', [{}])[0].get('linksetdbs', [])
    if not linksetdbs:
        print('Error: No GSM links found for this series.', file=sys.stderr)
        sys.exit(1)
    links = linksetdbs[0].get('links', [])
    # Exclude the series ID itself (same as GEO_ID)
    gsm_ids = [l for l in links if l != '$GEO_ID']
    if not gsm_ids:
        print('Error: GSM link list is empty.', file=sys.stderr)
        sys.exit(1)
    # Batch to first 500 to avoid URL length limits; RunTable fetch handles the rest
    print(','.join(gsm_ids[:500]))
except Exception as e:
    print(f'Error parsing GSM elink JSON: {e}', file=sys.stderr)
    sys.exit(1)
")

    SRA_IDS=$(curl -sg \
        "${BASE_URL}/elink.fcgi?dbfrom=gds&db=sra&id=${GSM_IDS}&retmode=json${API_PARAM}" \
        | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    all_links = set()
    for linkset in data.get('linksets', []):
        for ldb in linkset.get('linksetdbs', []):
            all_links.update(ldb.get('links', []))
    if not all_links:
        print('Error: No SRA IDs found via GSM fallback.', file=sys.stderr)
        sys.exit(1)
    print(','.join(sorted(all_links)))
except Exception as e:
    print(f'Error parsing GSM→SRA JSON: {e}', file=sys.stderr)
    sys.exit(1)
")
    echo "SRA IDs found via GSM fallback."
else
    echo "SRA IDs found via series-level link."
fi

# ── Step 3: SRA IDs → SraRunTable.csv ─────────────────────────────────────────
echo "Fetching run table..."
curl -fsg \
    "${BASE_URL}/efetch.fcgi?db=sra&id=${SRA_IDS}&rettype=runinfo&retmode=text${API_PARAM}" \
    -o "$RUN_TABLE"
echo "Saved: $RUN_TABLE"

# ── Step 4: GEO Series → SampleInfo.csv ───────────────────────────────────────
# FIX (Bug 2): some large/SuperSeries return samples=[] in esummary at the top
# series level. Original code only iterated record['samples'] — produced a
# header-only CSV for these accessions.
#
# Fix: if samples list is empty after parsing the series record, fall back to
# fetching GSM-level summaries via a targeted esummary on linked GSM IDs.

echo "Fetching sample info..."

ESUMMARY_RESPONSE=$(curl -sg \
    "${BASE_URL}/esummary.fcgi?db=gds&id=${GEO_ID}&retmode=json${API_PARAM}")

python3 - <<PYEOF > "$SAMPLE_INFO"
import sys, json, urllib.request, urllib.parse

base_url  = "${BASE_URL}"
geo_id    = "${GEO_ID}"
api_param = "${API_PARAM}"

raw = '''${ESUMMARY_RESPONSE}'''
try:
    data = json.loads(raw)
except Exception as e:
    print(f"Error parsing Step 4 JSON: {e}", file=sys.stderr)
    sys.exit(1)

print("LibraryName,LibraryDesc")

rows_written = 0
for uid, record in data.get("result", {}).items():
    if uid == "uids":
        continue
    for sample in record.get("samples", []):
        acc   = sample.get("accession", "").replace(",", " ")
        title = sample.get("title",     "").replace(",", " ")
        print(f"{acc},{title}")
        rows_written += 1

# Fallback: samples list was empty — fetch GSM records individually
if rows_written == 0:
    print("Sample list empty at series level — fetching via GSM elink...", file=sys.stderr)

    # Re-fetch the GSM IDs linked to this series
    url = (f"{base_url}/elink.fcgi?dbfrom=gds&db=gds"
           f"&id={geo_id}&retmode=json{api_param}")
    with urllib.request.urlopen(url) as resp:
        elink_data = json.load(resp)

    linksetdbs = elink_data.get("linksets", [{}])[0].get("linksetdbs", [])
    if not linksetdbs:
        print("Error: No GSM links found in fallback.", file=sys.stderr)
        sys.exit(1)

    gsm_ids = [l for l in linksetdbs[0].get("links", []) if l != geo_id]
    if not gsm_ids:
        print("Error: GSM link list empty in fallback.", file=sys.stderr)
        sys.exit(1)

    # Fetch in batches of 100 to stay within URL length limits
    batch_size = 100
    for i in range(0, len(gsm_ids), batch_size):
        batch = ",".join(gsm_ids[i:i + batch_size])
        url = (f"{base_url}/esummary.fcgi?db=gds"
               f"&id={batch}&retmode=json{api_param}")
        with urllib.request.urlopen(url) as resp:
            batch_data = json.load(resp)

        for uid, rec in batch_data.get("result", {}).items():
            if uid == "uids":
                continue
            # GSM-level records store title directly, no nested samples[]
            acc   = rec.get("accession", "").replace(",", " ")
            title = rec.get("title",     "").replace(",", " ")
            if acc.startswith("GSM"):
                print(f"{acc},{title}")
                rows_written += 1

if rows_written == 0:
    print("Warning: SampleInfo.csv contains only the header — no samples found.", file=sys.stderr)
PYEOF

echo "Saved: $SAMPLE_INFO"

# ── Summary preview ────────────────────────────────────────────────────────────
echo ""
echo "=== $RUN_TABLE (first 2 rows) ==="
head -2 "$RUN_TABLE"
echo ""
echo "=== $SAMPLE_INFO (first 2 rows) ==="
head -2 "$SAMPLE_INFO"
