#!/usr/bin/env bash
set -euo pipefail

BIOMART_URL="https://www.ensembl.org/biomart/martservice"
OUTPUT_DIR="/home/mahe9165/exon-structure-evolution/data/biomart/tab_data/eukaryotes21_paralogs_data_v115"
MAX_RETRIES=5
RETRY_WAIT=60

mkdir -p "$OUTPUT_DIR"

species=(
    hsapiens
    mmusculus
    clfamiliaris
    btaurus
    sscrofa
    mdomestica
    # oanatinus
    tguttata
    ggallus
    acarolinensis
    xtropicalis
    trubripes
    drerio
    cmilii
    eburgeri
    pmarinus
    lchalumnae
    cintestinalis
    dmelanogaster
    celegans
    scerevisiae
)

encode_query() {
    python3 -c "import urllib.parse; print(urllib.parse.quote(open('$1').read().strip()))"
}

fetch_paralogs() {
    local species=$1
    local output_file=$2

    # Write XML to temp file to avoid shell quoting issues
    local tmp_xml
    tmp_xml=$(mktemp)
    cat > "$tmp_xml" << XML
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="${species}_gene_ensembl" interface="default">
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="${species}_paralog_ensembl_gene"/>
        <Attribute name="${species}_paralog_subtype"/>
        <Attribute name="${species}_paralog_orthology_type"/>
        <Attribute name="${species}_paralog_perc_id"/>
        <Attribute name="${species}_paralog_perc_id_r1"/>
    </Dataset>
</Query>
XML

    local encoded
    encoded=$(encode_query "$tmp_xml")
    rm -f "$tmp_xml"

    local attempt=0
    while [ $attempt -lt $MAX_RETRIES ]; do
        attempt=$((attempt + 1))

        wget -q -O "$output_file" "${BIOMART_URL}?query=${encoded}"

        local first_line
        first_line=$(head -1 "$output_file" 2>/dev/null || echo "")
        if [[ "${first_line}" == *"Gene stable ID"* ]]; then
            local rows
            rows=$(tail -n +2 "$output_file" | wc -l)
            echo "  done - ${rows} rows (attempt ${attempt})"
            return 0
        else
            echo "  attempt ${attempt}/${MAX_RETRIES} failed"
            rm -f "$output_file"
            if [ $attempt -lt $MAX_RETRIES ]; then
                echo "  waiting ${RETRY_WAIT}s..."
                sleep $RETRY_WAIT
            fi
        fi
    done

    echo "  FAILED after ${MAX_RETRIES} attempts"
    return 1
}

for species in "${species[@]}"; do
    output_file="${OUTPUT_DIR}/${species}_paralogs.txt"

    if [[ -f "$output_file" ]]; then
        echo "Already downloaded: ${species}"
        continue
    fi

    echo "Querying: ${species}..."
    fetch_paralogs "$species" "$output_file" || true
    sleep 2
done

echo "Done."