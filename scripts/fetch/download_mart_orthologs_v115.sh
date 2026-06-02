#!/usr/bin/env bash
set -euo pipefail

BIOMART_URL="https://www.ensembl.org/biomart/martservice"
OUTPUT_DIR="/home/mahe9165/exonize/data/biomart/tab_data/eukaryotes18_orthologs_data_v115"
mkdir -p "$OUTPUT_DIR"

species=(
    hsapiens
    mmusculus
    clfamiliaris
    btaurus
    sscrofa
    mdomestica
    oanatinus
    tguttata
    ggallus
    acarolinensis
    xtropicalis
    trubripes
    drerio
    cmilii
    eburgeri
    dmelanogaster
    celegans
    scerevisiae
)

encode_query() {
    python3 -c "import urllib.parse, sys; print(urllib.parse.quote(sys.stdin.read()))"
}

for specie_i in "${species[@]}"; do
    for specie_j in "${species[@]}"; do
        if [ "$specie_i" == "$specie_j" ]; then
            continue
        fi

        output_file="${OUTPUT_DIR}/${specie_i}_to_${specie_j}_orthologs.txt"

        if [[ -f "$output_file" ]]; then
            echo "Already downloaded: ${specie_i} → ${specie_j}"
            continue
        fi

        echo "Querying: ${specie_i} → ${specie_j}"

        ENCODED_QUERY=$(cat <<EOF | tr -d '\n' | encode_query
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name="${specie_i}_gene_ensembl" interface="default">
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="${specie_j}_homolog_ensembl_gene"/>
        <Attribute name="${specie_j}_homolog_subtype"/>
        <Attribute name="${specie_j}_homolog_orthology_type"/>
        <Attribute name="${specie_j}_homolog_perc_id"/>
        <Attribute name="${specie_j}_homolog_perc_id_r1"/>
    </Dataset>
</Query>
EOF
)
        wget -q -O "$output_file" "${BIOMART_URL}?query=${ENCODED_QUERY}"

        # Check we got a valid TSV
        FIRST_LINE=$(head -1 "$output_file" 2>/dev/null || echo "")
        if [[ "${FIRST_LINE}" != *"Gene stable ID"* ]]; then
            echo "ERROR: unexpected response for ${specie_i} → ${specie_j}:"
            head -3 "$output_file"
            rm -f "$output_file"
        else
            ROWS=$(tail -n +2 "$output_file" | wc -l)
            echo "${ROWS} rows"
        fi

        sleep 1
    done
done

echo "Done."
