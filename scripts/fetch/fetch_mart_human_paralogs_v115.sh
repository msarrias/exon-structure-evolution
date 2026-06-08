#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="/home/mahe9165/exon-structure-evolution/data/biomart/tab_data"
OUTPUT="$OUTPUT_DIR/hsapiens_paralogs.txt"
mkdir -p "$OUTPUT_DIR"

XML='
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="hsapiens_gene_ensembl" interface="default">
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="hsapiens_paralog_ensembl_gene"/>
        <Attribute name="hsapiens_paralog_subtype"/>
        <Attribute name="hsapiens_paralog_orthology_type"/>
        <Attribute name="hsapiens_paralog_perc_id"/>
        <Attribute name="hsapiens_paralog_perc_id_r1"/>
    </Dataset>
</Query>
'

ENCODED=$(python3 -c "import urllib.parse, sys; print(urllib.parse.quote(sys.argv[1]))" "$XML")

curl -s -L \
    "https://www.ensembl.org/biomart/martservice?query=${ENCODED}" \
    -o "${OUTPUT}"

# Validate response
FIRST_LINE=$(head -1 "${OUTPUT}")
if [[ "${FIRST_LINE}" != *"Gene stable ID"* ]]; then
    echo "ERROR: unexpected response from BioMart:"
    head -3 "${OUTPUT}"
    exit 1
fi

ROWS=$(tail -n +2 "${OUTPUT}" | wc -l)
echo "Done — saved to ${OUTPUT} (${ROWS} rows)"
