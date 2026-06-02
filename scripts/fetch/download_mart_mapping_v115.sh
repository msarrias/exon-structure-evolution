#!/usr/bin/env bash
set -euo pipefail

OUTPUT="../tab_data/mart_export_canonical_transcripts_v115.txt"

XML='
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="uniprotswissprot"/>
    <Attribute name="transcript_is_canonical"/>
  </Dataset>
</Query>'

ENCODED=$(python3 -c "import urllib.parse, sys; print(urllib.parse.quote(sys.argv[1]))" "$XML")

echo "Downloading canonical transcript mapping from Ensembl 115 BioMart..."

curl -s -L \
    "https://www.ensembl.org/biomart/martservice?query=${ENCODED}" \
    -o "${OUTPUT}"

# Check we got a valid TSV
FIRST_LINE=$(head -1 "${OUTPUT}")
if [[ "${FIRST_LINE}" != *"Gene stable ID"* ]]; then
    echo "ERROR"
    head -3 "${OUTPUT}"
    exit 1
fi
