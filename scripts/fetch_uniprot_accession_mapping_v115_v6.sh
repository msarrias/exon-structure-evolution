#!/usr/bin/env bash
set -euo pipefail

OUTPUT="uniprotkb_proteome_UP000005640_reviewed_$(date +%Y_%m_%d).tsv"

echo "Downloading UniProt Swiss-Prot reviewed entries for UP000005640_9606_HUMAN_v6"

curl -s -L \
    "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640+AND+reviewed:true&format=tsv&fields=accession,id,length,xref_ensembl" \
    -o "${OUTPUT}"

ROWS=$(tail -n +2 "${OUTPUT}" | wc -l)
echo "Done — saved to ${OUTPUT} (${ROWS} entries)"