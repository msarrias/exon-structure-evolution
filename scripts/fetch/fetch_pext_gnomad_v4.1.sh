#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR=/home/mahe9165/exon-structure-evolution/data/pext
URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/pext/gnomad.pext.gtex_v10.base_level.tsv.gz"
OUTPUT="$OUTPUT_DIR/gnomad.pext.gtex_v10.base_level.tsv.gz"

mkdir -p "$OUTPUT_DIR"

echo "Downloading gnomAD v4.1 pext annotations (GTEx v10 base level)..."

wget -c \
    --progress=bar:force \
    -O "$OUTPUT" \
    "$URL"

echo "Done — saved to $OUTPUT"
echo "Size: $(du -sh "$OUTPUT" | cut -f1)"