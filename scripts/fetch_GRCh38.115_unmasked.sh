#!/bin/bash

TARGET_DIRECTORY='/home/mahe9165/data/exonize/ensembl/v115/data'

mkdir -p "$TARGET_DIRECTORY"

# Homo_sapiens
## annotations
curl -L -o "$TARGET_DIRECTORY/Homo_sapiens.GRCh38.115.chr.gff3.gz" \
  "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gff3.gz"

## genome
curl -L -o "$TARGET_DIRECTORY/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" \
  "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

for file in "$TARGET_DIRECTORY"/*.gff3.gz; do
    gunzip "$file"
done