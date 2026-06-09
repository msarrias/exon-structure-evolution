#!/bin/bash

BASE=/home/mahe9165/exon-structure-evolution/data/ensembl/v115

python3 build_protein_db.py \
    --gene-hierarchy-dictionary-path $BASE/GRCh38.115_gene_hierarchy.pkl \
    --protein-db-path $BASE/GRCh38.115_protein.db \
    --genome-path $BASE/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --timeout-db 30