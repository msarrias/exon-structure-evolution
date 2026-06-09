#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/mahe9165/exon-structure-evolution
DATA=$ROOT/data
SCRIPTS=$ROOT/scripts

PROTEIN_DB=$DATA/ensembl/v115/protein_db_v115.db
EXONIZE_DB=$DATA/human_structural_search/Homo_sapiens.GRCh38.115.80.90.80.struct_exonize/Homo_sapiens.GRCh38.115.80.90.80.struct_results.db
GENE_HIERARCHY=$DATA/ensembl/v115/GRCh38.115_gene_hierarchy.pkl

OUTPUT_DB=$DATA/aiupred/aiupred_GRCh38_v115.db
mkdir -p "$DATA/aiupred"

python $SCRIPTS/analysis/aiupred_idp_search.py \
    $PROTEIN_DB \
    $EXONIZE_DB \
    $GENE_HIERARCHY \
    --output-db-path $OUTPUT_DB \
    --anchor \
    --min-exon-length 30 \
    --ordered-score-threshold 0.5