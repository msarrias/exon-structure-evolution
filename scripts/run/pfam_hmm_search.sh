#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/mahe9165/exon-structure-evolution
DATA=$ROOT/data
SCRIPTS=$ROOT/scripts

# --- Input data ---
ECOD_DOMAINS=$DATA/ecod/ecod.latest.domains.txt
ECOD_PFAM_MAP=$DATA/ecod/ecod.latest.f_id_pfam_acc.txt
GENOME=$DATA/ensembl/v115/data/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
GENE_HIERARCHY=$DATA/ensembl/v115/GRCh38.115_gene_hierarchy.pkl
EXONIZE_DB=$DATA/human_structural_search/Homo_sapiens.GRCh38.115.80.90.80.struct_exonize/Homo_sapiens.GRCh38.115.80.90.80.struct_results.db

# --- Pfam ---
PFAM_HMM=$DATA/pfam/Pfam-A.hmm
HMMSCAN_CDS=$DATA/pfam/hmmscan_cdss.domtblout
HMMSCAN_MRNA=$DATA/pfam/hmmscan_mrna.domtblout

# --- Intermediate files ---
FASTA_CDS=$DATA/fasta/GRCh38_v115_cdss.fasta
FASTA_MRNA=$DATA/fasta/GRCh38_v115_mrna.fasta
mkdir -p "$DATA/fasta"

# --- Output databases ---
PROTEIN_DB=$DATA/ensembl/v115/protein_db_v115.db
ECOD_DB=$DATA/ecod/ecod_pfam_results_v115.db

python $SCRIPTS/analysis/pfam_hmm_search.py \
    --ecod-latest-domains-txt-path $ECOD_DOMAINS \
    --ecod-latest-domains-pfam-map-path $ECOD_PFAM_MAP \
    --protein-db-path $PROTEIN_DB \
    --ecod-db-path $ECOD_DB \
    --exonize-db-path $EXONIZE_DB \
    --genome-file-path $GENOME \
    --genome-hierarchy-pckl-path $GENE_HIERARCHY \
    --fasta-file-path-cdss $FASTA_CDS \
    --fasta-file-path-mrna $FASTA_MRNA \
    --hmm-profile-path $PFAM_HMM \
    --hmmsearch-cdss-path $HMMSCAN_CDS \
    --hmmsearch-mrna-path $HMMSCAN_MRNA
    