#!/bin/bash -l

DATA=/home/mahe9165/exon-structure-evolution/data
ENSEMBL=$DATA/ensembl/v115
ALPHAFOLD=$DATA/alphafold
OUTPUT=/home/mahe9165/exon-structure-evolution/human_structural_search

# --- Input data ---
GFF_FILE_PATH=$ENSEMBL/Homo_sapiens.GRCh38.115.gff3
GENOME_FILE_PATH=$ENSEMBL/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

# --- AlphaFold structures ---
PDB_STRUCTURES_PATH=$ALPHAFOLD/UP000005640_9606_HUMAN_v6
PDB_STRUCTURES_MAPPING_PATH=$ALPHAFOLD/af2_uniprot_accession_mapping_v115.db

# --- Output ---
SPECIE_ID=Homo_sapiens.GRCh38.115.80.90.80.struct

# --- Annotation features ---
CDS_ANNOT='CDS'
GENE_ANNOT='gene'
TRANS_ANNOT='mRNA'

# --- Sequence search parameters ---
MIN_SEQ_LENGTH=30
EVALUE=0.001
QUERY_COVERAGE=0.8
EXON_CLUSTERING_CUTOFF=0.9
TARGET_CLUSTERING_CUTOFF=0.8

# --- Structural search parameters ---
AF2_V=6
MIN_EXON_LEN_STRUCT=30
RMSD=2.0
TM_SCORE=0.6
PLDDT=70

CPUS=30

exonize $GFF_FILE_PATH \
  $GENOME_FILE_PATH \
  --gene-annot-feature $GENE_ANNOT \
  --cds-annot-feature $CDS_ANNOT \
  --transcript-annot-feature $TRANS_ANNOT \
  --output-prefix $SPECIE_ID \
  --output-directory-path $OUTPUT \
  --min-exon-length $MIN_SEQ_LENGTH \
  --evalue-threshold $EVALUE \
  --query-coverage-threshold $QUERY_COVERAGE \
  --exon-clustering-overlap-threshold $EXON_CLUSTERING_CUTOFF \
  --targets-clustering-overlap-threshold $TARGET_CLUSTERING_CUTOFF \
  --cpus $CPUS \
  --structural-search \
  --local-search \
  --global-search \
  --min-exon-length-structural-search $MIN_EXON_LEN_STRUCT \
  --RMSD-threshold $RMSD \
  --TM-score-threshold $TM_SCORE \
  --plddt-threshold $PLDDT \
  --pdb-structures-path $PDB_STRUCTURES_PATH \
  --pdb-ids-mapping-db-path $PDB_STRUCTURES_MAPPING_PATH \
  --af2-v $AF2_V