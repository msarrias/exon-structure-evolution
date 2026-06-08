# exon-structure-evolution

## Pipeline

### 1. Download data

1.1 Download Ensembl human genome and gff3 annotation v115:
```bash
bash scripts/fetch/fetch_GRCh38.115_unmasked.sh
```

1.2 Download UniProt accession mapping and AlphaFold structures:
```bash
bash scripts/fetch/fetch_uniprot_accession_mapping_v115_v6.sh
bash scripts/fetch/fetch_UP000005640_9606_HUMAN_v6.sh
```

1.3 Download human canonical transcript mapping from BioMart:
```bash
bash scripts/fetch/fetch_mart_human_canonical_transcripts_v115.sh
```

1.4 Download BioMart ortholog data (21 species):
```bash
bash scripts/fetch/fetch_mart_21_species_orthologs_v115.sh
```

1.5 Download BioMart paralog data (21 species):
```bash
bash scripts/fetch/fetch_mart_21_species_paralogs_v115.sh
```

1.6 Download pext score annotations:
```bash
bash scripts/fetch/fetch_pext_gnomad_v4.1.sh
```

### 2. Build databases

Build UniProt/AlphaFold mapping database:
```bash
python scripts/build/build_uniprot_accessions_mapping_db.py
```

Build BioMart ortholog database:
```bash
python scripts/build/build_biomart_orthologs_db.py
```

Build BioMart paralog database:
```bash
python scripts/build/build_biomart_human_paralogs_db.py
```

### 3. Run exonize

Constrained structural search:
```bash
bash scripts/exonize/exonize_GRCh38.116_80_90_80_const_struct.sh
```

Relaxed structural search:
```bash
bash scripts/exonize/exonize_GRCh38.116_80_90_80_relaxed_struct.sh
```