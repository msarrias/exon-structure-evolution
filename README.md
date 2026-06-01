## Setup
 
```bash
conda env create -f environment.yml
conda activate exon-structure-evolution
```
 
## Pipeline

### 1. Download Data

1.1 Download Ensembl genome and annotation:
```bash
bash scripts/fetch_GRCh38.115_unmasked.sh
```
 
1.2 Download UniProt accession mapping and AlphaFold structures:
```bash
bash scripts/fetch_uniprot_accession_mapping_v115_v6.sh
bash scripts/fetch_UP000005640_9606_HUMAN_v6.sh
```
 
1.3 Download BioMart mappings:
```bash
bash scripts/download_mart_mapping_v115.sh
bash scripts/download_mart_orthologs_v115.sh
```
 
### 2. Build AlphaFold database
 
```bash
python scripts/construct_uniprot_accessions_mapping_db.py
```
 
### 3. Run exonize
 
Constrained structural search:

```bash
bash scripts/exonize_GRCh38.116_80_90_80_const_struct.sh
```
 
Relaxed structural search:
```bash
bash scripts/exonize_GRCh38.116_80_90_80_relaxed_struct.sh
```
