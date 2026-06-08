#!/usr/bin/env bash
set -euo pipefail

DATA=/home/mahe9165/exon-structure-evolution/data
PFAM_DIR=$DATA/pfam
ECOD_DIR=$DATA/ecod
mkdir -p "$PFAM_DIR" "$ECOD_DIR"

#  Pfam v37.0 HMM profiles 
if [[ ! -f "$PFAM_DIR/Pfam-A.hmm" ]]; then
    echo "Downloading Pfam-A HMM profiles (v37.0)..."
    wget -c --progress=bar:force \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" \
        -O "$PFAM_DIR/Pfam-A.hmm.gz"
    gunzip -f "$PFAM_DIR/Pfam-A.hmm.gz"
else
    echo "Pfam-A.hmm already exists - skipping download"
fi

if [[ ! -f "$PFAM_DIR/Pfam-A.hmm.dat" ]]; then
    echo "Downloading Pfam-A metadata..."
    wget -c --progress=bar:force \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz" \
        -O "$PFAM_DIR/Pfam-A.hmm.dat.gz"
    gunzip -f "$PFAM_DIR/Pfam-A.hmm.dat.gz"
else
    echo "Pfam-A.hmm.dat already exists - skipping download"
fi

# hmmpress creates 4 index files
if [[ ! -f "$PFAM_DIR/Pfam-A.hmm.h3i" ]]; then
    echo "Pressing Pfam HMM database for hmmscan..."
    hmmpress "$PFAM_DIR/Pfam-A.hmm"
else
    echo "Pfam-A.hmm already pressed - skipping"
fi

#  ECOD F-group to Pfam mapping 
if [[ ! -f "$ECOD_DIR/ecod.latest.f_id_pfam_acc.txt" ]]; then
    echo "Downloading ECOD F-group to Pfam mapping..."
    wget -c --progress=bar:force \
        "http://prodata.swmed.edu/ecod/distributions/ecod.latest.f_id_pfam_acc.txt" \
        -O "$ECOD_DIR/ecod.latest.f_id_pfam_acc.txt"
else
    echo "ecod.latest.f_id_pfam_acc.txt already exists - skipping"
fi

#  ECOD domain definitions 
if [[ ! -f "$ECOD_DIR/ecod.latest.domains.txt" ]]; then
    echo "Downloading ECOD domain definitions..."
    wget -c --progress=bar:force \
        "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt" \
        -O "$ECOD_DIR/ecod.latest.domains.txt"
else
    echo "ecod.latest.domains.txt already exists - skipping"
fi

echo ""
echo "Done."
