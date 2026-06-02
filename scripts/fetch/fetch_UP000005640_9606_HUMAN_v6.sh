#!/bin/bash

URL="https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar"
TAR_FILE="UP000005640_9606_HUMAN_v6.tar"
TARGET_DIR="UP000005640_9606_HUMAN_v6"
mkdir -p "$TARGET_DIR"

wget -c "$URL" -O "$TAR_FILE"
tar -xf "$TAR_FILE" -C "$TARGET_DIR"

