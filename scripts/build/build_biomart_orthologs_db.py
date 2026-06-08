import csv
import logging
import sqlite3
from pathlib import Path


def parse_species(filename: str) -> tuple[str, str]:
    """Extract species_i and species_j from filename like 'hsapiens_to_mmusculus_orthologs.txt'."""
    parts = filename.split("_to_")
    species_i = parts[0]
    species_j = parts[1].replace("_orthologs.txt", "")
    return species_i, species_j


def build_table_name(species_i: str, species_j: str) -> str:
    return f"{species_i}_{species_j}_orthologs"


def load_tsv(path: Path) -> list[tuple]:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header
        rows = []
        for row in reader:
            if len(row) != len(COLUMNS):
                logger.warning(f"Skipping malformed row in {path.name}: {row}")
                continue
            # cast REAL columns
            parsed = []
            for val, (_, col_type) in zip(row, COLUMNS):
                if col_type == "REAL":
                    try:
                        parsed.append(float(val) if val.strip() else None)
                    except ValueError:
                        parsed.append(None)
                else:
                    parsed.append(val if val.strip() else None)
            rows.append(tuple(parsed))
    return rows


def main():
    biomart_files = sorted(
        f for f in BIOMART_DIR.iterdir()
        if f.is_file() and f.name.endswith(".txt")
    )
    if not biomart_files:
        logger.error(f"No .txt files found in {BIOMART_DIR}")
        return

    logger.info(f"Found {len(biomart_files)} BioMart files")
    logger.info(f"Output database: {DB_PATH}")
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)

    with sqlite3.connect(DB_PATH, timeout=60) as conn:
        cursor = conn.cursor()
        for idx, path in enumerate(biomart_files):
            try:
                species_i, species_j = parse_species(path.name)
            except (IndexError, ValueError) as e:
                logger.warning(f"Could not parse species from {path.name}: {e} — skipping")
                continue

            table = build_table_name(species_i, species_j)

            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS {table} (
                    {SCHEMA}
                )
            """)

            rows = load_tsv(path)
            if not rows:
                logger.warning(f"No rows loaded from {path.name}")
                continue

            cursor.executemany(
                f"INSERT INTO {table} ({', '.join(COL_NAMES)}) VALUES ({PLACEHOLDERS})",
                rows
            )
            conn.commit()
    logger.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("build_orthologs_db.log"),
    ]
    )
    logger = logging.getLogger(__name__)
    
    DB_PATH = Path("/home/mahe9165/exon-structure-evolution/data/biomart/eukaryotes21_orthologs.db")
    BIOMART_DIR = Path("/home/mahe9165/exon-structure-evolution/data/biomart/tab_data/eukaryotes21_orthologs_data_v115")
    
    COLUMNS = [
        ("gene_id", "TEXT"),
        ("homolog_gene_id", "TEXT"),
        ("homolog_subtype", "TEXT"),
        ("homology_type", "TEXT"),
        ("ident_query_to_target_gene", "REAL"),
        ("ident_target_to_query_gene", "REAL"),
    ]
    COL_NAMES  = [c for c, _ in COLUMNS]
    COL_TYPES  = [t for _, t in COLUMNS]
    SCHEMA     = ", ".join(f"{c} {t}" for c, t in COLUMNS)
    PLACEHOLDERS = ", ".join(["?"] * len(COLUMNS))
    main()
    