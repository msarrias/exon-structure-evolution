import csv
import logging
import sqlite3
from pathlib import Path


def create_database(db_path: Path):
    with sqlite3.connect(db_path) as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS Human_paralogs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_stable_id VARCHAR NOT NULL,
                paralogue_gene_stable_id VARCHAR,
                last_common_ancestor_with_human VARCHAR,
                paralogue_homology_type VARCHAR,
                paralogue_target_identity REAL,
                paralogue_query_identity REAL
            )
        """)
    logger.info("Database schema created")


def parse_tsv(path: Path) -> list[tuple]:
    rows = []
    skipped = 0
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                target_id = float(row["Paralogue %id. target Human gene identical to query gene"]) \
                    if row["Paralogue %id. target Human gene identical to query gene"].strip() else None
                query_id  = float(row["Paralogue %id. query gene identical to target Human gene"]) \
                    if row["Paralogue %id. query gene identical to target Human gene"].strip() else None
            except ValueError:
                skipped += 1
                continue
            rows.append((
                row["Gene stable ID"] or None,
                row["Human paralogue gene stable ID"] or None,
                row["Paralogue last common ancestor with Human"] or None,
                row["Human paralogue homology type"] or None,
                target_id,
                query_id,
            ))
    if skipped:
        logger.warning(f"Skipped {skipped:,} rows with unparseable identity values")
    return rows


def main():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Input:    {INPUT_TSV}")
    logger.info(f"Database: {DB_PATH}")

    create_database(DB_PATH)

    logger.info("Parsing TSV...")
    rows = parse_tsv(INPUT_TSV)
    logger.info(f"Parsed {len(rows):,} rows")

    logger.info("Inserting into database...")
    with sqlite3.connect(DB_PATH, timeout=60) as conn:
        conn.executemany("""
            INSERT INTO Human_paralogs (
                gene_stable_id,
                paralogue_gene_stable_id,
                last_common_ancestor_with_human,
                paralogue_homology_type,
                paralogue_target_identity,
                paralogue_query_identity
            ) VALUES (?, ?, ?, ?, ?, ?)
        """, rows)
        conn.commit()

    logger.info(f"Done — inserted {len(rows):,} rows into Human_paralogs")


if __name__ == "__main__":
    logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("build_paralogs_db.log"),
        ]
    )
    logger = logging.getLogger(__name__)
    
    INPUT_TSV = Path("/home/mahe9165/exon-structure-evolution/data/biomart/tab_data/mart_hsapiens_paralogs_data_v115.txt")
    DB_PATH = Path("/home/mahe9165/exon-structure-evolution/data/biomart/mart_hsapiens_paralogs_v115.db")

    main()