import sqlite3
from pathlib import Path


def check_if_empty_table(
        db_path: Path,
        table_name: str,
) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        return cursor.fetchone()[0] == 0


def check_if_empty_column(
        db_path: Path,
        table_name: str,
        column_name: str,
) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"""
            SELECT COUNT(*) FROM {table_name}
            WHERE {column_name} IS NOT NULL
        """)
        return cursor.fetchone()[0] == 0


def check_if_table_exists(
        db_path: Path,
        table_name: str,
        timeout_db: int = 60,
) -> bool:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        return any(table_name == col[0] for col in cursor.fetchall())


def add_column_to_table(
        db_path: Path,
        table_name: str,
        column_name: str,
        column_type: str,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(
            f"ALTER TABLE {table_name} ADD COLUMN {column_name} {column_type}"
        )
        db.commit()


def copy_table_between_data_bases(
        target_db_path: Path,
        source_db_path: Path,
        table_name: str,
) -> None:
    with sqlite3.connect(target_db_path) as target_conn, \
         sqlite3.connect(source_db_path) as source_conn:
        target_cursor = target_conn.cursor()
        source_cursor = source_conn.cursor()
        source_cursor.execute(f"""
            SELECT sql FROM sqlite_master
            WHERE type='table' AND name="{table_name}"
        """)
        create_table_sql = source_cursor.fetchone()[0]
        source_cursor.execute(f"SELECT * FROM {table_name}")
        rows = source_cursor.fetchall()
        target_cursor.execute(create_table_sql)
        for row in rows:
            placeholders = ", ".join("?" * len(row))
            target_cursor.execute(
                f"INSERT INTO {table_name} VALUES ({placeholders})", row
            )


def fetch_distinct_genes(
        cursor: sqlite3.Cursor,
        table_name: str,
        condition: str = None,
) -> set:
    query = f"SELECT DISTINCT GeneID FROM {table_name}"
    if condition:
        query += f" WHERE {condition}"
    cursor.execute(query)
    return {record[0] for record in cursor.fetchall()}

def query_genes(
        db_path: Path,
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            DISTINCT GeneID
        FROM Genes;
        """)
        # tailored for ensemble gff genome id format e.g. gene:ENSG00000186092
        return [i[0].rsplit(':')[1] for i in cursor.fetchall()]
