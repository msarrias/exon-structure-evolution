import sqlite3
from collections import defaultdict
import portion as P


def check_if_table_exists(db_path: str, table_name: str) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        return any(table_name == col[0] for col in cursor.fetchall())


def check_if_column_exists(db_path: str, table_name: str, column_name: str) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"PRAGMA table_info({table_name})")
        return any(column_name == col[1] for col in cursor.fetchall())


def check_if_empty_table(db_path: str, table_name: str) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        return cursor.fetchone()[0] == 0


def add_column_to_table(
        db_path: str,
        table_name: str,
        column_name: str,
        column_type: str,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {column_type}")


def copy_table_between_data_bases(
        target_db_path: str,
        source_db_path: str,
        table_name: str,
) -> None:
    with sqlite3.connect(target_db_path) as target_conn, sqlite3.connect(source_db_path) as source_conn:
        target_cursor = target_conn.cursor()
        source_cursor = source_conn.cursor()
        source_cursor.execute(
            f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}'"
        )
        create_table_sql = source_cursor.fetchone()[0]
        source_cursor.execute(f"SELECT * FROM {table_name}")
        rows = source_cursor.fetchall()
        target_cursor.execute(create_table_sql)
        for row in rows:
            placeholders = ", ".join("?" * len(row))
            target_cursor.execute(f"INSERT INTO {table_name} VALUES ({placeholders})", row)


def create_index_on_table(db_path: str, table_name: str, column_name: str) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute(f"CREATE INDEX IF NOT EXISTS {column_name}_index ON {table_name}({column_name})")


def create_aiupred_table(db_path: str, anchor: bool) -> None:
    with sqlite3.connect(db_path) as db:
        if anchor:
            db.execute("""
                CREATE TABLE IF NOT EXISTS aiupred_protein_sequences_scores (
                    GeneID TEXT NOT NULL,
                    TranscriptID VARCHAR(100) PRIMARY KEY,
                    IUPRED2 TEXT NOT NULL,
                    ANCHOR2 TEXT NOT NULL
                )
            """)
        else:
            db.execute("""
                CREATE TABLE IF NOT EXISTS aiupred_protein_sequences_scores (
                    GeneID TEXT NOT NULL,
                    TranscriptID VARCHAR(100) PRIMARY KEY,
                    IUPRED2 TEXT NOT NULL
                )
            """)


def create_representative_cds_human_table(db_path: str, anchor: bool) -> None:
    with sqlite3.connect(db_path) as db:
        if anchor:
            db.execute("""
                CREATE TABLE IF NOT EXISTS Representative_CDSs (
                    CdsID INTEGER PRIMARY KEY AUTOINCREMENT,
                    GeneID VARCHAR(100) NOT NULL,
                    CdsStart INTEGER NOT NULL,
                    CdsEnd INTEGER NOT NULL,
                    CdsDnaSequence TEXT NOT NULL,
                    IUPRED2Scores TEXT NOT NULL,
                    ANCHOR2Scores TEXT NOT NULL,
                    PercentageOrderedRegion REAL NOT NULL,
                    PercentageOrderedAnchor REAL NOT NULL
                )
            """)
        else:
            db.execute("""
                CREATE TABLE IF NOT EXISTS Representative_CDSs (
                    CdsID INTEGER PRIMARY KEY AUTOINCREMENT,
                    GeneID VARCHAR(100) NOT NULL,
                    CdsStart INTEGER NOT NULL,
                    CdsEnd INTEGER NOT NULL,
                    CdsDnaSequence TEXT NOT NULL,
                    IUPRED2Scores TEXT NOT NULL,
                    PercentageOrderedRegion REAL NOT NULL
                )
            """)


def insert_aiupred_scores(db_path: str, tuples_list: list, anchor: bool) -> None:
    if anchor:
        query = """
            INSERT INTO aiupred_protein_sequences_scores(
                GeneID,
                TranscriptID,
                IUPRED2,
                ANCHOR2
            )
            VALUES (?, ?, ?, ?)
        """
    else:
        query = """
            INSERT INTO aiupred_protein_sequences_scores(
                GeneID,
                TranscriptID,
                IUPRED2
            )
            VALUES (?, ?, ?)
        """
    with sqlite3.connect(db_path) as db:
        db.executemany(query, tuples_list)


def insert_representative_cds(db_path: str, anchor: bool, tuples_list: list) -> None:
    with sqlite3.connect(db_path) as db:
        if anchor:
            db.executemany("""
                INSERT INTO Representative_CDSs (
                    GeneID,
                    CdsStart,
                    CdsEnd,
                    CdsDnaSequence,
                    IUPRED2Scores,
                    ANCHOR2Scores,
                    PercentageOrderedRegion,
                    PercentageOrderedAnchor
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, tuples_list)
        else:
            db.executemany("""
                INSERT INTO Representative_CDSs (
                    GeneID,
                    CdsStart,
                    CdsEnd,
                    CdsDnaSequence,
                    IUPRED2Scores,
                    PercentageOrderedRegion
                )
                VALUES (?, ?, ?, ?, ?, ?)
            """, tuples_list)


def update_mapping_IUPRED2_peptides_columns(
        db_path: str,
        anchor: bool,
        tuples_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        if anchor:
            db.executemany("""
                UPDATE Gene_CDS_Proteins
                SET MappingIUPRED2Peptides=?,
                    MappingANCHOR2Peptides=?,
                    PercentageOrderedRegion=?,
                    PercentageOrderedAnchor=?
                WHERE CdsID=?
            """, tuples_list)
        else:
            db.executemany("""
                UPDATE Gene_CDS_Proteins
                SET MappingIUPRED2Peptides=?,
                    PercentageOrderedRegion=?
                WHERE CdsID=?
            """, tuples_list)


def update_cdss_table(anchor: bool, db_path: str, cds_table_name: str) -> None:
    """Add AIUPred score columns to Gene_CDS_Proteins table."""
    add_column_to_table(db_path, cds_table_name, 'MappingIUPRED2Peptides', 'TEXT')
    add_column_to_table(db_path, cds_table_name, 'PercentageOrderedRegion', 'REAL')
    if anchor:
        add_column_to_table(db_path, cds_table_name, 'MappingANCHOR2Peptides', 'TEXT')
        add_column_to_table(db_path, cds_table_name, 'PercentageOrderedAnchor', 'REAL')


def query_processed_transcripts(db_path: str) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("SELECT DISTINCT TranscriptID FROM aiupred_protein_sequences_scores")
        return [row[0] for row in cursor.fetchall()]


def query_proteins_table(db_path: str) -> defaultdict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("SELECT GeneID, TranscriptID, ProteinSequence FROM Gene_Transcript_Proteins")
        result = cursor.fetchall()
    proteins_dictionary = defaultdict(list)
    for gene_id, transcript, peptide_sequence in result:
        proteins_dictionary[transcript] = [gene_id, peptide_sequence]
    return proteins_dictionary


def query_protein_scores_table(db_path: str, anchor: bool) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        if anchor:
            cursor.execute("""
                SELECT 
                    GeneID,
                    TranscriptID,
                    IUPRED2,
                    ANCHOR2
                FROM aiupred_protein_sequences_scores
            """)
        else:
            cursor.execute("""
                SELECT 
                    GeneID,
                    TranscriptID,
                    IUPRED2
                FROM aiupred_protein_sequences_scores
            """)
        result = cursor.fetchall()

    scores = {}
    if anchor:
        for gene_id, transcript, iupred2, anchor2 in result:
            if gene_id not in scores:
                scores[gene_id] = {}
            scores[gene_id][transcript] = [iupred2, anchor2]
    else:
        for gene_id, transcript, iupred2 in result:
            if gene_id not in scores:
                scores[gene_id] = {}
            scores[gene_id][transcript] = [iupred2, None]
    return scores

def query_gene_cdss(db_path: str) -> defaultdict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT
                DISTINCT
                GeneID,
                TranscriptID,
                PeptideCdsStart,
                PeptideCdsEnd,
                CdsID
            FROM Gene_CDS_Proteins
        """)
        result = cursor.fetchall()
    gene_cdss = defaultdict(lambda: defaultdict(list))
    for gene_id, transcript_id, pep_start, pep_end, cds_id in result:
        gene_cdss[gene_id][transcript_id].append([P.open(pep_start, pep_end), cds_id])
    return gene_cdss


def query_protein_scores(db_path: str, anchor: bool) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        if anchor:
            cursor.execute("""
                           SELECT
                               TranscriptID,
                               IUPRED2,
                               ANCHOR2
                           FROM aiupred_protein_sequences_scores
                           """)
        else:
            cursor.execute("""
                           SELECT
                               TranscriptID,
                               IUPRED2
                           FROM aiupred_protein_sequences_scores
                           """)
        return cursor.fetchall()


def query_raw_cds_scores(path_db: str) -> dict:
    with sqlite3.connect(path_db) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT
                GeneID,
                CdsFrame,
                CdsStart,
                CdsEnd,
                CdsDnaSequence,
                MappingIUPRED2Peptides,
                MappingANCHOR2Peptides,
                MAX(PercentageOrderedRegion) AS PercentageOrderedRegion,
                PercentageOrderedAnchor
            FROM Gene_CDS_Proteins
            GROUP BY GeneID, CdsStart, CdsEnd
        """)
        records = cursor.fetchall()
    cdss_dictionary = defaultdict(lambda: defaultdict(list))
    for (gene_id, _, cds_start, cds_end, cds_dna_seq,
         mapping_iupred2, mapping_anchor2,
         perc_ordered, perc_ordered_anchor) in records:
        cdss_dictionary[gene_id][P.open(cds_start, cds_end)] = [
            cds_dna_seq, mapping_iupred2, mapping_anchor2,
            perc_ordered, perc_ordered_anchor
        ]
    return cdss_dictionary