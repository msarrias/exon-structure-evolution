import sqlite3
from pathlib import Path


def create_mrna_transcripts_cds_tables(db_path: Path) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS Gene_Transcript_Proteins (
                GeneID VARCHAR(100) NOT NULL,
                GeneChrom VARCHAR(100) NOT NULL,
                GeneStrand VARCHAR(1) NOT NULL,
                GeneStart INTEGER NOT NULL,
                GeneEnd INTEGER NOT NULL,
                TranscriptID VARCHAR(100) NOT NULL,
                TranscriptStart INTEGER NOT NULL,
                TranscriptEnd INTEGER NOT NULL,
                ProteinSequence VARCHAR NOT NULL,
                PRIMARY KEY (GeneID, TranscriptID)
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS Gene_CDS_Proteins (
                GeneID VARCHAR(100) NOT NULL,
                TranscriptID VARCHAR(100) NOT NULL,
                CdsID VARCHAR(100) NOT NULL,
                Rank INTEGER NOT NULL,
                CdsFrame INTEGER NOT NULL,
                CdsStart INTEGER NOT NULL,
                CdsEnd INTEGER NOT NULL,
                TranscriptCdsStart INTEGER NOT NULL,
                TranscriptCdsEnd INTEGER NOT NULL,
                PeptideCdsStart INTEGER NOT NULL,
                PeptideCdsEnd INTEGER NOT NULL,
                CdsDnaSequence VARCHAR NOT NULL,
                CdsProteinSequence VARCHAR NOT NULL,
                PRIMARY KEY (GeneID, TranscriptID, Rank, CdsID)
            )
        """)
        db.commit()


def insert_into_mrna_transcripts_table(
        db_path: Path,
        timeout_db: int,
        gene_args_tuple: list,
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.executemany("""
            INSERT INTO Gene_Transcript_Proteins (
                GeneID, GeneChrom, GeneStrand, GeneStart, GeneEnd,
                TranscriptID, TranscriptStart, TranscriptEnd, ProteinSequence
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, gene_args_tuple)
        db.commit()


def insert_into_cdss_table(
        db_path: Path,
        gene_args_tuple: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.executemany("""
            INSERT INTO Gene_CDS_Proteins (
                GeneID, TranscriptID, CdsID, Rank, CdsFrame,
                CdsStart, CdsEnd, TranscriptCdsStart, TranscriptCdsEnd,
                PeptideCdsStart, PeptideCdsEnd, CdsDnaSequence, CdsProteinSequence
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, gene_args_tuple)
        db.commit()


def query_for_mrna(db_path: Path, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT GeneID, TranscriptID, ProteinSequence
            FROM Gene_Transcript_Proteins
        """)
        return cursor.fetchall()


def query_for_cds(db_path: Path, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT GeneID, TranscriptID, CdsID, Rank, CdsProteinSequence
            FROM Gene_CDS_Proteins
        """)
        return cursor.fetchall()