import sqlite3
import portion as P
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
import os


def fetch_and_dump_data(
        file_path: Path,
        query_function: callable,
        key_format: str,
        protein_db_path: Path
) -> None:
    if not os.path.exists(protein_db_path):
        raise FileNotFoundError(
            f'Protein database not found: {protein_db_path},'
            f' use prot_sql_utils.py to generate it'
        )
    if not os.path.exists(file_path):
        results = query_function(
            protein_db_path
        )
        sequence_dictionary = {
            key_format.format(*sequence): (
                sequence[-1][:-1]
                if sequence[-1][-1] == '*'  # remove stop codon
                else sequence[-1]
            )
            for sequence in results if sequence[-1]
            # remove stop codon
        }
        dump_sequence_dictionary_into_fasta_file(
            out_file_path=file_path,
            sequence_dictionary=sequence_dictionary
        )


def query_xgroup_counts_ecod(
        db_path: Path,
) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            CAST(substr(FID, 1, instr(FID, '.') - 1) AS INTEGER) AS XGroupID,
            COUNT(*) AS count
        FROM ECOD_human_raw_domains
        GROUP BY XGroupID
        """)
        return {x_group_id: count for x_group_id, count in cursor.fetchall()}


def dump_sequence_dictionary_into_fasta_file(
        out_file_path: Path,
        sequence_dictionary: dict
) -> None:
    with open(out_file_path, "w") as handle:
        for annotation_id, annotation_sequence in sequence_dictionary.items():
            record = SeqRecord(
                Seq(annotation_sequence),
                id=str(annotation_id),
                description=''
            )
            SeqIO.write(record, handle, "fasta")


def query_for_cds(
        db_path: Path
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT
         GeneID,
         TranscriptID,
         CdsID,
         Rank,
         CdsProteinSequence
        FROM Gene_CDS_Proteins
        """)
        return cursor.fetchall()


def query_x_group_ids_names(
        db_path: str,
):
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        query = """
       SELECT DISTINCT
       CAST(substr(FID, 1, instr(f_id, '.') - 1) AS INTEGER) AS XGroupID,
       XGroupName 
       FROM 
       ECOD_latest_domains
       """
        cursor.execute(query)
        return cursor.fetchall()


def query_ecod_full_domain_duplication_records(
    db_path: str
) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        query = """
        SELECT
            SUBSTR(GeneID, INSTR(GeneID, ':') + 1) AS gene_id,
            EventStart,
            EventEnd
        FROM
            Expansions_full
        WHERE FullDomainExon > 0;
        """
        records_dict = defaultdict(list)
        for record in cursor.execute(query).fetchall():
            gene_id, start, end = record
            records_dict[gene_id].append(P.open(start, end))
        return records_dict
