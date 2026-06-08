import sqlite3
from collections import defaultdict
from pathlib import Path

import portion as P

from lib.shared_sql_utils import check_if_table_exists, fetch_distinct_genes


# Table creation

def create_ecod_cds_mrna_hmm_results_tables(
        db_path: Path,
        timeout_db: int,
) -> None:
    for table in ("ECOD_hmm_CDSs", "ECOD_hmm_mRNA"):
        if check_if_table_exists(db_path, table, timeout_db):
            with sqlite3.connect(db_path, timeout=timeout_db) as db:
                db.execute(f"DROP TABLE {table}")
                db.commit()

    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS ECOD_hmm_CDSs (
                HmmCdsID INTEGER PRIMARY KEY AUTOINCREMENT,
                GeneID VARCHAR(100) NOT NULL,
                TranscriptID VARCHAR(100) NOT NULL,
                CdsID VARCHAR(100) NOT NULL,
                TranscriptLen INTEGER NOT NULL,
                QueryName VARCHAR(100) NOT NULL,
                QueryAccession VARCHAR(100) NOT NULL,
                QueryLen INTEGER NOT NULL,
                SequenceEvalue FLOAT NOT NULL,
                SequenceScore FLOAT NOT NULL,
                SequenceBias FLOAT NOT NULL,
                DomainNum INTEGER NOT NULL,
                DomainTotal INTEGER NOT NULL,
                DomainCEvalue FLOAT NOT NULL,
                DomainIEvalue FLOAT NOT NULL,
                DomainScore FLOAT NOT NULL,
                DomainBias FLOAT NOT NULL,
                HmmFrom INTEGER NOT NULL,
                HmmTo INTEGER NOT NULL,
                AliFrom INTEGER NOT NULL,
                AliTo INTEGER NOT NULL,
                EnvFrom INTEGER NOT NULL,
                EnvTo INTEGER NOT NULL,
                Acc FLOAT NOT NULL,
                Description VARCHAR NOT NULL,
                QueryAlignPercentage FLOAT NOT NULL,
                TargetAlignPercentage FLOAT NOT NULL
            )
        """)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS ECOD_hmm_mRNA (
                HmmTranscriptID INTEGER PRIMARY KEY AUTOINCREMENT,
                GeneID VARCHAR(100) NOT NULL,
                TranscriptID VARCHAR(100) NOT NULL,
                TranscriptLen INTEGER NOT NULL,
                QueryName VARCHAR(100) NOT NULL,
                QueryAccession VARCHAR(100) NOT NULL,
                QueryLen INTEGER NOT NULL,
                SequenceEvalue FLOAT NOT NULL,
                SequenceScore FLOAT NOT NULL,
                SequenceBias FLOAT NOT NULL,
                DomainNum INTEGER NOT NULL,
                DomainTotal INTEGER NOT NULL,
                DomainCEvalue FLOAT NOT NULL,
                DomainIEvalue FLOAT NOT NULL,
                DomainScore FLOAT NOT NULL,
                DomainBias FLOAT NOT NULL,
                HmmFrom INTEGER NOT NULL,
                HmmTo INTEGER NOT NULL,
                AliFrom INTEGER NOT NULL,
                AliTo INTEGER NOT NULL,
                EnvFrom INTEGER NOT NULL,
                EnvTo INTEGER NOT NULL,
                Acc FLOAT NOT NULL,
                Description VARCHAR NOT NULL,
                QueryAlignPercentage FLOAT NOT NULL,
                TargetAlignPercentage FLOAT NOT NULL
            )
        """)
        db.commit()


def create_x_group_filtered_ecod_hmm_cdss_table(db_path: Path) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute("""
            CREATE TABLE IF NOT EXISTS x_group_filtered_ECOD_hmm_CDSs (
                HmmCdsID INTEGER PRIMARY KEY,
                GeneID VARCHAR(100) NOT NULL,
                CdsStart INTEGER NOT NULL,
                CDSEnd INTEGER NOT NULL,
                QueryName VARCHAR(100) NOT NULL,
                XGroupName VARCHAR(100) NOT NULL,
                XGroupID INTEGER NOT NULL,
                SymmetryScore FLOAT,
                HmmFrom INTEGER NOT NULL,
                HmmTo INTEGER NOT NULL,
                AliFrom INTEGER NOT NULL,
                AliTo INTEGER NOT NULL,
                SequenceEvalue FLOAT NOT NULL,
                FOREIGN KEY (HmmCdsID) REFERENCES ECOD_hmm_CDSs(HmmCdsID)
            )
        """)
        db.commit()


def create_x_group_filtered_ecod_hmm_mrna_table(db_path: Path) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute("""
            CREATE TABLE IF NOT EXISTS x_group_filtered_ECOD_hmm_mRNA (
                HmmTranscriptID INTEGER PRIMARY KEY,
                GeneID VARCHAR(100) NOT NULL,
                QueryName VARCHAR(100) NOT NULL,
                XGroupName VARCHAR(100) NOT NULL,
                XGroupID INTEGER NOT NULL,
                HmmFrom INTEGER NOT NULL,
                HmmTo INTEGER NOT NULL,
                AliFrom INTEGER NOT NULL,
                AliTo INTEGER NOT NULL,
                SequenceEvalue FLOAT NOT NULL,
                SymmetryScore FLOAT,
                MappingToGeneIntervals VARCHAR NOT NULL,
                FOREIGN KEY (HmmTranscriptID) REFERENCES ECOD_hmm_mRNA(HmmTranscriptID)
            )
        """)
        db.commit()


def create_event_domain_subdomain_matches_mapping_table(db_path: Path) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute("""
            CREATE TABLE IF NOT EXISTS Event_Domain_Subdomain_Matches_Mapping (
                MatchSource TEXT CHECK (MatchSource IN ('mRNA', 'CDS')),
                MatchCategory TEXT CHECK (MatchCategory IN (
                    'FullDomainExon', 'SubDomainExon', 'WithinDomainBoundary')),
                EventID INTEGER NOT NULL,
                GeneID VARCHAR(100) NOT NULL,
                EventStart INTEGER NOT NULL,
                EventEnd INTEGER NOT NULL,
                PRIMARY KEY (GeneID, EventID, EventStart, EventEnd)
            )
        """)
        db.commit()


def create_ecod_latest_domains_table(db_path: Path, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        db.execute("""
            CREATE TABLE IF NOT EXISTS ECOD_latest_domains (
                UID VARCHAR(100) PRIMARY KEY,
                EcodDomainID VARCHAR(100) NOT NULL,
                ManualRep VARCHAR(100) NOT NULL,
                FID VARCHAR(100) NOT NULL,
                PDB VARCHAR(100) NOT NULL,
                Chain VARCHAR(100) NOT NULL,
                PDBRange VARCHAR(100) NOT NULL,
                SeqIDRange VARCHAR(100) NOT NULL,
                UnpAcc VARCHAR(100) NOT NULL,
                ArchitectureName VARCHAR(100) NOT NULL,
                XGroupName VARCHAR(100) NOT NULL,
                HGroupName VARCHAR(100) NOT NULL,
                TGroupName VARCHAR(100) NOT NULL,
                FGroupName VARCHAR(100) NOT NULL,
                AssemblyStatus VARCHAR(100) NOT NULL,
                Ligand VARCHAR(100) NOT NULL
            )
        """)
        db.commit()


def create_ecod_human_raw_domains_table(db_path: Path) -> None:
    with sqlite3.connect(db_path) as db:
        db.execute("""
            CREATE TABLE IF NOT EXISTS ECOD_human_raw_domains (
                ID INTEGER PRIMARY KEY,
                UniprotID VARCHAR(100) NOT NULL,
                DpamDomainID VARCHAR(100) NOT NULL,
                DomainRange VARCHAR(100) NOT NULL,
                FID VARCHAR(100) NOT NULL,
                UID VARCHAR(100) NOT NULL,
                ConfidenceMeasure INTEGER NOT NULL,
                FOREIGN KEY (FID) REFERENCES ECOD_latest_domains(FID),
                FOREIGN KEY (UID) REFERENCES ECOD_latest_domains(UID)
            )
        """)
        db.commit()


# Inserts

def insert_into_ecod_hmm_cds_table(
        db_path: Path,
        timeout_db: int,
        matches_list: list,
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        db.executemany("""
            INSERT INTO ECOD_hmm_CDSs (
                GeneID, TranscriptID, CdsID, TranscriptLen,
                QueryName, QueryAccession, QueryLen,
                SequenceEvalue, SequenceScore, SequenceBias,
                DomainNum, DomainTotal, DomainCEvalue, DomainIEvalue,
                DomainScore, DomainBias, HmmFrom, HmmTo,
                AliFrom, AliTo, EnvFrom, EnvTo,
                Acc, Description, QueryAlignPercentage, TargetAlignPercentage
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, matches_list)
        db.commit()


def insert_into_ecod_hmm_mrna_table(
        db_path: Path,
        timeout_db: int,
        matches_list: list,
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        db.executemany("""
            INSERT INTO ECOD_hmm_mRNA (
                GeneID, TranscriptID, TranscriptLen,
                QueryName, QueryAccession, QueryLen,
                SequenceEvalue, SequenceScore, SequenceBias,
                DomainNum, DomainTotal, DomainCEvalue, DomainIEvalue,
                DomainScore, DomainBias, HmmFrom, HmmTo,
                AliFrom, AliTo, EnvFrom, EnvTo,
                Acc, Description, QueryAlignPercentage, TargetAlignPercentage
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, matches_list)
        db.commit()


def insert_into_x_group_filtered_ecod_hmm_cdss_table(
        db_path: Path,
        records: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            INSERT INTO x_group_filtered_ECOD_hmm_CDSs (
                HmmCdsID, GeneID, CdsStart, CdsEnd, QueryName,
                XGroupName, XGroupID, SymmetryScore,
                HmmFrom, HmmTo, AliFrom, AliTo, SequenceEvalue
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, records)
        db.commit()


def insert_into_x_group_filtered_ecod_hmm_mrna_table(
        db_path: Path,
        records: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            INSERT INTO x_group_filtered_ECOD_hmm_mRNA
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, records)
        db.commit()


def insert_into_event_domain_subdomain_matches_mapping_table(
        db_path: Path,
        tuples_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany(
            "INSERT INTO Event_Domain_Subdomain_Matches_Mapping VALUES (?,?,?,?,?,?)",
            tuples_list
        )
        db.commit()


def insert_in_ecod_latest_domains_table(
        db_path: Path,
        timeout_db: int,
        ecod_hmm_args_tuple: list,
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        db.executemany(
            "INSERT INTO ECOD_latest_domains VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ecod_hmm_args_tuple
        )
        db.commit()


def insert_ecod_human_raw_domains_table(
        db_path: Path,
        combined_domains_tuple_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            INSERT INTO ECOD_human_raw_domains (
                UniprotID, DpamDomainID, DomainRange,
                FID, UID, ConfidenceMeasure
            ) VALUES (?, ?, ?, ?, ?, ?)
        """, combined_domains_tuple_list)
        db.commit()


# Updates

def update_mapping_to_gene_intervals(
        db_path: Path,
        column_value_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            UPDATE x_group_filtered_ECOD_hmm_mRNA
            SET MappingToGeneIntervals=?
            WHERE HmmTranscriptID=?
        """, column_value_list)
        db.commit()


def update_x_group_filtered_ecod_hmm_mrna_match_nature_columns(
        db_path: Path,
        tuples_to_insert_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            UPDATE x_group_filtered_ECOD_hmm_mRNA
            SET WithinDomainBoundary=?, FullDomainExon=?, SubDomainExon=?
            WHERE HmmTranscriptID=?
        """, tuples_to_insert_list)
        db.commit()


def update_events_match_nature_columns(
        db_path: Path,
        tuples_to_insert_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany("""
            UPDATE Expansions_full
            SET WithinDomainBoundary=?, FullDomainExon=?, SubDomainExon=?
            WHERE GeneID=? AND EventStart=? AND EventEnd=?
        """, tuples_to_insert_list)
        db.commit()


def update_ecod_latest_domains_symmetry_score(
        db_path: Path,
        column_value_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        db.executemany(
            "UPDATE ECOD_latest_domains SET SymmetryScore=? WHERE UID=?",
            column_value_list
        )
        db.commit()


# Queries

def query_hmm_searches_matches(db_path: Path) -> tuple:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT GeneID, HmmTranscriptID, XGroupID,
                   MappingToGeneIntervals, SymmetryScore
            FROM x_group_filtered_ECOD_hmm_mRNA
        """)
        mrna_matches = cursor.fetchall()
        cursor.execute("""
            WITH x_groups_ids_symmetry AS (
                SELECT DISTINCT
                    CAST(substr(FID, 1, instr(FID, '.') - 1) AS INTEGER) AS XGroupID,
                    SymmetryScore
                FROM ECOD_latest_domains
            )
            SELECT cds.GeneID, cds.HmmCdsID, cds.XGroupID,
                   cds.CdsStart, cds.AliFrom, cds.AliTo, dom.SymmetryScore
            FROM x_group_filtered_ECOD_hmm_CDSs AS cds
            INNER JOIN x_groups_ids_symmetry AS dom ON cds.XGroupID = dom.XGroupID
        """)
        cds_matches = cursor.fetchall()
    return mrna_matches, cds_matches


def query_genes_dictionary(db_path: Path) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT GeneID, GeneStart, GeneEnd, GeneChrom, GeneStrand
            FROM Genes
        """)
        records = cursor.fetchall()
    return {
        gene_id: [P.open(start, end), chrom, strand]
        for gene_id, start, end, chrom, strand in records
    }


def query_events_dictionary(db_path: Path) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT e.GeneID,
                   e.EventStart - gene.GeneStart,
                   e.EventEnd   - gene.GeneStart,
                   e.ExpansionID,
                   e.Mode
            FROM Expansions_full AS e
            INNER JOIN Genes AS gene ON e.GeneID = gene.GeneID
        """)
        events = cursor.fetchall()
    events_dict = defaultdict(lambda: defaultdict(list))
    for gene_id, ref_start, ref_end, expansion_id, event_type in events:
        events_dict[gene_id][expansion_id].append(
            [P.open(ref_start, ref_end), event_type]
        )
    return events_dict


def query_ecod_mrna_events(db_path: Path) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT HmmTranscriptID, GeneID, QueryName, XGroupName, XGroupID,
                   AliFrom, AliTo, SymmetryScore, MappingToGeneIntervals
            FROM x_group_filtered_ECOD_hmm_mRNA
        """)
        return cursor.fetchall()


def query_xgroup_counts_ecod(db_path: Path) -> dict:
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


def query_filtered_ecod_hmm_mrna_records(db_path: Path) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT HmmTranscriptID, GeneID, TranscriptID, QueryName, AliFrom, AliTo
            FROM x_group_filtered_ECOD_hmm_mRNA
        """)
        return cursor.fetchall()


def query_ecod_subdomain_domain_counts(
        db_path: Path,
        hmm_mrna_table_name: str = "x_group_filtered_ECOD_hmm_mRNA",
        hmm_cds_table_name: str = "x_group_filtered_ECOD_hmm_CDSs",
        expansions_table_name: str = "Expansions_full",
) -> tuple:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        genes_with_full_domain_cds_hmm = fetch_distinct_genes(cursor, hmm_cds_table_name)
        genes_with_full_domain_mrna_hmm = fetch_distinct_genes(cursor, hmm_mrna_table_name, "FullDomainExon > 0")
        genes_with_exon_duplications = fetch_distinct_genes(cursor, expansions_table_name)
        genes_with_full_domain_dup = fetch_distinct_genes(cursor, expansions_table_name, "FullDomainExon > 0")
        genes_with_sub_domain_hmm = fetch_distinct_genes(cursor, hmm_mrna_table_name, "SubDomainExon > 0")
        genes_with_sub_domain_dup = fetch_distinct_genes(cursor, expansions_table_name, "SubDomainExon > 0")
        genes_within_domain_boundary_hmm = fetch_distinct_genes(cursor, hmm_mrna_table_name, "WithinDomainBoundary > 0")

    def strip_prefix(gene_set):
        return {g.rsplit(":")[1] if ":" in g else g for g in gene_set}

    return (
        strip_prefix(genes_with_exon_duplications),
        strip_prefix(genes_with_full_domain_dup),
        strip_prefix(genes_with_sub_domain_dup),
        genes_with_full_domain_mrna_hmm,
        genes_with_full_domain_cds_hmm,
        genes_with_sub_domain_hmm,
        genes_within_domain_boundary_hmm,
    )
