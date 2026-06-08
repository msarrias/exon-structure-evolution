import sqlite3
from collections import defaultdict
import portion as P
from pathlib import Path


def check_if_empty_table(
        db_path: Path,
        table_name: str,
) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"""
        SELECT COUNT(*) FROM {table_name}
        """)
        return cursor.fetchone()[0] == 0


def check_if_empty_column(
        db_path: Path,
        table_name: str,
        column_name: str,
) -> bool:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"""
        SELECT 
            COUNT(*) 
                FROM {table_name} 
                WHERE {column_name} IS NOT NULL
        """)
        return cursor.fetchone()[0] == 0


def check_if_table_exists(
        db_path: Path,
        table_name: str,
        timeout_db: int,
) -> bool:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            name
        FROM sqlite_master 
        WHERE type='table';
        """)
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
        target_db_path,
        source_db_path,
        table_name: str
):
    with sqlite3.connect(target_db_path) as target_conn, sqlite3.connect(source_db_path) as source_conn:
        target_cursor = target_conn.cursor()
        source_cursor = source_conn.cursor()
        source_cursor.execute(f"""
        SELECT 
            sql 
        FROM sqlite_master 
        WHERE type='table' AND name="{table_name}";
        """)
        create_table_sql = source_cursor.fetchone()[0]
        # Fetch all data from the source table
        source_cursor.execute(f"""SELECT * FROM {table_name}""")
        rows = source_cursor.fetchall()
        target_cursor.execute(create_table_sql)
        for row in rows:
            placeholders = ', '.join('?' * len(row))
            target_cursor.execute(
                f"""INSERT INTO {table_name} VALUES ({placeholders})""",
                row
            )


def create_mrna_transcripts_cds_tables(
        db_path: Path
) -> None:
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
            primary key (GeneID, TranscriptID)
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
            primary key (GeneID, TranscriptID, Rank, CdsID)
        )
        """)


def create_ecod_cds_mrna_hmm_results_tables(
        db_path: Path,
        timeout_db: int,
) -> None:
    if check_if_table_exists(
            db_path,
            'ECOD_hmm_CDSs',
            timeout_db
    ):
        with sqlite3.connect(db_path, timeout=timeout_db) as db:
            cursor = db.cursor()
            cursor.execute(
                """DROP TABLE ECOD_hmm_CDSs"""
            )
            db.commit()
    if check_if_table_exists(
            db_path,
            'ECOD_hmm_mRNA',
            timeout_db
    ):
        with sqlite3.connect(db_path, timeout=timeout_db) as db:
            cursor = db.cursor()
            cursor.execute(
                """DROP TABLE ECOD_hmm_mRNA"""
            )
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
                        );
        """)
        db.commit()


def create_x_group_filtered_ecod_hmm_cdss_table(
        db_path: Path,
):
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
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
            );
                """)
        db.commit()


def create_ecod_af2_human_domains_table(
        db_path: Path
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS ECOD_af2_human_domains (
            UID VARCHAR(100) PRIMARY KEY,
            ID VARCHAR(100) NOT NULL,
            TID VARCHAR(100) NOT NULL,
            DomainRange VARCHAR(100) NOT NULL,
            TGroupName VARCHAR(100) NOT NULL,	
            HGroupName VARCHAR(100) NOT NULL,
            XGroupName VARCHAR(100) NOT NULL,
            askDustin VARCHAR(100) NOT NULL
        )
        """)
        db.commit()


def create_ecod_latest_domains_table(
        db_path: Path,
        timeout_db: int,
) -> None:

    # Column 1: ECOD uid - internal domain unique identifier
    # Column 2: ECOD domain id - domain identifier
    # Column 3: ECOD representative status - manual (curated) or automated nonrep
    # Column 4: ECOD hierachy identifier - [X-group].[H-group].[T-group].[F-group]
    # Column 5: PDB identifier
    # Column 6: Chain identifier (note: case-sensitive)
    # Column 7: PDB residue number range
    # Column 8: seq_id number range (based on internal PDB indices)
    # Column 9: UniProt accession (if available, based on
    # available PDB/SIFTS mapping at time of release)
    # Column 10: Architecture name
    # Column 11: X-group name
    # Column 12: H-group name
    # Column 13: T-group name
    # Column 14: F-group name (F_UNCLASSIFIED denotes that
    #  domain has not been assigned to an F-group)
    # Column 15: Domain assembly status (if domain is member of assembly,
    #  partners' ecod domain ids listed)
    # Column 16: Comma-separated value list of non-polymer entities within
    #  4 A of at least one residue of domain
    # http://prodata.swmed.edu/ecod/distributions/README
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
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


def create_ecod_human_raw_domains_table(
        db_path: Path,
):
    # data comes from https://zenodo.org/records/6998803
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
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
            );
                """)
        db.commit()


def create_event_domain_subdomain_matches_mapping_table(
        db_path: Path,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Event_Domain_Subdomain_Matches_Mapping (
            MatchSource TEXT CHECK (MatchSource IN ('mRNA', 'CDS')),
            MatchCategory TEXT CHECK (MatchCategory IN ('FullDomainExon', 'SubDomainExon', 'WithinDomainBoundary')),
            EventID INTEGER NOT NULL,
            GeneID VARCHAR(100) NOT NULL,
            EventStart INTEGER NOT NULL,
            EventEnd INTEGER NOT NULL,
            PRIMARY KEY (GeneID, EventID, EventStart, EventEnd)
        )
        """)
        db.commit()


def create_x_group_filtered_ecod_hmm_mrna_table(
        db_path: Path,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
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


def create_filtered_ecod_hmm_mrna_records_table(
        db_path: Path,
        query_alignment_threshold: float,
        e_value_threshold: float,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS filtered_ECOD_hmm_mRNA AS
        WITH AmbiguousNames AS (
              SELECT 
                FGroupName
              FROM ECOD_latest_domains
              GROUP BY FGroupName
              HAVING COUNT(DISTINCT XGroupName) > 1
            )
        SELECT 
            DISTINCT
            hmm.HmmTranscriptID,
            hmm.GeneID,
            hmm.TranscriptID,
            hmm.QueryName,
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames)
                THEN "AMBIGUOUS" 
            ELSE ecod_dom.XGroupName
            END AS XGroupName, 
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames) THEN NULL
                ELSE ROUND(ecod_dom.SymmetryScore,3) 
            END AS SymmetryScore,
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames) THEN 0
                ELSE CAST(substr(ecod_dom.FID, 1, instr(ecod_dom.FID, '.') - 1) AS INTEGER)
            END AS XGroupID,
            hmm.HmmFrom,
            hmm.HmmTo,
            hmm.AliFrom,
            hmm.AliTo,
            hmm.SequenceEvalue
        FROM ECOD_hmm_mRNA AS hmm
        INNER JOIN ECOD_latest_domains AS ecod_dom 
            ON hmm.QueryName = ecod_dom.FGroupName
        WHERE 
            hmm.QueryAlignPercentage >= {query_alignment_threshold}
            AND hmm.SequenceEvalue <= {e_value_threshold}
            GROUP BY GeneID, XGroupName,  HmmFrom, HmmTo, AliFrom, AliTo;
        """)
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_hmm_mRNA_id ON filtered_ECOD_hmm_mRNA(HmmTranscriptID);
        """)
        db.commit()


def create_index_on_ecod_hmm_cdss(
        db_path: Path,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_gene_id ON ECOD_hmm_CDSs(GeneID);
                """)
        db.commit()


def create_index_on_ecod_latest_domains_table(
        db_path: Path,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS idx_f_group_name ON ECOD_latest_domains(FGroupName);
        """)
        db.commit()


def insert_into_event_domain_subdomain_matches_mapping_table(
        db_path: Path,
        tuples_list: list
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
                INSERT INTO Event_Domain_Subdomain_Matches_Mapping
                VALUES (?,?,?,?,?,?)
                """
        cursor.executemany(
            insert_gene_table_param,
            tuples_list
        )


def insert_in_x_group_filtered_ecod_hmm_mrna_table(
        db_path: Path,
        tuples_list: list
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO x_group_filtered_ECOD_hmm_mRNA
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """
        cursor.executemany(
            insert_gene_table_param,
            tuples_list
        )
        db.commit()


def insert_in_ecod_latest_domains_table(
        db_path: Path,
        timeout_db: int,
        ecod_hmm_args_tuple: list
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO ECOD_latest_domains
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_gene_table_param,
            ecod_hmm_args_tuple
        )


def insert_ecod_human_raw_domains_table(
        db_path: Path,
        combined_domains_tuple_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO ECOD_human_raw_domains(
        UniprotID,
        DpamDomainID,
        DomainRange,
        FID,
        UID,
        ConfidenceMeasure
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_gene_table_param,
            combined_domains_tuple_list
        )


def insert_ecod_af2_human_domains_table(
        db_path: Path,
        combined_domains_tuple_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO ECOD_af2_human_domains
        VALUES (?,?,?,?,?,?,?,?)
        """
        cursor.executemany(
            insert_gene_table_param,
            combined_domains_tuple_list
        )


def insert_in_x_group_filtered_ecod_hmm_cdss_table(
        db_path: Path,
        timeout_db: int,
        ecod_hmm_args_tuple: list
):
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_x_group_filtered_ECOD_hmm_CDSs = """  
        INSERT INTO x_group_filtered_ECOD_hmm_CDSs (
            HmmCdsID,
            GeneID,
            CdsStart,
            CdsEnd,
            QueryName,
            XGroupName,
            XGroupID,
            SymmetryScore,
            HmmFrom,
            HmmTo,
            AliFrom,
            AliTo,
            SequenceEvalue
            )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
    cursor.executemany(
        insert_x_group_filtered_ECOD_hmm_CDSs,
        ecod_hmm_args_tuple
    )
    db.commit()


def insert_in_ecod_hmm_cds_table(
        db_path: Path,
        timeout_db: int,
        ecod_hmm_args_tuple: list
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO ECOD_hmm_CDSs (
            GeneID,
            TranscriptID,
            CdsID,
            TranscriptLen,
            QueryName,
            QueryAccession,
            QueryLen,
            SequenceEvalue,
            SequenceScore,
            SequenceBias,
            DomainNum,
            DomainTotal,
            DomainCEvalue,
            DomainIEvalue,
            DomainScore,
            DomainBias,
            HmmFrom,
            HmmTo,
            AliFrom,
            AliTo,
            EnvFrom,
            EnvTo,
            Acc,
            Description,
            QueryAlignPercentage,
            TargetAlignPercentage
            )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_gene_table_param,
            ecod_hmm_args_tuple
        )
        db.commit()


def insert_in_ecod_hmm_mrna_table(
        db_path: Path,
        timeout_db: int,
        ecod_hmm_args_tuple: list
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_mRNA_table_param = """  
        INSERT INTO ECOD_hmm_mRNA (
            GeneID,
            TranscriptID,
            TranscriptLen, 
            QueryName,
            QueryAccession,
            QueryLen,
            SequenceEvalue,
            SequenceScore,
            SequenceBias,
            DomainNum,
            DomainTotal,
            DomainCEvalue,
            DomainIEvalue,
            DomainScore,
            DomainBias,
            HmmFrom,
            HmmTo,
            AliFrom,
            AliTo,
            EnvFrom,
            EnvTo,
            Acc,
            Description,
            QueryAlignPercentage,
            TargetAlignPercentage)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_mRNA_table_param,
            ecod_hmm_args_tuple
        )
        db.commit()


def insert_into_mrna_transcripts_table(
        db_path: Path,
        timeout_db: int,
        gene_args_tuple: list
) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO Gene_Transcript_Proteins (
            GeneID,
            GeneChrom,
            GeneStrand,
            GeneStart, 
            GeneEnd, 
            TranscriptID,
            TranscriptStart,
            TranscriptEnd,
            ProteinSequence) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_gene_table_param,
            gene_args_tuple
        )
        db.commit()


def insert_into_cdss_table(
        db_path: Path,
        gene_args_tuple: list
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO Gene_CDS_Proteins (
            GeneID,
            TranscriptID,
            CdsID,
            Rank,
            CdsFrame,
            CdsStart,
            CdsEnd, 
            TranscriptCdsStart,
            TranscriptCdsEnd,
            PeptideCdsStart,
            PeptideCdsEnd,
            CdsDnaSequence, 
            CdsProteinSequence
            ) 
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(
            insert_gene_table_param,
            gene_args_tuple
        )


def query_for_mrna(
        db_path: Path,
        timeout_db: int
) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            GeneID, 
            TranscriptID, 
            ProteinSequence
        FROM Gene_Transcript_Proteins
        """)
        return cursor.fetchall()


def query_for_cds(
        db_path: Path,
        timeout_db: int
) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


def query_exonize_gene_events(
        set_genes: set,
        db_path: Path
) -> list:
    records_list = list()
    for gene in set_genes:
        with sqlite3.connect(db_path) as db:
            cursor = db.cursor()
            query = f"""
            SELECT 
                * 
            FROM Expansions_full 
            WHERE GeneID="gene:{gene}"
            ORDER BY GeneID, ExpansionID
            """
            cursor.execute(query)
            records_list.extend(cursor.fetchall())
    return records_list


def query_filtered_ecod_hmm_cds_hits(
        db_path: Path,
        db_path2: str,
        query_coverage_threshold: float,
        target_coverage_threshold: float,
        evalue_threshold: float,
) -> dict:
    with sqlite3.connect(db_path) as db:
        db.execute(f"ATTACH DATABASE '{db_path2}' AS db2")
        cursor = db.cursor()
        cursor.execute(f"""
        WITH AmbiguousNames AS (
        SELECT 
            FGroupName
        FROM ECOD_latest_domains
        GROUP BY FGroupName
        HAVING COUNT(DISTINCT XGroupName) > 1
        )
        SELECT 
            DISTINCT
            hmm.HmmCdsID,
            hmm.GeneID, 
            cds.CdsStart - db2.mrna.GeneStart as CdsStart,
            cds.CdsEnd - db2.mrna.GeneStart as CdsEnd,
            hmm.QueryName,
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames)
                    THEN "AMBIGUOUS" 
                    ELSE ecod_dom.XGroupName
            END AS XGroupName,
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames) 
                    THEN 0
                    ELSE CAST(substr(ecod_dom.FID, 1, instr(ecod_dom.FID, '.') - 1) AS INTEGER)
            END AS XGroupID, 
            CASE 
                WHEN hmm.QueryName IN (SELECT FGroupName FROM AmbiguousNames) 
                    THEN NULL
                    ELSE ROUND(ecod_dom.SymmetryScore,3) 
            END AS SymmetryScore,
            hmm.HmmFrom, 
            hmm.HmmTo, 
            hmm.AliFrom, 
            hmm.AliTo, 
            MIN(hmm.SequenceEvalue) AS SequenceEvalue
        FROM ECOD_hmm_CDSs AS hmm
        INNER JOIN db2.Gene_CDS_Proteins AS cds ON hmm.CdsID = cds.CdsID
        INNER JOIN db2.Gene_Transcript_Proteins AS mrna ON hmm.GeneID = mrna.GeneID  
        INNER JOIN ECOD_latest_domains AS ecod_dom ON hmm.QueryName = ecod_dom.FGroupName
        WHERE
            QueryAlignPercentage >= {query_coverage_threshold}
            AND TargetAlignPercentage >= {target_coverage_threshold}
            AND SequenceEvalue <= {evalue_threshold}
        GROUP BY hmm.GeneID, CdsStart, CdsEnd, XGroupName,  hmm.HmmFrom, hmm.HmmTo, hmm.AliFrom, hmm.AliTo
        """)
        records = cursor.fetchall()
        genes_hits_dictionary = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for record in records:
            _, gene_id, cds_start, cds_end, query_name, x_group_name, *_ = record
            genes_hits_dictionary[gene_id][P.open(cds_start, cds_end)][x_group_name].append(record)
        return genes_hits_dictionary


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


def query_uids_for_given_f_id(
        db_path: Path,
        f_id: int,
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute(f"""
        SELECT
            UID,
            CAST(substr(FID, 1, instr(FID, '.') - 1) AS INTEGER) AS first_digit_before_dot
        FROM ECOD_latest_domains
        WHERE first_digit_before_dot={f_id}
        """)
        return [i[0] for i in cursor.fetchall()]


def query_filtered_ecod_hmm_mrna_with_mapping_to_gene_intervals(
        db_path: Path,
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            HmmTranscriptID,
            GeneID,
            QueryName,
            XGroupName,
            XGroupID,
            HmmFrom,
            HmmTo,
            AliFrom,
            AliTo,
            MIN(SequenceEvalue) as SequenceEvalue,
            SymmetryScore,
            MappingToGeneIntervals
        FROM filtered_ECOD_hmm_mRNA
        -- we group by x_group_name because "NO_X_NAME"
        -- can be associated to more than one x_group_id
        GROUP BY GeneID, XGroupName, MappingToGeneIntervals
        """)
        return cursor.fetchall()


def query_filtered_ecod_hmm_mrna_records(
        db_path: Path,
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            HmmTranscriptID,
            GeneID,
            TranscriptID,
            QueryName,
            AliFrom,
            AliTo
        FROM filtered_ECOD_hmm_mRNA;
        """)
        return cursor.fetchall()


def update_ecod_latest_domains_symmetry_score(
        db_path: Path,
        column_value_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.executemany(
            """
            UPDATE ECOD_latest_domains
            SET SymmetryScore=? 
            WHERE UID=?""", column_value_list)
        db.commit()


def update_mapping_to_gene_intervals(
        db_path: Path,
        column_value_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.executemany(
            """
            UPDATE filtered_ECOD_hmm_mRNA
            SET MappingToGeneIntervals=? 
            WHERE HmmTranscriptID=?""", column_value_list)
        db.commit()


def update_remove_quotes_from_string(
        db_path: Path,
        table_name: str,
        column_names_list: list,
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        for column_name in column_names_list:
            cursor.execute(f"""
            UPDATE {table_name}
            SET {column_name} = REPLACE({column_name}, '"', '')
            """)
        db.commit()


def update_x_group_filtered_ecod_hmm_mrna_match_nature_columns(
    db_path: Path,
    tuples_to_insert_list: list
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.executemany(
            "UPDATE x_group_filtered_ECOD_hmm_mRNA "
            "SET WithinDomainBoundary=?, FullDomainExon=?, SubDomainExon=? "
            "WHERE HmmTranscriptID=?", tuples_to_insert_list)


def update_events_match_nature_columns(
    db_path: Path,
    tuples_to_insert_list: list
) -> None:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.executemany(
            """
            UPDATE Expansions_full
            SET WithinDomainBoundary=?, FullDomainExon=?, SubDomainExon=?
            WHERE GeneID=? AND EventStart=? AND EventEnd=?
            """,
            tuples_to_insert_list)


def fetch_distinct_genes(
    cursor: sqlite3.Cursor,
    table_name: str,
    condition: str = None
) -> set:
    query = f"SELECT DISTINCT GeneID from {table_name}"
    if condition:
        query = f"SELECT DISTINCT GeneID from {table_name} WHERE {condition}"
    cursor.execute(query)
    return {record[0] for record in cursor.fetchall()}


def query_hmm_searches_matches(
        db_path: Path,
):
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        hmm_mrna_domain_matches = """
        SELECT
            GeneID,
            HmmTranscriptID,
            XGroupID,
            MappingToGeneIntervals,
            SymmetryScore
        FROM 
            x_group_filtered_ECOD_hmm_mRNA
        """
        cursor.execute(hmm_mrna_domain_matches)
        mrna_domain_matches_list = cursor.fetchall()
        hmm_cds_domain_matches = """
        WITH x_groups_ids_symmetry AS (
        SELECT DISTINCT
            CAST(substr(FID, 1, instr(FID, '.') - 1) AS INTEGER) AS XGroupID,
            SymmetryScore 
       FROM 
       ECOD_latest_domains
       )
        SELECT
            cds.GeneID,
            cds.HmmCdsID,
            cds.XGroupID,
            cds.CdsStart, 
            cds.AliFrom,
            cds.AliTo,
            dom.SymmetryScore

        FROM 
            x_group_filtered_ECOD_hmm_CDSs as cds
            INNER JOIN x_groups_ids_symmetry as dom
            ON cds.XGroupID = dom.XGroupID
        """
        cursor.execute(hmm_cds_domain_matches)
        cds_domain_matches_list = cursor.fetchall()
    return mrna_domain_matches_list, cds_domain_matches_list


def query_genes_dictionary(
        db_path: Path,
) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            GeneID,
            GeneStart,
            GeneEnd,
            GeneChrom,
            GeneStrand 
        FROM Genes
        """)
        genes_records = cursor.fetchall()
    genes_dict = {}
    for gene_record in genes_records:
        gene_id, start, end, chromosome, strand = gene_record
        genes_dict[gene_id] = [P.open(start, end), chromosome, strand]
    return genes_dict


def query_events_dictionary(
        db_path: Path,
) -> dict:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            e.GeneID,
            e.EventStart - gene.GeneStart,
            e.EventEnd - gene.GeneStart,
            e.ExpansionID,
            e.Mode
        FROM Expansions_full as e
        INNER JOIN Genes as gene
        ON e.GeneID = gene.GeneID
        """)
        events = cursor.fetchall()
    events_dict = defaultdict(lambda: defaultdict(list))
    for idx, event in enumerate(events):
        gene_id, ref_start, ref_end, expansion_id, event_type = event
        events_dict[gene_id][expansion_id].append(
            [P.open(ref_start, ref_end), event_type]
        )
    return events_dict


def query_ecod_mrna_events(
        db_path: Path
) -> list:
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        query = """
        SELECT
            HmmTranscriptID,
            GeneID,
            QueryName,
            XGroupName,
            XGroupID,
            AliFrom,
            AliTo,
            SymmetryScore,
            MappingToGeneIntervals
        FROM x_group_filtered_ECOD_hmm_mRNA
        """
        cursor.execute(query)
        return cursor.fetchall()


def query_ecod_subdomain_domain_counts(
    db_path: Path,
    hmm_mrna_table_name: str = "x_group_filtered_ECOD_hmm_mRNA",
    hmm_cds_table_name: str = "x_group_filtered_ECOD_hmm_CDSs",
    expansions_table_name: str = "Expansions_full"
) -> tuple[set, set, set, set, set, set, set]:
    with sqlite3.connect(db_path) as db:
        # Fetch genes with full domain cds search HMM
        cursor = db.cursor()
        genes_with_full_domain_cds_hmm = fetch_distinct_genes(
            cursor=cursor,
            table_name=hmm_cds_table_name,
            condition=None
        )
        # Fetch genes with full domain mrna search HMM
        genes_with_full_domain_mrna_hmm = fetch_distinct_genes(
            cursor=cursor,
            table_name=hmm_mrna_table_name,
            condition="FullDomainExon > 0"
        )
        # Fetch genes with full domain mrna exonize
        genes_with_exon_duplications = fetch_distinct_genes(
            cursor=cursor,
            table_name=expansions_table_name,
            condition=None
        )
        # Fetch genes with full domain mrna exonize
        genes_with_full_domain_exon_duplication_hmm_search = fetch_distinct_genes(
            cursor=cursor,
            table_name=expansions_table_name,
            condition="FullDomainExon > 0"
        )
        # Fetch genes with subdomain mrna search HMM
        genes_with_sub_domain_hmm_search = fetch_distinct_genes(
            cursor=cursor,
            table_name=hmm_mrna_table_name,
            condition="SubDomainExon > 0"
        )
        # Fetch genes with subdomain mrna search HMM
        genes_with_sub_domain_exon_duplication_hmm_search = fetch_distinct_genes(
            cursor=cursor,
            table_name=expansions_table_name,
            condition="SubDomainExon > 0"
        )
        # Fetch genes within domain boundary mrna search HMM
        genes_within_domain_boundary_hmm = fetch_distinct_genes(
            cursor=cursor,
            table_name=hmm_mrna_table_name,
            condition="WithinDomainBoundary > 0"
        )
    return (set([gene_id.rsplit(":")[1] if ":" in gene_id else gene_id
                 for gene_id in genes_with_exon_duplications]),
            set([gene_id.rsplit(":")[1] if ":" in gene_id else gene_id
                 for gene_id in genes_with_full_domain_exon_duplication_hmm_search]),
            set([gene_id.rsplit(":")[1] if ":" in gene_id else gene_id
                 for gene_id in genes_with_sub_domain_exon_duplication_hmm_search]),
            genes_with_full_domain_mrna_hmm,
            genes_with_full_domain_cds_hmm,
            genes_with_sub_domain_hmm_search,
            genes_within_domain_boundary_hmm
            )
