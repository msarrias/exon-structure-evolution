import logging
import sqlite3
import csv
import gzip
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("build_alphafold_db.log"),
    ]
)
logger = logging.getLogger(__name__)


def extract_alphafold_metrics_and_sequence(pdb_gz_path):
    parser = PDBParser(QUIET=True)
    with gzip.open(pdb_gz_path, "rt") as file:
        structure = parser.get_structure("model", file)
    plddt_scores = []
    chains = set()
    sequences = {}
    for model in structure:
        for chain in model:
            temp = {'plddt': '', 'aa': ''}
            chains.add(chain.id)
            sequence = []
            residue_plddt_by_position = []
            for residue in chain:
                residue_plddt = sum(atom.bfactor for atom in residue) / len(residue)
                sequence.append(seq1(residue.resname))
                plddt_scores.append(residue_plddt)
                residue_plddt_by_position.append(str(round(residue_plddt, 3)))
            temp['aa'] = "".join(sequence)
            temp['plddt'] = "_".join(residue_plddt_by_position)
            sequences[chain.id] = temp

    return {
        "num_chains": len(chains),
        "average_plddt": round(np.mean(plddt_scores), 3) if plddt_scores else None,
        "std_plddt": round(np.std(plddt_scores), 3) if plddt_scores else None,
        "sequences": sequences
    }


def main(
        alphafold_db_path: Path,
        alphafold_uniprot_des_path: Path,
        mart_ensembl_canonical_mapping: Path,
        database_path: Path
):
    logger.info("Starting AlphaFold database build")
    logger.info(f"AlphaFold structures: {alphafold_db_path}")
    logger.info(f"UniProt description:  {alphafold_uniprot_des_path}")
    logger.info(f"BioMart mapping:      {mart_ensembl_canonical_mapping}")
    logger.info(f"Output database:      {database_path}")

    logger.info("Scanning AlphaFold structure files...")
    proteome_pdb_tuples = {
        tuple(p.stem.rsplit('-model')[0].rsplit('AF-')[1].split('-'))
        for p in alphafold_db_path.iterdir()
        if p.is_file() and p.name.endswith('.pdb.gz')
    }
    pdb_structure_ids = {i for i, _ in proteome_pdb_tuples}
    proteome_pdb_paths = [
        p
        for p in alphafold_db_path.iterdir()
        if p.is_file() and p.name.endswith('.pdb.gz')
    ]
    logger.info(f"Found {len(proteome_pdb_paths)} PDB structure files")

    logger.info("Reading UniProt accession descriptions...")
    with open(alphafold_uniprot_des_path, mode="r", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        uniprot_accession_desc = {
            row['Entry']: {'trans': row['Ensembl'], 'length': row['Length']} for row in reader
        }
    logger.info(f"Loaded {len(uniprot_accession_desc)} UniProt entries")

    logger.info("Reading BioMart canonical transcript mapping...")
    with open(mart_ensembl_canonical_mapping, "r") as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        mart_canonical_mapping = [row for row in csv_reader]
    logger.info(f"Loaded {len(mart_canonical_mapping)} BioMart rows")

    reviewed_uniprot_structures_ids = set(uniprot_accession_desc.keys())

    logger.info("Creating database schema...")
    create_database(database_path)

    logger.info("Populating Ensembl-UniProt mapping tables...")
    populate_ensembl_uniprot_mapping(
        mart_canonical_mapping=mart_canonical_mapping,
        reviewed_uniprot_structures_ids=reviewed_uniprot_structures_ids,
        pdb_structure_ids=pdb_structure_ids,
        uniprot_accession_desc=uniprot_accession_desc,
        database_path=database_path
    )

    logger.info("Populating transcript-isoform mapping table...")
    populate_transcript_uniprot_isoform_mapping(
        uniprot_accession_desc=uniprot_accession_desc,
        database_path=database_path
    )

    logger.info("Parsing PDB structures and populating structure tables...")
    proteome_pdb_tables(
        proteome_pdb_paths=proteome_pdb_paths,
        database_path=database_path
    )

    logger.info("Database build complete.")


def populate_ensembl_uniprot_mapping(
        mart_canonical_mapping: list,
        reviewed_uniprot_structures_ids: set,
        pdb_structure_ids: set,
        uniprot_accession_desc: dict,
        database_path: Path,
):
    ensembl_genes_and_transcripts = set()
    ensembl_uniprot_mapping = set()
    for i in mart_canonical_mapping:
        if i['UniProtKB/Swiss-Prot ID'] in reviewed_uniprot_structures_ids:
            ensembl_genes_and_transcripts.add(
                (i['Gene stable ID'],
                 i['Transcript stable ID'],
                 1 if i['Ensembl Canonical'] else 0)
            )
            ensembl_uniprot_mapping.add((
                i['Gene stable ID'],
                i['UniProtKB/Swiss-Prot ID'],
                uniprot_accession_desc[i['UniProtKB/Swiss-Prot ID']]['length'],
                i['UniProtKB/Swiss-Prot ID'] in pdb_structure_ids,
                uniprot_accession_desc[i['UniProtKB/Swiss-Prot ID']]['trans']
            ))
    with sqlite3.connect(database_path) as db:
        cursor = db.cursor()
        cursor.executemany("""
        INSERT INTO ensembl_genes_and_transcripts (
            gene_stable_id,
            transcript_stable_id,
            canonical
            ) VALUES (?, ?, ?)
            """, list(ensembl_genes_and_transcripts))
        cursor.executemany("""
        INSERT INTO ensembl_uniprot_mapping (
            gene_stable_id,
            uniprot_swiss_prot_id,
            length,
            structure,
            uniprot_ensembl_cross_reference
            ) VALUES (?, ?, ?, ?, ?)
            """, list(ensembl_uniprot_mapping))
    logger.info(
        f"Inserted {len(ensembl_genes_and_transcripts)} gene-transcript pairs "
        f"and {len(ensembl_uniprot_mapping)} gene-UniProt mappings"
    )


def populate_transcript_uniprot_isoform_mapping(
        uniprot_accession_desc: dict,
        database_path: Path,
):
    tuples_to_insert = set()
    skipped = 0
    for key, value in uniprot_accession_desc.items():
        transcp = value['trans']
        if not isinstance(transcp, str) or not transcp.strip():
            skipped += 1
            continue
        for trans in [i.split() for i in transcp.split(';') if i]:
            if len(trans) > 1:
                transcript_id, uniprot_isoform = trans
                uniprot_id, isoform = uniprot_isoform.split('[')[1].split(']')[0].split('-')
                tuples_to_insert.add(
                    (transcript_id.split('.')[0], transcript_id, uniprot_id, f'F{isoform}')
                )
            else:
                tuples_to_insert.add(
                    (trans[0].split('.')[0], trans[0], key, 'F1')
                )
    with sqlite3.connect(database_path) as db:
        cursor = db.cursor()
        cursor.executemany("""
        INSERT INTO transcript_uniprot_isoform_mapping (
            transcript_id,
            transcript_id_version,
            uniprot_swiss_prot_id,
            isoform
        ) VALUES (?, ?, ?, ?)
        """, list(tuples_to_insert))
    logger.info(
        f"Inserted {len(tuples_to_insert)} transcript-isoform mappings "
        f"({skipped} entries skipped due to missing Ensembl cross-reference)"
    )


def create_database(database_path: Path):
    with sqlite3.connect(database_path) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS ensembl_genes_and_transcripts (
            gene_stable_id TEXT,
            transcript_stable_id TEXT,
            canonical BINARY(1) DEFAULT 0,
            PRIMARY KEY(gene_stable_id, transcript_stable_id)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS ensembl_uniprot_mapping (
            gene_stable_id TEXT,
            uniprot_swiss_prot_id TEXT,
            length INT,
            structure BINARY(1) DEFAULT 0,
            uniprot_ensembl_cross_reference TEXT,
            PRIMARY KEY(gene_stable_id, uniprot_swiss_prot_id)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS pdb_structures_description (
            uniprot_swiss_prot_id TEXT,
            isoform TEXT,
            n_chains INT,
            average_plddt REAL,
            std_plddt REAL,
            PRIMARY KEY(uniprot_swiss_prot_id, isoform)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS pdb_structures_chain_plddt (
            uniprot_swiss_prot_id TEXT,
            isoform TEXT,
            chain_id TEXT,
            residue_plddt_by_position TEXT,
            PRIMARY KEY(uniprot_swiss_prot_id, isoform, chain_id)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS pdb_structures_chain_sequence (
            uniprot_swiss_prot_id TEXT,
            isoform TEXT,
            chain_id VARCHAR,
            chain_amino_acid_sequence TEXT,
            PRIMARY KEY(uniprot_swiss_prot_id, isoform, chain_id)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS transcript_uniprot_isoform_mapping (
            transcript_id TEXT,
            transcript_id_version TEXT,
            uniprot_swiss_prot_id TEXT,
            isoform TEXT,
            PRIMARY KEY(transcript_id, uniprot_swiss_prot_id, isoform)
        )
        """)
    logger.info("Database schema created successfully")


def proteome_pdb_tables(
        proteome_pdb_paths: list,
        database_path: Path,
):
    total = len(proteome_pdb_paths)
    with sqlite3.connect(database_path, timeout=60) as db:
        cursor = db.cursor()
        for idx, pdb_path in enumerate(proteome_pdb_paths):
            if idx % 500 == 0:
                logger.info(f"Processing structures: {idx}/{total}")
            uniprot_id, isoform = pdb_path.name.split('-model')[0].split('AF-')[1].split('-')
            try:
                metrics = extract_alphafold_metrics_and_sequence(pdb_path)
            except Exception as e:
                logger.warning(f"Failed to parse {pdb_path.name}: {e}")
                continue
            cursor.execute("""
                INSERT INTO pdb_structures_description VALUES (?, ?, ?, ?, ?)
                """, (uniprot_id, isoform, metrics["num_chains"],
                      metrics["average_plddt"], metrics["std_plddt"])
            )
            for chain_id, chain_dict in metrics['sequences'].items():
                cursor.execute("""
                INSERT INTO pdb_structures_chain_plddt VALUES (?, ?, ?, ?)
                """, (uniprot_id, isoform, chain_id, chain_dict['plddt'])
                )
                cursor.execute("""
                INSERT INTO pdb_structures_chain_sequence VALUES (?, ?, ?, ?)
                """, (uniprot_id, isoform, chain_id, chain_dict['aa'])
                )
        db.commit()
    logger.info(f"Finished processing {total} PDB structures")


if __name__ == '__main__':
    root = Path('/home/mahe9165/exon-structure-evolution/data')
    alphafold_db_path = root / 'alphafold/UP000005640_9606_HUMAN_v6'
    alphafold_uniprot_des_path = root / 'alphafold/tab_data/uniprotkb_proteome_UP000005640_reviewed_2026_05_31.tsv'
    mart_ensembl_canonical_mapping = root / 'biomart/tab_data/mart_export_canonical_transcripts_v115.txt'
    database_path = root / 'alphafold/UP000005640_uniprot_accession_mapping_v115.db'

    main(
        alphafold_db_path=alphafold_db_path,
        alphafold_uniprot_des_path=alphafold_uniprot_des_path,
        mart_ensembl_canonical_mapping=mart_ensembl_canonical_mapping,
        database_path=database_path
    )