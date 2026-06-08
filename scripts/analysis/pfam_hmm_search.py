"""
pfam_hmm_search.py

Updated ECOD/domain analysis pipeline using Pfam HMMs instead of ECODf.

Key changes from original ecod_hmm_search.py:
  - Uses Pfam-A.hmm instead of ecodf.hmm
  - Uses hmmscan instead of hmmsearch (query sequences against Pfam db)
  - Maps Pfam accessions -> ECOD F-group -> X-group via ecod.latest.f_id_pfam_acc.txt
  - Replaces symmetry score with Pfam clan membership for subdomain classification
  - Updated to Ensembl v115

Pipeline steps:
  1. Build protein database from gene hierarchy (reuses prot_db.py / prot_utils.py)
  2. Dump CDS and mRNA FASTA files
  3. Run hmmscan against Pfam-A.hmm
  4. Parse results, map Pfam -> X-group, filter by coverage/evalue
  5. Map peptide coordinates back to genomic coordinates
  6. Classify exons as full-domain / subdomain / fragment
  7. Annotate exon duplication events with domain classification
"""

import logging
import os
import pickle
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path

import portion as P

from lib.pfam_hmm_search_utils import (
    build_pfam_xgroup_mapping,
    build_pfam_clan_mapping,
    collect_hmmscan_matches,
    filter_matches_by_coverage,
    deduplicate_overlapping_xgroup_matches,
    get_candidate_cds_coordinates,
    map_peptide_to_genomic_coords,
    get_coordinates_list,
    get_p_coordinates_list_union,
    min_perc_overlap,
    intervals_percentage_overlap,
)
from lib.prot_sql_utils import (
    check_if_table_exists,
    create_mrna_transcripts_cds_tables,
    create_ecod_cds_mrna_hmm_results_tables,
    create_x_group_filtered_ecod_hmm_cdss_table,
    insert_into_ecod_hmm_cds_table,
    insert_into_ecod_hmm_mrna_table,
    insert_into_x_group_filtered_ecod_hmm_cdss_table,
    insert_into_x_group_filtered_ecod_hmm_mrna_table,
    update_mapping_to_gene_intervals,
    update_x_group_filtered_ecod_hmm_mrna_match_nature_columns,
    update_events_match_nature_columns,
    query_genes_dictionary,
    query_events_dictionary,
    query_ecod_mrna_events,
)
from lib.prot_utils import populate_cds_and_mrna_tables, read_genome
from lib.ecod_x_group_sql_utils import fetch_and_dump_data, query_for_cds
from lib.pfam_hmm_search_parser import hmm_search_argument_parser

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("pfam_hmm_search.log"),
    ]
)
logger = logging.getLogger(__name__)


def run_hmmscan(
        out_hmmscan_path: str,
        pfam_hmm_path: str,
        fasta_file_path: str,
        evalue: float,
        cpus: int = 8,
) -> None:
    """Run hmmscan against Pfam-A.hmm. Note: hmmscan (not hmmsearch) for Pfam."""
    cmd = [
        "hmmscan",
        "--domtblout", out_hmmscan_path,
        "--cpu", str(cpus),
        "-E", str(evalue),
        "--domE", str(evalue),
        pfam_hmm_path,
        fasta_file_path,
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    logger.info(f"hmmscan complete: {out_hmmscan_path}")


def get_mrna_hmm_matches_dict(hmm_records_list):
    """Build nested dict: gene_id -> transcript_id -> pfam_acc -> list of matches."""
    hmm_mrna_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for record in hmm_records_list:
        hmm_match_id, gene_id, transcript_id, pfam_acc, align_start, align_end = record
        hmm_mrna_dict[gene_id][transcript_id][pfam_acc].append([
            hmm_match_id, P.open(align_start, align_end)
        ])
    return hmm_mrna_dict


def build_gene_hierarchy_from_pickle(pkl_path: Path) -> dict:
    with open(pkl_path, "rb") as f:
        return pickle.load(f)


def main():
    args = hmm_search_argument_parser()

    pfam_hmm_path         = args.hmm_profile_path
    protein_db_path       = Path(args.protein_db_path)
    ecod_db_path          = Path(args.ecod_db_path)
    exonize_db_path       = Path(args.exonize_db_path)
    fasta_cds_path        = args.fasta_file_path_cdss
    fasta_mrna_path       = args.fasta_file_path_mrna
    genome_hierarchy_path = args.genome_hierarchy_pckl_path
    hmmscan_cds_path      = args.hmmsearch_cdss_path
    hmmscan_mrna_path     = args.hmmsearch_mrna_path
    ecod_domains_path     = args.ecod_latest_domains_txt_path
    pfam_xgroup_map_path  = args.ecod_latest_domains_pfam_map_path  # new arg
    timeout_db            = args.timeout_db
    evalue_threshold      = args.evalue_threshold
    query_cov_threshold   = args.query_coverage_threshold
    target_cov_threshold  = args.target_coverage_threshold
    overlap_threshold     = args.overlapping_threshold
    hard_force            = args.hard_force
    soft_force            = args.soft_force

    # --- Step 1: Build Pfam -> X-group mapping ---
    logger.info("Building Pfam accession -> X-group mapping...")
    pfam_to_xgroup = build_pfam_xgroup_mapping(
        ecod_domains_path=Path(ecod_domains_path),
        pfam_map_path=Path(pfam_xgroup_map_path),
    )
    pfam_clan_map = build_pfam_clan_mapping(pfam_hmm_path=pfam_hmm_path)
    logger.info(f"  {len(pfam_to_xgroup):,} Pfam families mapped to X-groups")
    logger.info(f"  {len(pfam_clan_map):,} Pfam families with clan assignments")

    # --- Step 2: Build protein database ---
    if hard_force or not protein_db_path.exists():
        logger.info("Building protein database...")
        gene_hierarchy = build_gene_hierarchy_from_pickle(genome_hierarchy_path)
        genome = read_genome(Path(args.genome_file_path))
        create_mrna_transcripts_cds_tables(db_path=protein_db_path)
        populate_cds_and_mrna_tables(
            gene_hierarchy_dictionary=gene_hierarchy,
            genome_dictionary=genome,
            protein_db_path=protein_db_path,
            timeout_db=timeout_db,
        )
    else:
        logger.info(f"Using existing protein database: {protein_db_path}")
        gene_hierarchy = build_gene_hierarchy_from_pickle(genome_hierarchy_path)

    # --- Step 3: Dump FASTA files ---
    logger.info("Dumping CDS FASTA file...")
    fetch_and_dump_data(
        file_path=Path(fasta_cds_path),
        query_function=query_for_cds,
        key_format="{0}-{1}-{2}",
        protein_db_path=protein_db_path,
        timeout_db=timeout_db,
    )
    logger.info("Dumping mRNA FASTA file...")
    fetch_and_dump_data(
        file_path=Path(fasta_mrna_path),
        query_function=lambda db, _: _query_mrna(db),
        key_format="{0}-{1}",
        protein_db_path=protein_db_path,
        timeout_db=timeout_db,
    )

    # --- Step 4: Run hmmscan ---
    if soft_force or hard_force or not os.path.exists(hmmscan_cds_path):
        logger.info("Running hmmscan on CDS sequences...")
        run_hmmscan(
            out_hmmscan_path=hmmscan_cds_path,
            pfam_hmm_path=pfam_hmm_path,
            fasta_file_path=fasta_cds_path,
            evalue=evalue_threshold,
        )
    if soft_force or hard_force or not os.path.exists(hmmscan_mrna_path):
        logger.info("Running hmmscan on mRNA sequences...")
        run_hmmscan(
            out_hmmscan_path=hmmscan_mrna_path,
            pfam_hmm_path=pfam_hmm_path,
            fasta_file_path=fasta_mrna_path,
            evalue=evalue_threshold,
        )

    # --- Step 5: Parse and store raw hmmscan results ---
    if soft_force or hard_force or not check_if_table_exists(
            ecod_db_path, "ECOD_hmm_CDSs", timeout_db
    ):
        logger.info("Parsing hmmscan results and storing to database...")
        create_ecod_cds_mrna_hmm_results_tables(
            db_path=ecod_db_path, timeout_db=timeout_db
        )
        cds_matches  = collect_hmmscan_matches(hmmscan_cds_path)
        mrna_matches = collect_hmmscan_matches(hmmscan_mrna_path)
        logger.info(f"  CDS raw matches:  {len(cds_matches):,}")
        logger.info(f"  mRNA raw matches: {len(mrna_matches):,}")

        insert_into_ecod_hmm_cds_table(
            db_path=ecod_db_path, timeout_db=timeout_db,
            matches_list=cds_matches
        )
        insert_into_ecod_hmm_mrna_table(
            db_path=ecod_db_path, timeout_db=timeout_db,
            matches_list=mrna_matches
        )

    # --- Step 6: Filter by coverage + evalue, map to X-groups ---
    if soft_force or hard_force or not check_if_table_exists(
            ecod_db_path, "x_group_filtered_ECOD_hmm_CDSs", timeout_db
    ):
        logger.info("Filtering matches and mapping to X-groups...")
        create_x_group_filtered_ecod_hmm_cdss_table(db_path=ecod_db_path)

        # CDS: filter + map
        cds_filtered = filter_matches_by_coverage(
            db_path=ecod_db_path,
            table="ECOD_hmm_CDSs",
            pfam_to_xgroup=pfam_to_xgroup,
            pfam_clan_map=pfam_clan_map,
            evalue_threshold=evalue_threshold,
            query_cov_threshold=query_cov_threshold,
            target_cov_threshold=target_cov_threshold,
            is_cds=True,
        )
        cds_dedup = deduplicate_overlapping_xgroup_matches(
            matches=cds_filtered,
            overlap_threshold=overlap_threshold,
        )
        insert_into_x_group_filtered_ecod_hmm_cdss_table(
            db_path=ecod_db_path, records=cds_dedup
        )
        logger.info(f"  CDS filtered matches: {len(cds_dedup):,}")

        # mRNA: filter + map
        mrna_filtered = filter_matches_by_coverage(
            db_path=ecod_db_path,
            table="ECOD_hmm_mRNA",
            pfam_to_xgroup=pfam_to_xgroup,
            pfam_clan_map=pfam_clan_map,
            evalue_threshold=evalue_threshold,
            query_cov_threshold=query_cov_threshold,
            target_cov_threshold=target_cov_threshold,
            is_cds=False,
        )
        mrna_dedup = deduplicate_overlapping_xgroup_matches(
            matches=mrna_filtered,
            overlap_threshold=overlap_threshold,
        )
        insert_into_x_group_filtered_ecod_hmm_mrna_table(
            db_path=ecod_db_path, records=mrna_dedup
        )
        logger.info(f"  mRNA filtered matches: {len(mrna_dedup):,}")

    # --- Step 7: Map peptide coordinates back to genomic coordinates ---
    logger.info("Mapping peptide coordinates to genomic coordinates...")
    mrna_events = query_ecod_mrna_events(db_path=ecod_db_path)
    mapping_tuples = map_peptide_to_genomic_coords(
        mrna_events=mrna_events,
        gene_hierarchy=gene_hierarchy,
    )
    update_mapping_to_gene_intervals(
        db_path=ecod_db_path,
        column_value_list=mapping_tuples,
    )

    # --- Step 8: Classify exon/event domain categories ---
    logger.info("Classifying exon domain categories...")
    genes_dict   = query_genes_dictionary(db_path=ecod_db_path)
    events_dict  = query_events_dictionary(db_path=ecod_db_path)

    from pfam_hmm_search_utils import (
        get_mrna_matches_domain_subdomain_categories,
        get_events_domain_subdomain_categories,
        get_hmm_mrna_cds_search_matches_dictionary,
    )
    from prot_sql_utils import query_hmm_searches_matches

    mrna_matches_list, cds_matches_list = query_hmm_searches_matches(
        db_path=ecod_db_path
    )
    hmm_mrna_cds_dict = get_hmm_mrna_cds_search_matches_dictionary(
        hmm_matches_mrna_list=mrna_matches_list,
        hmm_matches_cds_list=cds_matches_list,
        overlap_threshold=overlap_threshold,
    )

    # Classify mRNA matches
    mrna_tuples, _ = get_mrna_matches_domain_subdomain_categories(
        query_mrna_hmm_matches_dictionary=hmm_mrna_cds_dict,
        genes_dictionary=genes_dict,
        exonize_events_dictionary=events_dict,
        gene_hierarchy_dictionary=gene_hierarchy,
        overlapping_threshold=overlap_threshold,
        clan_as_subdomain=True,   # use Pfam clan instead of symmetry score
    )
    update_x_group_filtered_ecod_hmm_mrna_match_nature_columns(
        db_path=ecod_db_path,
        tuples_to_insert_list=mrna_tuples,
    )

    # Classify expansion events
    event_tuples, event_match_tuples, _ = get_events_domain_subdomain_categories(
        exonize_events_dictionary=events_dict,
        hmm_mrna_cds_matches_dict=hmm_mrna_cds_dict,
        gene_hierarchy_dictionary=gene_hierarchy,
        overlapping_threshold=overlap_threshold,
        clan_as_subdomain=True,
    )
    update_events_match_nature_columns(
        db_path=ecod_db_path,
        tuples_to_insert_list=event_tuples,
    )
    logger.info("Done.")


def _query_mrna(db_path):
    """Query full mRNA protein sequences."""
    with sqlite3.connect(db_path) as db:
        cursor = db.cursor()
        cursor.execute("""
            SELECT GeneID, TranscriptID, ProteinSequence
            FROM Gene_Transcript_Proteins
        """)
        return cursor.fetchall()


if __name__ == "__main__":
    main()
