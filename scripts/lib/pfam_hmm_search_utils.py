"""
pfam_hmm_search_utils.py

Updated utility functions for Pfam-based domain analysis.

Key differences from ecod_hmm_search_utils.py:
  - build_pfam_xgroup_mapping: maps Pfam accessions -> ECOD X-group IDs
    using ecod.latest.f_id_pfam_acc.txt + ecod.latest.domains.txt
  - build_pfam_clan_mapping: reads clan assignments from Pfam-A.hmm.dat
    (replaces symmetry score as subdomain indicator)
  - collect_hmmscan_matches: same logic as collect_hmm_matches but handles
    hmmscan output format (query/target swapped vs hmmsearch)
  - filter_matches_by_coverage: queries raw hmmscan table and maps
    Pfam accession -> X-group, applies coverage thresholds
  - get_mrna_matches_domain_subdomain_categories: clan_as_subdomain flag
    replaces symmetry_score_threshold
"""

import json
import os
import re
import subprocess
import pickle
from collections import defaultdict
from functools import reduce
from pathlib import Path
from typing import Union

import portion as P
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ============================================================
# Pfam -> X-group mapping
# ============================================================

def build_pfam_xgroup_mapping(
        ecod_domains_path: Path,
        pfam_map_path: Path,
) -> dict:
    """
    Build a mapping from Pfam accession -> (x_group_id, x_group_name, h_group_name).

    Uses two files:
      - ecod.latest.f_id_pfam_acc.txt: maps F-group ID -> Pfam accession(s)
      - ecod.latest.domains.txt: maps F-group ID -> X-group hierarchy
    """
    # Step 1: parse ecod.latest.domains.txt -> f_id -> (x_group_id, x_name, h_name)
    f_id_to_xgroup = {}
    with open(ecod_domains_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            try:
                f_id      = parts[3]   # e.g. "1.1.1"
                x_name    = parts[9].strip('"')
                h_name    = parts[10].strip('"')
                x_group_id = int(f_id.split(".")[0])
                if f_id not in f_id_to_xgroup:
                    f_id_to_xgroup[f_id] = (x_group_id, x_name, h_name)
            except (IndexError, ValueError):
                continue

    # Step 2: parse ecod.latest.f_id_pfam_acc.txt -> pfam_acc -> f_id
    # Format: f_id <tab> pfam_accession(s) (comma-separated)
    pfam_to_xgroup = {}
    with open(pfam_map_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            f_id = parts[0].strip()
            pfam_accs = [p.strip() for p in parts[1].split(",") if p.strip()]
            xgroup_info = f_id_to_xgroup.get(f_id)
            if xgroup_info:
                for pfam_acc in pfam_accs:
                    # Handle ambiguous (multiple X-groups for one F-group) -> X-group 0
                    if pfam_acc in pfam_to_xgroup:
                        existing = pfam_to_xgroup[pfam_acc]
                        if existing[0] != xgroup_info[0]:
                            pfam_to_xgroup[pfam_acc] = (0, "AMBIGUOUS", "")
                    else:
                        pfam_to_xgroup[pfam_acc] = xgroup_info

    return pfam_to_xgroup


def build_pfam_clan_mapping(pfam_hmm_path: str) -> dict:
    """
    Parse Pfam-A.hmm.dat to build Pfam accession -> clan accession mapping.
    Clan membership replaces symmetry score as a subdomain indicator:
    if two exons share the same Pfam clan, they may be evolutionarily related
    subdomains of a larger repeat unit.

    Returns dict: {pfam_acc: clan_acc or None}
    """
    dat_path = pfam_hmm_path.replace(".hmm", ".hmm.dat")
    if not os.path.exists(dat_path):
        # Try adjacent .dat file
        dat_path = pfam_hmm_path + ".dat"
    if not os.path.exists(dat_path):
        return {}

    clan_map = {}
    current_acc = None
    current_clan = None
    with open(dat_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("ACC"):
                current_acc = line.split()[1].split(".")[0]  # strip version
                current_clan = None
            elif line.startswith("CL"):
                current_clan = line.split()[1]
            elif line == "//":
                if current_acc:
                    clan_map[current_acc] = current_clan
                current_acc = None
                current_clan = None
    return clan_map


# ============================================================
# hmmscan output parsing
# ============================================================

def collect_hmmscan_matches(hmmscan_file_path: str) -> list:
    """
    Parse hmmscan --domtblout output.

    Note: in hmmscan output (vs hmmsearch), query and target are swapped:
      - target = Pfam profile (col 0-1)
      - query  = your sequence (col 2-3)
    Column layout is otherwise identical to hmmsearch domtblout.
    """
    matches_list = []
    with open(hmmscan_file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 23:
                continue

            # In hmmscan: col0=profile_name, col1=profile_acc, col2=seq_name, col3=seq_acc
            target_name = fields[0]   # Pfam profile name
            target_acc  = fields[1]   # Pfam accession (PFxxxxx.N)
            query_name  = fields[2]   # your sequence ID (gene-transcript-cds or gene-transcript)
            qlen        = fields[5]   # query sequence length
            t_len       = fields[2]   # target (Pfam) model length -- col 5 is qlen

            # Re-extract fields matching the original collect_hmm_matches layout
            (
                _target_name, _target_acc, _query_name, _query_acc,
                seq_evalue, seq_score, seq_bias,
                domain_num, domain_total,
                domain_c_evalue, domain_i_evalue,
                domain_score, domain_bias,
                hmm_from, hmm_to,
                ali_from, ali_to,
                env_from, env_to,
                acc,
            ) = fields[:20]

            t_len = int(fields[2+1] if len(fields) > 3 else fields[5])  # profile length
            qlen  = int(fields[5])

            # HMMER coords are 1-based inclusive → convert to 0-based
            ali_from = int(ali_from) - 1
            ali_to   = int(ali_to)
            hmm_from = int(hmm_from) - 1
            hmm_to   = int(hmm_to)

            # target_name for hmmscan is "gene-transcript[-cds]"
            id_parts = query_name.rsplit("-")

            pfam_acc_clean = target_acc.split(".")[0]  # strip version number

            matches_list.append([
                *id_parts,
                qlen,                                         # sequence length
                target_name,                                  # Pfam profile name
                pfam_acc_clean,                               # Pfam accession
                t_len,                                        # HMM model length
                seq_evalue, seq_score, seq_bias,
                domain_num, domain_total,
                domain_c_evalue, domain_i_evalue,
                domain_score, domain_bias,
                hmm_from, hmm_to,
                ali_from, ali_to,
                env_from, env_to,
                acc,
                fields[22] if len(fields) > 22 else "",       # description
                round((hmm_to - hmm_from) / int(t_len), 3),  # query (HMM) coverage
                round((ali_to - ali_from) / int(qlen), 3),   # target (seq) coverage
            ])

    return matches_list


# ============================================================
# Coverage filtering and X-group deduplication
# ============================================================

def filter_matches_by_coverage(
        db_path: Path,
        table: str,
        pfam_to_xgroup: dict,
        pfam_clan_map: dict,
        evalue_threshold: float,
        query_cov_threshold: float,
        target_cov_threshold: float,
        is_cds: bool,
) -> list:
    """
    Query raw hmmscan results table, apply coverage/evalue filters,
    map Pfam accession -> X-group, attach clan info (replaces symmetry score).

    For CDS: requires both query and target coverage >= threshold (full domain).
    For mRNA: requires only query (HMM) coverage >= threshold.
    """
    import sqlite3
    id_col = "CdsID" if is_cds else None
    cov_condition = (
        f"QueryAlignPercentage >= {query_cov_threshold} "
        f"AND TargetAlignPercentage >= {target_cov_threshold}"
        if is_cds else
        f"QueryAlignPercentage >= {query_cov_threshold}"
    )

    select_cols = (
        "HmmCdsID, GeneID, TranscriptID, CdsID, CdsStart, QueryName, QueryAccession, "
        "QueryLen, HmmFrom, HmmTo, AliFrom, AliTo, DomainIEvalue"
        if is_cds else
        "HmmTranscriptID, GeneID, TranscriptID, QueryName, QueryAccession, "
        "QueryLen, HmmFrom, HmmTo, AliFrom, AliTo, DomainIEvalue"
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(f"""
            SELECT {select_cols}
            FROM {table}
            WHERE DomainIEvalue <= {evalue_threshold}
              AND {cov_condition}
        """)
        rows = cursor.fetchall()

    filtered = []
    for row in rows:
        pfam_acc = row[6] if is_cds else row[5]   # QueryAccession
        xgroup_info = pfam_to_xgroup.get(pfam_acc)
        if xgroup_info is None:
            xgroup_info = (0, "NO_X_NAME", "")
        x_group_id, x_name, h_name = xgroup_info
        clan = pfam_clan_map.get(pfam_acc)
        filtered.append((*row, x_group_id, x_name, h_name, clan))

    return filtered


def deduplicate_overlapping_xgroup_matches(
        matches: list,
        overlap_threshold: float,
) -> list:
    """
    For overlapping alignment regions with the same X-group,
    keep only the match with the lowest i-Evalue.
    Mirrors the original logic from ecod_hmm_search_utils.py.
    """
    # Group by (gene_id, x_group_id)
    grouped = defaultdict(list)
    for match in matches:
        gene_id    = match[1]
        x_group_id = match[-3]  # appended x_group_id
        grouped[(gene_id, x_group_id)].append(match)

    deduplicated = []
    for (gene_id, x_group_id), group in grouped.items():
        # Sort by evalue ascending
        group = sorted(group, key=lambda m: float(m[-6]) if m[-6] != "." else 1.0)
        kept = []
        for match in group:
            ali_from = match[-8] if len(match) > 10 else match[9]
            ali_to   = match[-7] if len(match) > 10 else match[10]
            coord    = P.open(ali_from, ali_to)
            overlap  = any(
                min_perc_overlap(coord, P.open(k[-8], k[-7])) >= overlap_threshold
                for k in kept
            )
            if not overlap:
                kept.append(match)
        deduplicated.extend(kept)

    return deduplicated


# ============================================================
# Coordinate mapping (unchanged from original)
# ============================================================

def get_coordinates_list(gene_mapping_coords_str: str) -> list:
    """Parse stored genomic coordinate string back into P.Interval list."""
    return [
        P.open(int(coord.split(",")[0].replace("(", "")),
               int(coord.split(",")[1].replace(")", "")))
        for coord in gene_mapping_coords_str.split("|")
        if coord.strip()
    ]


def get_p_coordinates_list_union(p_coordinates_list: list) -> P.Interval:
    return reduce(lambda a, b: a | b, p_coordinates_list)


def min_perc_overlap(intv_i: P.Interval, intv_j: P.Interval) -> float:
    if intv_i.empty or intv_j.empty:
        return 0.0
    intersection = intv_i & intv_j
    if intersection.empty:
        return 0.0
    len_i = intv_i.upper - intv_i.lower
    len_j = intv_j.upper - intv_j.lower
    if len_i == 0 or len_j == 0:
        return 0.0
    intersection_len = intersection.upper - intersection.lower
    return intersection_len / min(len_i, len_j)


def intervals_percentage_overlap(intv_i: P.Interval, intv_j: P.Interval) -> float:
    if intv_i.empty or intv_j.empty:
        return 0.0
    intersection = intv_i & intv_j
    if intersection.empty:
        return 0.0
    len_i = intv_i.upper - intv_i.lower
    if len_i == 0:
        return 0.0
    return (intersection.upper - intersection.lower) / len_i


def get_candidate_cds_coordinates(
        gene_id: str,
        gene_hierarchy_dictionary: dict,
        min_exon_length: int = 0,
) -> dict:
    gene_dict = gene_hierarchy_dictionary[gene_id]
    candidates = set()
    for transcript_dict in gene_dict["mRNAs"].values():
        for cds in transcript_dict["structure"]:
            if cds["type"] == "CDS":
                coord = cds["coordinate"]
                if coord.upper - coord.lower >= min_exon_length:
                    candidates.add(coord)
    return {"candidates_cds_coordinates": candidates}


def map_peptide_to_genomic_coords(
        mrna_events: list,
        gene_hierarchy: dict,
) -> list:
    """
    Map peptide alignment coordinates back to genomic coordinates.
    Replicates the logic from prot_utils.peptide_nucleotide_mapping.
    Returns list of (mapping_string, hmm_transcript_id) for database update.
    """
    from prot_utils import get_cds_coordinates_list, peptide_nucleotide_mapping, extract_ranges

    mapping_tuples = []
    for record in mrna_events:
        (hmm_transcript_id, gene_id, query_name,
         x_group_name, x_group_id, ali_from, ali_to,
         symmetry_score, _) = record

        ensembl_gene_id = f"gene:{gene_id}" if ":" not in gene_id else gene_id
        if ensembl_gene_id not in gene_hierarchy:
            continue

        gene_dict = gene_hierarchy[ensembl_gene_id]
        transcript_id = f"transcript:{gene_id.split('-')[1]}" if "-" in gene_id else None
        if transcript_id not in gene_dict["mRNAs"]:
            continue

        transcript_dict = gene_dict["mRNAs"][transcript_id]
        cds_coords = get_cds_coordinates_list(
            strand=gene_dict["strand"],
            transcript_dictionary=transcript_dict
        )
        pep_nuc_map = peptide_nucleotide_mapping(
            cds_coordinates_list=cds_coords,
            gene_start=gene_dict["coordinate"].lower,
        )
        codon_indices = pep_nuc_map[ali_from * 3: ali_to * 3]
        if not codon_indices:
            continue
        ranges = extract_ranges(codon_indices)
        mapping_str = "|".join(
            f"({r.lower},{r.upper})" for r in ranges
        )
        mapping_tuples.append((mapping_str, hmm_transcript_id))

    return mapping_tuples


def verify_f_group_matches_are_non_overlapping(
        hmm_mrna_matches_dict: dict,
        overlapping_threshold: float,
) -> bool:
    for gene_id, transcript_dict in hmm_mrna_matches_dict.items():
        for transcript_id, pfam_dict in transcript_dict.items():
            for pfam_acc, matches in pfam_dict.items():
                coords = [coord for _, coord in matches]
                for i, ci in enumerate(coords):
                    for cj in coords[i + 1:]:
                        if min_perc_overlap(ci, cj) >= overlapping_threshold:
                            return True
    return False


# ============================================================
# Domain/subdomain classification (updated for Pfam clans)
# ============================================================

def get_hmm_mrna_cds_search_matches_dictionary(
        hmm_matches_mrna_list: list,
        hmm_matches_cds_list: list,
        overlap_threshold: float,
) -> dict:
    """Combines mRNA and CDS search results non-redundantly. Unchanged from original."""
    gene_matches_dictionary = defaultdict(lambda: defaultdict(list))

    for gene_id, hmm_mrna_id, x_group_id, mapping_to_gene_intervals, clan in hmm_matches_mrna_list:
        gene_matches_dictionary[gene_id][x_group_id].append((
            "mRNA",
            hmm_mrna_id,
            get_p_coordinates_list_union(
                get_coordinates_list(mapping_to_gene_intervals)
            ),
            clan   # replaces symmetry_score
        ))

    for gene_id, hmm_cds_id, x_group_id, cds_start, ali_from, ali_to, clan in hmm_matches_cds_list:
        align_start = cds_start + (ali_from * 3)
        align_end   = cds_start + (ali_from * 3) + (ali_to * 3)
        gene_coord  = P.open(align_start, align_end)

        if gene_matches_dictionary[gene_id][x_group_id]:
            if not any(
                    min_perc_overlap(gene_coord, mrna_coord) >= overlap_threshold
                    for _, _, mrna_coord, _ in gene_matches_dictionary[gene_id][x_group_id]
            ):
                gene_matches_dictionary[gene_id][x_group_id].append(
                    ("CDS", hmm_cds_id, gene_coord, clan)
                )
        else:
            gene_matches_dictionary[gene_id][x_group_id].append(
                ("CDS", hmm_cds_id, gene_coord, clan)
            )

    return gene_matches_dictionary


def get_mrna_matches_domain_subdomain_categories(
        query_mrna_hmm_matches_dictionary: dict,
        genes_dictionary: dict,
        exonize_events_dictionary: dict,
        gene_hierarchy_dictionary: dict,
        overlapping_threshold: float = 0.9,
        clan_as_subdomain: bool = True,
) -> tuple:
    """
    Classify mRNA HMM matches as full_domain / subdomain / within_domain_boundary.

    clan_as_subdomain=True: uses Pfam clan membership instead of symmetry score.
    A match is classified as SubDomainExon if it partially overlaps a domain
    AND belongs to a Pfam clan (indicating it's part of a repeat/subdomain family).
    """
    groups_missing_clan = set()
    tuples_to_insert = []

    for gene_id, mrna_hmm_matches in query_mrna_hmm_matches_dictionary.items():
        ensembl_gene_id = f"gene:{gene_id}"
        coords, chromosome, strand = genes_dictionary[ensembl_gene_id]
        cds_coords_set = set([
            P.open(c.lower - coords.lower, c.upper - coords.lower)
            for c in get_candidate_cds_coordinates(
                gene_id=ensembl_gene_id,
                gene_hierarchy_dictionary=gene_hierarchy_dictionary,
                min_exon_length=0,
            )["candidates_cds_coordinates"]
        ])

        gene_events = {
            coord
            for events_list in exonize_events_dictionary.get(gene_id, {}).values()
            for coord, *_ in events_list
        }
        cds_coords_set = cds_coords_set.union(gene_events)

        for full_domain_match in mrna_hmm_matches:
            mrna_match_id, *_, x_group_id, clan, _, gene_mapping_coord = full_domain_match
            is_within_domain = is_full_domain = is_subdomain = 0

            cds_domain_matches = [
                c for c in cds_coords_set
                if min_perc_overlap(gene_mapping_coord, c) >= overlapping_threshold
            ]

            if cds_domain_matches:
                is_full_domain += len(set(cds_domain_matches))
            else:
                cds_subdomain_matches = [
                    c for c in cds_coords_set
                    if intervals_percentage_overlap(gene_mapping_coord, c) >= overlapping_threshold
                ]
                if cds_subdomain_matches:
                    n = len(set(cds_subdomain_matches))
                    is_within_domain += n
                    if clan_as_subdomain:
                        if clan is not None:
                            is_within_domain -= n
                            is_subdomain += n
                        else:
                            groups_missing_clan.add(x_group_id)
                    # else: use symmetry_score (legacy mode) — not implemented here

            tuples_to_insert.append((
                is_within_domain, is_full_domain, is_subdomain, mrna_match_id
            ))

    return tuples_to_insert, groups_missing_clan


def get_events_domain_subdomain_categories(
        exonize_events_dictionary: dict,
        hmm_mrna_cds_matches_dict: dict,
        gene_hierarchy_dictionary: dict,
        overlapping_threshold: float = 0.9,
        clan_as_subdomain: bool = True,
) -> tuple:
    """Classify expansion events as full_domain/subdomain. Updated for Pfam clans."""
    x_group_ids_missing_clan = set()
    records_to_insert = []
    records_to_insert_matches = []

    for gene_id, expansion_dict in exonize_events_dictionary.items():
        trimmed_gene_id = gene_id.rsplit(":")[1] if ":" in gene_id else gene_id
        gene_start = gene_hierarchy_dictionary[gene_id]["coordinate"].lower

        mrna_hmm_matches = hmm_mrna_cds_matches_dict.get(trimmed_gene_id, {})

        for expansion_id, events_list in expansion_dict.items():
            for event in events_list:
                event_coord, event_type = event
                is_within = is_full = is_sub = 0

                full_domain_hits = [
                    (match_id, match_type, gene_coord, event_coord)
                    for x_group_id, coords_list in mrna_hmm_matches.items()
                    for match_type, match_id, gene_coord, _ in coords_list
                    if min_perc_overlap(gene_coord, event_coord) >= overlapping_threshold
                ]

                if full_domain_hits:
                    is_full = len(full_domain_hits)
                    records_to_insert_matches.extend([
                        (match_type, "FullDomainExon", match_id,
                         trimmed_gene_id, ev.lower, ev.upper)
                        for match_id, match_type, _, ev in full_domain_hits
                    ])
                else:
                    subdomain_hits = [
                        (x_group_id, match_id, match_type, gene_coord, event_coord, clan)
                        for x_group_id, coords_list in mrna_hmm_matches.items()
                        for match_type, match_id, gene_coord, clan in coords_list
                        if intervals_percentage_overlap(
                            gene_coord, event_coord
                        ) >= overlapping_threshold
                    ]
                    for (x_group_id, match_id, match_type,
                         gene_coord, ev_coord, clan) in subdomain_hits:
                        is_within += 1
                        if clan_as_subdomain:
                            if clan is not None:
                                records_to_insert_matches.append((
                                    match_type, "SubDomainExon", match_id,
                                    trimmed_gene_id, ev_coord.lower, ev_coord.upper
                                ))
                                is_sub += 1
                                is_within -= 1
                            else:
                                records_to_insert_matches.append((
                                    match_type, "WithinDomainBoundary", match_id,
                                    trimmed_gene_id, ev_coord.lower, ev_coord.upper
                                ))
                                x_group_ids_missing_clan.add(x_group_id)

                records_to_insert.append((
                    is_within, is_full, is_sub,
                    gene_id,
                    event_coord.lower + gene_start,
                    event_coord.upper + gene_start,
                ))

    return records_to_insert, records_to_insert_matches, x_group_ids_missing_clan
