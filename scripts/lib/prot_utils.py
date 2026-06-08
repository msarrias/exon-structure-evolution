import pickle
import portion as P
from Bio.Seq import Seq
from typing import Union
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import gzip
import sys
from prot_sql_utils import (insert_into_cdss_table,
                            insert_into_mrna_transcripts_table,
                            query_filtered_ecod_hmm_mrna_records)
from ecod_hmm_search_utils import verify_f_group_matches_are_non_overlapping


def read_genome(
        genome_file_path: Path
) -> dict:
    """
    read_genome is a function that reads a FASTA file and stores
    the masked/unmasked genome sequence in a dictionary.
    The dictionary has the following structure: {chromosome: sequence}
    """
    genome_dictionary = {}
    try:
        if genome_file_path.suffix == '.gz':
            with gzip.open(genome_file_path, mode='rt') as genome_file:  # 'rt' for textmode
                parsed_genome = SeqIO.parse(genome_file, 'fasta')
                for record in parsed_genome:
                    genome_dictionary[record.id] = str(record.seq)
        else:
            with open(genome_file_path, mode='r') as genome_file:
                parsed_genome = SeqIO.parse(genome_file, 'fasta')
                for record in parsed_genome:
                    genome_dictionary[record.id] = str(record.seq)
    except (ValueError, FileNotFoundError) as e:
        print(f"Incorrect genome file path: {e}")
        sys.exit()
    except OSError as e:
        print(f"There was an error reading the genome file: {e}")
        sys.exit()
    return genome_dictionary


def read_pkl_file(
        file_path: Path,
) -> dict:
    with open(file_path, 'rb') as handle:
        read_file = pickle.load(handle)
    return read_file


def construct_mrna_sequence(
        chromosome: str,
        strand: str,
        genome_dictionary: dict,
        coordinates_list: list
) -> str:
    mrna_sequence = ''
    for cds_coordinate_idx, coordinate in enumerate(coordinates_list):
        start, end = coordinate['coordinate'].lower, coordinate['coordinate'].upper
        exon_sequence = genome_dictionary[chromosome][start:end]
        if strand == '-':
            exon_sequence = str(Seq(exon_sequence).reverse_complement())
        mrna_sequence += exon_sequence
    return mrna_sequence


def get_cds_coordinates_list(
        strand: str,
        transcript_dictionary: dict
) -> list:
    cds_coordinates_list = [
        annotation for annotation in transcript_dictionary['structure']
        if annotation['type'] == 'CDS'
    ]
    if strand == '-':
        # if the gene is in the negative strand the direction of
        # transcription and translation is opposite to
        # the direction the DNA sequence is represented meaning that
        # translation starts from the  last CDS
        cds_coordinates_list = sorted(
            cds_coordinates_list,
            key=lambda x: (x['coordinate'].lower, x['coordinate'].upper),
            reverse=True
        )
    return cds_coordinates_list


def extract_ranges(
    nucleotides_indices_list: list
) -> list:
    coordinates_list = []
    start = nucleotides_indices_list[0]
    for nucleotide_idx, nucleotide in enumerate(nucleotides_indices_list[:-1]):
        next_nucleotide = nucleotides_indices_list[nucleotide_idx + 1]
        diff = next_nucleotide - nucleotide
        if diff != 1:
            coordinates_list.append(P.open(start, nucleotide + 1))
            start = next_nucleotide
        if nucleotide_idx == len(nucleotides_indices_list) - 2:
            coordinates_list.append(P.open(start, next_nucleotide + 1))
    return coordinates_list


def populate_cds_and_mrna_tables(
        gene_hierarchy_dictionary: dict,
        genome_dictionary: dict,
        protein_db_path: Path,
        timeout_db: int
) -> None:
    for gene_id, gene_dictionary in gene_hierarchy_dictionary.items():
        chromosome = gene_dictionary['chrom']
        strand = gene_dictionary['strand']
        mrna_transcripts_tuples_list = []
        for transcript_id, transcript_dictionary in gene_dictionary['mRNAs'].items():
            cds_coordinates_list = get_cds_coordinates_list(
                strand=strand,
                transcript_dictionary=transcript_dictionary
            )
            mrna_dna_sequence = construct_mrna_sequence(
                chromosome=chromosome,
                strand=strand,
                genome_dictionary=genome_dictionary,
                coordinates_list=cds_coordinates_list
            )
            peptide_sequence, cds_dna_peptide_sequences_tuples = construct_protein_sequence(
                gene_id=gene_id,
                transcript_id=transcript_id,
                mrna_sequence=mrna_dna_sequence,
                coordinates_list=cds_coordinates_list
            )
            mrna_transcripts_tuples_list.append(
                (
                    gene_id.rsplit(':')[1],
                    chromosome,
                    strand,
                    gene_dictionary['coordinate'].lower,
                    gene_dictionary['coordinate'].upper,
                    transcript_id.rsplit(':')[1],
                    transcript_dictionary['coordinate'].lower,
                    transcript_dictionary['coordinate'].upper,
                    peptide_sequence
                )
            )
            insert_into_cdss_table(
                db_path=protein_db_path,
                gene_args_tuple=cds_dna_peptide_sequences_tuples
            )
        insert_into_mrna_transcripts_table(
            db_path=protein_db_path,
            timeout_db=timeout_db,
            gene_args_tuple=mrna_transcripts_tuples_list
        )


def peptide_nucleotide_mapping(
    cds_coordinates_list: list,
    gene_start: int,
) -> list[int]:
    base_pair_indices_list = [
        base_pair_pos - gene_start
        for coordinate in cds_coordinates_list
        for base_pair_pos in range(coordinate['coordinate'].lower, coordinate['coordinate'].upper)
    ]
    start_coord = 0
    in_frame_codon_indices = []
    for cds_coordinate_idx, cds_coordinate in enumerate(cds_coordinates_list):
        frame = int(cds_coordinate['frame'])
        cds_start = cds_coordinate['coordinate'].lower
        cds_end = cds_coordinate['coordinate'].upper
        cds_length = cds_end - cds_start
        next_frame = 0
        end_coord = cds_length + start_coord
        if cds_coordinate_idx < len(cds_coordinates_list) - 1:
            next_frame = int(cds_coordinates_list[cds_coordinate_idx + 1]['frame'])
        base_pairs_cds_list = check_for_overhangs(
            sequence=base_pair_indices_list[start_coord + frame: end_coord + next_frame],
            cds_idx=cds_coordinate_idx,
             numb_cdss=len(cds_coordinates_list)
        )
        in_frame_codon_indices.extend(base_pairs_cds_list)
        start_coord, frame = end_coord, next_frame
    return in_frame_codon_indices


def get_hmm_mrna_matches_dict(
        hmm_records_list: list
) -> dict:
    hmm_mrna_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for record in hmm_records_list:
        hmm_match_id, gene_id, transcript_id, f_group, align_start, align_end = record
        hmm_mrna_dict[gene_id][transcript_id][f_group].append([
            hmm_match_id,
            P.open(align_start, align_end)
        ])
    return hmm_mrna_dict


def get_dna_peptide_mapping_tuples(
        db_path: Path,
        gene_hierarchy_dictionary: dict
) -> list:
    filtered_ecod_hmm_mrna_records = query_filtered_ecod_hmm_mrna_records(
        db_path=db_path
    )
    hmm_mrna_dict = get_hmm_mrna_matches_dict(
        hmm_records_list=filtered_ecod_hmm_mrna_records
    )
    if verify_f_group_matches_are_non_overlapping(
            hmm_mrna_matches_dict=hmm_mrna_dict,
            overlapping_threshold=0.9
    ):
        raise ValueError(
            'Overlapping matches found for f_group'
        )
    insert_dna_mapping_tuples_list = []
    for gene_id, trancript_matches_dictionary in hmm_mrna_dict.items():
        gene_dictionary = gene_hierarchy_dictionary[f'gene:{gene_id}']
        for transcript_id, f_groups_dictionary in trancript_matches_dictionary.items():
            transcript_dictionary = gene_dictionary['mRNAs'][f'transcript:{transcript_id}']
            cds_coordinates_list = get_cds_coordinates_list(
                strand=gene_dictionary['strand'],
                transcript_dictionary=transcript_dictionary
            )
            peptide_nucleotide_mapping_list = peptide_nucleotide_mapping(
                cds_coordinates_list=cds_coordinates_list,
                gene_start=gene_hierarchy_dictionary[f'gene:{gene_id}']['coordinate'].lower,
            )
            for f_group, matches_list in f_groups_dictionary.items():
                for match in matches_list:
                    hmm_match_id, align_coords = match
                    codon_ind = peptide_nucleotide_mapping_list[align_coords.lower * 3:align_coords.upper * 3]
                    peptide_nucleotide_mapping_ranges = extract_ranges(
                        nucleotides_indices_list=codon_ind
                    )
                    insert_dna_mapping_tuples_list.append(
                        ('_'.join([str(gene_range) for gene_range in peptide_nucleotide_mapping_ranges]),
                         hmm_match_id)
                    )
    return insert_dna_mapping_tuples_list


def construct_protein_sequence(
        gene_id: str,
        transcript_id: str,
        mrna_sequence: str,
        coordinates_list: list
) -> tuple[str, list]:

    cds_annotations_list = []
    mrna_peptide_sequence = ''
    overhangs_list = []
    start_coord = 0

    n_CDSs = len(coordinates_list)
    for cds_coordinate_idx, cds_coordinate in enumerate(coordinates_list):
        frame = int(cds_coordinate['frame'])
        cds_start = cds_coordinate['coordinate'].lower
        cds_end = cds_coordinate['coordinate'].upper
        cds_length = cds_end - cds_start
        frame_next = 0
        end_coord = cds_length + start_coord
        if cds_coordinate_idx != n_CDSs - 1:
            frame_next = int(coordinates_list[cds_coordinate_idx + 1]['frame'])
        exon_dna_sequence = check_for_overhangs(
            sequence=mrna_sequence[start_coord + frame: end_coord + frame_next],
            cds_idx=cds_coordinate_idx,
             numb_cdss=n_CDSs
        )
        exon_peptide_sequence = str(
            Seq(exon_dna_sequence).translate()
        )
        mrna_cds_start = start_coord + frame
        mrna_cds_end = end_coord + frame_next
        protein_cds_start = len(mrna_peptide_sequence)
        mrna_peptide_sequence += exon_peptide_sequence
        protein_cds_end = len(mrna_peptide_sequence)
        overhangs_list.append(len(exon_dna_sequence) % 3)
        start_coord, frame = end_coord, frame_next
        cds_annotations_list.append(
            (
                gene_id.rsplit(':')[1],
                transcript_id.rsplit(':')[1],
                cds_coordinate['id'].rsplit(':')[1],
                cds_coordinate_idx,
                frame,
                cds_start,
                cds_end,
                mrna_cds_start,
                mrna_cds_end,
                protein_cds_start,
                protein_cds_end,
                exon_dna_sequence,
                exon_peptide_sequence
            )
        )
    if max(overhangs_list) > 0:
        if not (overhangs_list[-1] != 0) and all(overhang == 0 for overhang in overhangs_list[:-1]):
            raise ValueError(
                f'Overhang in the middle of the sequence: {gene_id}, '
                f'{transcript_id}'
            )
    return mrna_peptide_sequence, cds_annotations_list


def check_for_overhangs(
        sequence: Union[str, list[int]],
        cds_idx: int,
         numb_cdss: int
) -> str | list[int]:
    overhang = len(sequence) % 3
    if overhang:
        if cds_idx == ( numb_cdss - 1):
            sequence = sequence[:-overhang]
        else:
            raise ValueError(
                'Overhang in the middle of the sequence: '
            )
    return sequence
