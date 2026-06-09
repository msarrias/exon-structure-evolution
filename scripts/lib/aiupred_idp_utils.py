import argparse
import pickle
from collections.abc import Iterator
from typing import Union
import portion as P
from lib.aiupred_idp_sqlite_utils import insert_aiupred_scores
import aiupred


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'protein_db_path',
        type=str,
        help='protein database path',
    )
    parser.add_argument(
        'exonize_db_path',
        type=str,
        help='exonize database path',
    )
    parser.add_argument(
        'gene_hierarchy_path',
        type=str,
        help='exonize database path',
    )
    parser.add_argument(
        'iupred3_path',
        type=str,
        help='path to iupred3.py',
    )
    parser.add_argument(
        '--output-db-path',
        type=str,
        default='iupred3_protein_scores.db',
        help='output database path'
    )
    parser.add_argument(
        '--anchor',
        action='store_true',
        default=False,
        help='enable ANCHOR2 prediction',
    )
    parser.add_argument(
        '--min-exon-length',
        type=int,
        default=20,
        help='minimum exon length for CDS coordinates',
    )
    parser.add_argument(
        '--smoothing',
        type=str,
        default="medium",
        help='smoothing type: "no", "medium" or "strong". Default is "medium"',
    )
    parser.add_argument(
        '--iupred-type',
        type=str,
        default='long',
        help='Analysis type: "long", "short" or "glob"',
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=5,
        help='batch size for processing the sequences',
    )
    parser.add_argument(
        '--ordered-score-threshold',
        type=float,
        default=0.5,
        help='threshold for considering a protein as ordered',
    )
    args = parser.parse_args()
    if args.iupred_type not in ['long', 'short', 'glob']:
        raise ValueError(
            'Analysis type (--iupred-type) must be either long, short or glob!'
        )
    if args.smoothing not in ['no', 'medium', 'strong']:
        raise ValueError(
            'Smoothing (--smoothing) must be either no, medium or strong!'
        )
    return args


def read_pkl_file(
        file_path: str,
) -> dict:
    with open(file_path, 'rb') as handle:
        read_file = pickle.load(handle)
    return read_file


def dump_sequence_dictionary_into_fasta_file(
        out_file_path: str,
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


def run_iupred3(
        iupred3_path: str,
        input_fasta_path: str,
        iupred_type: str,
        anchor: bool,
        smoothing: str,
):
    iupred_call = [
        'python3',
        iupred3_path,
        input_fasta_path,
        iupred_type,
        '-a' if anchor else '',
        '-s', smoothing
    ]
    result = run(
        iupred_call,
        capture_output=True,
        text=True
    )
    return result


def parse_iupred3_output(
        subprocess_result,
        anchor: bool
) -> list:
    scores_list = []
    anchor_list = []
    for line in subprocess_result.stdout.split('\n'):
        if not line.startswith('#') and line.strip():
            parts = line.split('\t')
            if len(parts) > 1:
                if anchor:
                    anchor_list.append(parts[-1])
                    scores_list.append(parts[-2])
                else:
                    scores_list.append(parts[-1])
    if anchor:
        return ['_'.join(scores_list), '_'.join(anchor_list)]
    return ['_'.join(scores_list)]


def even_batches(
            data: list,
            number_of_batches: int = 1,
) -> Iterator[list]:
    """
    Given a list and a number of batches, returns 'number_of_batches'
     consecutive subsets elements of 'data' of even size each, except
     for the last one whose size is the remainder of the division of
    the length of 'data' by 'number_of_batches'.
    """
    # We round up to the upper integer value to guarantee that there
    # will be 'number_of_batches' batches
    even_batch_size = (len(data) // number_of_batches) + 1
    for batch_number in range(number_of_batches):
        batch_start_index = batch_number * even_batch_size
        batch_end_index = min((batch_number + 1) * even_batch_size, len(data))
        yield data[batch_start_index:batch_end_index]


def process_gene_transcripts_and_update_database(
        list_transcripts_batch,
        gene_transcripts_dictionary: dict,
        iupred3_path: str,
        iupred_type: str,
        anchor: bool,
        smoothing: str,
        output_db_path: str
) -> None:
    """
    Processes gene transcripts in batches, runs iupred3 on each,
     and updates the database with the scores.
    :param list_transcripts_batch:
    :param gene_transcripts_dictionary: Dictionary of gene transcripts.
    :param iupred3_path: Path to the iupred3 executable.
    :param iupred_type: The type parameter for iupred3.
    :param anchor: The anchor parameter for iupred3.
    :param smoothing: The smoothing parameter for iupred3.
    :param output_db_path: Path to the database where the results will be stored.
    """
    for transcript in list_transcripts_batch:
        tuples_to_insert = []
        gene_id, peptide_sequence = gene_transcripts_dictionary[transcript]
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta_path = os.path.join(temporary_directory, 'seq.fa')
            dump_sequence_dictionary_into_fasta_file(
                out_file_path=fasta_path,
                sequence_dictionary={transcript: peptide_sequence}
            )
            iupred3_out = run_iupred3(
                iupred3_path=iupred3_path,
                input_fasta_path=fasta_path,
                iupred_type=iupred_type,
                anchor=anchor,
                smoothing=smoothing
            )
            iupred3_scores = parse_iupred3_output(
                subprocess_result=iupred3_out,
                anchor=anchor
            )
            if not any(seq == '' for seq in iupred3_scores):
                tuples_to_insert.append(
                    (gene_id, transcript, *iupred3_scores)
                )
        insert_iupred3_scores(
            db_path=output_db_path,
            tuples_list=tuples_to_insert,
            anchor=anchor
        )


def parse_scores_string(
        probability_string: str,
) -> list:
    probabilities = [float(prob) for prob in probability_string.split('_') if prob]
    return probabilities


def scores_to_string(scores_list):
    return '_'.join([str(score) for score in scores_list])


def get_ordered_scores_percentage(
        scores_list: list,
        threshold: float
) -> float:
    return round(len([score for score in scores_list if score <= threshold])/len(scores_list), 3)


def update_scores_cds_mapping(
        gene_cdss_dictionary: dict,
        protein_scores_dictionary: dict,
        anchor: bool,
        ordered_score_threshold: float
):
    tuples_to_insert = []
    for gene_id, transcript_dictionary in gene_cdss_dictionary.items():
        for transcript_id, cdss_list in transcript_dictionary.items():
            iupred_scores, anchor_scores = protein_scores_dictionary[gene_id][transcript_id]
            scores_list = parse_scores_string(
                probability_string=iupred_scores
            )
            if anchor:
                scores_list_anchor = parse_scores_string(
                    probability_string=anchor_scores
                )
                tuples_to_insert.extend([
                    (scores_to_string(scores_list[pep_cds_coord.lower: pep_cds_coord.upper]),
                     scores_to_string(scores_list_anchor[pep_cds_coord.lower: pep_cds_coord.upper]),
                     get_ordered_scores_percentage(
                         scores_list=scores_list[pep_cds_coord.lower: pep_cds_coord.upper],
                         threshold=ordered_score_threshold
                     ),
                     get_ordered_scores_percentage(
                         scores_list=scores_list_anchor[pep_cds_coord.lower: pep_cds_coord.upper],
                         threshold=ordered_score_threshold
                     ),
                     cds_id)
                    for pep_cds_coord, cds_id in cdss_list
                    if pep_cds_coord
                ])
            else:
                tuples_to_insert.extend([
                    (scores_to_string(scores_list[pep_cds_coord.lower: pep_cds_coord.upper]),
                     get_ordered_scores_percentage(
                         scores_list=scores_list[pep_cds_coord.lower: pep_cds_coord.upper],
                         threshold=ordered_score_threshold
                     ),
                     cds_id
                     )
                    for pep_cds_coord, cds_id in cdss_list
                    if pep_cds_coord
                ])
    return tuples_to_insert


def get_candidate_cds_coordinates(
        gene_id: str,
        gene_hierarchy_dictionary: dict,
        min_exon_length: int,
) -> dict:
    """
    get_candidate_cds_coordinates is a function that given a gene_id,
    collects all the CDS coordinates with a length greater than the
    minimum exon length across all transcript.
    If there are overlapping CDS coordinates, they are resolved
    according to the following criteria:

    * If they overlap by more than the overlapping threshold
      of both CDS lengths the one with the highest overlapping
      percentage is selected.
    * Both CDSs are selected otherwise

    The rationale behind choosing the interval with the higher percentage
    of overlap is that it will favor the selection of shorter intervals
    reducing the exclusion of short duplications due to coverage thresholds.

    :param gene_id: gene identifier
    :param gene_hierarchy_dictionary: gene hierarchy dictionary
    :param min_exon_length: minimum exon length
    :return: list of representative CDS coordinates across all transcripts
    """
    # collect all CDS coordinates and frames across all transcripts
    # we are interested in the frames to account for the unlikely event
    # that two CDS with same coordinates in different transcripts
    # have different frames
    cds_coordinates_and_frames: list[tuple[P.Interval, str]] = list(
        set(
            (coordinate, annotation_structure['frame'])
            for mrna_annotation in gene_hierarchy_dictionary[gene_id]['mRNAs'].values()
            for annotation_structure in mrna_annotation['structure']
            for coordinate in (annotation_structure['coordinate'],)
            if (
                    annotation_structure['type'] == 'CDS' and
                    (coordinate.upper - coordinate.lower) >= min_exon_length
            )
        )
    )
    if cds_coordinates_and_frames:
        representative_cds_frame_dictionary = dict()
        for cds_coordinate, frame in cds_coordinates_and_frames:
            if cds_coordinate in representative_cds_frame_dictionary:
                representative_cds_frame_dictionary[cds_coordinate] += f'_{str(frame)}'
            else:
                representative_cds_frame_dictionary[cds_coordinate] = str(frame)
        sorted_cds_coordinates_list: list[P.Interval] = sorted(
            [coordinate for coordinate, _ in cds_coordinates_and_frames],
            key=lambda coordinate: (coordinate.lower, coordinate.upper),
        )
        return dict(
            candidates_cds_coordinates=resolve_overlaps_between_coordinates(
                sorted_cds_coordinates=sorted_cds_coordinates_list
            ),
            cds_frame_dict=representative_cds_frame_dictionary,
        )
    return dict()


def resolve_overlaps_between_coordinates(
        sorted_cds_coordinates: list[P.Interval],
        overlapping_threshold: float = 0.9,
) -> list[P.Interval]:
    """
    resolve_overlaps_between_coordinates is a recursive function that given
    a list of coordinates, resolves overlaps between them.
    :param sorted_cds_coordinates: list of unique and sorted coordinates.
    :param overlapping_threshold: percentage of overlap
    :return: list of coordinates without "overlaps"
    """
    new_list = list()
    # We only process the first pair of overlapping intervals since
    # the resolved overlap could also overlap with the next interval in the list.
    intv_i, intv_j = get_first_overlapping_intervals(
        sorted_intervals=sorted_cds_coordinates,
        overlapping_threshold=overlapping_threshold
    )
    if all((intv_i, intv_j)):
        candidate = get_single_candidate_cds_coordinate(
            intv_i=intv_i,
            intv_j=intv_j
        )
        # Note: new_list is populated through enumeration
        # of 'sorted_cds_coordinates', so it will also be sorted
        # according to the same criteria used for 'sorted_cds_coordinates'.
        for idx, cds_coordinate in enumerate(sorted_cds_coordinates):
            if cds_coordinate != intv_i:
                new_list.append(cds_coordinate)
            else:
                new_list.append(candidate)
                new_list.extend(sorted_cds_coordinates[idx + 2:])
                return resolve_overlaps_between_coordinates(
                    sorted_cds_coordinates=new_list
                )
    return sorted_cds_coordinates


def min_perc_overlap(
        intv_i: P.Interval,
        intv_j: P.Interval,
) -> float:
    """
    Given two intervals, the function returns the percentage of the overlapping
    region relative to the longest interval. The percentage overlap of the shortest
    interval will always be greater or equal than that of the longest interval.
    :param intv_i:
    :param intv_j:
    :return:
    """

    def get_interval_length(
            interval: P.Interval,
    ):
        return sum(intv.upper - intv.lower for intv in interval)
    if intv_i.overlaps(intv_j):
        intersection_span = get_interval_length(intv_i.intersection(intv_j))
        longest_length = max(get_interval_length(intv_i), get_interval_length(intv_j))
        return round(intersection_span / longest_length, 3)
    return 0.0


def get_first_overlapping_intervals(
        sorted_intervals: list[P.Interval],
        overlapping_threshold: float,
) -> Union[tuple[P.Interval, P.Interval], tuple[None, None]]:
    """
    Given a list of intervals, get_first_overlapping_intervals returns
    the first consecutive interval tuples with overlapping pairs.
    :param sorted_intervals: list of intervals
    :param overlapping_threshold: percentage of overlap
    :return: overlapping intervals pair
    """
    first_overlap_index = 0
    while first_overlap_index < len(sorted_intervals) - 1:
        current_interval = sorted_intervals[first_overlap_index]
        next_interval = sorted_intervals[first_overlap_index + 1]
        if min_perc_overlap(
                    intv_i=current_interval,
                    intv_j=next_interval) >= overlapping_threshold:
            return current_interval, next_interval
        first_overlap_index += 1
    return None, None


def get_single_candidate_cds_coordinate(
        intv_i: P.Interval,
        intv_j: P.Interval,
) -> P.Interval:
    """
    Given two sorted intervals this function returns the interval with
    the largerst percentage of overlap.
    Assumption: intv_i =/= intv_j and they are sorted such that
    intv_i.lower < intv_j.lower.
    param intv_i: the first interval
    param intv_j: the second interval
    Returns: P.interval
    """
    if (
            get_overlap_percentage(
                intv_i=intv_i,
                intv_j=intv_j
            ) >
            get_overlap_percentage(
                intv_i=intv_j,
                intv_j=intv_i
            )
    ):
        return intv_j
    return intv_i


def get_overlap_percentage(
        intv_i: P.Interval,
        intv_j: P.Interval,
) -> float:
    """
    Given two intervals, the function get_overlap_percentage returns the percentage
    of the overlapping region relative to an interval j.
    """
    intersection = intv_i & intv_j
    if intersection:
        return (intersection.upper - intersection.lower) / (intv_j.upper - intv_j.lower)
    return 0.0


def get_representative_cds_tuples(
        gene_hierarchy_dictionary,
        raw_scores_dict,
        min_exon_length,
):
    tuples_to_insert = []
    gene_ids_list = list(gene_hierarchy_dictionary.keys())
    for gene_id in gene_ids_list:
        representative_cdss_dictionary = get_candidate_cds_coordinates(
            gene_id=gene_id,
            gene_hierarchy_dictionary=gene_hierarchy_dictionary,
            min_exon_length=min_exon_length
        )
        if representative_cdss_dictionary:
            gene_id = gene_id.rsplit('gene:')[1]
            for cds in representative_cdss_dictionary['candidates_cds_coordinates']:
                if not raw_scores_dict[gene_id][cds]:
                    print(gene_id, cds)
                tuples_to_insert.append(
                    (gene_id, cds.lower, cds.upper, *raw_scores_dict[gene_id][cds]))
    return tuples_to_insert
