import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import defaultdict
import portion as P
from pathlib import Path
from ecod_x_group_sql_utils import query_x_group_ids_names
from ecod_hmm_search_utils import (
    min_perc_overlap,
    get_coordinates_list,
    get_p_coordinates_list_union
)


def dump_query_into_csv_file(
        output_file: Path,
        column_names_list: list,
        rows_list: list) -> None:
    with open(output_file, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(column_names_list)  # Write the column headers
        writer.writerows(rows_list)


def plot_side_by_side_piechart(
        pie_chart_hmm_data: list,
        pie_chart_ecod_data: list,
        x_names_dictionary: list,
        fig_path: str,
):
    np.random.seed(42)
    pie_chart_hmm_y = np.array([perc for _, perc in pie_chart_hmm_data])
    pie_chart_hmm_labels = [str(x_group_id) for x_group_id, _ in pie_chart_hmm_data]
    pie_chart_ecod_y = np.array([perc for _, perc in pie_chart_ecod_data])
    pie_chart_ecod_labels = [str(x_group_id) for x_group_id, _ in pie_chart_ecod_data]
    common_labels = list(set(pie_chart_hmm_labels + pie_chart_ecod_labels))

    # Get the initial set of colors from 'tab20c'
    colors_list = plt.cm.tab20c(np.linspace(0, 1, 20))

    # If more than 20 unique colors are needed, extend with 'tab20b'
    if len(common_labels) > 20:
        extended_colors = plt.cm.tab20b(np.linspace(0, 1, len(common_labels) - 20))
        colors_list = np.concatenate((colors_list, extended_colors))

    colors_list_hex = [
        '#%02x%02x%02x' % (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
        for color in colors_list
    ]

    colors_dict = {
        lable: color
        for lable, color in zip(common_labels, colors_list_hex)
    }
    plt.figure(figsize=(15, 10), dpi=100)
    # HMM search
    plt.subplot(1, 2, 1)
    pie_wedges_hmm = plt.pie(
        pie_chart_hmm_y,
        labels=pie_chart_hmm_labels,
        colors=[colors_dict[label] for label in pie_chart_hmm_labels],
        labeldistance=1.1  # Position the labels outside the slices
    )
    plt.title('HMM - Human Domains', fontweight='bold', fontsize=16)
    # ECOD human raw distribution
    plt.subplot(1, 2, 2)
    pie_wedges_ecod = plt.pie(
        pie_chart_ecod_y,
        colors=[colors_dict[label] for label in pie_chart_ecod_labels],
        labels=pie_chart_ecod_labels,
        labeldistance=1.1  # Position the labels outside the slices
    )
    plt.title('ECOD - Human Raw Domains', fontweight='bold', fontsize=16)
    wedges = pie_wedges_hmm[0] + pie_wedges_ecod[0]
    labels = pie_chart_hmm_labels + pie_chart_ecod_labels
    unique_labels = list(dict.fromkeys(labels))
    # Filter out wedges whose labels are already in the legend
    unique_wedges = []
    added_labels = set()
    for wedge, label in zip(wedges, labels):
        if label not in added_labels:
            unique_wedges.append(wedge)
            added_labels.add(label)
    plt.figlegend(
        unique_wedges,
        [f"{xid}-{x_names_dictionary[xid]}" for xid in unique_labels],
        loc='center',
        fontsize=14,
        bbox_to_anchor=(1.15, 0.5),
        labelspacing=0.15,
        frameon=False
    )
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.show()


def get_distribution_clustering_small_values(
        list_records: list,
        cluster_percentage_threshold: float = 0.01
) -> list:
    new_list_records = []
    less_percentage_cumm_sum = 0
    for record in list_records:
        support, percentage = record
        if percentage >= cluster_percentage_threshold:
            new_list_records.append((support, percentage))
        else:
            less_percentage_cumm_sum += percentage
    new_list_records.append((
        f"<{cluster_percentage_threshold}",
        less_percentage_cumm_sum
    ))
    return new_list_records


def get_hmm_mrna_cds_search_matches_dictionary(
        hmm_matches_mrna_list: list,
        hmm_matches_cds_list: list,
        overlap_threshold: float,
):
    # WE FIRST GO THROUGH THE MRNA SEARCH MATCHES
    gene_matches_dictionary = defaultdict(lambda: defaultdict(list))
    for gene_id, hmm_mrna_id, x_group_id, mapping_to_gene_intervals, symmetry_score in hmm_matches_mrna_list:
        gene_matches_dictionary[gene_id][x_group_id].append((
            'mRNA',
            hmm_mrna_id,
            get_p_coordinates_list_union(
                p_coordinates_list=get_coordinates_list(
                    gene_mapping_coords_str=mapping_to_gene_intervals
                )
            ),
            symmetry_score
        ))
    # WE COMPLEMENT THE PREVIOUS MATCHES WITH THE ONES FOUND IN THE CDS SEARCH MATCHES
    # THAT WERE NOT FOUND IN THE MRNA ONE - IF THERE ARE OVERLAPPING MATCHES WE KEEP THE
    # ONES FOUND IN THE MRNA SEARCH
    for gene_id, hmm_CDS_id, x_group_id, cds_start, ali_from, ali_to, symmetry_score in hmm_matches_cds_list:
        # the alignment coordinates map to a protein, we want gene coordinates
        align_start_gene_coordinate = cds_start + (ali_from * 3)
        align_end_gene_coordinate = cds_start + (ali_from * 3) + (ali_to * 3)
        gene_alignment_coordinate = P.open(align_start_gene_coordinate, align_end_gene_coordinate)
        if gene_matches_dictionary[gene_id][x_group_id]:
            if not any(
                    min_perc_overlap(
                        intv_i=gene_alignment_coordinate,
                        intv_j=mrna_match_coordinate,
                    ) >= overlap_threshold
                    for _, hmm_mrna_id, mrna_match_coordinate, symmetry_score in gene_matches_dictionary[gene_id][x_group_id]
            ):
                gene_matches_dictionary[gene_id][x_group_id].append(
                    ('CDS', hmm_CDS_id, gene_alignment_coordinate, symmetry_score)
                )
        else:
            gene_matches_dictionary[gene_id][x_group_id].append(
                ('CDS', hmm_CDS_id, gene_alignment_coordinate, symmetry_score)
            )
    return gene_matches_dictionary


def get_distribution_list(
        counts_dictionary: dict
) -> list:
    n_instances = sum(counts_dictionary.values())
    distribution_tuples = [
        (support, round(count / n_instances, 3))
        for support, count in counts_dictionary.items()
    ]
    return sorted(distribution_tuples, key=lambda x: -x[1])


def dump_venn_diagram_csv_input(
    genes_with_full_domain_duplication,
    genes_with_sub_domain_duplication,
    genes_with_exon_duplication,
    genes_with_full_domain_hmm,
    genes_with_sub_domain_hmm,
    output_csv_path='data/csvs/genes_sets2.csv'
) -> None:
    #
    mrna_match_mapping_to_sub_doman_domain = [gene for gene in genes_with_sub_domain_duplication
                                              if gene not in genes_with_sub_domain_hmm]
    if mrna_match_mapping_to_sub_doman_domain:
        genes_with_sub_domain_duplication = [gene for gene in genes_with_sub_domain_duplication
                                             if gene not in mrna_match_mapping_to_sub_doman_domain]

    # Find the maximum length
    max_length = max(
        len(genes_with_full_domain_duplication),
        len(genes_with_sub_domain_duplication),
        len(genes_with_exon_duplication),
        len(genes_with_full_domain_hmm),
        len(genes_with_sub_domain_hmm)
    )

    # Pad shorter lists with empty string
    genes_with_full_domain_duplication += [''] * (max_length - len(genes_with_full_domain_duplication))
    genes_with_sub_domain_duplication += [''] * (max_length - len(genes_with_sub_domain_duplication))
    genes_with_exon_duplication += [''] * (max_length - len(genes_with_exon_duplication))
    genes_with_full_domain_hmm += [''] * (max_length - len(genes_with_full_domain_hmm))
    genes_with_sub_domain_hmm += [''] * (max_length - len(genes_with_sub_domain_hmm))
    with open(output_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'Full Domain Exon Duplication',
            'Sub Domain Exon Duplication',
            'Genes with Exon Duplication',
            'Genes with Full Domain Exons',
            'Genes with Sub Domain Exons'
        ])
        for i in range(max_length):
            writer.writerow([
                genes_with_full_domain_duplication[i],
                genes_with_sub_domain_duplication[i],
                genes_with_exon_duplication[i],
                genes_with_full_domain_hmm[i],
                genes_with_sub_domain_hmm[i]
            ])


def get_hmm_matches_and_events_distributions(
        x_group_id_dictionary: dict,
        ecod_x_groups_distribution_dictionary: dict,
        clustering_percentage_threshold: float
) -> tuple:
    hmm_x_groups_distribution = get_distribution_clustering_small_values(
        list_records=get_distribution_list(
            counts_dictionary=x_group_id_dictionary
        ),
        cluster_percentage_threshold=clustering_percentage_threshold
    )
    ecod_x_groups_distribution = get_distribution_clustering_small_values(
        list_records=get_distribution_list(
            counts_dictionary=ecod_x_groups_distribution_dictionary
        ),
        cluster_percentage_threshold=clustering_percentage_threshold
    )
    return hmm_x_groups_distribution, ecod_x_groups_distribution


def get_x_groups_distribution(
    gene_matches_dictionary: dict,
    full_domain_exon_duplication_records: dict,
    overlap_threshold: float,
    db_path: str,
    clustering_percentage_threshold: float
):
    # X GROUP IDS GLOBAL DISTRIBUTION
    x_group_id_dictionary = defaultdict(int)
    for gene_id, x_group_dictionary in gene_matches_dictionary.items():
        for x_group_id, x_group_matches_coordinates_list in x_group_dictionary.items():
            x_group_id_dictionary[x_group_id] += len(x_group_matches_coordinates_list)
    # X GROUP IDS GLOBAL DISTRIBUTION FOR THE FOLLOWING CATEGORIES
    # 1. FULL DOMAIN EXON IN GENE WITHOUT DUPLICATION
    # 2. FULL DOMAIN EXON IN GENE WITH DUPLICATION
    # 3. FULL DOMAIN EXON DUPLICATION
    full_domain_without_duplication = defaultdict(int)
    full_domain_with_duplication = defaultdict(int)
    full_domain_duplication = defaultdict(int)
    for gene_id, x_group_dictionary in gene_matches_dictionary.items():
        # IS IT A GENE WITH DUPLICATIONS?
        # since it's a default dictionary
        if gene_id in full_domain_exon_duplication_records:
            for x_group_id, x_group_matches_list in x_group_dictionary.items():
                for *_, coordinate_match, _ in x_group_matches_list:
                    if any(min_perc_overlap(event_coordinate, coordinate_match)
                           for event_coordinate in full_domain_exon_duplication_records[gene_id]
                           ) >= overlap_threshold:
                        full_domain_duplication[x_group_id] += 1
                    else:
                        full_domain_with_duplication[x_group_id] += 1
        else:
            for x_group_id, x_group_matches_list in x_group_dictionary.items():
                full_domain_without_duplication[x_group_id] += len(x_group_matches_list)

    lables = get_x_groups_dictionary_count(
        db_path=db_path,
        percent_threshold=clustering_percentage_threshold
    )
    distribution_list = [
        [f"{lables[str(xgroup_id)][0]} ({global_count})",
         round(full_domain_duplication[xgroup_id] / global_count, 3) if xgroup_id in full_domain_duplication else 0,
         round(full_domain_with_duplication[xgroup_id] / global_count,
               3) if xgroup_id in full_domain_with_duplication else 0,
         round(full_domain_without_duplication[xgroup_id] / global_count,
               3) if xgroup_id in full_domain_without_duplication else 0
         ]
        for xgroup_id, global_count in x_group_id_dictionary.items()
        if (full_domain_with_duplication[xgroup_id] > 0 or full_domain_duplication[xgroup_id] > 0)
    ]
    # sort by the number of genes with full domain exon duplication
    distribution_list = sorted(distribution_list, key=lambda x: (-x[1], -x[2], -x[3]))
    column_list = [
        'x group ID',
        'Full domain exon duplication',
        'Full domain exon in gene with dup.',
        'Full domain exon in gene without dup.'
    ]
    return column_list, distribution_list


def get_x_groups_dictionary_count(
        db_path: str,
        percent_threshold: float,
):
    x_groups_list = query_x_group_ids_names(db_path=db_path)
    x_names_dict = defaultdict(list)
    x_names_dict[str(0)] = ['AMBIGUOUS']
    for x_group_id, x_group_name in x_groups_list:
        if x_group_name == "NO_X_NAME":
            x_group_name = f"X{x_group_id}"
        else:
            x_group_name = f"X{x_group_id} - {x_group_name}"
        x_names_dict[str(x_group_id)].append(x_group_name)
    if percent_threshold > 0:
        x_names_dict[f'<{percent_threshold}'] = ['']
    return x_names_dict


def get_counts_tuples_hmm_vs_ecod(
    ecod_x_groups_distribution_dictionary: dict,
    x_group_id_dictionary: dict
):
    x_groups_ids = sorted(
        list(set(ecod_x_groups_distribution_dictionary.keys()).union(set(x_group_id_dictionary.keys())))
    )
    x_groups_coords_list = []
    for x_group_id in x_groups_ids:
        if x_group_id in x_group_id_dictionary:
            hmm_x_coord = x_group_id_dictionary[x_group_id]
        else:
            continue
        if x_group_id in ecod_x_groups_distribution_dictionary:
            ecod_y_coord = ecod_x_groups_distribution_dictionary[x_group_id]
        else:
            continue
        x_groups_coords_list.append([hmm_x_coord, ecod_y_coord, x_group_id])
    return x_groups_coords_list


def plot_counts_comparison(
        x_groups_coords_list: list,
        fig_file_path: str = 'ecod-hmm-counts-x-groups.png',
):
    plt.figure(figsize=(10, 6))
    x_values = [x for x, y, _ in x_groups_coords_list]
    y_values = [y for x, y, _ in x_groups_coords_list]

    corr_coefficient = np.corrcoef(x_values, y_values)[0, 1]
    plt.scatter(x_values, y_values, color='blue', marker='x', label="X group")

    # Setting x and y axes to log scale
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    plt.xlabel('Hmm - Human Domains Counts')
    plt.ylabel('ECOD - Human Raw Domain Counts')
    plt.title('')
    text_x = max(x_values) * 0.6
    text_y = min(y_values) * 1.1
    plt.text(
        text_x,
        text_y,
        f'Pearson r = {corr_coefficient:.2f}',
        fontsize=12,
        ha='right'
    )
    plt.title(
        "Comparison of Human Protein Domain Counts: HMM search versus ECOD"
    )
    plt.savefig(fig_file_path)
