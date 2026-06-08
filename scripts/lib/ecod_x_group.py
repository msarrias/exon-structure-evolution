import os.path
from datetime import date
from prot_sql_utils import (query_hmm_searches_matches)
from collections import defaultdict
from pathlib import Path
from ecod_x_group_utils import (
    dump_query_into_csv_file,
    get_x_groups_distribution,
    get_x_groups_dictionary_count,
    get_hmm_matches_and_events_distributions,
    get_counts_tuples_hmm_vs_ecod,
    plot_side_by_side_piechart,
    dump_venn_diagram_csv_input,
    get_hmm_mrna_cds_search_matches_dictionary
)
from prot_sql_utils import (
    query_xgroup_counts_ecod,
    query_ecod_subdomain_domain_counts,
)
from ecod_x_group_sql_utils import query_ecod_full_domain_duplication_records
import logging
from ecod_hmm_search_parser import x_group_analysis_argument_parser


def main():
    args = x_group_analysis_argument_parser()
    ecod_db_path = args.ecod_db_path
    output_directory = args.output_directory
    cluster_xgroups_pie_chart_threshold = args.cluster_xgroups_pie_chart_threshold
    overlap_threshold = args.overlap_threshold
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(
        logging.Formatter('[%(levelname)s]: %(message)s')
    )
    logging.basicConfig(
        level=logging.INFO, handlers=[console_handler]
    )

    output_directory_name = str(date.today()).replace('-', '_')
    output_directory_path = os.path.join(output_directory, f'{output_directory_name}')
    os.makedirs(os.path.join(output_directory_path), exist_ok=True)
    logging.info(
        f"Output files will be saved in: {output_directory_path}"
    )
    logging.info(
        "The following items will be created: "
    )
    mrna_domain_matches_list, cds_domain_matches_list = query_hmm_searches_matches(
        db_path=ecod_db_path,
    )
    gene_matches_dictionary = get_hmm_mrna_cds_search_matches_dictionary(
        hmm_matches_mrna_list=mrna_domain_matches_list,
        hmm_matches_cds_list=cds_domain_matches_list,
        overlap_threshold=overlap_threshold
    )
    # 1. Get Venn diagram - Relationship Between Gene Sets Featuring Full Domain,
    # Sub Domain Exons, and Exon Duplication - the diagram is generated in Origin 2023b
    logging.info(
        "- Venn diagram csv file - Relationship Between Gene Sets"
    )
    (genes_with_exon_duplications,
     genes_with_full_domain_exon_duplication_mrna_cds_hmm_search,
     genes_with_sub_domain_exon_duplication_mrna_hmm_search,
     genes_with_full_domain_mrna_hmm,
     genes_with_full_domain_cds_hmm,
     genes_with_sub_domain_mrna_hmm,
     genes_within_domain_boundary_hmm
     ) = query_ecod_subdomain_domain_counts(
        db_path=ecod_db_path,
        hmm_mrna_table_name='x_group_filtered_ECOD_hmm_mRNA',
        hmm_cds_table_name='x_group_filtered_ECOD_hmm_CDSs',
        expansions_table_name='Expansions_full',
    )
    dump_venn_diagram_csv_input(
        genes_with_full_domain_duplication=list(genes_with_full_domain_exon_duplication_mrna_cds_hmm_search),
        genes_with_sub_domain_duplication=list(genes_with_sub_domain_exon_duplication_mrna_hmm_search),
        genes_with_exon_duplication=list(genes_with_exon_duplications),
        genes_with_full_domain_hmm=list(genes_with_full_domain_mrna_hmm.union(genes_with_full_domain_cds_hmm)),
        genes_with_sub_domain_hmm=list(genes_with_sub_domain_mrna_hmm),
        output_csv_path=os.path.join(output_directory_path, 'venn_diagram.csv')
    )
    logging.info(
        "- Pie chart - Global Distribution of x-groups based on ECOD distribution and HMM search matches"
    )
    #  2. Get global distribution of x-groups in ECOD and HMM
    x_group_id_dictionary = defaultdict(int)
    for gene_id, x_group_dictionary in gene_matches_dictionary.items():
        for x_group_id, x_group_matches_coordinates_list in x_group_dictionary.items():
            if x_group_id != 0:
                x_group_id_dictionary[x_group_id] += len(x_group_matches_coordinates_list)
    ecod_x_groups_distribution_dictionary = query_xgroup_counts_ecod(
        db_path=ecod_db_path,
    )
    x_groups_coords_list = get_counts_tuples_hmm_vs_ecod(
        x_group_id_dictionary=x_group_id_dictionary,
        ecod_x_groups_distribution_dictionary=ecod_x_groups_distribution_dictionary
    )
    # 2.1. Dump the distribution of x-groups in a csv file for plotting in Origin 2023b
    dump_query_into_csv_file(
        output_file=Path(output_directory_path, 'distributions_scatter_plot.csv'),
        column_names_list=['Hmm - Human Domains Counts',
                           'ECOD - Human Raw Domain Counts',
                           'X group'],
        rows_list=x_groups_coords_list
    )
    # for plotting using matplotlib
    # plot_counts_comparison(
    #     hmm_x_groups_distribution=x_group_id_dictionary,
    #     ecod_x_groups_distribution=ecod_x_groups_distribution_dictionary,
    #     fig_path=os.path.join(output_directory_path, 'x_groups_counts_comparison.png')
    # )
    hmm_x_groups_distribution, ecod_x_groups_distribution = get_hmm_matches_and_events_distributions(
        x_group_id_dictionary=x_group_id_dictionary,
        ecod_x_groups_distribution_dictionary=ecod_x_groups_distribution_dictionary,
        clustering_percentage_threshold=cluster_xgroups_pie_chart_threshold
    )
    plot_side_by_side_piechart(
        pie_chart_hmm_data=hmm_x_groups_distribution,
        pie_chart_ecod_data=ecod_x_groups_distribution,
        x_names_dictionary=get_x_groups_dictionary_count(
            db_path=ecod_db_path,
            percent_threshold=cluster_xgroups_pie_chart_threshold
        ),
        fig_path=os.path.join(output_directory_path, 'pie_chart_x_groups_distribution.png')
    )
    logging.info(
        "- X groups csv file distribution in genes with full domain exon duplications"
    )
    #  3. Get distribution of x-groups in genes with full domain exon duplications
    # the figure is generated in Origin 2023b

    full_domain_exon_duplication_records_dictionary = query_ecod_full_domain_duplication_records(
        db_path=ecod_db_path
    )
    column_list, new_records = get_x_groups_distribution(
        gene_matches_dictionary=gene_matches_dictionary,
        full_domain_exon_duplication_records=full_domain_exon_duplication_records_dictionary,
        db_path=ecod_db_path,
        overlap_threshold=overlap_threshold,
        clustering_percentage_threshold=cluster_xgroups_pie_chart_threshold
    )
    dump_query_into_csv_file(
        output_file=Path(output_directory_path, 'xgroups_counts.csv'),
        column_names_list=column_list,
        rows_list=new_records,
    )


if __name__ == '__main__':
    main()
# ecod_db_path = "/Users/marinaherrerasarrias/Desktop/PhD/01-Projects/exonize_private/ECOD_hmm_analysis/data/protein_db_human_unmasked_genome.db"
# ecod_latest_domains_txt_path = "/Users/marinaherrerasarrias/Desktop/PhD/01-Projects/exonize_private/ECOD_hmm_analysis/ecod.latest.domains.txt"
# timeout_db = 60