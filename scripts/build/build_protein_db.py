import argparse
from pathlib import Path
from lib.prot_utils import populate_cds_and_mrna_tables, read_pkl_file, read_genome
from lib.protein_db_sql_utils import create_mrna_transcripts_cds_tables


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--gene-hierarchy-dictionary-path',
        type=str,
        help='Path to gene hierarchy dictionary.',
    )
    parser.add_argument(
        '--genome-path',
        type=str,
        help='Path to genome.',
    )
    parser.add_argument(
        '--protein-db-path',
        type=str,
        help='Path to protein database.',
    )
    parser.add_argument(
        '--timeout-db',
        type=int,
        help='Timeout for database.',
    )
    return parser.parse_args()



def main():
    args = args_parser()
    gene_hierarchy_dictionary_path = Path(args.gene_hierarchy_dictionary_path)
    genome_path = Path(args.genome_path)
    protein_db_path = Path(args.protein_db_path)
    timeout_db = args.timeout_db
    gene_hierarchy_dictionary = read_pkl_file(
        file_path=Path(gene_hierarchy_dictionary_path)
    )
    genome_dictionary = read_genome(
        genome_file_path=genome_path
    )
    create_mrna_transcripts_cds_tables(
        db_path=protein_db_path
    )
    populate_cds_and_mrna_tables(
        gene_hierarchy_dictionary=gene_hierarchy_dictionary,
        genome_dictionary=genome_dictionary,
        protein_db_path=protein_db_path,
        timeout_db=timeout_db
    )


if __name__ == '__main__':
    main()
