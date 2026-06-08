from lib.iupred3_idp_utils import *
import logging
# IUpred scores are probabilities of amino acids being in
# an Intrinsically Disordered Region (IDR) of a protein.
# scores are calculated for each amino acid residue.
# close to 1, likely to be disordered
# close to 0, likely to be structured
# IUPred3 choices: "long", "short" or "glob"
# ANCHOR method: predicts disordered regions that become
# structured upon binding to other proteins.


def main():
    args = argument_parser()
    protein_db_path = args.protein_db_path
    anchor = args.anchor
    smoothing = args.smoothing
    iupred_type = args.iupred_type
    iupred3_path = args.iupred3_path
    output_db_path = args.output_db_path
    exonize_db_path = args.exonize_db_path
    gene_hierarchy_path = args.gene_hierarchy_path
    min_exon_length = args.min_exon_length
    ordered_score_threshold = args.ordered_score_threshold
    BATCH_SIZE = 100

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(
        logging.Formatter('[%(levelname)s]: %(message)s')
    )
    logging.basicConfig(
        level=logging.INFO, handlers=[console_handler]
    )
    logging.info(
        "Checking file requirements and database integrity."
    )
    if not os.path.exists(protein_db_path):
        raise FileNotFoundError(
            f'Protein database not found at:'
            f'{protein_db_path}'
        )
    if not os.path.exists(output_db_path):
        create_iupred3_table(
            db_path=output_db_path,
            anchor=anchor
        )
    # COPY TABLES FROM THE PROTEIN DB
    if not check_if_table_exists(
            db_path=output_db_path,
            table_name='Expansions'
    ):
        logging.info(
            "Copying the 'Expansions' table from the Exonize database."
        )
        copy_table_between_data_bases(
            target_db_path=output_db_path,
            source_db_path=exonize_db_path,
            table_name='Expansions'
        )
    if not check_if_table_exists(
            db_path=output_db_path,
            table_name='Gene_CDS_Proteins'
    ):
        logging.info(
            "Copying the 'Gene_CDS_Proteins' table from the Protein database."
        )
        copy_table_between_data_bases(
            target_db_path=output_db_path,
            source_db_path=protein_db_path,
            table_name='Gene_CDS_Proteins'
        )
        create_index_on_table(
            db_path=output_db_path,
            table_name='Gene_CDS_Proteins',
            column_name="CdsID"
        )
        update_cdss_table(
            db_path=output_db_path,
            anchor=anchor,
            cds_table_name='Gene_CDS_Proteins'
        )
    # POPULATE THE iupred3_protein_sequences_scores TABLE IF EMPTY
    gene_transcripts_dictionary = query_proteins_table(
        db_path=protein_db_path
    )
    processed_transcripts_list = query_processed_transcripts(
        db_path=output_db_path
    )
    gene_transcripts_ids = list(gene_transcripts_dictionary.keys())
    transcript_ids_to_process = list(set(gene_transcripts_ids) - set(processed_transcripts_list))

    transactions_pks: set[int]
    status: int
    code: int
    forks: int = 0
    FORKS_NUMBER = 5 #os.cpu_count()
    if transcript_ids_to_process:
        logging.info(
            f"Starting the IUPred3 predictions for {len(transcript_ids_to_process)} "
            f"transcripts. This may take a while..."
        )
        for balanced_genes_batch in even_batches(
                data=list(transcript_ids_to_process),
                number_of_batches=FORKS_NUMBER,
        ):
            if os.fork():
                forks += 1
                if forks >= FORKS_NUMBER:
                    _, status = os.wait()
                    code = os.waitstatus_to_exitcode(status)
                    assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                    assert code != os.EX_SOFTWARE
                    forks -= 1
            else:
                status = os.EX_OK
                try:
                    process_gene_transcripts_and_update_database(
                        list_transcripts_batch=balanced_genes_batch,
                        gene_transcripts_dictionary=gene_transcripts_dictionary,
                        iupred3_path=iupred3_path,
                        iupred_type=iupred_type,
                        anchor=anchor,
                        smoothing=smoothing,
                        output_db_path=output_db_path
                    )
                except Exception as exception:
                    print(exception)
                    status = os.EX_SOFTWARE
                finally:
                    # This prevents the child process forked above to keep along the for loop upon completion
                    # of the try/except block. If this was not present, it would resume were it left off, and
                    # fork in turn its own children, duplicating the work done, and creating a huge mess.
                    # We do not want that, so we gracefully exit the process when it is done.
                    os._exit(status)

    gene_cdss_dictionary = query_gene_cdss(
        db_path=output_db_path
    )
    protein_scores_dictionary = query_protein_scores_table(
        db_path=output_db_path,
        anchor=anchor
    )
    tuples_to_insert = update_scores_cds_mapping(
        gene_cdss_dictionary=gene_cdss_dictionary,
        protein_scores_dictionary=protein_scores_dictionary,
        anchor=anchor,
        ordered_score_threshold=ordered_score_threshold
    )
    logging.info(
        f"Mapping IUPRED2 scores to CDSs"
    )
    if tuples_to_insert:
        for batch in even_batches(
            data=tuples_to_insert,
            number_of_batches=BATCH_SIZE,
        ):
            update_mapping_IUPRED2_peptides_columns(
                db_path=output_db_path,
                anchor=anchor,
                tuples_list=batch
            )

    if not check_if_table_exists(
            db_path=output_db_path,
            table_name='Representative_CDSs'
    ):
        logging.info(
            f"Fetching representative CDSs"
        )
        create_representative_cds_human_table(
            db_path=output_db_path,
            anchor=anchor
        )
    if not check_if_empty_table(
            db_path=output_db_path,
            table_name='Representative_CDSs'
    ):
        logging.info(
            f"Populating the 'Representative_CDSs' table"
        )
        gene_hierarchy_dictionary = read_pkl_file(
            file_path=gene_hierarchy_path
        )
        raw_scores_dict = query_raw_cds_scores(
            path_db=iupred_db_path
        )
        tuples_to_insert = get_representative_cds_tuples(
            gene_hierarchy_dictionary=gene_hierarchy_dictionary,
            raw_scores_dict=raw_scores_dict,
            min_exon_length=min_exon_length
        )
        for batch in even_batches(
                data=tuples_to_insert,
                number_of_batches=BATCH_SIZE
        ):
            insert_representative_cds(
                db_path=iupred_db_path,
                anchor=anchor,
                tuples_list=batch
            )
    logging.info("Done!")


if __name__ == '__main__':
    main()
