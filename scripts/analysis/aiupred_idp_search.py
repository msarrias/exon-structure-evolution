import logging
import os

from lib.aiupred_idp_utils import (
    argument_parser,
    even_batches,
    process_gene_transcripts_and_update_database,
    update_scores_cds_mapping,
    get_representative_cds_tuples,
    read_pkl_file,
)
from lib.aiupred_idp_sqlite_utils import (
    check_if_table_exists,
    check_if_empty_table,
    create_aiupred_table,
    create_representative_cds_human_table,
    create_index_on_table,
    copy_table_between_data_bases,
    update_cdss_table,
    update_mapping_IUPRED2_peptides_columns,
    insert_representative_cds,
    query_proteins_table,
    query_processed_transcripts,
    query_gene_cdss,
    query_protein_scores_table,
    query_raw_cds_scores,
)

# AIUPred scores are probabilities of amino acids being in
# an Intrinsically Disordered Region (IDR) of a protein.
# close to 1, likely disordered
# close to 0, likely structured
# ANCHOR method: predicts disordered regions that become
# structured upon binding to other proteins.

FORKS_NUMBER = 5
BATCH_SIZE   = 100


def main():
    args = argument_parser()
    protein_db_path = args.protein_db_path
    exonize_db_path = args.exonize_db_path
    gene_hierarchy_path = args.gene_hierarchy_path
    output_db_path = args.output_db_path
    anchor = args.anchor
    min_exon_length = args.min_exon_length
    ordered_score_threshold = args.ordered_score_threshold

    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s]: %(message)s",
        handlers=[logging.StreamHandler()]
    )
    logging.info("Checking file requirements and database integrity.")

    if not os.path.exists(protein_db_path):
        raise FileNotFoundError(f"Protein database not found: {protein_db_path}")

    if not os.path.exists(output_db_path):
        create_aiupred_table(db_path=output_db_path, anchor=anchor)

    # # --- Copy Expansions_full table from exonize db ---
    # if not check_if_table_exists(db_path=output_db_path, table_name="Expansions_full"):
    #     logging.info("Copying 'Expansions_full' table from exonize database.")
    #     copy_table_between_data_bases(
    #         target_db_path=output_db_path,
    #         source_db_path=exonize_db_path,
    #         table_name="Expansions_full",
    #     )

    if not check_if_table_exists(
            db_path=output_db_path,
            table_name="Gene_CDS_Proteins"
    ):
        logging.info("Copying 'Gene_CDS_Proteins' table from protein database.")
        copy_table_between_data_bases(
            target_db_path=output_db_path,
            source_db_path=protein_db_path,
            table_name="Gene_CDS_Proteins",
        )
        create_index_on_table(
            db_path=output_db_path,
            table_name="Gene_CDS_Proteins",
            column_name="CdsID",
        )
        update_cdss_table(
            db_path=output_db_path,
            anchor=anchor,
            cds_table_name="Gene_CDS_Proteins",
        )
    # Run AIUPred predictions
    gene_transcripts_dictionary = query_proteins_table(
        db_path=protein_db_path
    )
    processed_transcripts = query_processed_transcripts(
        db_path=output_db_path
    )
    transcript_ids_to_process = list(
        set(gene_transcripts_dictionary.keys()) - set(processed_transcripts)
    )

    if transcript_ids_to_process:
        logging.info(
            f"Starting AIUPred predictions for "
            f"{len(transcript_ids_to_process):,} transcripts"
        )
        forks = 0
        for batch in even_batches(
                data=transcript_ids_to_process,
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
                        list_transcripts_batch=batch,
                        gene_transcripts_dictionary=gene_transcripts_dictionary,
                        anchor=anchor,
                        output_db_path=output_db_path,
                    )
                except Exception as e:
                    print(e)
                    status = os.EX_SOFTWARE
                finally:
                    os._exit(status)

    # Map AIUPred scores back to CDS coordinates
    logging.info("Mapping AIUPred scores to CDS coordinates")
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
        ordered_score_threshold=ordered_score_threshold,
    )
    if tuples_to_insert:
        for batch in even_batches(
                data=tuples_to_insert,
                number_of_batches=BATCH_SIZE
        ):
            update_mapping_IUPRED2_peptides_columns(
                db_path=output_db_path,
                anchor=anchor,
                tuples_list=batch
            )

    # Build representative CDS table
    if not check_if_table_exists(
            db_path=output_db_path,
            table_name="Representative_CDSs"
    ):
        logging.info(
            "Creating 'Representative_CDSs' table"
        )
        create_representative_cds_human_table(
            db_path=output_db_path,
            anchor=anchor
        )
    if check_if_empty_table(
            db_path=output_db_path,
            table_name="Representative_CDSs"
    ):
        logging.info("Populating 'Representative_CDSs' table")
        gene_hierarchy_dictionary = read_pkl_file(file_path=gene_hierarchy_path)
        raw_scores_dict = query_raw_cds_scores(path_db=output_db_path)
        tuples_to_insert = get_representative_cds_tuples(
            gene_hierarchy_dictionary=gene_hierarchy_dictionary,
            raw_scores_dict=raw_scores_dict,
            min_exon_length=min_exon_length,
        )
        for batch in even_batches(
                data=tuples_to_insert,
                number_of_batches=BATCH_SIZE
        ):
            insert_representative_cds(
                db_path=output_db_path,
                anchor=anchor,
                tuples_list=batch
            )
    logging.info("Done!")


if __name__ == "__main__":
    main()
