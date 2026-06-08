import argparse


def hmm_search_argument_parser():
    parser = argparse.ArgumentParser(
        description="Pfam-based domain analysis (updated from ECODf)"
    )
    parser.add_argument("--ecod-latest-domains-txt-path", type=str, required=True,
                        help="Path to ecod.latest.domains.txt")
    parser.add_argument("--ecod-latest-domains-pfam-map-path", type=str, required=True,
                        help="Path to ecod.latest.f_id_pfam_acc.txt")
    parser.add_argument("--protein-db-path", type=str, required=True,
                        help="Path to the protein database (will be created if absent)")
    parser.add_argument("--ecod-db-path", type=str, required=True,
                        help="Path to the ECOD results database (will be created)")
    parser.add_argument("--exonize-db-path", type=str, required=True,
                        help="Path to the exonize results database")
    parser.add_argument("--genome-file-path", type=str,
                        help="Path to genome FASTA (.fa.gz)")
    parser.add_argument("--genome-hierarchy-pckl-path", type=str,
                        help="Path to gene hierarchy pickle file")
    parser.add_argument("--fasta-file-path-cdss", type=str, required=True,
                        help="Path to output FASTA file for CDS sequences")
    parser.add_argument("--fasta-file-path-mrna", type=str, required=True,
                        help="Path to output FASTA file for mRNA sequences")
    parser.add_argument("--hmm-profile-path", type=str, required=True,
                        help="Path to Pfam-A.hmm (pressed with hmmpress)")
    parser.add_argument("--hmmsearch-cdss-path", type=str, required=True,
                        help="Path to hmmscan output for CDS (created if absent)")
    parser.add_argument("--hmmsearch-mrna-path", type=str, required=True,
                        help="Path to hmmscan output for mRNA (created if absent)")
    parser.add_argument("--timeout-db", type=int, default=160,
                        help="SQLite connection timeout (seconds)")
    parser.add_argument("--evalue-threshold", type=float, default=0.01,
                        help="E-value threshold for hmmscan")
    parser.add_argument("--query-coverage-threshold", type=float, default=0.9,
                        help="Minimum HMM profile coverage (query)")
    parser.add_argument("--target-coverage-threshold", type=float, default=0.9,
                        help="Minimum sequence coverage (target, CDS search only)")
    parser.add_argument("--overlapping-threshold", type=float, default=0.9,
                        help="Overlap threshold for deduplication and classification")
    parser.add_argument("--hard-force", action="store_true", default=False,
                        help="Remove and recreate all files from scratch")
    parser.add_argument("--soft-force", action="store_true", default=False,
                        help="Rerun hmmscan and rebuild ECOD database only")
    parser.add_argument("--out-filename", type=str, default="odds_ratio.txt",
                        help="Path to odds ratio output file")
    return parser.parse_args()


def x_group_analysis_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ecod-db-path", type=str, required=True)
    parser.add_argument("--ecod-latest-domains-txt-path", type=str, required=True)
    parser.add_argument("--cluster-xgroups-pie-chart-threshold", type=float, default=0.01)
    parser.add_argument("--output-directory", type=str, required=True)
    parser.add_argument("--overlap-threshold", type=float, default=0.9)
    return parser.parse_args()
