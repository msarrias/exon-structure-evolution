import pickle
import pandas as pd
from pathlib import Path

PROJECT_PATH = Path("/home/mahe9165/exon-structure-evolution/data/")
PEXT_PATH = PROJECT_PATH / Path( "pext/gnomad.pext.gtex_v10.base_level.tsv.gz")
OUTPUT_PATH = PROJECT_PATH / Path("pext/GRCh38_v115_human_pext_scores.pkl")
HIERARCHY_PATH = PROJECT_PATH / Path("ensembl/v115/GRCh38.115_gene_hierarchy.pkl")

# --- Load pext data ---
# Note: some base scores are missing, others are set to NaN:
#   - NaN entries: genes not expressed in GTEx tissues
#   - None entries: genes with transcripts/annotations for which no expression
#     data is available (e.g. different annotations)
df = pd.read_csv(PEXT_PATH, sep="\t", compression="gzip")
df[["chrom", "pos"]] = df["locus"].str.split(":", expand=True)
df["pos"] = df["pos"].astype(int)
df_grouped = {
    gene_id: sub_df
    for gene_id, sub_df in df.groupby("gene_id")
}

# --- Compute mean pext score per CDS ---
pext_gene_cds_dict = {}
missing_genes = set()

with open(HIERARCHY_PATH, 'rb') as handle:
    gene_hierarchy_dictionary = pickle.load(handle)

for gene, gene_dict in gene_hierarchy_dictionary.items():
    _, gene_id = gene.rsplit(":")
    if gene_id not in df_grouped:
        missing_genes.add(gene_id)
        continue

    set_cds = {
        cds["coordinate"]
        for transcript_dict in gene_dict["mRNAs"].values()
        for cds in transcript_dict["structure"]
    }

    gene_df = df_grouped[gene_id]
    pext_gene_cds_dict[gene] = {}

    for cds in set_cds:
        mask = (gene_df["pos"] >= cds.lower) & (gene_df["pos"] <= cds.upper)
        base_pext_scores = gene_df.loc[mask, "exp_prop_mean"]
        pext_gene_cds_dict[gene][cds] = (
            base_pext_scores.mean() if not base_pext_scores.empty else None
        )

print(f"Processed {len(pext_gene_cds_dict)} genes ({len(missing_genes)} missing from pext data)")

with open(OUTPUT_PATH, "wb") as handle:
    pickle.dump(pext_gene_cds_dict, handle)
