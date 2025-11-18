import gseapy as gp

# -------------------------------------------------------------------
# Load gene list
# -------------------------------------------------------------------
with open("MYC_genes.txt") as f:
    genes = [g.strip() for g in f if g.strip()]

print(f"Loaded {len(genes)} gene names")

if len(genes) == 0:
    raise ValueError("ERROR: Gene list is empty. Check MYC_genes.txt.")

# -------------------------------------------------------------------
# Run Enrichr with updated, valid gene sets
# -------------------------------------------------------------------
print("Running enrichment analysis...")

go_res = gp.enrichr(
    gene_list=genes,
    gene_sets=[
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'KEGG_2021_Human'
    ],
    organism='Human',
    outdir='enrichment_results',
    cutoff=0.05  # FDR threshold
)

print("\nEnrichment Completed!")
print("Results saved to: enrichment_results/")

