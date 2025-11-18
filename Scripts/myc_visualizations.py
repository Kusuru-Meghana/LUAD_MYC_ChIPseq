import os
import pandas as pd
import matplotlib.pyplot as plt

# 1) Load annotated MYC peaks
peaks = pd.read_csv("MYC_peaks_annotated.bed", sep="\t", header=None)

# Rename only the columns we care about.
# Based on your bedtools output, these indices are:
# 0: chr
# 1: start
# 2: end
# 3: peak ID (MYC_combined_peak_1, etc.)
# 6: fold enrichment (from MACS2)
# 13: gene symbol (HES4, SOX2, etc.)
# 16: distance to TSS (bp, negative = upstream)
peaks = peaks.rename(
    columns={
        0: "chrom",
        1: "start",
        2: "end",
        3: "peak_id",
        6: "fold_enrichment",
        13: "gene",
        16: "dist_to_tss"
    }
)

print(f"Loaded {len(peaks)} MYC binding peaks")
print(f"Unique genes with MYC peaks: {peaks['gene'].nunique()}")

os.makedirs("myc_plots", exist_ok=True)

# 2) Genome-wide MYC binding strength
plt.figure(figsize=(6, 4))
peaks["fold_enrichment"].hist(bins=60)
plt.xlabel("MYC ChIP-seq fold enrichment")
plt.ylabel("Number of peaks")
plt.title("Genome-wide MYC binding strength (LUAD MYC ChIP)")
plt.tight_layout()
plt.savefig("myc_plots/MYC_genome_wide_peaks.png", dpi=150)
plt.close()
print("Saved: myc_plots/MYC_genome_wide_peaks.png")

# 3) Helper to plot MYC binding around a specific gene
def plot_gene(gene_name: str):
    sub = peaks[peaks["gene"] == gene_name].copy()
    if sub.empty:
        print(f"{gene_name} not found in annotated peaks.")
        return

    # Sort peaks by distance to TSS
    sub = sub.sort_values("dist_to_tss")

    plt.figure(figsize=(6, 4))
    plt.scatter(sub["dist_to_tss"], sub["fold_enrichment"])
    plt.axvline(0, color="black", linestyle="--", linewidth=1)
    plt.xlabel("Distance to TSS (bp; negative = upstream)")
    plt.ylabel("Fold enrichment")
    plt.title(f"MYC binding around {gene_name} TSS")
    plt.tight_layout()

    out_path = f"myc_plots/MYC_{gene_name}_binding.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")

# 4) Make gene-specific plots
genes_of_interest = ["HES4", "SOX2", "MYC", "FOSL1", "EGFR", "TP53"]

for g in genes_of_interest:
    plot_gene(g)

print("DONE.")

