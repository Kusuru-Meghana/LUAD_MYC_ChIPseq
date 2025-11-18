#!/usr/bin/env bash
# ============================================================
# EXACT COMMANDS MEGHANA RAN FOR LUAD MYC ChIP-seq PROJECT
# ============================================================

# ------------------------------------------------------------
# 0. Set working directory
# ------------------------------------------------------------
cd /home/meghana-20/projects

mkdir LUAD_MYC_ChIP
cd LUAD_MYC_ChIP

# ------------------------------------------------------------
# 1. FASTQ FILES 
# ------------------------------------------------------------
# Files:
# SRR568131.fastq.gz   (MYC rep1)
# SRR568132.fastq.gz   (MYC rep2)
# SRR568133.fastq.gz   (IgG rep1)
# SRR568134.fastq.gz   (IgG rep2)

# Renamed them:
# MYC_rep1.fastq.gz
# MYC_rep2.fastq.gz
# IgG_rep1.fastq.gz
# IgG_rep2.fastq.gz

# ------------------------------------------------------------
# 2. QC 
# ------------------------------------------------------------

# Installed fastqc and multiqc manually.
# Then created the qc folder and ran:

mkdir qc
fastqc *.fastq.gz -o qc

multiqc qc -o qc/multiqc

# ------------------------------------------------------------
# 3. BOWTIE2 
# ------------------------------------------------------------

# downloaded the hg38 FASTA:
# hg38.fa.gz -> unzipped -> hg38.fa

# Generated Bowtie2 index 
bowtie2-build hg38.fa hg38

mkdir align

# Alignment commands 

bowtie2 -x hg38 -U MYC_rep1.fastq.gz -S align/MYC_rep1.sam
bowtie2 -x hg38 -U MYC_rep2.fastq.gz -S align/MYC_rep2.sam
bowtie2 -x hg38 -U IgG_rep1.fastq.gz -S align/IgG_rep1.sam
bowtie2 -x hg38 -U IgG_rep2.fastq.gz -S align/IgG_rep2.sam

# ------------------------------------------------------------
# 4. SAMTOOLS 
# ------------------------------------------------------------

cd align

# Convert SAM → BAM
samtools view -bS MYC_rep1.sam > MYC_rep1.bam
samtools view -bS MYC_rep2.sam > MYC_rep2.bam
samtools view -bS IgG_rep1.sam > IgG_rep1.bam
samtools view -bS IgG_rep2.sam > IgG_rep2.bam

# Sort BAM
samtools sort -o MYC_rep1.sorted.bam MYC_rep1.bam
samtools sort -o MYC_rep2.sorted.bam MYC_rep2.bam
samtools sort -o IgG_rep1.sorted.bam IgG_rep1.bam
samtools sort -o IgG_rep2.sorted.bam IgG_rep2.bam

# Index BAM
samtools index MYC_rep1.sorted.bam
samtools index MYC_rep2.sorted.bam
samtools index IgG_rep1.sorted.bam
samtools index IgG_rep2.sorted.bam

# ------------------------------------------------------------
# 5. MACS2
# ------------------------------------------------------------

cd ..

mkdir peaks

macs2 callpeak \
  -t align/MYC_rep1.sorted.bam align/MYC_rep2.sorted.bam \
  -c align/IgG_rep1.sorted.bam align/IgG_rep2.sorted.bam \
  -g hs \
  -f BAM \
  -n MYC_combined \
  --outdir peaks \
  --call-summits

# Output:
# MYC_combined_peaks.narrowPeak
# MYC_combined_summits.bed
# MYC_combined_peaks.xls

# ------------------------------------------------------------
# 6. CHIPSEEKER 
# ------------------------------------------------------------

# R script 
# Rscript annotation_script.R


cat <<EOF > annotation_script.R
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

peak <- readPeakFile("peaks/MYC_combined_peaks.narrowPeak")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(peak, TxDb=txdb, annoDb="org.Hs.eg.db")
anno <- as.data.frame(peakAnno)
write.table(anno, "results/MYC_peaks_annotated.bed",
            sep="\t", quote=FALSE, row.names=FALSE)

genes <- unique(anno$SYMBOL)
writeLines(genes, "results/MYC_genes.txt")
EOF

Rscript annotation_script.R

# ------------------------------------------------------------
# 7. GSEAPY
# ------------------------------------------------------------

cat <<EOF > enrichment.py
import gseapy as gp

with open("results/MYC_genes.txt") as f:
    genes = [x.strip() for x in f]

gp.enrichr(gene_list=genes,
           gene_sets=["GO_Biological_Process_2023",
                      "GO_Molecular_Function_2023",
                      "GO_Cellular_Component_2023",
                      "KEGG_2021_Human"],
           outdir="results/enrichment_results",
           cutoff=0.5)
EOF

python3 enrichment.py

# ------------------------------------------------------------
# 8. PLOT 
# ------------------------------------------------------------

cat <<EOF > plot_peaks.py
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table("peaks/MYC_combined_peaks.xls", comment="#")
plt.hist(df["fold_enrichment"], bins=40)
plt.savefig("results/figures/MYC_genome_wide_peaks.png", dpi=300)
EOF

python3 plot_peaks.py

# ------------------------------------------------------------
# 9. HOMER 
# ------------------------------------------------------------

findMotifsGenome.pl peaks/MYC_combined_summits.bed hg38 results/homer -size 200 -mask



