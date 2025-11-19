#!/bin/bash
set -e

### -------------------------
### 1. Create directories
### -------------------------
mkdir -p data fastqc multiqc alignment peaks annotation enrichment motif

### -------------------------
### 2. Download FASTQ files
### -------------------------
cd data

# MYC replicates
prefetch SRR568131
prefetch SRR568132

# IgG controls
prefetch SRR568133
prefetch SRR568134

# Convert to FASTQ
fasterq-dump SRR568131 -O .
fasterq-dump SRR568132 -O .
fasterq-dump SRR568133 -O .
fasterq-dump SRR568134 -O .

# Rename
mv SRR568131.fastq MYC_rep1.fastq
mv SRR568132.fastq MYC_rep2.fastq
mv SRR568133.fastq IgG_rep1.fastq
mv SRR568134.fastq IgG_rep2.fastq

cd ..

### -------------------------
### 3. FastQC + MultiQC
### -------------------------
fastqc data/*.fastq -o fastqc/
multiqc fastqc/ -o multiqc/

### -------------------------
### 4. Alignment (Bowtie2)
### -------------------------
# (Assumes hg38 index is present)
bowtie2 -x hg38 -U data/MYC_rep1.fastq -S alignment/MYC_rep1.sam
bowtie2 -x hg38 -U data/MYC_rep2.fastq -S alignment/MYC_rep2.sam
bowtie2 -x hg38 -U data/IgG_rep1.fastq -S alignment/IgG_rep1.sam
bowtie2 -x hg38 -U data/IgG_rep2.fastq -S alignment/IgG_rep2.sam

### -------------------------
### 5. Convert SAM → BAM, sort & index
### -------------------------
for sample in MYC_rep1 MYC_rep2 IgG_rep1 IgG_rep2
do
    samtools view -S -b alignment/${sample}.sam > alignment/${sample}.bam
    samtools sort alignment/${sample}.bam -o alignment/${sample}_sorted.bam
    samtools index alignment/${sample}_sorted.bam
done

### -------------------------
### 6. Peak calling (MACS2)
### -------------------------
macs2 callpeak -t alignment/MYC_rep1_sorted.bam alignment/MYC_rep2_sorted.bam \
               -c alignment/IgG_rep1_sorted.bam alignment/IgG_rep2_sorted.bam \
               -f BAM -g hs -n MYC --outdir peaks/

### -------------------------
### 7. Annotation (ChIPseeker in R)
### -------------------------
Rscript scripts/annotate_peaks.R

### -------------------------
### 8. GO/KEGG enrichment
### -------------------------
python3 scripts/enrichment_py.py

### -------------------------
### 9. Motif analysis (HOMER)
### -------------------------
findMotifsGenome.pl peaks/MYC_peaks.narrowPeak hg38 motif/ -size given

echo "Pipeline complete!"
