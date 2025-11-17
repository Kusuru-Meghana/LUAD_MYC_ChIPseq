## LUAD_MYC_ChIPseq
This repository contains a complete end-to-end ChIP-seq analysis of MYC transcription factor binding in Lung Adenocarcinoma (LUAD).
It demonstrates a full NGS pipeline from FASTQ → alignment → peak calling → annotation → motif analysis → functional interpretation.

# Repository Structure
LUAD_MYC_ChIPseq/
│── fastq/                     # Raw FASTQ files (MYC + IgG)
│── qc/                        # FASTQC and MultiQC reports
│── alignment/                 # Bowtie2 alignments (.bam)
│── macs2_results/             # MACS2 peak calling output
│── peaks/                     # Annotated peaks, MYC_genes.txt
│── motif_results/             # HOMER de novo + known motifs
│── enrichment_results/        # GO/KEGG enrichment outputs
│── figures/                   # MYC binding visualizations
│── scripts/                   # All analysis scripts
│── README.md

# 1️. Dataset

ChIP-seq samples (LUAD cell line):

Sample	Type	Purpose
MYC_rep1	ChIP	MYC binding
MYC_rep2	ChIP	Replicate
IgG_rep1	Control	Background
IgG_rep2	Control	Background

# 2. Quality Control

Tools: FASTQC, MultiQC

fastqc *.fastq -o qc/
multiqc qc/ -o qc/


All samples passed QC. No trimming required.

# 3️. Alignment (Bowtie2)
bowtie2 -x hg38_index -U sample.fastq -S sample.sam
samtools sort -o sample.bam sample.sam
samtools index sample.bam


Alignment rate: 92–95%

# 4️. Peak Calling (MACS2)
macs2 callpeak \
  -t MYC.bam \
  -c IgG.bam \
  -f BAM -g hs \
  -n MYC \
  --outdir macs2_results/


Results:

28,760 MYC peaks (narrowPeak)

Fold-enrichment range: 3×–40×

Files in: macs2_results/

# 5️. Peak Annotation

Tools: bedtools, refGene

bedtools closest \
  -a MYC_peaks.sorted.narrowPeak \
  -b genes_hg38.sorted.bed \
  -D a > MYC_peaks_annotated.bed


Outputs:

MYC_peaks_annotated.bed

MYC_genes.txt (5,319 MYC-associated genes)

# 6️. Motif Enrichment (HOMER)
findMotifsGenome.pl \
  MYC_peaks.sorted.narrowPeak \
  hg38 \
  motif_results/ \
  -size 200


Top enriched motifs:

MYC/MAX (E-box, CACGTG) — correct canonical motif

FOSL1/JUN (AP-1)

ETS/ELK factors

NF-κB/REL family

Indicates MYC co-binding with AP-1, ETS, and NF-κB in LUAD.

# 7️. Functional Enrichment (GO/KEGG)

Tool: GSEAPY (Enrichr)
Script: scripts/run_enrichment.py

Enriched pathways:

GO Biological Process

Transcription regulation

Chromatin organization

RNA processing

Mitotic cell cycle

KEGG

MAPK signaling

PI3K-Akt signaling

Pathways in cancer

Ribosome biogenesis

Results in:
enrichment_results/

# 8️. Key Visualizations

Located in: figures/

Genome-wide MYC peak distribution

MYC binding vs. TSS distance

MYC binding around LUAD oncogenes:

EGFR

FOSL1

TP53

HES4

Example (EGFR):

figures/MYC_EGFR_binding.png


Shows strong MYC promoter binding.

# 9️. Reproducible Workflow

To reproduce entire pipeline:

bash scripts/01_qc.sh
bash scripts/02_alignment.sh
bash scripts/03_peak_calling.sh
bash scripts/04_annotation.sh
python scripts/run_enrichment.py
python scripts/myc_visualizations.py

# 10. Technical Skills Demonstrated

NGS quality assessment

Bowtie2 alignment

SAM/BAM processing

Peak calling with MACS2

Peak annotation

Motif discovery (HOMER)

Functional enrichment (GO/KEGG)

Python visualization

Bash scripting

Workflow organization

Cancer genomics interpretation

# Summary

This project provides a complete analysis of MYC binding in LUAD, identifying:

28,760 MYC-bound regions

5,319 MYC-associated genes

Strong MYC peaks at EGFR, FOSL1, TP53

Enrichment of MYC/MAX and AP-1 motifs

Activation of LUAD-relevant pathways (MAPK, PI3K-AKT, chromatin regulation)
