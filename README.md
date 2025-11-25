# LUAD MYC ChIP-seq - Genome-wide mapping of MYC binding sites in A549 lung adenocarcinoma cells using ChIP-Seq

## Overview

This repository contains a complete ChIP-seq analysis of the MYC transcription factor in A549 lung adenocarcinoma cells. The goal of this project is to define the genome-wide MYC binding landscape (cistrome) and characterize the biological pathways, regulatory programs, and transcriptional processes controlled by MYC in lung cancer.

This analysis constitutes Layer 1 of a planned multi-omics regulatory map (ChIP-seq → RNA-seq → ATAC-seq).

## Objective
- Identify high-confidence MYC binding sites across the genome.

- Annotate peaks to genes and genomic features.

- Characterize biological processes and pathways enriched among MYC-bound genes.

- Establish the role of MYC in regulating transcription, proliferation, and oncogenic signaling in A549 cells.

### What is MYC?

MYC (c-Myc) is a master regulatory transcription factor controlling:

- Ribosome biogenesis

- Metabolism

- Cell-cycle progression

- Nucleotide synthesis

In cancer (including LUAD, lung adenocarcinoma):

- MYC becomes amplified or overactive
- It drives continuous proliferation

- Rewires metabolism

- Blocks differentiation

- Cooperates with KRAS and other oncogenic drivers

Understanding where MYC binds DNA in cancer cells reveals:

- The gene networks it controls

- Oncogenic programs it activates

- Potential drug targets

## Why ChIP-Seq for MYC?
hIP-Seq identifies genomic regions bound by transcription factors.
For MYC in A549 cells, it answers - "Which genes does MYC directly regulate in lung cancer?”

ChIP-Seq provides:

- Promoter/enhancer binding sites

- Peak summits

- Motif enrichment

- Pathway activation

- Regulatory circuits

- Cancer-specific MYC targets

## Project Structure

```
LUAD_MYC_ChIPseq/
│
├── raw_fastq/                     # Raw FASTQ files (not included due to size)
├── qc/                            # FastQC and MultiQC reports
├── aligned_bam/                   # Sorted BAM + index files (not included)
├── peaks/                         # MACS2 peak outputs (.narrowPeak, .xls, summits)
│
├── Annotation/
│   ├── MYC_A549_peak_annotations.csv
│   └── MYC_A549_peak_annotation_plots.pdf
│
├── Functional_enrichment/
│   ├── GO_results.csv
│   ├── KEGG_results.csv
│   └── Hallmark_results.csv
│
├── scripts/
│   ├── chipseq_pipeline.sh        # Complete ChIP-seq workflow
│   └── enrichment_analysis.R      # GO/KEGG/Hallmark enrichment
│
├── MYC_ChIPseq_Friendly_Summary.pdf
└── README.md
```



## Biological Samples

Four SRA datasets were used (2 MYC ChIP, 2 IgG controls).
This matches ENCODE guidelines for reproducible TF ChIP-Seq.

| Condition       | Replicate | SRA ID        | File Name            | Source Link                                                                              |
| --------------- | --------- | ------------- | -------------------- | ---------------------------------------------------------------------------------------- |
| **MYC ChIP**    | Rep 1     | **SRR568131** | `SRR568131.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568131](https://www.ncbi.nlm.nih.gov/sra/SRR568131) |
| **MYC ChIP**    | Rep 2     | **SRR568132** | `SRR568132.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568132](https://www.ncbi.nlm.nih.gov/sra/SRR568132) |
| **IgG Control** | Rep 1     | **SRR568133** | `SRR568133.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568133](https://www.ncbi.nlm.nih.gov/sra/SRR568133) |
| **IgG Control** | Rep 2     | **SRR568134** | `SRR568134.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568134](https://www.ncbi.nlm.nih.gov/sra/SRR568134) |


All datasets originate from a single ENCODE ChIP-seq experiment.

## ChIP-Seq Workflow

## 1. Quality Control

Tools used:

- FastQC

- MultiQC

QC Results:

- Very low duplication rates (4–18%)

- Perfect base quality (Q35–40)

- Correct GC distribution (40–43%)

- No adapter contamination

- No N-content issues

These indicate high-quality ENCODE-grade datasets.

## 2. Alignment (Bowtie2 → hg38)

Steps:

- Download hg38 reference genome

- Build Bowtie2 index

- Align MYC and IgG reads

- Convert SAM → BAM

- Sort and index BAM

Bowtie2 was selected because it is optimised for:

- Short reads (36 bp)

- TF ChIP-Seq

- Large genomes (human hg38)


## 3. Peak Calling (MACS2)

Peak calling was performed with

```bash
macs2 callpeak -t MYC_rep1.bam MYC_rep2.bam \
               -c IgG_rep1.bam IgG_rep2.bam \
               -f BAM -g hs --outdir peak_calling/
```
MACS2 outputs produced:

- narrowPeak files (primary peaks)

- summits.bed (motif sites)

- .xls statistics

- bigWig / bedGraph files (optional)

MACS2 effectively finds locations where MYC signal exceeds IgG background.

## 4. Annotation of Peaks

Tools used:

- ChIPseeker

- HOMER (motifs)

Outputs:

- Promoter/enhancer assignments

- Distance to TSS

- Gene names for MYC-bound regions

- Motif enrichment (E-box CACGTG expected)

## 5. Functional Enrichment

Performed using GSEApy:

- GO Biological Process

- GO Molecular Function

- GO Cellular Component

- KEGG pathways

Enrichment for MYC targets included:

- Ribosome biogenesis

- Cell-cycle regulation

- RNA metabolism

- Oncogenic signalling pathways

All results saved in:
[Enrichment Results Folder](results/enrichment_results/)



## 6. Visualization

Figures include:
### **Genome-wide MYC binding strength**
![Genome-wide MYC peaks](MYC_genome_wide_peaks.png)

### **MYC binding around EGFR TSS**
![EGFR binding](MYC_EGFR_binding.png)

### **MYC binding around FOSL1 TSS**
![FOSL1 binding](MYC_FOSL1_binding.png)

### **MYC binding around TP53 TSS**
![TP53 binding](MYC_TP53_binding.png)

### **MYC binding upstream of HES4**
![HES4 binding](MYC_HES4_binding.png)


Located in:
[results/figures/](results/figures)


## Key Findings (summary)

- High-confidence MYC peaks were identified across the genome.

- MYC binds promoters of major cancer-related genes (EGFR, FOSL1, etc.).

- Motif analysis reveals strong enrichment of canonical E-box motifs.

- Functional enrichment highlights MYC’s role in ribosome biogenesis, cell proliferation, and metabolic regulation.

- Results align with known MYC oncogenic behaviour in LUAD.

## How to Run This Pipeline
Clone the repo

```bash
git clone https://github.com/Kusuru-Meghana/LUAD_MYC_ChIPseq.git
cd LUAD_MYC_ChIPseq
```

Run the entire pipeline
```bash
Scripts/run_chipseq_pipeline.sh
```

Run enrichment analysis
```
python3 Scripts/run_enrichment.py
```


## Software Requirements

- FastQC

- MultiQC

- Bowtie2

- Samtools

- MACS2

- ChIPseeker (R)

- HOMER

- GSEApy (Python)

- Environment files (conda/Docker) can be added on request.

# Author

Meghana Kusuru
LUAD MYC ChIP-Seq Analysis
2025


