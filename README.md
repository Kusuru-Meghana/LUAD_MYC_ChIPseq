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

1. Quality Control

- Performed using FastQC and MultiQC.

- All samples demonstrated strong base quality, uniform GC content, and low adapter contamination.

2. Alignment

- Reads aligned to the hg38 reference genome using Bowtie2.

- Alignment rates ranged between 85–98%, indicating high-quality ChIP enrichment.

- BAM files were sorted, indexed, and quality-checked with samtools.

3. Peak Calling

Peaks were called using MACS2, comparing MYC ChIP replicates to IgG controls.

Outputs include:

- *.narrowPeak → High-confidence MYC binding regions

- *_summits.bed → Peak summits (highest enrichment)

- *.xls → Peak statistics

These represent the MYC cistrome in A549 cells.

4. Peak Annotation

Peak annotation was performed using ChIPseeker.

Analyses included:

- Genomic feature distribution

- Peak-to-TSS distance

- Promoter/enhancer annotation

- Peak-to-gene mapping

Key findings:

- ~60% of peaks were located in promoter regions (<3 kb from TSS)

- Strong enrichment around transcription start sites

- Additional binding within introns and distal regulatory elements

5. Functional Enrichment Analysis

Gene-level functional interpretation was performed using:

- clusterProfiler

- org.Hs.eg.db

- msigdbr (MSigDB Hallmark gene sets)

Analyses included:

- Gene Ontology (Biological Processes)

- KEGG Pathway Enrichment

- Hallmark Pathway Enrichment

Enriched biological programs include:

- Cell cycle progression

- DNA replication

- RNA processing and transcription

- Ribosomal biogenesis and translation

- Metabolic reprogramming

- Oxidative phosphorylation

- MYC targets (Hallmark V1/V2)

- E2F target activation

- G2/M checkpoint regulation

These results represent canonical MYC-driven oncogenic pathways.

Key Biological Insights

The analysis reveals that MYC:

1. Binds predominantly at promoters

High promoter occupancy indicates MYC directly regulates transcription initiation in A549 cells.

2. Drives proliferative and metabolic gene programs

Enrichment of MYC hallmark signatures, E2F targets, and cell cycle pathways confirms MYC's central role in cancer proliferation.

3. Activates transcriptional and translational machinery

GO and KEGG analyses highlight MYC regulation of ribosome production, RNA metabolism, and DNA replication.

4. Supports oncogenic signaling in lung adenocarcinoma

Enrichment of G2M checkpoint and oxidative phosphorylation pathways suggests MYC drives both growth and metabolic rewiring.

Overall, this dataset captures a high-fidelity MYC regulatory landscape typical of MYC-driven tumors.


## Future Work

This project is Layer 1 of a multi-omics regulatory analysis.
Upcoming analyses will include:

Layer 2: RNA-seq Integration

- Differential gene expression

- Identification of transcriptionally active MYC targets

- MYC-activated vs. MYC-repressed gene classification

Layer 3: ATAC-seq Integration

- Chromatin accessibility at MYC binding sites

- Open vs. closed chromatin characterization

- Integration with MYC-driven transcriptional activity

## Author

Meghana Kusuru

Bioinformatics & Computational Biology

University of Delaware
