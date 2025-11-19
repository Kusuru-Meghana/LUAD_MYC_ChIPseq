# LUAD MYC ChIP-seq - Genome-wide mapping of MYC binding sites in A549 lung adenocarcinoma cells using ChIP-Seq

## Overview

This repository contains a complete, reproducible ChIP-Seq workflow to identify genome-wide MYC binding sites in the A549 lung adenocarcinoma cell line.
The project includes:

- Data acquisition from SRA

- Quality control (FastQC + MultiQC)

- Alignment to hg38 using Bowtie2

- Peak calling with MACS2

- Annotation of MYC peaks

- GO/KEGG functional enrichment

- Visualization of MYC binding profiles

- Figures, logs, and processed results

This project serves as a full end-to-end example of transcription factor ChIP-Seq analysis, following ENCODE standards.

## Biological Motivation

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

2. Alignment (Bowtie2 → hg38)

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


3. Peak Calling (MACS2)

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

4. Annotation of Peaks

Tools used:

- ChIPseeker

- HOMER (motifs)

Outputs:

- Promoter/enhancer assignments

- Distance to TSS

- Gene names for MYC-bound regions

- Motif enrichment (E-box CACGTG expected)

5. Functional Enrichment

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



6. Visualization

Figures include:

- Genome-wide MYC peak distribution

- Binding tracks (EGFR, FOSL1, HES4, TP53, etc.)

- Peak annotation plots

- Enrichment dotplots

Located in:
results/figures/







1. QC (FastQC + MultiQC)

- Q-scores: 35–40 (excellent)

- GC content: ~41–43%

- Duplicate rate: low (4–18%)

- No adapters → no trimming needed

2. Alignment (Bowtie2 → hg38)

- 26–31 million reads per sample

- High alignment rate

- Sorted + indexed BAMs used for peak calling

3. Peak Calling (MACS2)

```bash
macs2 callpeak -t MYC.bam -c IgG.bam -g hs -f BAM -n MYC
```

Results:

- ~5,000+ MYC peaks

- Enrichment up to 40×

- Summits used for motif analysis

4. Peak Annotation (ChIPseeker)

- ~50% promoter peaks

- ~935 MYC-bound genes

5. Enrichment Analysis (GSEApy + Enrichr)

- GO: Biological Process

- GO: Molecular Function

- GO: Cellular Component

- KEGG Pathways

6. Motif Analysis (HOMER)

Top enriched motif: MYC E-box (CACGTG)

Example Figures
- Genome-wide MYC Peak Distribution

- MYC Binding at EGFR Promoter

## Key Results
1. Strong MYC binding at LUAD oncogenes
   
| Gene             | Fold Enrichment | Interpretation                             |
| ---------------- | --------------- | ------------------------------------------ |
| **EGFR**         | 18.2×           | Direct promoter binding → growth signaling |
| **FOSL1 (AP-1)** | 9×              | Oncogenic transcription factor             |
| **TP53**         | 8.6×            | Stress & checkpoint regulation             |
| **HES4**         | 5.7×            | Notch pathway regulator                    |

2. Genome-wide MYC activity

- Peaks cluster around 4–6× enrichment

- Long high-confidence tail up to 40×

- Classic transcription factor ChIP-seq profile

3. GO Biological Process — MYC drives biosynthesis

- Top enriched biological processes:

- Ribosome biogenesis

- Translation

- rRNA processing

- Macromolecule biosynthesis

- Peptide biosynthetic process

→ MYC activates biosynthetic and proliferative programs.

4. GO Molecular Function

- RNA binding

- mRNA 5′-UTR binding

- rRNA binding

- snoRNA binding

- Cadherin binding

→ MYC regulates post-transcriptional control & adhesion.

5. GO Cellular Component

- Nucleolus

- Ribosomal subunits

- Small-subunit processome

- Focal adhesions

→ MYC activates nucleolar, translational, and migration modules.

6. KEGG Pathway Enrichment

- Most enriched pathways:

- Ribosome (84 genes)

- Ribosome biogenesis

- RNA transport

- Purine metabolism

- PI3K/AKT & insulin signaling

- Ferroptosis

- Spliceosome

→ MYC controls metabolic, proliferative, and stress-survival pathways.

## Final Biological Interpretation

In A549 LUAD cells, MYC directly binds and activates a multi-layer regulatory program that drives tumor progression.
MYC controls:

- Growth signaling (EGFR, STAT3, FOSL1)

- Cell-cycle progression (CDK4, CCND1, WEE1)

- Metabolic rewiring (HK2, LDHA, FASN)

- Ribosome & translation machinery (hundreds of RPL/RPS/EIF genes)

- RNA processing & splicing

- Ferroptosis regulators

- Stress and DNA-repair programs

This represents a canonical, MYC-driven LUAD regulatory state.

## Repository Structure


```
LUAD_MYC_ChIPseq/
│
├── LICENSE
├── README.md
│
├── Scripts/
│   ├── .gitkeep
│   ├── QC/
│   │   └── 01_qc.sh
│   ├── alignment/
│   │   └── bowtie2_alignment.sh
│   ├── annotation/
│   │   └── annotate_peaks.sh
│   ├── myc_visualizations.py
│   ├── peak_calling/
│   │   └── macs2_callpeaks.sh
│   ├── run_chipseq_pipeline.sh
│   └── run_enrichment.py
│
├── results/
│   ├── MYC_genes.txt/
│   │   ├── .gitkeep
│   │   └── MYC_genes.txt
│   │
│   ├── annotated_peaks/
│   │   ├── .gitkeep
│   │   └── MYC_peaks_annotated.bed
│   │
│   ├── enrichment_results/
│   │   ├── .gitkeep
│   │   ├── GO_Biological_Process_2023.Human.enrichr.reports.pdf
│   │   ├── GO_Cellular_Component_2023.Human.enrichr.reports.pdf
│   │   ├── GO_Molecular_Function_2023.Human.enrichr.reports.pdf
│   │   ├── KEGG_2021_Human.Human.enrichr.reports.pdf
│   │   ├── gseapy.enrichr.127968526904224.log
│   │   ├── gseapy.enrichr.129469226776096.log
│   │   └── gseapy.enrichr.139063081760288.log
│   │
│   └── figures/
│       ├── .gitkeep
│       ├── MYC_EGFR_binding.png
│       ├── MYC_FOSL1_binding.png
│       ├── MYC_HES4_binding.png
│       ├── MYC_TP53_binding.png
│       └── MYC_genome_wide_peaks.png
│
└── (root directory)

```

*(Large FASTQ/BAM files not included.)*


## Tools Used

- FastQC

- MultiQC

- Bowtie2

- samtools

- MACS2

- ChIPseeker (R)

- HOMER

- GSEApy


## License

Released under the MIT License.
