# LUAD MYC ChIP-seq - Mapping MYC Binding Sites in Lung Adenocarcinoma (A549)

## Overview

This project performs a complete ChIP-seq pipeline to identify genome-wide binding sites of the oncogenic transcription factor MYC in A549 lung adenocarcinoma (LUAD) cells.
Using ENCODE-matched replicates (2 MYC + 2 IgG controls), the analysis recovers:

- High-confidence MYC binding peaks

- Direct MYC-regulated oncogenes

- Key LUAD pathways controlled by MYC

- Peak annotation, motif enrichment, and pathway biology


## Biological Motivation

MYC is a master regulator of tumor growth. In LUAD, MYC drives:

- Ribosome biogenesis

- Protein synthesis

- Metabolic rewiring

- Cell-cycle progression

- Stress and survival signaling

## Goal

Identify where MYC binds the genome in A549 cells, and what biological programs MYC directly regulates.

## Data Summary

| Condition       | Replicate | SRA ID        | File Name            | Source Link                                                                              |
| --------------- | --------- | ------------- | -------------------- | ---------------------------------------------------------------------------------------- |
| **MYC ChIP**    | Rep 1     | **SRR568131** | `SRR568131.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568131](https://www.ncbi.nlm.nih.gov/sra/SRR568131) |
| **MYC ChIP**    | Rep 2     | **SRR568132** | `SRR568132.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568132](https://www.ncbi.nlm.nih.gov/sra/SRR568132) |
| **IgG Control** | Rep 1     | **SRR568133** | `SRR568133.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568133](https://www.ncbi.nlm.nih.gov/sra/SRR568133) |
| **IgG Control** | Rep 2     | **SRR568134** | `SRR568134.fastq.gz` | [https://www.ncbi.nlm.nih.gov/sra/SRR568134](https://www.ncbi.nlm.nih.gov/sra/SRR568134) |


All datasets originate from a single ENCODE ChIP-seq experiment.

## Pipeline
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
├── peaks/
├── results/
│   ├── annotated_peaks/
│   ├── enrichment_results/
│   ├── figures/
├── scripts/
└── README.md
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
