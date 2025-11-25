# Analysis Scripts

This folder contains reproducible workflows for ChIP-seq processing and enrichment.

### Scripts
- **chipseq_pipeline.sh**  
  Bash pipeline for:
  - alignment (bowtie2)
  - sorting/indexing
  - MACS2 peak calling
  - QC

- **enrichment_analysis.R**  
  R script for functional analysis:
  - GO, KEGG, and Hallmark enrichment
  - visualization

All steps can be reproduced by running these scripts.
