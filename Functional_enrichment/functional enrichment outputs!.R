setwd("C:/Users/megha/projects/MYC_Project/LAYER1_ChIPseq/")





BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("clusterProfiler")
install.packages("msigdbr")



library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(DOSE)




peak_file <- "peaks/MYC_A549_peaks.narrowPeak"
peak <- readPeakFile(peak_file)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")


genes <- unique(peak_anno@anno$geneId)


ego <- enrichGO(genes, OrgDb="org.Hs.eg.db", keyType="ENTREZID")



ekegg <- enrichKEGG(genes, organism="hsa")


m_df <- msigdbr(species="Homo sapiens", category="H")
hallmark_list <- split(m_df$entrez_gene, m_df$gs_name)
ehall <- enricher(genes, TERM2GENE=hallmark_list)



library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Extract HUMAN Hallmark gene sets
m_df <- msigdbr(species = "Homo sapiens", category = "H")

# Convert msigdbr into TERM2GENE for clusterProfiler
hallmark_list <- m_df[, c("gs_name", "entrez_gene")]

head(hallmark_list)


ehall <- enricher(genes, TERM2GENE = hallmark_list)


head(ehall)



write.csv(as.data.frame(ego),   "GO_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_results.csv")
write.csv(as.data.frame(ehall), "Hallmark_results.csv")







