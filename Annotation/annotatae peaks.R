setwd("C:/Users/megha/projects/MYC_Project/LAYER1_ChIPseq")



a


setwd("C:/Users/megha/projects/MYC_Project/LAYER1_ChIPseq/")
getwd()

a

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


#####peak annotation
setwd("C:/Users/megha/projects/MYC_Project/LAYER1_ChIPseq/")

peak_file <- "peaks/MYC_A549_peaks.narrowPeak"
peak <- readPeakFile(peak_file)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_anno <- annotatePeak(
  peak,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)



##view the annotated peaks table
annot_df <- as.data.frame(peak_anno)
head(annot_df)

write.csv(annot_df, "MYC_A549_peak_annotations.csv", row.names = FALSE)


##Peak Annotation Pie Chart
plotAnnoPie(peak_anno)


##TSS Enrichment Plot (Distance to TSS)
plotDistToTSS(peak_anno,
              title = "MYC ChIP-seq Peak Distance to TSS")

##Genomic Feature Distribution
plotAnnoBar(peak_anno)


#Peak Coverage Over Chromosomes
plotChrDistribution(peak_anno)





pdf("MYC_A549_peak_annotation_plots.pdf")
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)
plotAnnoBar(peak_anno)
dev.off()




