#JRT 22May2020
#Input files downloaded from Omicsoft result folders for each project
#P-20170104-0001_Scleroderma_RNA-Seq.QCMetrics.txt
#P-20170104-0001_Lung_RNA-Seq.QCMetrics.txt
#Summarize QC stats 
#was read length, insert size distribution and other technical details about the RNA-seq
#
## ReadLength = 75bp PE 
## Library Pre  = TruSeq Stranded Total RNA Ribo Zero-Gold
## Genome = Mouse.B38
## Gene Model = Ensembl.R86
## 

library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(matrixStats)

qc_skin <- read_delim("P-20170104-0001_Scleroderma_RNA-Seq.QCMetrics.txt", delim="\t")
qc_lung <- read_delim("P-20170104-0001_Lung_RNA-Seq.QCMetrics.txt", delim="\t")
#print available metrics
qc$Metric

metricNames= c("Alignment_All", "Alignment_MappedRate", "Alignment_PairedRate",
               "Alignment_UniquelyMappedRate","InsertSize_Median", "InsertSize_Mean",
               "Profile_ExonRate", "Profile_IntronRate", "Profile_InterGeneRate",
               "Profile_InterGene_FPK", "Source_antisenseRate", "Source_MitochondrialRate")

skinPlots <- QCplots(qc_skin,
        sampleNames = as.character(1:(ncol(qc_skin)-1)),
        metricNames= metricNames)
for (metric in names(skinPlots)){
  ggsave(str_c("skin_", metric, ".png"), skinPlots[[metric]])
}

qcstats_skin <- qc_skin
qcstats_skin$mean <- rowMeans(qc_skin[,2:ncol(qc_skin)])
qcstats_skin$sd <- rowSds(as.matrix(qc_skin[,2:ncol(qc_skin)]))
qcstats_skin$min <- rowMins(as.matrix(qc_skin[,2:ncol(qc_skin)]))
qcstats_skin$max <- rowMaxs(as.matrix(qc_skin[,2:ncol(qc_skin)]))
qcstats_skin[,2:ncol(qc_skin)] <- NULL
#just the metrics selected
idx <- qcstats_skin$Metric %in% metricNames
qcstats_skin <- qcstats_skin[idx,]
write_delim(qcstats_skin, path="qcstats_skin.txt", delim = "\t")

###  Now same for lung
lungPlots <- QCplots(qc_lung,
        sampleNames = as.character(1:(ncol(qc_lung)-1)),
        metricNames= metricNames)

for (metric in names(lungPlots)){
  ggsave(str_c("lung_", metric, ".png"), skinPlots[[metric]])
}

qcstats_lung <- qc_lung
qcstats_lung$mean <- rowMeans(qc_lung[,2:ncol(qc_lung)])
qcstats_lung$sd <- rowSds(as.matrix(qc_lung[,2:ncol(qc_lung)]))
qcstats_lung$min <- rowMins(as.matrix(qc_lung[,2:ncol(qc_lung)]))
qcstats_lung$max <- rowMaxs(as.matrix(qc_lung[,2:ncol(qc_lung)]))
qcstats_lung[,2:ncol(qc_lung)] <- NULL
#just the metrics selected
idx <- qcstats_lung$Metric %in% metricNames
qcstats_lung <- qcstats_lung[idx,]
write_delim(qcstats_lung, path="qcstats_lung.txt", delim = "\t")
