library(Rsamtools)
library(ggplot2)
library(magrittr)
library(dplyr)
library(DESeq2)
library(maptools)
library(tidyr)

counts_erato <- read.table("MACS2/dem_hyd.MACS2.peaks.all.merged.sort.counts", h=F)
counts_melp <- read.table("MACS2/melp_ros.MACS2.peaks.all.merged.sort.counts", h=F)

sampleInfo_erato <- read.table("MACS2/header_samples_dem_hyd.txt", sep = '\t', h=F)
sampleInfo_melp <- read.table("MACS2/header_samples_melp_ros.txt", sep = '\t', h=F)


# add headers
colnames(sampleInfo_erato) <- c('id','stage','tissue','species')
head(sampleInfo_erato)

header <- c('scaffold','start','end',as.character(sampleInfo_erato$id))

colnames(counts_erato) <- header
head(counts_erato)

counts_erato1 <- counts_erato[,c(1:3)]
counts_erato2 <- counts_erato[,c(4:ncol(counts_erato))]

counts_erato <- cbind(paste(counts_erato1$scaf,counts_erato1$start,counts_erato1$end, sep='_'), counts_erato2)
colnames(counts_erato) <- c('peak', as.character(sampleInfo_erato$id))
head(counts_erato)

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)


# add headers
colnames(sampleInfo_melp) <- c('id','stage','tissue','species')
head(sampleInfo_melp)

header <- c('scaffold','start','end',as.character(sampleInfo_melp$id))

colnames(counts_melp) <- header
head(counts_melp)

counts_melp1 <- counts_melp[,c(1:3)]
counts_melp2 <- counts_melp[,c(4:ncol(counts_melp))]

counts_melp <- cbind(paste(counts_melp1$scaf,counts_melp1$start,counts_melp1$end, sep='_'), counts_melp2)
colnames(counts_melp) <- c('peak', as.character(sampleInfo_melp$id))
head(counts_melp)

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)



sampleInfo_sub <- subset(sampleInfo_melp, sampleInfo_melp$id != 'LI15_melpomene_HW')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI22_melpomene_FP')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI20_demophoon_HP')



counts_melpF2 <- select(counts_melpF, c(as.character(sampleInfo_sub$id)))



dds_erato <- DESeqDataSetFromMatrix(countData = counts_eratoF,
                              colData = sampleInfo_erato,
                              design = ~tissue)

dds_melp <- DESeqDataSetFromMatrix(countData = counts_melpF2,
                                    colData = sampleInfo_sub,
                                    design = ~tissue)


sizes_erato <- estimateSizeFactors(dds_erato)
sizes_melp <- estimateSizeFactors(dds_melp)

write.table(cbind(as.character(sizes_erato$id), 1/sizes_erato$sizeFactor), file="erato_sizes.txt", row.names = F, quote = F, col.names = F)
write.table(cbind(as.character(sizes_melp$id), 1/sizes_melp$sizeFactor), file="melp_sizes.txt", row.names = F, quote = F, col.names = F)

