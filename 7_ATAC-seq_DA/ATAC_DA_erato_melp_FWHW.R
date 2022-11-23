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

peaks_all_erato <- read.table("MACS2/erato_peaks_start_end_pan.txt", h=F)
peaks_all_melp <- read.table("MACS2/melp_peaks_start_end_pan.txt", h=F)

peaks_all_erato_sort <- as.data.frame(t(apply(peaks_all_erato[,c(2,3)],1,sort)))
peaks_all_melp_sort <- as.data.frame(t(apply(peaks_all_melp[,c(2,3)],1,sort)))

peaks_all_erato <- cbind(peaks_all_erato[,1],peaks_all_erato_sort)
peaks_all_melp <- cbind(peaks_all_melp[,1],peaks_all_melp_sort)

# add headers
colnames(sampleInfo_erato) <- c('id','stage','tissue','race')
colnames(sampleInfo_melp) <- c('id','stage','tissue','race')

sampleInfo_erato$species <- 'erato'
sampleInfo_melp$species <- 'melp'

sampleInfo <- rbind(sampleInfo_erato, sampleInfo_melp)

header_erato <- c('scaffold', 'start', 'end', as.character(sampleInfo_erato$id))
header_melp <- c('scaffold', 'start', 'end', as.character(sampleInfo_melp$id))

colnames(counts_erato) <- header_erato
colnames(counts_melp) <- header_melp

colnames(peaks_all_erato) <- c('genome', 'pan_start', 'pan_end')
colnames(peaks_all_melp) <- c('genome', 'pan_start', 'pan_end')

nrow(counts_erato)
nrow(counts_melp)

nrow(peaks_all_erato)
nrow(peaks_all_melp)

# combine counts with pan genome positions
counts_erato_pan <- cbind(peaks_all_erato, counts_erato)
counts_melp_pan <- cbind(peaks_all_melp, counts_melp)

# load intersect files
match_ME <- read.table("MACS2/match_melp_to_erato.0.5.txt", h=F)
match_EM <- read.table("MACS2/match_erato_to_melp.0.5.txt", h=F)

colnames(match_ME) <- c('genome1', 'pan_start_E', 'pan_end_E', 'genome2', 'pan_start_M', 'pan_end_M', 'length')
colnames(match_EM) <- c('genome1', 'pan_start_M', 'pan_end_M', 'genome2', 'pan_start_E', 'pan_end_E', 'length')

match_ME_A <- subset(match_ME, match_ME$genome2 != '.')
match_EM_A <- subset(match_EM, match_EM$genome2 != '.')

nrow(match_ME_A)
nrow(match_EM_A)

# merge intersect files with peak counts
counts_erato_pan_match <- merge(counts_erato_pan, match_ME, by.x = c('pan_start','pan_end'), by.y =c('pan_start_E', 'pan_end_E'))
# counts_melp_pan_match <- merge(counts_melp_pan, match_EM, by.x = c('pan_start','pan_end'), by.y =c('pan_start_M', 'pan_end_M'))

counts_erato_melp_pan_match <- merge(counts_erato_pan_match, counts_melp_pan, by.x = c('pan_start_M', 'pan_end_M'), by.y = c('pan_start','pan_end'))
counts_erato_melp_pan_match <- merge(counts_erato_pan_match, counts_melp_pan, by.x = c('pan_start_M', 'pan_end_M'), by.y = c('pan_start','pan_end'), all.x = TRUE, all.y = TRUE)

head(counts_erato_melp_pan_match)

# cts1 <- counts_erato_melp_pan_match[,c(5,3,4,1,2)] 5,3,4,1,2,6,7,8
cts1 <- counts_erato_melp_pan_match[,c('genome.x', 'pan_start', 'pan_end', 'pan_start_M', 'pan_end_M', 'scaffold.x', 'start.x', 'end.x', 'scaffold.y','start.y', 'end.y')]
cts2 <- counts_erato_melp_pan_match[,c(as.character(sampleInfo_erato$id), as.character(sampleInfo_melp$id))]

# cts <- cbind(paste(cts1$genome.x,cts1$pan_start,cts1$pan_end,cts1$pan_start_M,cts1$pan_end_M, sep='_'), cts2)
cts <- cbind(paste(cts1$genome.x,cts1$pan_start,cts1$pan_end,cts1$pan_start_M,cts1$pan_end_M, cts1$scaffold.x,cts1$start.x,cts1$end.x, cts1$scaffold.y, cts1$start.y, cts1$end.y, sep='_'), cts2)
colnames(cts) <- c('peak', c(as.character(sampleInfo_erato$id), as.character(sampleInfo_melp$id)))
head(cts)

ctsF <- cts[,-1]
rownames(ctsF) <- cts[,1]
head(ctsF)

ctsF[is.na(ctsF)] <- 0

sampleInfo_sub <- subset(sampleInfo, sampleInfo$id != 'LI15_melpomene_HW')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI22_melpomene_FP')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI20_demophoon_HP')

ctsF2 <- select(ctsF, c(as.character(sampleInfo_sub$id)))

ctsF2A <- ctsF2[-grep("NA", rownames(ctsF2)), ]
######################
# 5th FW vs HW
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato' &
                                             sampleInfo_sub$stage == '5thinstar'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp' &
                                             sampleInfo_sub$stage == '5thinstar'))

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue == 'forewing') &
                                             # sampleInfo_sub$species == 'erato' &
                                             sampleInfo_sub$stage == '5thinstar'))

sampleInfo_subAHW <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue == 'hindwing') &
                                             # sampleInfo_sub$species == 'erato' &
                                             sampleInfo_sub$stage == '5thinstar'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))
ctsF2AFW <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))
ctsF2AHW <- select(ctsF2A, c(as.character(sampleInfo_subAHW$id)))
###
# setup model
###

dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                              colData = sampleInfo_sub,
                              design = ~tissue)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)

ddsAFW <- DESeqDataSetFromMatrix(countData = ctsF2AFW,
                              colData = sampleInfo_subA,
                              design = ~species)

ddsAHW <- DESeqDataSetFromMatrix(countData = ctsF2AHW,
                               colData = sampleInfo_subAHW,
                               design = ~species)
###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
atacDDS_AFW <- DESeq(ddsA)
atacDDS_AHW <- DESeq(ddsAHW)

res_E_5thFWHW <- results(atacDDS_E)
res_M_5thFWHW <- results(atacDDS_M)
res_A_5thFW <- results(atacDDS_AFW)
res_A_5thHW <- results(atacDDS_AHW)

save(res_E_5thFWHW, file='sign_E_5th_FWHW_p05fc1_ALL.rda')
save(res_M_5thFWHW, file='sign_M_5th_FWHW_p05fc1_ALL.rda')

sign_E_5th_FWHW_HWup <- subset(res_E, res_E$padj < 0.05 & res_E$log2FoldChange > 1)
sign_M_5th_FWHW_HWup <- subset(res_M, res_M$padj < 0.05 & res_M$log2FoldChange > 1)

sign_E_5th_FWHW_FWup <- subset(res_E, res_E$padj < 0.05 & res_E$log2FoldChange < -1)
sign_M_5th_FWHW_FWup <- subset(res_M, res_M$padj < 0.05 & res_M$log2FoldChange < -1)

sign_A_5th_FW <- subset(res_A_5thFW, res_A_5thFW$padj < 0.05 & abs(res_A_5thFW$log2FoldChange > 1))
sign_A_5th_HW <- subset(res_A_5thHW, res_A_5thHW$padj < 0.05 & abs(res_A_5thHW$log2FoldChange > 1))

save(sign_E_5th_FWHW_HWup, file='sign_E_5th_FWHW_p05fc1_HWup.rda')
save(sign_M_5th_FWHW_HWup, file='sign_M_5th_FWHW_p05fc1_HWup.rda')

save(sign_E_5th_FWHW_FWup, file='sign_E_5th_FWHW_p05fc1_FWup.rda')
save(sign_M_5th_FWHW_FWup, file='sign_M_5th_FWHW_p05fc1_FWup.rda')

nrow(sign_A_5th_FW)
nrow(merge(as.data.frame(sign_A_5th_FW), as.data.frame(sign_E_5th_FWHW_FWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_5th_FW), as.data.frame(sign_M_5th_FWHW_FWup), by = 'row.names'))

nrow(sign_A_5th_HW)
nrow(merge(as.data.frame(sign_A_5th_HW), as.data.frame(sign_E_5th_FWHW_HWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_5th_HW), as.data.frame(sign_M_5th_FWHW_HWup), by = 'row.names'))

######################
# D1 vs 5th/D2
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato' &
                                             sampleInfo_sub$stage == 'pupalDAY1'))

sampleInfo_subE$comp <- 'forewing'
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] == 'HA' | sampleInfo_subE$tissue[e] == 'HP'){sampleInfo_subE$comp[e] <- 'hindwing' }}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp' &
                                             sampleInfo_sub$stage == 'pupalDAY1'))

sampleInfo_subM$comp <- 'forewing'
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] == 'HA' | sampleInfo_subM$tissue[e] == 'HP'){sampleInfo_subM$comp[e] <- 'hindwing' }}

sampleInfo_sub$comp <- 'forewing'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$tissue[e] == 'HA' | sampleInfo_sub$tissue[e] == 'HP'){sampleInfo_sub$comp[e] <- 'hindwing' }}

sampleInfo_subAFW <- subset(sampleInfo_sub, ((sampleInfo_sub$comp == 'forewing') &
                                               # sampleInfo_sub$species == 'erato' &
                                               sampleInfo_sub$stage == 'pupalDAY1'))

sampleInfo_subAHW <- subset(sampleInfo_sub, ((sampleInfo_sub$comp == 'hindwing') &
                                               # sampleInfo_sub$species == 'erato' &
                                               sampleInfo_sub$stage == 'pupalDAY1'))


ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ctsF2AFW <- select(ctsF2A, c(as.character(sampleInfo_subAFW$id)))
ctsF2AHW <- select(ctsF2A, c(as.character(sampleInfo_subAHW$id)))


###
# setup model
###

dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                              colData = sampleInfo_sub,
                              design = ~tissue)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~comp)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~comp)

ddsAFW <- DESeqDataSetFromMatrix(countData = ctsF2AFW,
                               colData = sampleInfo_subAFW,
                               design = ~species)

ddsAHW <- DESeqDataSetFromMatrix(countData = ctsF2AHW,
                                 colData = sampleInfo_subAHW,
                                 design = ~species)

###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E_D1FWHW <- results(atacDDS_E)
res_M_D1FWHW <- results(atacDDS_M)

atacDDS_AFW <- DESeq(ddsAFW)
atacDDS_AHW <- DESeq(ddsAHW)
res_E_D1AFW <- results(atacDDS_AFW)
res_M_D1AHW <- results(atacDDS_AHW)

save(res_E_D1FWHW, file='sign_E_D1_FWHW_p05fc1_ALL.rda')
save(res_M_D1FWHW, file='sign_M_D1_FWHW_p05fc1_ALL.rda')

sign_E_D1_FWHW_HWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_FWHW_HWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D1_FWHW_FWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_FWHW_FWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

sign_A_D1_FW <- subset(res_E_D1AFW, res_E_D1AFW$padj < 0.05 & abs(res_E_D1AFW$log2FoldChange > 1))
sign_A_D1_HW <- subset(res_M_D1AHW, res_M_D1AHW$padj < 0.05 & abs(res_M_D1AHW$log2FoldChange > 1))

save(sign_E_D1_FWHW_HWup, file='sign_E_D1_FWHW_p05fc1_HWup.rda')
save(sign_M_D1_FWHW_HWup, file='sign_M_D1_FWHW_p05fc1_HWup.rda')

save(sign_E_D1_FWHW_FWup, file='sign_E_D1_FWHW_p05fc1_FWup.rda')
save(sign_M_D1_FWHW_FWup, file='sign_M_D1_FWHW_p05fc1_FWup.rda')


nrow(sign_A_D1_FW)
nrow(merge(as.data.frame(sign_A_D1_FW), as.data.frame(sign_E_D1_FWHW_FWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FW), as.data.frame(sign_M_D1_FWHW_FWup), by = 'row.names'))

nrow(sign_A_D1_HW)
nrow(merge(as.data.frame(sign_A_D1_HW), as.data.frame(sign_E_D1_FWHW_HWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_HW), as.data.frame(sign_M_D1_FWHW_HWup), by = 'row.names'))

######################
# D2 vs 5th/D1
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato'&
                                             sampleInfo_sub$stage == 'pupalDAY2'))

sampleInfo_subE$comp <- 'forewing'
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] == 'HA' | sampleInfo_subE$tissue[e] == 'HP'){sampleInfo_subE$comp[e] <- 'hindwing' }}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp' &
                                             sampleInfo_sub$stage == 'pupalDAY2'))

sampleInfo_subM$comp <- 'forewing'
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] == 'HA' | sampleInfo_subM$tissue[e] == 'HP'){sampleInfo_subM$comp[e] <- 'hindwing' }}


sampleInfo_subAFW <- subset(sampleInfo_sub, ((sampleInfo_sub$comp == 'forewing') &
                                               # sampleInfo_sub$species == 'erato' &
                                               sampleInfo_sub$stage == 'pupalDAY2'))

sampleInfo_subAHW <- subset(sampleInfo_sub, ((sampleInfo_sub$comp == 'hindwing') &
                                               # sampleInfo_sub$species == 'erato' &
                                               sampleInfo_sub$stage == 'pupalDAY2'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ctsF2AFW <- select(ctsF2A, c(as.character(sampleInfo_subAFW$id)))
ctsF2AHW <- select(ctsF2A, c(as.character(sampleInfo_subAHW$id)))

###
# setup model
###

dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                              colData = sampleInfo_sub,
                              design = ~tissue)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~comp)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~comp)

ddsAFW <- DESeqDataSetFromMatrix(countData = ctsF2AFW,
                                 colData = sampleInfo_subAFW,
                                 design = ~species)

ddsAHW <- DESeqDataSetFromMatrix(countData = ctsF2AHW,
                                 colData = sampleInfo_subAHW,
                                 design = ~species)
###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E_D2FWHW <- results(atacDDS_E)
res_M_D2FWHW <- results(atacDDS_M)

atacDDS_AFW <- DESeq(ddsAFW)
atacDDS_AHW <- DESeq(ddsAHW)
res_E_D2AFW <- results(atacDDS_AFW)
res_M_D2AHW <- results(atacDDS_AHW)

save(res_E_D2FWHW, file='sign_E_D2_FWHW_p05fc1_ALL.rda')
save(res_M_D2FWHW, file='sign_M_D2_FWHW_p05fc1_ALL.rda')

sign_E_D2_FWHW_HWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FWHW_HWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FWHW_FWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FWHW_FWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

sign_A_D2_FW <- subset(res_E_D2AFW, res_E_D2AFW$padj < 0.05 & abs(res_E_D2AFW$log2FoldChange > 1))
sign_A_D2_HW <- subset(res_M_D2AHW, res_M_D2AHW$padj < 0.05 & abs(res_M_D2AHW$log2FoldChange > 1))

save(sign_E_D2_FWHW_HWup, file='sign_E_D2_FWHW_p05fc1_HWup.rda')
save(sign_M_D2_FWHW_HWup, file='sign_M_D2_FWHW_p05fc1_HWup.rda')

save(sign_E_D2_FWHW_FWup, file='sign_E_D2_FWHW_p05fc1_FWup.rda')
save(sign_M_D2_FWHW_FWup, file='sign_M_D2_FWHW_p05fc1_FWup.rda')

nrow(sign_A_D2_FW)
nrow(merge(as.data.frame(sign_A_D2_FW), as.data.frame(sign_E_D2_FWHW_FWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D2_FW), as.data.frame(sign_M_D2_FWHW_FWup), by = 'row.names'))

nrow(sign_A_D2_HW)
nrow(merge(as.data.frame(sign_A_D2_HW), as.data.frame(sign_E_D2_FWHW_HWup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D2_HW), as.data.frame(sign_M_D2_FWHW_HWup), by = 'row.names'))

### load

load("sign_E_5th_FWHW_p05fc1_HWup.rda")
load("sign_M_5th_FWHW_p05fc1_HWup.rda")

load("sign_E_D1_FWHW_p05fc1_HWup.rda")
load("sign_M_D1_FWHW_p05fc1_HWup.rda")

load("sign_E_D2_FWHW_p05fc1_HWup.rda")
load("sign_M_D2_FWHW_p05fc1_HWup.rda")

load("sign_E_5th_FWHW_p05fc1_FWup.rda")
load("sign_M_5th_FWHW_p05fc1_FWup.rda")

load("sign_E_D1_FWHW_p05fc1_FWup.rda")
load("sign_M_D1_FWHW_p05fc1_FWup.rda")

load("sign_E_D2_FWHW_p05fc1_FWup.rda")
load("sign_M_D2_FWHW_p05fc1_FWup.rda")


load("sign_E_5th_FWHW_p05fc1.rda")
load("sign_M_5th_FWHW_p05fc1.rda")

load("sign_E_D1_FWHW_p05fc1.rda")
load("sign_M_D1_FWHW_p05fc1.rda")

load("sign_E_D2_FWHW_p05fc1.rda")
load("sign_M_D2_FWHW_p05fc1.rda")


sign_E_5th_FWHW

###
# PCA
###
sampleInfo_sub$col <- NA 
sampleInfo_sub$shape <- 19 
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$stage[e] == '5thinstar') { sampleInfo_sub$col[e] = '#ff0000b5'}
  if(sampleInfo_sub$stage[e] == 'pupalDAY1') { sampleInfo_sub$col[e] = '#009e73ff'}
  if(sampleInfo_sub$stage[e] == 'pupalDAY2') { sampleInfo_sub$col[e] = '#0072b2ff'}
  if(sampleInfo_sub$tissue[e] == 'brain') { sampleInfo_sub$col[e] = '#f0e442ff'}
  
  if(sampleInfo_sub$tissue[e] == 'brain' & sampleInfo_sub$species[e] == 'erato') { sampleInfo_sub$shape[e] = 19}
  if(sampleInfo_sub$tissue[e] == 'brain' & sampleInfo_sub$species[e] == 'melp') { sampleInfo_sub$shape[e] = 1}
  
  if(sampleInfo_sub$stage[e] == '5thinstar' & sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'forewing') { sampleInfo_sub$shape[e] = 19}
  if(sampleInfo_sub$stage[e] == '5thinstar' & sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'hindwing') { sampleInfo_sub$shape[e] = 17}
  if(sampleInfo_sub$stage[e] == '5thinstar' & sampleInfo_sub$species[e] == 'melp'  & sampleInfo_sub$tissue[e] == 'forewing') { sampleInfo_sub$shape[e] = 1}
  if(sampleInfo_sub$stage[e] == '5thinstar' & sampleInfo_sub$species[e] == 'melp'  & sampleInfo_sub$tissue[e] == 'hindwing') { sampleInfo_sub$shape[e] = 2}
  
  if(sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'FP') { sampleInfo_sub$shape[e] = 19}
  if(sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'FM') { sampleInfo_sub$shape[e] = 17}
  if(sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'FD') { sampleInfo_sub$shape[e] = 15}
  if(sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'HA') { sampleInfo_sub$shape[e] = 23}
  if(sampleInfo_sub$species[e] == 'erato' & sampleInfo_sub$tissue[e] == 'HP') { sampleInfo_sub$shape[e] = 25}
  
  if(sampleInfo_sub$species[e] == 'melp' & sampleInfo_sub$tissue[e] == 'FP') { sampleInfo_sub$shape[e] = 1}
  if(sampleInfo_sub$species[e] == 'melp' & sampleInfo_sub$tissue[e] == 'FM') { sampleInfo_sub$shape[e] = 2}
  if(sampleInfo_sub$species[e] == 'melp' & sampleInfo_sub$tissue[e] == 'FD') { sampleInfo_sub$shape[e] = 0}
  if(sampleInfo_sub$species[e] == 'melp' & sampleInfo_sub$tissue[e] == 'HA') { sampleInfo_sub$shape[e] = 5}
  if(sampleInfo_sub$species[e] == 'melp' & sampleInfo_sub$tissue[e] == 'HP') { sampleInfo_sub$shape[e] = 6}
}
atac_Rlog <- vst(dds)

pcaData <- plotPCA(atac_Rlog, intgroup = c("stage", "species", "tissue"), ntop = nrow(atac_Rlog), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=3, color=sampleInfo_sub$col, shape = sampleInfo_sub$shape, stroke=2, fill = sampleInfo_sub$col) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() 

###
# genome plot
###
library(berryFunctions)

chromTable <- read.table("chromosome_dem_start_end.txt", h=T)

sign_list_E <- list(sign_E_5th_FWHW, sign_E_D1_FWHW, sign_E_D2_FWHW)
sign_list_M <- list(sign_M_5th_FWHW, sign_M_D1_FWHW, sign_M_D2_FWHW)

colList <- c('#ff0000b5', '#009e73ff', '#0072b2ff')

sign_list_E_df <- list()
sign_list_E_df_all <- list()
for(i in 1:length(sign_list_E)){
  sign_list_E_df[[i]] <- as.data.frame(rownames(sign_list_E[[i]]))
  sign_list_E_df[[i]] <- separate(data = sign_list_E_df[[i]], col = 'rownames(sign_list_E[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE"), sep = "_")
  sign_list_E_df[[i]][,c(2:5,7,8)] <- sapply(sign_list_E_df[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_E_df_all[[i]] <- cbind(sign_list_E[[i]], sign_list_E_df[[i]])
}

sign_list_M_df <- list()
sign_list_M_df_all <- list()
for(i in 1:length(sign_list_M)){
  sign_list_M_df[[i]] <- as.data.frame(rownames(sign_list_M[[i]]))
  sign_list_M_df[[i]] <- separate(data = sign_list_M_df[[i]], col = 'rownames(sign_list_M[[i]])', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE"), sep = "_")
  sign_list_M_df[[i]][,c(2:5,7,8)] <- sapply(sign_list_M_df[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_M_df_all[[i]] <- cbind(sign_list_M[[i]], sign_list_M_df[[i]])
} 

####
# Plot fold change
####

layout(matrix(c(1:6), 2, 3, byrow = TRUE))
par(mar=c(4,2,2,2), oma=c(1,1,1,1))
# layout.show(6)
shared_E_5th <- subset(sign_list_E_df_all[[1]], sign_list_E_df_all[[1]]$start2 != -1)
unique_E_5th <- subset(sign_list_E_df_all[[1]], sign_list_E_df_all[[1]]$start2 == -1)

d_shared_E_5th <- density(abs(shared_E_5th$log2FoldChange))
d_unique_E_5th <- density(abs(unique_E_5th$log2FoldChange))

plot(d_shared_E_5th, xlim=c(1,6), ylim=c(0,2), main='5th erato', xlab='')
polygon(d_shared_E_5th, col=adjustcolor("red", alpha=0.2), border="red")
par(new=T)
plot(d_unique_E_5th, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='')


shared_E_D1 <- subset(sign_list_E_df_all[[2]], sign_list_E_df_all[[2]]$start2 != -1)
unique_E_D1 <- subset(sign_list_E_df_all[[2]], sign_list_E_df_all[[2]]$start2 == -1)

d_shared_E_D1 <- density(abs(shared_E_D1$log2FoldChange))
d_unique_E_D1 <- density(abs(unique_E_D1$log2FoldChange))

plot(d_shared_E_D1, xlim=c(1,6), ylim=c(0,2),  main='D1 erato', xlab='')
polygon(d_shared_E_D1, col=adjustcolor("blue", alpha=0.2), border="blue")
par(new=T)
plot(d_unique_E_D1, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='')


shared_E_D2 <- subset(sign_list_E_df_all[[3]], sign_list_E_df_all[[3]]$start2 != -1)
unique_E_D2 <- subset(sign_list_E_df_all[[3]], sign_list_E_df_all[[3]]$start2 == -1)

d_shared_E_D2 <- density(abs(shared_E_D2$log2FoldChange))
d_unique_E_D2 <- density(abs(unique_E_D2$log2FoldChange))

plot(d_shared_E_D2, xlim=c(1,6), ylim=c(0,2),  main='D2 erato', xlab='')
polygon(d_shared_E_D2, col=adjustcolor("darkgreen", alpha=0.2), border="darkgreen")
par(new=T)
plot(d_unique_E_D2, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='')



shared_M_5th <- subset(sign_list_M_df_all[[1]], is.na(sign_list_M_df_all[[1]]$start) == F)
unique_M_5th <- subset(sign_list_M_df_all[[1]], is.na(sign_list_M_df_all[[1]]$start) == T)

d_shared_M_5th <- density(abs(shared_M_5th$log2FoldChange))
d_unique_M_5th <- density(abs(unique_M_5th$log2FoldChange))

plot(d_shared_M_5th, xlim=c(1,6), ylim=c(0,2), main='5th melp', xlab='')
polygon(d_shared_M_5th, col=adjustcolor("red", alpha=0.2), border="red")
par(new=T)
plot(d_unique_M_5th, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='log2FC')


shared_M_D1 <- subset(sign_list_M_df_all[[2]], is.na(sign_list_M_df_all[[2]]$start) == F)
unique_M_D1 <- subset(sign_list_M_df_all[[2]], is.na(sign_list_M_df_all[[2]]$start) == T)

d_shared_M_D1 <- density(abs(shared_M_D1$log2FoldChange))
d_unique_M_D1 <- density(abs(unique_M_D1$log2FoldChange))

plot(d_shared_M_D1, xlim=c(1,6), ylim=c(0,2),  main='D1 melp', xlab='')
polygon(d_shared_M_D1, col=adjustcolor("blue", alpha=0.2), border="blue")
par(new=T)
plot(d_unique_M_D1, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='log2FC')


shared_M_D2 <- subset(sign_list_M_df_all[[3]], is.na(sign_list_M_df_all[[3]]$start) == F)
unique_M_D2 <- subset(sign_list_M_df_all[[3]], is.na(sign_list_M_df_all[[3]]$start) == T)

d_shared_M_D2 <- density(abs(shared_M_D2$log2FoldChange))
d_unique_M_D2 <- density(abs(unique_M_D2$log2FoldChange))

plot(d_shared_M_D2, xlim=c(1,6), ylim=c(0,2),  main='D2 melp', xlab='')
polygon(d_shared_M_D2, col=adjustcolor("darkgreen", alpha=0.2), border="darkgreen")
par(new=T)
plot(d_unique_M_D2, xlim=c(1,6), ylim=c(0,2), main='', lwd=2, xlab='log2FC')
################

write.table(shared_E_5th[,c("scafE","startE","endE")], file = "shared_FWHW_E_5th_p05fc1.txt", col.names = F, row.names = F, quote = F)
write.table(shared_E_D1[,c("scafE","startE","endE")], file = "shared_FWHW_E_D1_p05fc1.txt", col.names = F, row.names = F, quote = F)
write.table(shared_E_D2[,c("scafE","startE","endE")], file = "shared_FWHW_E_D2_p05fc1.txt", col.names = F, row.names = F, quote = F)

write.table(shared_M_5th[,c("scafE","startE","endE")], file = "shared_FWHW_M_5th_p05fc1.txt", col.names = F, row.names = F, quote = F)
write.table(shared_M_D1[,c("scafE","startE","endE")], file = "shared_FWHW_M_D1_p05fc1.txt", col.names = F, row.names = F, quote = F)
write.table(shared_M_D2[,c("scafE","startE","endE")], file = "shared_FWHW_M_D2_p05fc1.txt", col.names = F, row.names = F, quote = F)


# png(filename = 'chrom_5th.png', height = 800, width = 2400)

par(mar=c(1,2,1,1),oma=c(2,2,1,1))
plot(NULL, xlim = c(-1000000,max(chromTable$end_pan-chromTable$start_pan)), ylim = c(0,21), axes=F, ylab = '', xlab = '')
axis(seq(1,30,1))

for(e in 1:21){
  bot=e-0.7
  top=e-0.3
  roundedRect(0, bot, chromTable$end_pan[e]-chromTable$start_pan[e], top, rounding = 0.5, col = 'gray', border = ' gray')
  text(-1000000, e-0.5, e)
  
  for(p in 1:length(sign_list_E_df)){
    # for(p in 3){
    signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e])
    if(nrow(signs_E)>0){
      for(i in 1:nrow(signs_E)){
        rect(as.numeric(signs_E$start[i])-chromTable$start_pan[e], bot-0.15, as.numeric(signs_E$end[i])-chromTable$start_pan[e], top-0.2, col = colList[p], border = colList[p])
      }
    }
    
    signs_M <- subset(sign_list_M_df[[p]], sign_list_M_df[[p]]$start2 > chromTable$start_pan[e] & sign_list_M_df[[p]]$end2 < chromTable$end_pan[e])
    signs_M <- subset(signs_M, abs(signs_M$end2 - signs_M$start2) < 100000) 
    if(nrow(signs_M)>0){
      for(i in 1:nrow(signs_M)){
        rect(as.numeric(signs_M$start2[i])-chromTable$start_pan[e], bot+0.2, as.numeric(signs_M$end2[i])-chromTable$start_pan[e], top+0.15, col = colList[p], border = colList[p])
      }
    }
  }
}
mtext('chromosome',2)



########################
blocks <- 100000
for(n in 1:3){
  
  pos <- blocks
  pattern_E <- c()
  
  sign_file <- sign_list_E_df[[n]][order(sign_list_E_df[[n]]$start),]
  sign_file <- subset(sign_file, sign_file$end - sign_file$start < 100000) 
  
  while(pos <= 483353880){
    sign_file_sub <- subset(sign_file, sign_file$start >= pos - blocks & sign_file$end <= pos)
    count <- nrow(sign_file_sub)
    
    pattern_E <- rbind(pattern_E, c(pos - blocks, pos, count))
    
    pos = pos + blocks
  }
  
  print(nrow(pattern_E))
  if(n == 1){
    pattern_E_merge <- as.data.frame(pattern_E)
  }
  else{
    pattern_E_merge <- merge(pattern_E_merge, as.data.frame(pattern_E), by=c('V1','V2'), all.x=T, all.y=T)
  }
}
pattern_E_merge[is.na(pattern_E_merge)] <- 0
pattern_E_merge <- pattern_E_merge[order(pattern_E_merge$V1),]

for(n in 1:3){
  
  pos <- blocks
  pattern_M <- c()
  
  sign_file <- sign_list_M_df[[n]][order(sign_list_M_df[[n]]$start2),]
  sign_file <- subset(sign_file, abs(sign_file$end2 - sign_file$start2) < 100000) 
  
  while(pos <= 483353880){
    sign_file_sub <- subset(sign_file, sign_file$start2 >= pos - blocks & sign_file$end2 <= pos)
    count <- nrow(sign_file_sub)
    
    pattern_M <- rbind(pattern_M, c(pos - blocks, pos, count))
    
    pos = pos + blocks
  }
  
  print(nrow(pattern_M))
  if(n == 1){
    pattern_M_merge <- as.data.frame(pattern_M)
  }
  else{
    pattern_M_merge <- merge(pattern_M_merge, as.data.frame(pattern_M), by=c('V1','V2'), all.x=T, all.y=T)
  }
}
pattern_M_merge[is.na(pattern_M_merge)] <- 0
pattern_M_merge <- pattern_M_merge[order(pattern_M_merge$V1),]

####
####
colList <- c('hotpink2', 'orange', 'dodgerblue', 'gray')

layout(matrix(c(46:1,49,rep(48,21),47), 23, 3, byrow = FALSE), width=c(1,1,5))
layout.show(n=49)

par(mar=c(0.5,0,0.5,1),oma=c(0,0,0,0))

#all peaks
# plot(NULL, xlim=c(400,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
# axis(1,seq(400,0,by=-100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)
# for(e in 1:21){
#   
#   countpeaks_E <- c()
#   countpeaks_M <- c()
#   for(p in 1:length(sign_list_E_df)){
# 
#     signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e])
#     signs_M <- subset(sign_list_M_df[[p]], sign_list_M_df[[p]]$start2 > chromTable$start_pan[e] & sign_list_M_df[[p]]$end2 < chromTable$end_pan[e])
#     signs_M <- subset(signs_M, abs(signs_M$end2 - signs_M$start2) < 100000) 
#     
#     countpeaks_E <- c(countpeaks_E, nrow(signs_E))
#     countpeaks_M <- c(countpeaks_M, nrow(signs_M))
#   }
#   countpeaks <- as.table(cbind(countpeaks_E,countpeaks_M))
#   barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(400,0), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
# }

#shared peaks
plot(NULL, xlim=c(0,100), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(0,100,by=20),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

ptest <- c()
for(e in 1:21){
  
  countpeaks_E <- c()
  countpeaks_M <- c()
  for(p in 1:length(sign_list_E_df)){
    
    signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e] & sign_list_E_df[[p]]$start2 != -1)
    signs_M <- subset(sign_list_M_df[[p]], sign_list_M_df[[p]]$start2 > chromTable$start_pan[e] & sign_list_M_df[[p]]$end2 < chromTable$end_pan[e] & is.na(sign_list_M_df[[p]]$start) == FALSE)
    signs_M <- subset(signs_M, abs(signs_M$end2 - signs_M$start2) < 100000) 
    
    countpeaks_E <- c(countpeaks_E, nrow(signs_E))
    countpeaks_M <- c(countpeaks_M, nrow(signs_M))
  }
  countpeaks <- as.table(cbind(countpeaks_E,countpeaks_M))
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(0,100), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  # pval <- c()
  # for(k in 1:3) {
  #   pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
  #   pval <- c(pval, pstat$p.value)
  # }
  # ptest <- rbind(ptest, pval)
}

plot(NULL, xlim=c(0,100), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "shared")

#unique peaks
plot(NULL, xlim=c(100,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(100,0,by=-20),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

ptest <- c()
for(e in 1:21){
  
  countpeaks_E <- c()
  countpeaks_M <- c()
  for(p in 1:length(sign_list_E_df)){
    
    signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e] & sign_list_E_df[[p]]$start2 == -1)
    signs_M <- subset(sign_list_M_df[[p]], sign_list_M_df[[p]]$start2 > chromTable$start_pan[e] & sign_list_M_df[[p]]$end2 < chromTable$end_pan[e] & is.na(sign_list_M_df[[p]]$start) == TRUE)
    signs_M <- subset(signs_M, abs(signs_M$end2 - signs_M$start2) < 100000) 
    
    countpeaks_E <- c(countpeaks_E, nrow(signs_E))
    countpeaks_M <- c(countpeaks_M, nrow(signs_M))
  }
  countpeaks <- as.table(cbind(countpeaks_E,countpeaks_M))
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(100,0), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  # pval <- c()
  # for(k in 1:3) {
  #   pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
  #   pval <- c(pval, pstat$p.value)
  # }
  # ptest <- rbind(ptest, pval)
}

plot(NULL, xlim=c(100,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "unique")



plot(NULL, xlim = c(-2000000,max(chromTable$end_pan-chromTable$start_pan)+1000000), ylim = c(0,1), axes=F, ylab = '', xlab = '', xaxs="i", yaxs="i")
axis(1,seq(0,max(chromTable$end_pan-chromTable$start_pan),by=3000000),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)


par(mar=c(0,0,0,1))

plot(NULL, xlim = c(-2000000,max(chromTable$end_pan-chromTable$start_pan)+1000000), ylim = c(0,21), axes=F, ylab = '', xlab = '', xaxs="i", yaxs="i", xaxt='none', yaxt='none')

for(e in 1:21){
  bot=e-0.8
  top=e-0.2
  
  for(p in 1:length(sign_list_E_df)){
    # for(p in 3){
    signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start2 != -1 & sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e])
    if(nrow(signs_E)>0){
      for(i in 1:nrow(signs_E)){
        rect(as.numeric(signs_E$start[i])-chromTable$start_pan[e], top-0.35, as.numeric(signs_E$end[i])-chromTable$start_pan[e], bot+0.35, col = 'black', border = 'black')
      }
    }
  }
}

for(e in 1:21){
  bot=e-0.8
  top=e-0.2
  roundedRect(-100000, bot, chromTable$end_pan[e]-chromTable$start_pan[e]+100000, top-0.35, rounding = 0.5, col = 'gray', border = ' gray')
  roundedRect(-100000, bot+0.35, chromTable$end_pan[e]-chromTable$start_pan[e]+100000, top, rounding = 0.5, col = 'gray', border = ' gray')
  
  text(-1000000, e-0.5, e)
  text(-500000, e-0.75, 'e', cex=0.5)
  text(-500000, e-0.25, 'm', cex=0.5)
  
  signs_E <- subset(pattern_E_merge, pattern_E_merge$V1 >= chromTable$start_pan[e]-1 & pattern_E_merge$V2 < chromTable$end_pan[e])
  signs_M <- subset(pattern_M_merge, pattern_M_merge$V1 >= chromTable$start_pan[e]-1 & pattern_M_merge$V2 < chromTable$end_pan[e])
  
  for(i in 1:nrow(signs_E)){
    m <- as.matrix(signs_E[i,][,c(3:5)])
    if(nrow(which(m == max(m), arr.ind=T))>1){
      indx <- 4
      alphaV <- 1
    }
    else{
      indx <- which(m == max(m), arr.ind=T)[2]
      if(indx ==1) alphaV <- m[indx]/max(signs_E[,3])
      if(indx ==2) alphaV <- m[indx]/max(signs_E[,4])
      if(indx ==3) alphaV <- m[indx]/max(signs_E[,5])
    }
    
    rect(as.numeric(signs_E$V1[i])-chromTable$start_pan[e], bot, as.numeric(signs_E$V2[i])-chromTable$start_pan[e], top-0.35, col = adjustcolor(colList[indx], alphaV), border = adjustcolor(colList[indx], alphaV))
  }
  for(i in 1:nrow(signs_M)){
    m <- as.matrix(signs_M[i,][,c(3:5)])
    if(nrow(which(m == max(m), arr.ind=T))>1){
      indx <- 4
      alphaV <- 1
    }
    else{
      indx <- which(m == max(m), arr.ind=T)[2]
      if(indx ==1) alphaV <- m[indx]/max(signs_E[,3])
      if(indx ==2) alphaV <- m[indx]/max(signs_E[,4])
      if(indx ==3) alphaV <- m[indx]/max(signs_E[,5])
    }
    
    rect(as.numeric(signs_M$V1[i])-chromTable$start_pan[e], bot+0.35, as.numeric(signs_M$V2[i])-chromTable$start_pan[e], top, col = adjustcolor(colList[indx], alphaV), border = adjustcolor(colList[indx], alphaV))
  }
}












###
# shared melp/erato diff
###

shared_5th <- merge(sign_list_E_df[[1]],sign_list_M_df[[1]], by=c('start','end'))
shared_D1 <- merge(sign_list_E_df[[2]],sign_list_M_df[[2]], by=c('start','end'))
shared_D2 <- merge(sign_list_E_df[[3]],sign_list_M_df[[3]], by=c('start','end'))

nrow(shared_5th)
nrow(shared_D1)
nrow(shared_D2)

write.table(shared_5th[c(3,1,2)], file='shared_erato_5th_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(shared_D1[c(3,1,2)], file='shared_erato_D1_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(shared_D2[c(3,1,2)], file='shared_erato_D2_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')

write.table(shared_5th[c(3,4,5)], file='shared_melp_5th_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(shared_D1[c(3,4,5)], file='shared_melp_D1_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')
write.table(shared_D2[c(3,4,5)], file='shared_melp_D2_up.bed', quote = F, row.names = F, col.names = F, sep = '\t')



############ create tables

sign_list_E_dir <- list(sign_E_5th_FWHW_HWup, sign_E_D1_FWHW_HWup, sign_E_D2_FWHW_HWup, sign_E_5th_FWHW_FWup, sign_E_D1_FWHW_FWup, sign_E_D2_FWHW_FWup)
sign_list_M_dir <- list(sign_M_5th_FWHW_HWup, sign_M_D1_FWHW_HWup, sign_M_D2_FWHW_HWup, sign_M_5th_FWHW_FWup, sign_M_D1_FWHW_FWup, sign_M_D2_FWHW_FWup)

colList <- c('#ff0000b5', '#009e73ff', '#0072b2ff')

sign_list_E_df_dir <- list()
sign_list_E_df_all_dir <- list()
for(i in 1:length(sign_list_E_dir)){
  sign_list_E_df_dir[[i]] <- as.data.frame(rownames(sign_list_E_dir[[i]]))
  sign_list_E_df_dir[[i]] <- separate(data = sign_list_E_df_dir[[i]], col = 'rownames(sign_list_E_dir[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_E_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_E_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_E_df_all_dir[[i]] <- cbind(sign_list_E_dir[[i]], sign_list_E_df_dir[[i]])
}

sign_list_M_df_dir <- list()
sign_list_M_df_all_dir <- list()
for(i in 1:length(sign_list_M_dir)){
  sign_list_M_df_dir[[i]] <- as.data.frame(rownames(sign_list_M_dir[[i]]))
  sign_list_M_df_dir[[i]] <- separate(data = sign_list_M_df_dir[[i]], col = 'rownames(sign_list_M_dir[[i]])', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_M_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_M_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_M_df_all_dir[[i]] <- cbind(sign_list_M_dir[[i]], sign_list_M_df_dir[[i]])
} 

write.table(sign_list_E_df_all_dir[[1]][,c("genome","start","end")], file = '5th_p05FC1_erato_HWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[4]][,c("genome","start","end")], file = '5th_p05FC1_erato_FWup_pan.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("genome","start","end")], file = '5th_p05FC1_melp_HWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[4]][,c("genome","start","end")], file = '5th_p05FC1_melp_FWup_pan.bed', col.names = F, row.names = F, quote = F)




write.table(sign_list_E_df_dir[[1]][,c("genome","start","end")], file = "FWHW_5th_p05FC1_erato_HWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_dir[[2]][,c("genome","start","end")], file = "FWHW_D1_p05FC1_erato_HWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_dir[[3]][,c("genome","start","end")], file = "FWHW_D2_p05FC1_erato_HWup_PAN.txt", col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_dir[[1]][,c("genome","start","end")], file = "FWHW_5th_p05FC1_melp_HWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_dir[[2]][,c("genome","start","end")], file = "FWHW_D1_p05FC1_melp_HWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_dir[[3]][,c("genome","start","end")], file = "FWHW_D2_p05FC1_melp_HWup_PAN.txt", col.names = F, row.names = F, quote = F)

write.table(sign_list_E_df_dir[[4]][,c("genome","start","end")], file = "FWHW_5th_p05FC1_erato_FWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_dir[[5]][,c("genome","start","end")], file = "FWHW_D1_p05FC1_erato_FWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_dir[[6]][,c("genome","start","end")], file = "FWHW_D2_p05FC1_erato_FWup_PAN.txt", col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_dir[[4]][,c("genome","start","end")], file = "FWHW_5th_p05FC1_melp_FWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_dir[[5]][,c("genome","start","end")], file = "FWHW_D1_p05FC1_melp_FWup_PAN.txt", col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_dir[[6]][,c("genome","start","end")], file = "FWHW_D2_p05FC1_melp_FWup_PAN.txt", col.names = F, row.names = F, quote = F)

HW_D1_shared <- merge(sign_list_E_df_dir[[2]],sign_list_M_df_dir[[2]], by=c('start','end'))
HW_D2_shared <- merge(sign_list_E_df_dir[[3]],sign_list_M_df_dir[[3]], by=c('start','end'))

FW_5th_shared <- merge(sign_list_E_df_dir[[4]],sign_list_M_df_dir[[4]], by=c('start','end'))
FW_D1_shared <- merge(sign_list_E_df_dir[[5]],sign_list_M_df_dir[[5]], by=c('start','end'))
FW_D2_shared <- merge(sign_list_E_df_dir[[6]],sign_list_M_df_dir[[6]], by=c('start','end'))


write.table(HW_5th_shared[,c(6,7,8)], file = "FWHW_5th_p05FC1_shared_HWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D1_shared[,c(6,7,8)], file = "FWHW_D1_p05FC1_shared_HWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D2_shared[,c(6,7,8)], file = "FWHW_D2_p05FC1_shared_HWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)

write.table(FW_5th_shared[,c(6,7,8)], file = "FWHW_5th_p05FC1_shared_FWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D1_shared[,c(6,7,8)], file = "FWHW_D1_p05FC1_shared_FWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D2_shared[,c(6,7,8)], file = "FWHW_D2_p05FC1_shared_FWup_eratoCoords.bed", col.names = F, row.names = F, quote = F)


write.table(HW_5th_shared[,c(3,4,5)], file = "FWHW_5th_p05FC1_shared_HWup_panCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D1_shared[,c(3,4,5)], file = "FWHW_D1_p05FC1_shared_HWup_panCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D2_shared[,c(3,4,5)], file = "FWHW_D2_p05FC1_shared_HWup_panCoords.bed", col.names = F, row.names = F, quote = F)

write.table(FW_5th_shared[,c(3,4,5)], file = "FWHW_5th_p05FC1_shared_FWup_panCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D1_shared[,c(3,4,5)], file = "FWHW_D1_p05FC1_shared_FWup_panCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D2_shared[,c(3,4,5)], file = "FWHW_D2_p05FC1_shared_FWup_panCoords.bed", col.names = F, row.names = F, quote = F)


write.table(HW_5th_shared[,c(9,10,11)], file = "FWHW_5th_p05FC1_shared_HWup_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D1_shared[,c(9,10,11)], file = "FWHW_D1_p05FC1_shared_HWup_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(HW_D2_shared[,c(9,10,11)], file = "FWHW_D2_p05FC1_shared_HWup_melpCoords.bed", col.names = F, row.names = F, quote = F)

write.table(FW_5th_shared[,c(9,10,11)], file = "FWHW_5th_p05FC1_shared_FWup_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D1_shared[,c(9,10,11)], file = "FWHW_D1_p05FC1_shared_FWup_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(FW_D2_shared[,c(9,10,11)], file = "FWHW_D2_p05FC1_shared_FWup_melpCoords.bed", col.names = F, row.names = F, quote = F)



write.table(subset(sign_list_E_df_dir[[1]], sign_list_E_df_dir[[1]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_5th_p05FC1_eratoUnique_HWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[2]], sign_list_E_df_dir[[2]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_D1_p05FC1_eratoUnique_HWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[3]], sign_list_E_df_dir[[3]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_D2_p05FC1_eratoUnique_HWup.txt", col.names = F, row.names = F, quote = F)

write.table(subset(sign_list_M_df_dir[[1]], sign_list_M_df_dir[[1]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_5th_p05FC1_melpUnique_HWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_dir[[2]], sign_list_M_df_dir[[2]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_D1_p05FC1_melpUnique_HWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_dir[[3]], sign_list_M_df_dir[[3]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_D2_p05FC1_melpUnique_HWup.txt", col.names = F, row.names = F, quote = F)

write.table(subset(sign_list_E_df_dir[[4]], sign_list_E_df_dir[[4]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_5th_p05FC1_eratoUnique_FWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[5]], sign_list_E_df_dir[[5]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_D1_p05FC1_eratoUnique_FWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[6]], sign_list_E_df_dir[[6]]$scafM == "NA")[,c("scafE","startE","endE")], file = "FWHW_D2_p05FC1_eratoUnique_FWup.txt", col.names = F, row.names = F, quote = F)

write.table(subset(sign_list_M_df_dir[[4]], sign_list_M_df_dir[[4]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_5th_p05FC1_melpUnique_FWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_dir[[5]], sign_list_M_df_dir[[5]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_D1_p05FC1_melpUnique_FWup.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_dir[[6]], sign_list_M_df_dir[[6]]$scafE == "NA")[,c("scafM","startM","endM")], file = "FWHW_D2_p05FC1_melpUnique_FWup.txt", col.names = F, row.names = F, quote = F)

## not strict shared (heterochrony)
write.table(subset(sign_list_E_df_dir[[1]], sign_list_E_df_dir[[1]]$scafM != "NA")[,c("genome","start","end")], file = "FWHW_5th_p05FC1_shared_HWup_notStrict.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[2]], sign_list_E_df_dir[[2]]$scafM != "NA")[,c("genome","start","end")], file = "FWHW_D1_p05FC1_shared_HWup_notStrict.txt", col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_dir[[3]], sign_list_E_df_dir[[3]]$scafM != "NA")[,c("genome","start","end")], file = "FWHW_D2_p05FC1_shared_HWup_notStrict.txt", col.names = F, row.names = F, quote = F)





##############################
# FC correlations


conserved_E <- read.table('ALL_ERA_50_pan.merged.sort.clean.bed', h=F)
conserved_M <- read.table('ALL_MELP_50_pan.merged.sort.clean.bed', h=F)

colnames(conserved_E) <- c('chrom','start','end')
colnames(conserved_M) <- c('chrom','start','end')

conserved_S <- bed_intersect(conserved_E,conserved_M)
conserved_S <- subset(conserved_S, conserved_S$.overlap/(conserved_S$end.x-conserved_S$start.x) >= 0.5 & conserved_S$.overlap/(conserved_S$end.y-conserved_S$start.y) >= 0.5)

conserved_S <- as.data.frame(conserved_S[,c(1:3)])
colnames(conserved_S) <- c('chrom','start','end')

library(plotrix)
library(valr)


load('sign_E_5th_FWHW_p05fc1_ALL.rda')
load('sign_M_5th_FWHW_p05fc1_ALL.rda')

load('sign_E_D1_FWHW_p05fc1_ALL.rda')
load('sign_M_D1_FWHW_p05fc1_ALL.rda')

load('sign_E_D2_FWHW_p05fc1_ALL.rda')
load('sign_M_D2_FWHW_p05fc1_ALL.rda')

res_E_5thFWHW
res_M_5thFWHW

sign_E_5th_FWHW_FWup <- subset(res_E_5thFWHW, res_E_5thFWHW$padj < 0.05 & res_E_5thFWHW$log2FoldChange < -1)
sign_M_5th_FWHW_FWup <- subset(res_M_5thFWHW, res_M_5thFWHW$padj < 0.05 & res_M_5thFWHW$log2FoldChange < -1)

sign_E_5th_FWHW_HWup <- subset(res_E_5thFWHW, res_E_5thFWHW$padj < 0.05 & res_E_5thFWHW$log2FoldChange > 1)
sign_M_5th_FWHW_HWup <- subset(res_M_5thFWHW, res_M_5thFWHW$padj < 0.05 & res_M_5thFWHW$log2FoldChange > 1)


sign_E_D1_FWHW_FWup <- subset(res_E_D1FWHW, res_E_D1FWHW$padj < 0.05 & res_E_D1FWHW$log2FoldChange < -1)
sign_M_D1_FWHW_FWup <- subset(res_M_D1FWHW, res_M_D1FWHW$padj < 0.05 & res_M_D1FWHW$log2FoldChange < -1)

sign_E_D1_FWHW_HWup <- subset(res_E_D1FWHW, res_E_D1FWHW$padj < 0.05 & res_E_D1FWHW$log2FoldChange > 1)
sign_M_D1_FWHW_HWup <- subset(res_M_D1FWHW, res_M_D1FWHW$padj < 0.05 & res_M_D1FWHW$log2FoldChange > 1)


sign_E_D2_FWHW_FWup <- subset(res_E_D2FWHW, res_E_D2FWHW$padj < 0.05 & res_E_D2FWHW$log2FoldChange < -1)
sign_M_D2_FWHW_FWup <- subset(res_M_D2FWHW, res_M_D2FWHW$padj < 0.05 & res_M_D2FWHW$log2FoldChange < -1)

sign_E_D2_FWHW_HWup <- subset(res_E_D2FWHW, res_E_D2FWHW$padj < 0.05 & res_E_D2FWHW$log2FoldChange > 1)
sign_M_D2_FWHW_HWup <- subset(res_M_D2FWHW, res_M_D2FWHW$padj < 0.05 & res_M_D2FWHW$log2FoldChange > 1)

# 5th FWup

merged_EM_5thFWup <- merge(as.data.frame(res_E_5thFWHW), as.data.frame(sign_M_5th_FWHW_FWup), by = 'row.names')
merged_ME_5thFWup <- merge(as.data.frame(res_M_5thFWHW), as.data.frame(sign_E_5th_FWHW_FWup), by = 'row.names')

nrow(sign_E_5th_FWHW_FWup)-nrow(merged_ME_5thFWup[grep("NA", merged_ME_5thFWup$Row.names), ])
nrow(sign_M_5th_FWHW_FWup)-nrow(merged_EM_5thFWup[grep("NA", merged_EM_5thFWup$Row.names), ])

merged_shared_5thFWup <- merge(as.data.frame(sign_E_5th_FWHW_FWup), as.data.frame(sign_M_5th_FWHW_FWup), by = 'row.names')

merged_EM_5thFWup$log2FoldChange.x <- merged_EM_5thFWup$log2FoldChange.x*-1
merged_EM_5thFWup$log2FoldChange.y <- merged_EM_5thFWup$log2FoldChange.y*-1
merged_ME_5thFWup$log2FoldChange.x <- merged_ME_5thFWup$log2FoldChange.x*-1
merged_ME_5thFWup$log2FoldChange.y <- merged_ME_5thFWup$log2FoldChange.y*-1
merged_shared_5thFWup$log2FoldChange.x <- merged_shared_5thFWup$log2FoldChange.x*-1
merged_shared_5thFWup$log2FoldChange.y <- merged_shared_5thFWup$log2FoldChange.y*-1

merged_EM_5thFWup_unique <- merged_EM_5thFWup[!merged_EM_5thFWup$Row.names %in% merged_shared_5thFWup$Row.names, ]
merged_ME_5thFWup_unique <- merged_ME_5thFWup[!merged_ME_5thFWup$Row.names %in% merged_shared_5thFWup$Row.names, ]

par(mfrow=c(3,2))

plot(merged_EM_5thFWup_unique$log2FoldChange.x, merged_EM_5thFWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "5th instar FW up")
par(new=T)
plot(merged_ME_5thFWup_unique$log2FoldChange.x, merged_ME_5thFWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_5thFWup$log2FoldChange.x, merged_shared_5thFWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')


ablineclip(lm(merged_EM_5thFWup_unique$log2FoldChange.y ~ merged_EM_5thFWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_5thFWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_5thFWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_5thFWup_unique$log2FoldChange.y ~ merged_ME_5thFWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_5thFWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_5thFWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_5thFWup$log2FoldChange.y ~ merged_shared_5thFWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_5thFWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_5thFWup$log2FoldChange.x)))

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")

######################## % completely missing

sign_E_5th_FWHW_FWup$Row.names <- rownames(sign_E_5th_FWHW_FWup)
sign_E_5th_FWHW_FWup <- separate(data = as.data.frame(sign_E_5th_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_5th_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_E_5th_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_E_5th_FWHW_FWup$Genome <- 'pan'

write.table(sign_E_5th_FWHW_FWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

sign_M_5th_FWHW_FWup$Row.names <- rownames(sign_M_5th_FWHW_FWup)
sign_M_5th_FWHW_FWup <- separate(data = as.data.frame(sign_M_5th_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_M_5th_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_M_5th_FWHW_FWup$Genome <- 'pan'


sign_M_5th_FWHW_FWup_write <- sign_M_5th_FWHW_FWup[,c(18,10,11)]
sign_M_5th_FWHW_FWup_write <- subset(sign_M_5th_FWHW_FWup_write, as.numeric(sign_M_5th_FWHW_FWup_write$end2)-as.numeric(sign_M_5th_FWHW_FWup_write$start2) < 10000)

write.table(sign_M_5th_FWHW_FWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_5thFWup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_FWHW_FWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_5th_FWHW_FWup_bed <- sign_E_5th_FWHW_FWup[,c(18,8,9)]
colnames(sign_E_5th_FWHW_FWup_bed) <- c('chrom','start','end')
sign_E_5th_FWHW_FWup_bed$start <- as.numeric(sign_E_5th_FWHW_FWup_bed$start)
sign_E_5th_FWHW_FWup_bed$end <- as.numeric(sign_E_5th_FWHW_FWup_bed$end)

sign_M_5th_FWHW_FWup_bed <- sign_M_5th_FWHW_FWup_write
colnames(sign_M_5th_FWHW_FWup_bed) <- c('chrom','start','end')
sign_M_5th_FWHW_FWup_bed$start <- as.numeric(sign_M_5th_FWHW_FWup_bed$start)
sign_M_5th_FWHW_FWup_bed$end <- as.numeric(sign_M_5th_FWHW_FWup_bed$end)

merged_shared_5thFWup_bed <- merged_shared_5thFWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_5thFWup_bed) <- c('chrom','start','end')
merged_shared_5thFWup_bed$start <- as.numeric(merged_shared_5thFWup_bed$start)
merged_shared_5thFWup_bed$end <- as.numeric(merged_shared_5thFWup_bed$end)

sign_E_5th_FWHW_FWup_bed_match <- bed_intersect(sign_E_5th_FWHW_FWup_bed,conserved_E)
sign_M_5th_FWHW_FWup_bed_match <- bed_intersect(sign_M_5th_FWHW_FWup_bed,conserved_M)
merged_shared_5thFWup_bed_match <- bed_intersect(merged_shared_5thFWup_bed,conserved_S)


sign_E_5th_FWHW_FWup_cons <- merge(sign_E_5th_FWHW_FWup, sign_E_5th_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_5th_FWHW_FWup_cons <- merge(sign_M_5th_FWHW_FWup, sign_M_5th_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_5thFWup_cons <- merge(merged_ME_5thFWup, sign_E_5th_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_5thFWup_cons <- merge(merged_EM_5thFWup, sign_M_5th_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_5th_FWHW_FWup_cons)-nrow(merged_ME_5thFWup_cons[grep("NA", merged_ME_5thFWup_cons$Row.names), ])
nrow(sign_M_5th_FWHW_FWup_cons)-nrow(merged_EM_5thFWup_cons[grep("NA", merged_EM_5thFWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_5th_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_5th_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_5thFWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_5th_FWHW_FWup_bed_match)
nrow(sign_M_5th_FWHW_FWup_bed_match)
nrow(merged_shared_5thFWup_bed_match)

nrow(merged_ME_5thFWup_cons[grep("NA", merged_ME_5thFWup_cons$Row.names), ])
nrow(merged_EM_5thFWup_cons[grep("NA", merged_EM_5thFWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################



# 5th HWup

merged_EM_5thHWup <- merge(as.data.frame(res_E_5thFWHW), as.data.frame(sign_M_5th_FWHW_HWup), by = 'row.names')
merged_ME_5thHWup <- merge(as.data.frame(res_M_5thFWHW), as.data.frame(sign_E_5th_FWHW_HWup), by = 'row.names')

nrow(sign_E_5th_FWHW_HWup)-nrow(merged_ME_5thHWup[grep("NA", merged_ME_5thHWup$Row.names), ])
nrow(sign_M_5th_FWHW_HWup)-nrow(merged_EM_5thHWup[grep("NA", merged_EM_5thHWup$Row.names), ])

merged_shared_5thHWup <- merge(as.data.frame(sign_E_5th_FWHW_HWup), as.data.frame(sign_M_5th_FWHW_HWup), by = 'row.names')

merged_EM_5thHWup_unique <- merged_EM_5thFWup[!merged_EM_5thHWup$Row.names %in% merged_shared_5thHWup$Row.names, ]
merged_ME_5thHWup_unique <- merged_ME_5thFWup[!merged_ME_5thHWup$Row.names %in% merged_shared_5thHWup$Row.names, ]

plot(merged_EM_5thHWup_unique$log2FoldChange.x, merged_EM_5thHWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "5th instar HW up")
par(new=T)
plot(merged_ME_5thHWup_unique$log2FoldChange.x, merged_ME_5thHWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_5thHWup$log2FoldChange.x, merged_shared_5thHWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                    "significant in both (shared)"),
      col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")


ablineclip(lm(merged_EM_5thHWup_unique$log2FoldChange.y ~ merged_EM_5thHWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_5thHWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_5thHWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_5thHWup_unique$log2FoldChange.y ~ merged_ME_5thHWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_5thHWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_5thHWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_5thHWup$log2FoldChange.y ~ merged_shared_5thHWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_5thHWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_5thHWup$log2FoldChange.x)))

######################## % completely missing

sign_E_5th_FWHW_HWup$Row.names <- rownames(sign_E_5th_FWHW_HWup)
sign_E_5th_FWHW_HWup <- separate(data = as.data.frame(sign_E_5th_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_5th_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_E_5th_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_E_5th_FWHW_HWup$Genome <- 'pan'

write.table(sign_E_5th_FWHW_HWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

sign_M_5th_FWHW_HWup$Row.names <- rownames(sign_M_5th_FWHW_HWup)
sign_M_5th_FWHW_HWup <- separate(data = as.data.frame(sign_M_5th_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_M_5th_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_M_5th_FWHW_HWup$Genome <- 'pan'


sign_M_5th_FWHW_HWup_write <- sign_M_5th_FWHW_HWup[,c(18,10,11)]
sign_M_5th_FWHW_HWup_write <- subset(sign_M_5th_FWHW_HWup_write, as.numeric(sign_M_5th_FWHW_HWup_write$end2)-as.numeric(sign_M_5th_FWHW_HWup_write$start2) < 10000)

write.table(sign_M_5th_FWHW_HWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_5thHWup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_FWHW_HWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_5th_FWHW_HWup_bed <- sign_E_5th_FWHW_HWup[,c(18,8,9)]
colnames(sign_E_5th_FWHW_HWup_bed) <- c('chrom','start','end')
sign_E_5th_FWHW_HWup_bed$start <- as.numeric(sign_E_5th_FWHW_HWup_bed$start)
sign_E_5th_FWHW_HWup_bed$end <- as.numeric(sign_E_5th_FWHW_HWup_bed$end)

sign_M_5th_FWHW_HWup_bed <- sign_M_5th_FWHW_HWup_write
colnames(sign_M_5th_FWHW_HWup_bed) <- c('chrom','start','end')
sign_M_5th_FWHW_HWup_bed$start <- as.numeric(sign_M_5th_FWHW_HWup_bed$start)
sign_M_5th_FWHW_HWup_bed$end <- as.numeric(sign_M_5th_FWHW_HWup_bed$end)

merged_shared_5thHWup_bed <- merged_shared_5thHWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_5thHWup_bed) <- c('chrom','start','end')
merged_shared_5thHWup_bed$start <- as.numeric(merged_shared_5thHWup_bed$start)
merged_shared_5thHWup_bed$end <- as.numeric(merged_shared_5thHWup_bed$end)

sign_E_5th_FWHW_HWup_bed_match <- bed_intersect(sign_E_5th_FWHW_HWup_bed,conserved_E)
sign_M_5th_FWHW_HWup_bed_match <- bed_intersect(sign_M_5th_FWHW_HWup_bed,conserved_M)
merged_shared_5thHWup_bed_match <- bed_intersect(merged_shared_5thHWup_bed,conserved_S)


sign_E_5th_FWHW_HWup_cons <- merge(sign_E_5th_FWHW_HWup, sign_E_5th_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_5th_FWHW_HWup_cons <- merge(sign_M_5th_FWHW_HWup, sign_M_5th_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_5thHWup_cons <- merge(merged_ME_5thHWup, sign_E_5th_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_5thHWup_cons <- merge(merged_EM_5thHWup, sign_M_5th_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_5th_FWHW_HWup_cons)-nrow(merged_ME_5thHWup_cons[grep("NA", merged_ME_5thHWup_cons$Row.names), ])
nrow(sign_M_5th_FWHW_HWup_cons)-nrow(merged_EM_5thHWup_cons[grep("NA", merged_EM_5thHWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_5th_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_5th_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_5thHWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_5th_FWHW_HWup_bed_match)
nrow(sign_M_5th_FWHW_HWup_bed_match)
nrow(merged_shared_5thHWup_bed_match)

nrow(merged_ME_5thHWup_cons[grep("NA", merged_ME_5thHWup_cons$Row.names), ])
nrow(merged_EM_5thHWup_cons[grep("NA", merged_EM_5thHWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################     

# D1 FWup

merged_EM_D1FWup <- merge(as.data.frame(res_E_D1FWHW), as.data.frame(sign_M_D1_FWHW_FWup), by = 'row.names')
merged_ME_D1FWup <- merge(as.data.frame(res_M_D1FWHW), as.data.frame(sign_E_D1_FWHW_FWup), by = 'row.names')

nrow(sign_E_D1_FWHW_FWup)-nrow(merged_ME_D1FWup[grep("NA", merged_ME_D1FWup$Row.names), ])
nrow(sign_M_D1_FWHW_FWup)-nrow(merged_EM_D1FWup[grep("NA", merged_EM_D1FWup$Row.names), ])

merged_shared_D1FWup <- merge(as.data.frame(sign_E_D1_FWHW_FWup), as.data.frame(sign_M_D1_FWHW_FWup), by = 'row.names')

merged_EM_D1FWup$log2FoldChange.x <- merged_EM_D1FWup$log2FoldChange.x*-1
merged_EM_D1FWup$log2FoldChange.y <- merged_EM_D1FWup$log2FoldChange.y*-1
merged_ME_D1FWup$log2FoldChange.x <- merged_ME_D1FWup$log2FoldChange.x*-1
merged_ME_D1FWup$log2FoldChange.y <- merged_ME_D1FWup$log2FoldChange.y*-1
merged_shared_D1FWup$log2FoldChange.x <- merged_shared_D1FWup$log2FoldChange.x*-1
merged_shared_D1FWup$log2FoldChange.y <- merged_shared_D1FWup$log2FoldChange.y*-1

merged_EM_D1FWup_unique <- merged_EM_D1FWup[!merged_EM_D1FWup$Row.names %in% merged_shared_D1FWup$Row.names, ]
merged_ME_D1FWup_unique <- merged_ME_D1FWup[!merged_ME_D1FWup$Row.names %in% merged_shared_D1FWup$Row.names, ]


plot(merged_EM_D1FWup_unique$log2FoldChange.x, merged_EM_D1FWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h pupa FW up")
par(new=T)
plot(merged_ME_D1FWup_unique$log2FoldChange.x, merged_ME_D1FWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D1FWup$log2FoldChange.x, merged_shared_D1FWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')


ablineclip(lm(merged_EM_D1FWup_unique$log2FoldChange.y ~ merged_EM_D1FWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D1FWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D1FWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D1FWup_unique$log2FoldChange.y ~ merged_ME_D1FWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D1FWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D1FWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D1FWup$log2FoldChange.y ~ merged_shared_D1FWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D1FWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D1FWup$log2FoldChange.x)))

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")


######################## % completely missing

sign_E_D1_FWHW_FWup$Row.names <- rownames(sign_E_D1_FWHW_FWup)
sign_E_D1_FWHW_FWup <- separate(data = as.data.frame(sign_E_D1_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_E_D1_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FWHW_FWup$Genome <- 'pan'

write.table(sign_E_D1_FWHW_FWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FWHW_FWup$Row.names <- rownames(sign_M_D1_FWHW_FWup)
sign_M_D1_FWHW_FWup <- separate(data = as.data.frame(sign_M_D1_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_M_D1_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FWHW_FWup$Genome <- 'pan'

write.table(sign_M_D1_FWHW_FWup[,c(18,10,11)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FWHW_FWup_write <- sign_M_D1_FWHW_FWup[,c(18,10,11)]
sign_M_D1_FWHW_FWup_write <- subset(sign_M_D1_FWHW_FWup_write, as.numeric(sign_M_D1_FWHW_FWup_write$end2)-as.numeric(sign_M_D1_FWHW_FWup_write$start2) < 10000)

write.table(sign_M_D1_FWHW_FWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_FWHW_FWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FWHW_FWup_bed <- sign_E_D1_FWHW_FWup[,c(18,8,9)]
colnames(sign_E_D1_FWHW_FWup_bed) <- c('chrom','start','end')
sign_E_D1_FWHW_FWup_bed$start <- as.numeric(sign_E_D1_FWHW_FWup_bed$start)
sign_E_D1_FWHW_FWup_bed$end <- as.numeric(sign_E_D1_FWHW_FWup_bed$end)

sign_M_D1_FWHW_FWup_bed <- sign_M_D1_FWHW_FWup_write
colnames(sign_M_D1_FWHW_FWup_bed) <- c('chrom','start','end')
sign_M_D1_FWHW_FWup_bed$start <- as.numeric(sign_M_D1_FWHW_FWup_bed$start)
sign_M_D1_FWHW_FWup_bed$end <- as.numeric(sign_M_D1_FWHW_FWup_bed$end)

merged_shared_D1FWup_bed <- merged_shared_D1FWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D1FWup_bed) <- c('chrom','start','end')
merged_shared_D1FWup_bed$start <- as.numeric(merged_shared_D1FWup_bed$start)
merged_shared_D1FWup_bed$end <- as.numeric(merged_shared_D1FWup_bed$end)

sign_E_D1_FWHW_FWup_bed_match <- bed_intersect(sign_E_D1_FWHW_FWup_bed,conserved_E)
sign_M_D1_FWHW_FWup_bed_match <- bed_intersect(sign_M_D1_FWHW_FWup_bed,conserved_M)
merged_shared_D1FWup_bed_match <- bed_intersect(merged_shared_D1FWup_bed,conserved_S)


sign_E_D1_FWHW_FWup_cons <- merge(sign_E_D1_FWHW_FWup, sign_E_D1_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FWHW_FWup_cons <- merge(sign_M_D1_FWHW_FWup, sign_M_D1_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D1FWup_cons <- merge(merged_ME_D1FWup, sign_E_D1_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D1FWup_cons <- merge(merged_EM_D1FWup, sign_M_D1_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D1FWup_cons[grep("NA", merged_ME_D1FWup_cons$Row.names), ])
nrow(merged_EM_D1FWup_cons[grep("NA", merged_EM_D1FWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D1FWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FWHW_FWup_bed_match)
nrow(sign_M_D1_FWHW_FWup_bed_match)
nrow(merged_shared_D1FWup_bed_match)

nrow(merged_ME_D1FWup_cons[grep("NA", merged_ME_D1FWup_cons$Row.names), ])
nrow(merged_EM_D1FWup_cons[grep("NA", merged_EM_D1FWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################   

# D1 HWup

merged_EM_D1HWup <- merge(as.data.frame(res_E_D1FWHW), as.data.frame(sign_M_D1_FWHW_HWup), by = 'row.names')
merged_ME_D1HWup <- merge(as.data.frame(res_M_D1FWHW), as.data.frame(sign_E_D1_FWHW_HWup), by = 'row.names')

nrow(sign_E_D1_FWHW_HWup)-nrow(merged_ME_D1HWup[grep("NA", merged_ME_D1HWup$Row.names), ])
nrow(sign_M_D1_FWHW_HWup)-nrow(merged_EM_D1HWup[grep("NA", merged_EM_D1HWup$Row.names), ])

merged_shared_D1HWup <- merge(as.data.frame(sign_E_D1_FWHW_HWup), as.data.frame(sign_M_D1_FWHW_HWup), by = 'row.names')

merged_EM_D1HWup_unique <- merged_EM_D1FWup[!merged_EM_D1HWup$Row.names %in% merged_shared_D1HWup$Row.names, ]
merged_ME_D1HWup_unique <- merged_ME_D1FWup[!merged_ME_D1HWup$Row.names %in% merged_shared_D1HWup$Row.names, ]

plot(merged_EM_D1HWup_unique$log2FoldChange.x, merged_EM_D1HWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h pupa HW up")
par(new=T)
plot(merged_ME_D1HWup_unique$log2FoldChange.x, merged_ME_D1HWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D1HWup$log2FoldChange.x, merged_shared_D1HWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")


ablineclip(lm(merged_EM_D1HWup_unique$log2FoldChange.y ~ merged_EM_D1HWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D1HWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D1HWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D1HWup_unique$log2FoldChange.y ~ merged_ME_D1HWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D1HWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D1HWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D1HWup$log2FoldChange.y ~ merged_shared_D1HWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D1HWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D1HWup$log2FoldChange.x)))


######################## % completely missing

sign_E_D1_FWHW_HWup$Row.names <- rownames(sign_E_D1_FWHW_HWup)
sign_E_D1_FWHW_HWup <- separate(data = as.data.frame(sign_E_D1_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_E_D1_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FWHW_HWup$Genome <- 'pan'

write.table(sign_E_D1_FWHW_HWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FWHW_HWup$Row.names <- rownames(sign_M_D1_FWHW_HWup)
sign_M_D1_FWHW_HWup <- separate(data = as.data.frame(sign_M_D1_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_M_D1_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FWHW_HWup$Genome <- 'pan'

write.table(sign_M_D1_FWHW_HWup[,c(18,10,11)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FWHW_HWup_write <- sign_M_D1_FWHW_HWup[,c(18,10,11)]
sign_M_D1_FWHW_HWup_write <- subset(sign_M_D1_FWHW_HWup_write, as.numeric(sign_M_D1_FWHW_HWup_write$end2)-as.numeric(sign_M_D1_FWHW_HWup_write$start2) < 10000)


write.table(sign_M_D1_FWHW_HWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_FWHW_HWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FWHW_HWup_bed <- sign_E_D1_FWHW_HWup[,c(18,8,9)]
colnames(sign_E_D1_FWHW_HWup_bed) <- c('chrom','start','end')
sign_E_D1_FWHW_HWup_bed$start <- as.numeric(sign_E_D1_FWHW_HWup_bed$start)
sign_E_D1_FWHW_HWup_bed$end <- as.numeric(sign_E_D1_FWHW_HWup_bed$end)

sign_M_D1_FWHW_HWup_bed <- sign_M_D1_FWHW_HWup_write
colnames(sign_M_D1_FWHW_HWup_bed) <- c('chrom','start','end')
sign_M_D1_FWHW_HWup_bed$start <- as.numeric(sign_M_D1_FWHW_HWup_bed$start)
sign_M_D1_FWHW_HWup_bed$end <- as.numeric(sign_M_D1_FWHW_HWup_bed$end)

merged_shared_D1HWup_bed <- merged_shared_D1HWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D1HWup_bed) <- c('chrom','start','end')
merged_shared_D1HWup_bed$start <- as.numeric(merged_shared_D1HWup_bed$start)
merged_shared_D1HWup_bed$end <- as.numeric(merged_shared_D1HWup_bed$end)

sign_E_D1_FWHW_HWup_bed_match <- bed_intersect(sign_E_D1_FWHW_HWup_bed,conserved_E)
sign_M_D1_FWHW_HWup_bed_match <- bed_intersect(sign_M_D1_FWHW_HWup_bed,conserved_M)
merged_shared_D1HWup_bed_match <- bed_intersect(merged_shared_D1HWup_bed,conserved_S)


sign_E_D1_FWHW_HWup_cons <- merge(sign_E_D1_FWHW_HWup, sign_E_D1_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FWHW_HWup_cons <- merge(sign_M_D1_FWHW_HWup, sign_M_D1_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D1HWup_cons <- merge(merged_ME_D1HWup, sign_E_D1_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D1HWup_cons <- merge(merged_EM_D1HWup, sign_M_D1_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D1HWup_cons[grep("NA", merged_ME_D1HWup_cons$Row.names), ])
nrow(merged_EM_D1HWup_cons[grep("NA", merged_EM_D1HWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D1HWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FWHW_HWup_bed_match)
nrow(sign_M_D1_FWHW_HWup_bed_match)
nrow(merged_shared_D1HWup_bed_match)

nrow(merged_ME_D1HWup_cons[grep("NA", merged_ME_D1HWup_cons$Row.names), ])
nrow(merged_EM_D1HWup_cons[grep("NA", merged_EM_D1HWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################   

# D2 FWup

merged_EM_D2FWup <- merge(as.data.frame(res_E_D2FWHW), as.data.frame(sign_M_D2_FWHW_FWup), by = 'row.names')
merged_ME_D2FWup <- merge(as.data.frame(res_M_D2FWHW), as.data.frame(sign_E_D2_FWHW_FWup), by = 'row.names')

nrow(sign_E_D2_FWHW_FWup)-nrow(merged_ME_D2FWup[grep("NA", merged_ME_D2FWup$Row.names), ])
nrow(sign_M_D2_FWHW_FWup)-nrow(merged_EM_D2FWup[grep("NA", merged_EM_D2FWup$Row.names), ])

merged_shared_D2FWup <- merge(as.data.frame(sign_E_D2_FWHW_FWup), as.data.frame(sign_M_D2_FWHW_FWup), by = 'row.names')

merged_EM_D2FWup$log2FoldChange.x <- merged_EM_D2FWup$log2FoldChange.x*-1
merged_EM_D2FWup$log2FoldChange.y <- merged_EM_D2FWup$log2FoldChange.y*-1
merged_ME_D2FWup$log2FoldChange.x <- merged_ME_D2FWup$log2FoldChange.x*-1
merged_ME_D2FWup$log2FoldChange.y <- merged_ME_D2FWup$log2FoldChange.y*-1
merged_shared_D2FWup$log2FoldChange.x <- merged_shared_D2FWup$log2FoldChange.x*-1
merged_shared_D2FWup$log2FoldChange.y <- merged_shared_D2FWup$log2FoldChange.y*-1

merged_EM_D2FWup_unique <- merged_EM_D2FWup[!merged_EM_D2FWup$Row.names %in% merged_shared_D2FWup$Row.names, ]
merged_ME_D2FWup_unique <- merged_ME_D2FWup[!merged_ME_D2FWup$Row.names %in% merged_shared_D2FWup$Row.names, ]



plot(merged_EM_D2FWup_unique$log2FoldChange.x, merged_EM_D2FWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "60h pupa FW up")
par(new=T)
plot(merged_ME_D2FWup_unique$log2FoldChange.x, merged_ME_D2FWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D2FWup$log2FoldChange.x, merged_shared_D2FWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')


ablineclip(lm(merged_EM_D2FWup_unique$log2FoldChange.y ~ merged_EM_D2FWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D2FWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D2FWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D2FWup_unique$log2FoldChange.y ~ merged_ME_D2FWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D2FWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D2FWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D2FWup$log2FoldChange.y ~ merged_shared_D2FWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D2FWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D2FWup$log2FoldChange.x)))

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")

######################## % completely missing

sign_E_D2_FWHW_FWup$Row.names <- rownames(sign_E_D2_FWHW_FWup)
sign_E_D2_FWHW_FWup <- separate(data = as.data.frame(sign_E_D2_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_E_D2_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_E_D2_FWHW_FWup$Genome <- 'pan'

write.table(sign_E_D2_FWHW_FWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_FWHW_FWup$Row.names <- rownames(sign_M_D2_FWHW_FWup)
sign_M_D2_FWHW_FWup <- separate(data = as.data.frame(sign_M_D2_FWHW_FWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_FWHW_FWup[,c(2:5,7,8)] <- sapply(sign_M_D2_FWHW_FWup[,c(2:5,7,8)],as.numeric) 
sign_M_D2_FWHW_FWup$Genome <- 'pan'

write.table(sign_M_D2_FWHW_FWup[,c(18,10,11)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_FWHW_FWup_write <- sign_M_D1_FWHW_FWup[,c(18,10,11)]
sign_M_D2_FWHW_FWup_write <- subset(sign_M_D2_FWHW_FWup_write, as.numeric(sign_M_D2_FWHW_FWup_write$end2)-as.numeric(sign_M_D2_FWHW_FWup_write$start2) < 10000)


write.table(sign_M_D2_FWHW_FWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_FWHW_FWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_FWHW_FWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_FWHW_FWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D2_FWHW_FWup_bed <- sign_E_D2_FWHW_FWup[,c(18,8,9)]
colnames(sign_E_D2_FWHW_FWup_bed) <- c('chrom','start','end')
sign_E_D2_FWHW_FWup_bed$start <- as.numeric(sign_E_D2_FWHW_FWup_bed$start)
sign_E_D2_FWHW_FWup_bed$end <- as.numeric(sign_E_D2_FWHW_FWup_bed$end)

sign_M_D2_FWHW_FWup_bed <- sign_M_D2_FWHW_FWup_write
colnames(sign_M_D2_FWHW_FWup_bed) <- c('chrom','start','end')
sign_M_D2_FWHW_FWup_bed$start <- as.numeric(sign_M_D2_FWHW_FWup_bed$start)
sign_M_D2_FWHW_FWup_bed$end <- as.numeric(sign_M_D2_FWHW_FWup_bed$end)

merged_shared_D2FWup_bed <- merged_shared_D2FWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D2FWup_bed) <- c('chrom','start','end')
merged_shared_D2FWup_bed$start <- as.numeric(merged_shared_D2FWup_bed$start)
merged_shared_D2FWup_bed$end <- as.numeric(merged_shared_D2FWup_bed$end)

sign_E_D2_FWHW_FWup_bed_match <- bed_intersect(sign_E_D2_FWHW_FWup_bed,conserved_E)
sign_M_D2_FWHW_FWup_bed_match <- bed_intersect(sign_M_D2_FWHW_FWup_bed,conserved_M)
merged_shared_D2FWup_bed_match <- bed_intersect(merged_shared_D2FWup_bed,conserved_S)


sign_E_D2_FWHW_FWup_cons <- merge(sign_E_D2_FWHW_FWup, sign_E_D2_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D2_FWHW_FWup_cons <- merge(sign_M_D2_FWHW_FWup, sign_M_D2_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D2FWup_cons <- merge(merged_ME_D2FWup, sign_E_D2_FWHW_FWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D2FWup_cons <- merge(merged_EM_D2FWup, sign_M_D2_FWHW_FWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D2FWup_cons[grep("NA", merged_ME_D2FWup_cons$Row.names), ])
nrow(merged_EM_D2FWup_cons[grep("NA", merged_EM_D2FWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D2_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D2_FWHW_FWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D2FWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D2_FWHW_FWup_bed_match)
nrow(sign_M_D2_FWHW_FWup_bed_match)
nrow(merged_shared_D2FWup_bed_match)

nrow(merged_ME_D2FWup_cons[grep("NA", merged_ME_D2FWup_cons$Row.names), ])
nrow(merged_EM_D2FWup_cons[grep("NA", merged_EM_D2FWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))


# D2 HWup

merged_EM_D2HWup <- merge(as.data.frame(res_E_D2FWHW), as.data.frame(sign_M_D2_FWHW_HWup), by = 'row.names')
merged_ME_D2HWup <- merge(as.data.frame(res_M_D2FWHW), as.data.frame(sign_E_D2_FWHW_HWup), by = 'row.names')

nrow(sign_E_D2_FWHW_HWup)-nrow(merged_ME_D2HWup[grep("NA", merged_ME_D2HWup$Row.names), ])
nrow(sign_M_D2_FWHW_HWup)-nrow(merged_EM_D2HWup[grep("NA", merged_EM_D2HWup$Row.names), ])

merged_shared_D2HWup <- merge(as.data.frame(sign_E_D2_FWHW_HWup), as.data.frame(sign_M_D2_FWHW_HWup), by = 'row.names')

merged_EM_D2HWup_unique <- merged_EM_D2FWup[!merged_EM_D2HWup$Row.names %in% merged_shared_D2HWup$Row.names, ]
merged_ME_D2HWup_unique <- merged_ME_D2FWup[!merged_ME_D2HWup$Row.names %in% merged_shared_D2HWup$Row.names, ]

plot(merged_EM_D2HWup_unique$log2FoldChange.x, merged_EM_D2HWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "blue", xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "60h pupa HW up")
par(new=T)
plot(merged_ME_D2HWup_unique$log2FoldChange.x, merged_ME_D2HWup_unique$log2FoldChange.y, pch = 21, xlim = c(-1,4), ylim = c(0,4), col = "orange", xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D2HWup$log2FoldChange.x, merged_shared_D2HWup$log2FoldChange.y, pch = 19, xlim = c(-1,4), ylim = c(0,4), col = "black", xlab = '', ylab = '')

legend(-1, 4, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(21,21,19), cex=1, bty = "n")


ablineclip(lm(merged_EM_D2HWup_unique$log2FoldChange.y ~ merged_EM_D2HWup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D2HWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D2HWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D2HWup_unique$log2FoldChange.y ~ merged_ME_D2HWup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D2HWup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D2HWup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D2HWup$log2FoldChange.y ~ merged_shared_D2HWup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D2HWup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D2HWup$log2FoldChange.x)))

######################## % completely missing

sign_E_D2_FWHW_HWup$Row.names <- rownames(sign_E_D2_FWHW_HWup)
sign_E_D2_FWHW_HWup <- separate(data = as.data.frame(sign_E_D2_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_E_D2_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_E_D2_FWHW_HWup$Genome <- 'pan'

write.table(sign_E_D2_FWHW_HWup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_FWHW_HWup$Row.names <- rownames(sign_M_D2_FWHW_HWup)
sign_M_D2_FWHW_HWup <- separate(data = as.data.frame(sign_M_D2_FWHW_HWup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_FWHW_HWup[,c(2:5,7,8)] <- sapply(sign_M_D2_FWHW_HWup[,c(2:5,7,8)],as.numeric) 
sign_M_D2_FWHW_HWup$Genome <- 'pan'

write.table(sign_M_D2_FWHW_HWup[,c(18,10,11)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_FWHW_HWup_write <- sign_M_D1_FWHW_HWup[,c(18,10,11)]
sign_M_D2_FWHW_HWup_write <- subset(sign_M_D2_FWHW_HWup_write, as.numeric(sign_M_D2_FWHW_HWup_write$end2)-as.numeric(sign_M_D2_FWHW_HWup_write$start2) < 10000)

write.table(sign_M_D2_FWHW_HWup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_FWHW_HWup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_FWHW_HWup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_FWHW_HWup_PAN_IDY.txt', h=T)

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D2_FWHW_HWup_bed <- sign_E_D2_FWHW_HWup[,c(18,8,9)]
colnames(sign_E_D2_FWHW_HWup_bed) <- c('chrom','start','end')
sign_E_D2_FWHW_HWup_bed$start <- as.numeric(sign_E_D2_FWHW_HWup_bed$start)
sign_E_D2_FWHW_HWup_bed$end <- as.numeric(sign_E_D2_FWHW_HWup_bed$end)

sign_M_D2_FWHW_HWup_bed <- sign_M_D2_FWHW_HWup_write
colnames(sign_M_D2_FWHW_HWup_bed) <- c('chrom','start','end')
sign_M_D2_FWHW_HWup_bed$start <- as.numeric(sign_M_D2_FWHW_HWup_bed$start)
sign_M_D2_FWHW_HWup_bed$end <- as.numeric(sign_M_D2_FWHW_HWup_bed$end)

merged_shared_D2HWup_bed <- merged_shared_D2HWup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D2HWup_bed) <- c('chrom','start','end')
merged_shared_D2HWup_bed$start <- as.numeric(merged_shared_D2HWup_bed$start)
merged_shared_D2HWup_bed$end <- as.numeric(merged_shared_D2HWup_bed$end)

sign_E_D2_FWHW_HWup_bed_match <- bed_intersect(sign_E_D2_FWHW_HWup_bed,conserved_E)
sign_M_D2_FWHW_HWup_bed_match <- bed_intersect(sign_M_D2_FWHW_HWup_bed,conserved_M)
merged_shared_D2HWup_bed_match <- bed_intersect(merged_shared_D2HWup_bed,conserved_S)


sign_E_D2_FWHW_HWup_cons <- merge(sign_E_D2_FWHW_HWup, sign_E_D2_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D2_FWHW_HWup_cons <- merge(sign_M_D2_FWHW_HWup, sign_M_D2_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D2HWup_cons <- merge(merged_ME_D2HWup, sign_E_D2_FWHW_HWup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D2HWup_cons <- merge(merged_EM_D2HWup, sign_M_D2_FWHW_HWup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D2HWup_cons[grep("NA", merged_ME_D2HWup_cons$Row.names), ])
nrow(merged_EM_D2HWup_cons[grep("NA", merged_EM_D2HWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D2_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D2_FWHW_HWup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D2HWup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D2_FWHW_HWup_bed_match)
nrow(sign_M_D2_FWHW_HWup_bed_match)
nrow(merged_shared_D2HWup_bed_match)

nrow(merged_ME_D2HWup_cons[grep("NA", merged_ME_D2HWup_cons$Row.names), ])
nrow(merged_EM_D2HWup_cons[grep("NA", merged_EM_D2HWup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))


#########################   