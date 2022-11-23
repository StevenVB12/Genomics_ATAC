library(Rsamtools)
library(ggplot2)
library(magrittr)
library(dplyr)
library(DESeq2)
library(maptools)
library(tidyr)

counts_erato <- read.table("MACS2/dem_hyd.MACS2.peaks.all.merged.sort.counts", h=F)
counts_melp <- read.table("MACS2/melp_ros.MACS2.peaks.all.merged.sort.counts", h=F)

optix_intervals <- subset(counts_erato, (counts_erato$V1 == 'Herato1801') & 
              (counts_erato$V2 > 1239943-2000) & 
              (counts_erato$V3 < 1428307) & ((counts_erato$V3-counts_erato$V2) > 100))[,c(1:3)]

write.table(optix_intervals, "optix_ATAC_peaks_Herato1801_1237943_1428307.bed", quote = F, col.names = F, row.names = F)

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

nrow(counts_erato_pan)
nrow(counts_melp_pan)

nrow(match_ME)
nrow(match_EM)

nrow(match_ME_A)
nrow(match_EM_A)

# merge intersect files with peak counts
counts_erato_pan_match <- merge(counts_erato_pan, match_ME, by.x = c('pan_start','pan_end'), by.y =c('pan_start_E', 'pan_end_E'), all.x = TRUE, all.y = TRUE)
# counts_melp_pan_match <- merge(counts_melp_pan, match_EM, by.x = c('pan_start','pan_end'), by.y =c('pan_start_M', 'pan_end_M'))
head(counts_erato_pan_match)
nrow(counts_erato_pan_match)
# nrow(counts_melp_pan_match)
counts_erato_melp_pan_match <- merge(counts_erato_pan_match, counts_melp_pan, by.x = c('pan_start_M', 'pan_end_M'), by.y = c('pan_start','pan_end'))
counts_erato_melp_pan_match <- merge(counts_erato_pan_match, counts_melp_pan, by.x = c('pan_start_M', 'pan_end_M'), by.y = c('pan_start','pan_end'), all.x = TRUE, all.y = TRUE)

nrow(counts_erato_melp_pan_match)

cts1 <- counts_erato_melp_pan_match[,c('genome.x', 'pan_start', 'pan_end', 'pan_start_M', 'pan_end_M', 'scaffold.x', 'start.x', 'end.x', 'scaffold.y','start.y', 'end.y')]
cts2 <- counts_erato_melp_pan_match[,c(as.character(sampleInfo_erato$id), as.character(sampleInfo_melp$id))]

cts <- cbind(paste(cts1$genome.x,cts1$pan_start,cts1$pan_end,cts1$pan_start_M,cts1$pan_end_M, cts1$scaffold.x,cts1$start.x,cts1$end.x, cts1$scaffold.y, cts1$start.y, cts1$end.y, sep='_'), cts2)
colnames(cts) <- c('peak', c(as.character(sampleInfo_erato$id), as.character(sampleInfo_melp$id)))
head(cts)

ctsF <- cts[,-1]
rownames(ctsF) <- cts[,1]
head(ctsF)

ctsF2A <- ctsF2[-grep("NA", rownames(ctsF2)), ]
# 

sampleInfo_sub <- subset(sampleInfo, sampleInfo$id != 'LI15_melpomene_HW')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI22_melpomene_FP')
sampleInfo_sub <- subset(sampleInfo_sub, sampleInfo_sub$id != 'LI20_demophoon_HP')

ctsF2 <- select(ctsF, c(as.character(sampleInfo_sub$id)))


nrow(ctsF2[is.na(ctsF2$BR10_Demophoon_Brain),])
nrow(ctsF2[!is.na(ctsF2$LI6_melpomene_HP),])

nrow(counts_erato[!is.na(counts_erato$BR10_Demophoon_Brain),])

ctsF[is.na(ctsF)] <- 0
ctsF2[is.na(ctsF2)] <- 0

ctsF2A <- ctsF2[-grep("NA", rownames(ctsF2)), ]
######################
# 5th vs D1/D2
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato'))

sampleInfo_subE$comp <- 'pupa'
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$stage[e] == '5thinstar'){sampleInfo_subE$comp[e] <- '5thinstar' }}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp'))

sampleInfo_subM$comp <- 'pupa'
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$stage[e] == '5thinstar'){sampleInfo_subM$comp[e] <- '5thinstar' }}

sampleInfo_sub$comp <- 'pupa'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$stage[e] == '5thinstar'){sampleInfo_sub$comp[e] <- '5thinstar' }}

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == '5thinstar') &
                                             sampleInfo_sub$tissue != 'brain'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ctsF2A_5th <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))
###
# setup model
###
# 
dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                               colData = sampleInfo_sub,
                               design = ~tissue)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                              colData = sampleInfo_subE,
                              design = ~comp)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~comp)

ddsA <- DESeqDataSetFromMatrix(countData = ctsF2A_5th,
                               colData = sampleInfo_subA,
                               design = ~species)
###
# run model
###
# dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E_5th_D1D2 <- results(atacDDS_E)
res_M_5th_D1D2 <- results(atacDDS_M)

atacDDS_A <- DESeq(ddsA)
res_A_5th_D1D2 <- results(atacDDS_A)
sign_A_5th_D1D2 <- subset(res_A_5th_D1D2, res_A_5th_D1D2$padj < 0.05 & abs(res_A_5th_D1D2$log2FoldChange > 1))

nrow(sign_A_5th_D1D2)
nrow(merge(as.data.frame(sign_A_5th_D1D2), as.data.frame(sign_E_5th_D1D2_5thup), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_5th_D1D2), as.data.frame(sign_M_5th_D1D2_5thup), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_5th_D1D2), as.data.frame(sign_E_5th_D1D2_5thdown), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_5th_D1D2), as.data.frame(sign_M_5th_D1D2_5thdown), by = 'row.names'))

save(res_E_5th_D1D2, file='sign_E_5th_D1D2_p05fc1_ALL.rda')
save(res_M_5th_D1D2, file='sign_M_5th_D1D2_p05fc1_ALL.rda')

sign_E_5th_D1D2_5thup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_5th_D1D2_5thup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

sign_E_5th_D1D2_5thdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_5th_D1D2_5thdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

nrow(sign_E_5th_D1D2_5thup)
nrow(sign_M_5th_D1D2_5thup)
nrow(sign_E_5th_D1D2_5thdown)
nrow(sign_M_5th_D1D2_5thdown)

# save(sign_E_5th_D1D2_5thup, file='sign_E_5th_D1D2_p05fc1_5thup.rda')
# save(sign_E_5th_D1D2_5thdown, file='sign_E_5th_D1D2_p05fc1_5thdown.rda')
# save(sign_M_5th_D1D2_5thup, file='sign_M_5th_D1D2_p05fc1_5thup.rda')
# save(sign_M_5th_D1D2_5thdown, file='sign_M_5th_D1D2_p05fc1_5thdown.rda')
# 
# load('sign_E_5th_D1D2_p05fc1_5thup.rda')
# load('sign_E_5th_D1D2_p05fc1_5thdown.rda')
# load('sign_M_5th_D1D2_p05fc1_5thup.rda')
# load('sign_M_5th_D1D2_p05fc1_5thdown.rda')

sign_E_5th_D1D2_5thup_rowname <- as.data.frame(rownames(sign_E_5th_D1D2_5thup))
sign_M_5th_D1D2_5thup_rowname <- as.data.frame(rownames(sign_M_5th_D1D2_5thup))
sign_E_5th_D1D2_5thdown_rowname <- as.data.frame(rownames(sign_E_5th_D1D2_5thdown))
sign_M_5th_D1D2_5thdown_rowname <- as.data.frame(rownames(sign_M_5th_D1D2_5thdown))

sign_E_5th_D1D2_5thup_rowname <- separate(data = sign_E_5th_D1D2_5thup_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_D1D2_5thup_rowname <- separate(data = sign_M_5th_D1D2_5thup_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_5th_D1D2_5thdown_rowname <- separate(data = sign_E_5th_D1D2_5thdown_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_D1D2_5thdown_rowname <- separate(data = sign_M_5th_D1D2_5thdown_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

sign_E_5th_D1D2_5thup_rowname$genome <- 'pan'
sign_M_5th_D1D2_5thup_rowname$genome <- 'pan'
sign_E_5th_D1D2_5thdown_rowname$genome <- 'pan'
sign_M_5th_D1D2_5thdown_rowname$genome <- 'pan'

head(subset(sign_E_5th_D1D2_5thup_rowname, sign_E_5th_D1D2_5thup_rowname$scafM == 'NA'))
head(subset(sign_M_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname$scafE == 'NA'))
head(subset(sign_E_5th_D1D2_5thup_rowname, sign_E_5th_D1D2_5thup_rowname$scafM != 'NA'), 10)
head(subset(sign_M_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname$scafE != 'NA'), 10)

write.table(subset(sign_E_5th_D1D2_5thup_rowname, sign_E_5th_D1D2_5thup_rowname$scafM == 'NA')[,c(1:3)], file='5thup_eratoUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(sign_M_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname$scafE == 'NA' & (as.numeric(sign_M_5th_D1D2_5thup_rowname$end2) - as.numeric(sign_M_5th_D1D2_5thup_rowname$start2) < 10000))[,c(1,4,5)], file='5thup_melpUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

sign_shared_5th <- unique.data.frame(rbind(subset(sign_E_5th_D1D2_5thup_rowname, sign_E_5th_D1D2_5thup_rowname$scafM != 'NA'), subset(sign_M_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname$scafE != 'NA')))
sign_shared_5th$start_min <- apply(sign_shared_5th[,c(2,4)], 1, FUN = min)
sign_shared_5th$start_max <- apply(sign_shared_5th[,c(3,5)], 1, FUN = min)

write.table(sign_shared_5th[,c(1,12,13)], file='5thup_shared_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

head(sign_shared_5th)

nrow(subset(sign_E_5th_D1D2_5thup_rowname, sign_E_5th_D1D2_5thup_rowname$scafM != 'NA'))
nrow(subset(sign_M_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname$scafE != 'NA'))
nrow(subset(sign_E_5th_D1D2_5thdown_rowname, sign_E_5th_D1D2_5thdown_rowname$scafM != 'NA'))
nrow(subset(sign_M_5th_D1D2_5thdown_rowname, sign_M_5th_D1D2_5thdown_rowname$scafE != 'NA'))


nrow(merge(sign_E_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname, by = "start"))

shared_5thup_strict <- unique(merge(sign_E_5th_D1D2_5thup_rowname, sign_M_5th_D1D2_5thup_rowname, by = c("start"), all = F)[,c(2,1,3)])
colnames(shared_5thup_strict) <- c("genome", "start", "end")

unique_erato_5th_strict <- sign_E_5th_D1D2_5thup_rowname[!sign_E_5th_D1D2_5thup_rowname$start %in% shared_5thup_strict$start,][,c(1:3)]
unique_melp_5th_strict <- sign_M_5th_D1D2_5thup_rowname[!sign_M_5th_D1D2_5thup_rowname$start %in% shared_5thup_strict$start,][,c(1,4,5)]


write.table(shared_5thup_strict, file='5thup_shared_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(unique_erato_5th_strict, file='5thup_eratoUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(unique_melp_5th_strict, (as.numeric(unique_melp_5th_strict$end) - as.numeric(unique_melp_5th_strict$start) < 10000)), file='5thup_melpUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')

nrow(shared_5thup_strict)

######################
# D1 vs 5th/D2
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato'))

sampleInfo_subE$comp <- 'class2'
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$stage[e] == 'pupalDAY1'){sampleInfo_subE$comp[e] <- 'pupalDAY1' }}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp'))

sampleInfo_subM$comp <- 'class2'
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$stage[e] == 'pupalDAY1'){sampleInfo_subM$comp[e] <- 'pupalDAY1' }}

sampleInfo_sub$comp <- 'class2'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$stage[e] == 'pupalDAY1'){sampleInfo_sub$comp[e] <- 'pupalDAY1' }}

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$tissue != 'brain'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ctsF2A_D1 <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))
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

ddsA <- DESeqDataSetFromMatrix(countData = ctsF2A_D1,
                               colData = sampleInfo_subA,
                               design = ~species)
###
# run model
###
# dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E_D1_5thD2 <- results(atacDDS_E)
res_M_D1_5thD2 <- results(atacDDS_M)

atacDDS_A <- DESeq(ddsA)
res_A_D1_5thD2 <- results(atacDDS_A)
sign_A_D1_5thD2 <- subset(res_A_D1_5thD2, res_A_D1_5thD2$padj < 0.05 & abs(res_A_D1_5thD2$log2FoldChange > 1))

nrow(sign_A_D1_5thD2)
nrow(merge(as.data.frame(sign_A_D1_5thD2), as.data.frame(sign_E_D1_5thD2_D1up), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_5thD2), as.data.frame(sign_M_D1_5thD2_D1up), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_D1_5thD2), as.data.frame(sign_E_D1_5thD2_D1down), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_5thD2), as.data.frame(sign_M_D1_5thD2_D1down), by = 'row.names'))

save(res_E_D1_5thD2, file='sign_E_D1_5thD2_p05fc1_ALL.rda')
save(res_M_D1_5thD2, file='sign_M_D1_5thD2_p05fc1_ALL.rda')


sign_E_D1_5thD2_D1up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_5thD2_D1up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D1_5thD2_D1down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_5thD2_D1down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

nrow(sign_E_D1_5thD2_D1up)
nrow(sign_M_D1_5thD2_D1up)
nrow(sign_E_D1_5thD2_D1down)
nrow(sign_M_D1_5thD2_D1down)

# save(sign_E_D1_5thD2_D1up, file='sign_E_D1_5thD2_p05fc1_D1up.rda')
# save(sign_M_D1_5thD2_D1up, file='sign_M_D1_5thD2_p05fc1_D1up.rda')
# save(sign_E_D1_5thD2_D1down, file='sign_E_D1_5thD2_p05fc1_D1down.rda')
# save(sign_M_D1_5thD2_D1down, file='sign_M_D1_5thD2_p05fc1_D1down.rda')



sign_E_D1_5thD2_D1up_rowname <- as.data.frame(rownames(sign_E_D1_5thD2_D1up))
sign_M_D1_5thD2_D1up_rowname <- as.data.frame(rownames(sign_M_D1_5thD2_D1up))
sign_E_D1_5thD2_D1down_rowname <- as.data.frame(rownames(sign_E_D1_5thD2_D1down))
sign_M_D1_5thD2_D1down_rowname <- as.data.frame(rownames(sign_M_D1_5thD2_D1down))

sign_E_D1_5thD2_D1up_rowname <- separate(data = sign_E_D1_5thD2_D1up_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_5thD2_D1up_rowname <- separate(data = sign_M_D1_5thD2_D1up_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_5thD2_D1down_rowname <- separate(data = sign_E_D1_5thD2_D1down_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_5thD2_D1down_rowname <- separate(data = sign_M_D1_5thD2_D1down_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

sign_E_D1_5thD2_D1up_rowname$genome <- 'pan'
sign_M_D1_5thD2_D1up_rowname$genome <- 'pan'
sign_E_D1_5thD2_D1down_rowname$genome <- 'pan'
sign_M_D1_5thD2_D1down_rowname$genome <- 'pan'

head(subset(sign_E_D1_5thD2_D1up_rowname, sign_E_D1_5thD2_D1up_rowname$scafM == 'NA'))
head(subset(sign_M_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname$scafE == 'NA'))
head(subset(sign_E_D1_5thD2_D1up_rowname, sign_E_D1_5thD2_D1up_rowname$scafM != 'NA'), 10)
head(subset(sign_M_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname$scafE != 'NA'), 10)

write.table(subset(sign_E_D1_5thD2_D1up_rowname, sign_E_D1_5thD2_D1up_rowname$scafM == 'NA')[,c(1:3)], file='D1up_eratoUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(sign_M_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname$scafE == 'NA' & (as.numeric(sign_M_D1_5thD2_D1up_rowname$end2) - as.numeric(sign_M_D1_5thD2_D1up_rowname$start2) < 10000))[,c(1,4,5)], file='D1up_melpUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

sign_shared_D1 <- unique.data.frame(rbind(subset(sign_E_D1_5thD2_D1up_rowname, sign_E_D1_5thD2_D1up_rowname$scafM != 'NA'), subset(sign_M_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname$scafE != 'NA')))
sign_shared_D1$start_min <- apply(sign_shared_D1[,c(2,4)], 1, FUN = min)
sign_shared_D1$start_max <- apply(sign_shared_D1[,c(3,5)], 1, FUN = min)

write.table(sign_shared_D1[,c(1,12,13)], file='D1up_shared_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

head(sign_shared_D1)

nrow(subset(sign_E_D1_5thD2_D1up_rowname, sign_E_D1_5thD2_D1up_rowname$scafM != 'NA'))
nrow(subset(sign_M_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname$scafE != 'NA'))
nrow(subset(sign_E_D1_5thD2_D1down_rowname, sign_E_D1_5thD2_D1down_rowname$scafM != 'NA'))
nrow(subset(sign_M_D1_5thD2_D1down_rowname, sign_M_D1_5thD2_D1down_rowname$scafE != 'NA'))


shared_D1up_strict <- unique(merge(sign_E_D1_5thD2_D1up_rowname, sign_M_D1_5thD2_D1up_rowname, by = c("start"), all = F)[,c(2,1,3)])
colnames(shared_D1up_strict) <- c("genome", "start", "end")

unique_erato_D1_strict <- sign_E_D1_5thD2_D1up_rowname[!sign_E_D1_5thD2_D1up_rowname$start %in% shared_D1up_strict$start,][,c(1:3)]
unique_melp_D1_strict <- sign_M_D1_5thD2_D1up_rowname[!sign_M_D1_5thD2_D1up_rowname$start %in% shared_D1up_strict$start,][,c(1,4,5)]


write.table(shared_D1up_strict, file='D1up_shared_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(unique_erato_D1_strict, file='D1up_eratoUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(unique_melp_D1_strict, (as.numeric(unique_melp_D1_strict$end) - as.numeric(unique_melp_D1_strict$start) < 10000)), file='D1up_melpUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')



######################
# D2 vs 5th/D1
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'erato'))

sampleInfo_subE$comp <- 'class2'
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$stage[e] == 'pupalDAY2'){sampleInfo_subE$comp[e] <- 'pupalDAY2' }}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue != 'brain') &
                                             sampleInfo_sub$species == 'melp'))

sampleInfo_subM$comp <- 'class2'
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$stage[e] == 'pupalDAY2'){sampleInfo_subM$comp[e] <- 'pupalDAY2' }}

sampleInfo_sub$comp <- 'class2'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$stage[e] == 'pupalDAY2'){sampleInfo_sub$comp[e] <- 'pupalDAY2' }}

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$tissue != 'brain'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ctsF2A_D2 <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))

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

ddsA <- DESeqDataSetFromMatrix(countData = ctsF2A_D2,
                               colData = sampleInfo_subA,
                               design = ~species)
###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E_D2_5thD1 <- results(atacDDS_E)
res_M_D2_5thD1 <- results(atacDDS_M)


atacDDS_A <- DESeq(ddsA)
res_A_D2_5thD1 <- results(atacDDS_A)
sign_A_D2_5thD1 <- subset(res_A_D2_5thD1, res_A_D2_5thD1$padj < 0.05 & abs(res_A_D2_5thD1$log2FoldChange > 1))

nrow(sign_A_D2_5thD1)
nrow(merge(as.data.frame(sign_A_D2_5thD1), as.data.frame(sign_E_D2_5thD1_D2up), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D2_5thD1), as.data.frame(sign_M_D2_5thD1_D2up), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_D2_5thD1), as.data.frame(sign_E_D2_5thD1_D2down), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D2_5thD1), as.data.frame(sign_M_D2_5thD1_D2down), by = 'row.names'))

save(res_E_D2_5thD1, file='sign_E_D2_5thD1_p05fc1_ALL.rda')
save(res_M_D2_5thD1, file='sign_M_D2_5thD1_p05fc1_ALL.rda')



sign_E_D2_5thD1_D2up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_5thD1_D2up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_5thD1_D2down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_5thD1_D2down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

nrow(sign_E_D2_5thD1_D2up)
nrow(sign_M_D2_5thD1_D2up)
nrow(sign_E_D2_5thD1_D2down)
nrow(sign_M_D2_5thD1_D2down)

# save(sign_E_D2_5thD1_D2up, file='sign_E_D2_5thD1_p05fc1_D2up.rda')
# save(sign_M_D2_5thD1_D2up, file='sign_M_D2_5thD1_p05fc1_D2up.rda')
# save(sign_E_D2_5thD1_D2down, file='sign_E_D2_5thD1_p05fc1_D2down.rda')
# save(sign_M_D2_5thD1_D2down, file='sign_M_D2_5thD1_p05fc1_D2down.rda')



sign_E_D2_5thD1_D2up_rowname <- as.data.frame(rownames(sign_E_D2_5thD1_D2up))
sign_M_D2_5thD1_D2up_rowname <- as.data.frame(rownames(sign_M_D2_5thD1_D2up))
sign_E_D2_5thD1_D2down_rowname <- as.data.frame(rownames(sign_E_D2_5thD1_D2down))
sign_M_D2_5thD1_D2down_rowname <- as.data.frame(rownames(sign_M_D2_5thD1_D2down))

sign_E_D2_5thD1_D2up_rowname <- separate(data = sign_E_D2_5thD1_D2up_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_5thD1_D2up_rowname <- separate(data = sign_M_D2_5thD1_D2up_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_5thD1_D2down_rowname <- separate(data = sign_E_D2_5thD1_D2down_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_5thD1_D2down_rowname <- separate(data = sign_M_D2_5thD1_D2down_rowname, col = 1, into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

sign_E_D2_5thD1_D2up_rowname$genome <- 'pan'
sign_M_D2_5thD1_D2up_rowname$genome <- 'pan'
sign_E_D2_5thD1_D2down_rowname$genome <- 'pan'
sign_M_D2_5thD1_D2down_rowname$genome <- 'pan'

head(subset(sign_E_D2_5thD1_D2up_rowname, sign_E_D2_5thD1_D2up_rowname$scafM == 'NA'))
head(subset(sign_M_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname$scafE == 'NA'))
head(subset(sign_E_D2_5thD1_D2up_rowname, sign_E_D2_5thD1_D2up_rowname$scafM != 'NA'), 10)
head(subset(sign_M_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname$scafE != 'NA'), 10)

write.table(subset(sign_E_D2_5thD1_D2up_rowname, sign_E_D2_5thD1_D2up_rowname$scafM == 'NA')[,c(1:3)], file='D2up_eratoUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(sign_M_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname$scafE == 'NA' & (as.numeric(sign_M_D2_5thD1_D2up_rowname$end2) - as.numeric(sign_M_D2_5thD1_D2up_rowname$start2) < 10000))[,c(1,4,5)], file='D2up_melpUnique_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

sign_shared_D2 <- unique.data.frame(rbind(subset(sign_E_D2_5thD1_D2up_rowname, sign_E_D2_5thD1_D2up_rowname$scafM != 'NA'), subset(sign_M_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname$scafE != 'NA')))
sign_shared_D2$start_min <- apply(sign_shared_D2[,c(2,4)], 1, FUN = min)
sign_shared_D2$start_max <- apply(sign_shared_D2[,c(3,5)], 1, FUN = min)

write.table(sign_shared_D2[,c(1,12,13)], file='D2up_shared_PANpos_50perccOverlap.bed', col.names = F, quote = F, row.names = F, sep = '\t')

head(sign_shared_D2)

nrow(subset(sign_E_D2_5thD1_D2up_rowname, sign_E_D2_5thD1_D2up_rowname$scafM != 'NA'))
nrow(subset(sign_M_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname$scafE != 'NA'))
nrow(subset(sign_E_D2_5thD1_D2down_rowname, sign_E_D2_5thD1_D2down_rowname$scafM != 'NA'))
nrow(subset(sign_M_D2_5thD1_D2down_rowname, sign_M_D2_5thD1_D2down_rowname$scafE != 'NA'))

shared_D2up_strict <- unique(merge(sign_E_D2_5thD1_D2up_rowname, sign_M_D2_5thD1_D2up_rowname, by = c("start"), all = F)[,c(2,1,3)])
colnames(shared_D2up_strict) <- c("genome", "start", "end")

unique_erato_D2_strict <- sign_E_D2_5thD1_D2up_rowname[!sign_E_D2_5thD1_D2up_rowname$start %in% shared_D2up_strict$start,][,c(1:3)]
unique_melp_D2_strict <- sign_M_D2_5thD1_D2up_rowname[!sign_M_D2_5thD1_D2up_rowname$start %in% shared_D2up_strict$start,][,c(1,4,5)]


write.table(shared_D2up_strict, file='D2up_shared_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(unique_erato_D2_strict, file='D2up_eratoUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')
write.table(subset(unique_melp_D2_strict, (as.numeric(unique_melp_D2_strict$end) - as.numeric(unique_melp_D2_strict$start) < 10000)), file='D2up_melpUnique_PANpos_50perccOverlap_strict.bed', col.names = F, quote = F, row.names = F, sep = '\t')


### load

load("sign_E_5th_D1D2_p05fc1.rda")
load("sign_M_5th_D1D2_p05fc1.rda")

load("sign_E_D1_5thD2_p05fc1.rda")
load("sign_M_D1_5thD2_p05fc1.rda")

load("sign_E_D2_5thD1_p05fc1.rda")
load("sign_M_D2_5thD1_p05fc1.rda")


load("sign_E_5th_D1D2_p05fc1_5thup.rda")
load("sign_M_5th_D1D2_p05fc1_5thup.rda")

load("sign_E_D1_5thD2_p05fc1_D1up.rda")
load("sign_M_D1_5thD2_p05fc1_D1up.rda")

load("sign_E_D2_5thD1_p05fc1_D2up.rda")
load("sign_M_D2_5thD1_p05fc1_D2up.rda")

load("sign_E_5th_D1D2_p05fc1_5thdown.rda")
load("sign_M_5th_D1D2_p05fc1_5tdown.rda")

load("sign_E_D1_5thD2_p05fc1_D1down.rda")
load("sign_M_D1_5thD2_p05fc1_D1down.rda")

load("sign_E_D2_5thD1_p05fc1_D2down.rda")
load("sign_M_D2_5thD1_p05fc1_D2down.rda")

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

sign_list_E_dir <- list(sign_E_5th_D1D2_5thup, sign_E_D1_5thD2_D1up, sign_E_D2_5thD1_D2up)
sign_list_M_dir <- list(sign_M_5th_D1D2_5thup, sign_M_D1_5thD2_D1up, sign_M_D2_5thD1_D2up)

colList <- c('#ff0000b5', '#009e73ff', '#0072b2ff')

sign_list_E_df_dir <- list()
sign_list_E_df_all <- list()
for(i in 1:length(sign_list_E_dir)){
  sign_list_E_df_dir[[i]] <- as.data.frame(rownames(sign_list_E_dir[[i]]))
  sign_list_E_df_dir[[i]] <- separate(data = sign_list_E_df_dir[[i]], col = 'rownames(sign_list_E_dir[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_E_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_E_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_E_df_all[[i]] <- cbind(sign_list_E_dir[[i]], sign_list_E_df_dir[[i]])
}

sign_list_M_df_dir <- list()
sign_list_M_df_all <- list()
for(i in 1:length(sign_list_M_dir)){
  sign_list_M_df_dir[[i]] <- as.data.frame(rownames(sign_list_M_dir[[i]]))
  sign_list_M_df_dir[[i]] <- separate(data = sign_list_M_df_dir[[i]], col = 'rownames(sign_list_M_dir[[i]])', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_M_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_M_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_M_df_all[[i]] <- cbind(sign_list_M_dir[[i]], sign_list_M_df_dir[[i]])
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

plot(d_shared_E_5th, xlim=c(1,6), ylim=c(0,2.5), main='5th erato', xlab='', lwd=2)
polygon(d_shared_E_5th, col=adjustcolor("red", alpha=0.2), border="red")
par(new=T)
plot(d_unique_E_5th, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='')


shared_E_D1 <- subset(sign_list_E_df_all[[2]], sign_list_E_df_all[[2]]$start2 != -1)
unique_E_D1 <- subset(sign_list_E_df_all[[2]], sign_list_E_df_all[[2]]$start2 == -1)

d_shared_E_D1 <- density(abs(shared_E_D1$log2FoldChange))
d_unique_E_D1 <- density(abs(unique_E_D1$log2FoldChange))

plot(d_shared_E_D1, xlim=c(1,6), ylim=c(0,2.5),  main='D1 erato', xlab='', lwd=2)
polygon(d_shared_E_D1, col=adjustcolor("blue", alpha=0.2), border="blue")
par(new=T)
plot(d_unique_E_D1, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='')


shared_E_D2 <- subset(sign_list_E_df_all[[3]], sign_list_E_df_all[[3]]$start2 != -1)
unique_E_D2 <- subset(sign_list_E_df_all[[3]], sign_list_E_df_all[[3]]$start2 == -1)

d_shared_E_D2 <- density(abs(shared_E_D2$log2FoldChange))
d_unique_E_D2 <- density(abs(unique_E_D2$log2FoldChange))

plot(d_shared_E_D2, xlim=c(1,6), ylim=c(0,2.5),  main='D2 erato', xlab='', lwd=2)
polygon(d_shared_E_D2, col=adjustcolor("darkgreen", alpha=0.2), border="darkgreen")
par(new=T)
plot(d_unique_E_D2, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='')



shared_M_5th <- subset(sign_list_M_df_all[[1]], is.na(sign_list_M_df_all[[1]]$start) == F)
unique_M_5th <- subset(sign_list_M_df_all[[1]], is.na(sign_list_M_df_all[[1]]$start) == T)

d_shared_M_5th <- density(abs(shared_M_5th$log2FoldChange))
d_unique_M_5th <- density(abs(unique_M_5th$log2FoldChange))

plot(d_shared_M_5th, xlim=c(1,6), ylim=c(0,2.5), main='5th melp', xlab='', lwd=2)
polygon(d_shared_M_5th, col=adjustcolor("red", alpha=0.2), border="red")
par(new=T)
plot(d_unique_M_5th, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='log2FC')


shared_M_D1 <- subset(sign_list_M_df_all[[2]], is.na(sign_list_M_df_all[[2]]$start) == F)
unique_M_D1 <- subset(sign_list_M_df_all[[2]], is.na(sign_list_M_df_all[[2]]$start) == T)

d_shared_M_D1 <- density(abs(shared_M_D1$log2FoldChange))
d_unique_M_D1 <- density(abs(unique_M_D1$log2FoldChange))

plot(d_shared_M_D1, xlim=c(1,6), ylim=c(0,2.5),  main='D1 melp', xlab='', lwd=2)
polygon(d_shared_M_D1, col=adjustcolor("blue", alpha=0.2), border="blue")
par(new=T)
plot(d_unique_M_D1, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='log2FC')


shared_M_D2 <- subset(sign_list_M_df_all[[3]], is.na(sign_list_M_df_all[[3]]$start) == F)
unique_M_D2 <- subset(sign_list_M_df_all[[3]], is.na(sign_list_M_df_all[[3]]$start) == T)

d_shared_M_D2 <- density(abs(shared_M_D2$log2FoldChange))
d_unique_M_D2 <- density(abs(unique_M_D2$log2FoldChange))

plot(d_shared_M_D2, xlim=c(1,6), ylim=c(0,2.5),  main='D2 melp', xlab='', lwd=2)
polygon(d_shared_M_D2, col=adjustcolor("darkgreen", alpha=0.2), border="darkgreen")
par(new=T)
plot(d_unique_M_D2, xlim=c(1,6), ylim=c(0,2.5), main='', lwd=2, xlab='log2FC')
################






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

# shared_erato_5th_strict <- sign_E_5th_D1D2_5thup_rowname[!sign_E_5th_D1D2_5thup_rowname$start %in% shared_5thup_strict$start,][,c(1:3)]
# 
# shared_5thup_strict
# unique_erato_5th_strict
# unique_melp_5th_strict 

shared_5thup_strict$start <- as.numeric(shared_5thup_strict$start)
shared_D1up_strict$start <- as.numeric(shared_D1up_strict$start)
shared_D2up_strict$start <- as.numeric(shared_D2up_strict$start)

shared_5thup_strict$end <- as.numeric(shared_5thup_strict$end)
shared_D1up_strict$end <- as.numeric(shared_D1up_strict$end)
shared_D2up_strict$end <- as.numeric(shared_D2up_strict$end)

sign_list_shared_dir <- list(shared_5thup_strict, shared_D1up_strict, shared_D2up_strict)


#shared peaks
plot(NULL, xlim=c(0,1600), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(0,1600,by=100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

ptest <- c()
for(e in 1:21){
  
  countpeaks_S <- c()
  for(p in 1:length(sign_list_shared_dir)){
    
    signs_S <- subset(sign_list_shared_dir[[p]], sign_list_shared_dir[[p]]$start > chromTable$start_pan[e] & sign_list_shared_dir[[p]]$end < chromTable$end_pan[e])

    countpeaks_S <- c(countpeaks_S, nrow(signs_S))
  }
  countpeaks <- as.table(cbind(countpeaks_S,0))
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = T, xlim=c(0,1500), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  # pval <- c()
  # for(k in 1:3) {
  #   pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
  #   pval <- c(pval, pstat$p.value)
  # }
  # ptest <- rbind(ptest, pval)
}

# #shared peaks
# plot(NULL, xlim=c(0,1000), ylim = c(0,1), axes=FALSE, ann=FALSE)
# axis(1,seq(0,1000,by=100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)
# 
# ptest <- c()
# for(e in 1:21){
#   
#   countpeaks_E <- c()
#   countpeaks_M <- c()
#   for(p in 1:length(sign_list_E_df)){
#     
#     signs_E <- subset(sign_list_E_df[[p]], sign_list_E_df[[p]]$start > chromTable$start_pan[e] & sign_list_E_df[[p]]$end < chromTable$end_pan[e] & sign_list_E_df[[p]]$start2 != -1)
#     signs_M <- subset(sign_list_M_df[[p]], sign_list_M_df[[p]]$start2 > chromTable$start_pan[e] & sign_list_M_df[[p]]$end2 < chromTable$end_pan[e] & is.na(sign_list_M_df[[p]]$start) == FALSE)
#     signs_M <- subset(signs_M, abs(signs_M$end2 - signs_M$start2) < 100000) 
#     
#     countpeaks_E <- c(countpeaks_E, nrow(signs_E))
#     countpeaks_M <- c(countpeaks_M, nrow(signs_M))
#   }
#   countpeaks <- as.table(cbind(countpeaks_E,countpeaks_M))
#   barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(0,1000), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
#   
#   pval <- c()
#   for(k in 1:3) {
#     pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
#     pval <- c(pval, pstat$p.value)
#   }
#   ptest <- rbind(ptest, pval)
# }

plot(NULL, xlim=c(0,1600), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "shared")

unique_erato_5th_strict$start <- as.numeric(unique_erato_5th_strict$start)
unique_erato_D1_strict$start <- as.numeric(unique_erato_D1_strict$start)
unique_erato_D2_strict$start <- as.numeric(unique_erato_D2_strict$start)

unique_erato_5th_strict$end <- as.numeric(unique_erato_5th_strict$end)
unique_erato_D1_strict$end <- as.numeric(unique_erato_D1_strict$end)
unique_erato_D2_strict$end <- as.numeric(unique_erato_D2_strict$end)

unique_melp_5th_strict$start <- as.numeric(unique_melp_5th_strict$start)
unique_melp_D1_strict$start <- as.numeric(unique_melp_D1_strict$start)
unique_melp_D2_strict$start <- as.numeric(unique_melp_D2_strict$start)

unique_melp_5th_strict$end <- as.numeric(unique_melp_5th_strict$end)
unique_melp_D1_strict$end <- as.numeric(unique_melp_D1_strict$end)
unique_melp_D2_strict$end <- as.numeric(unique_melp_D2_strict$end)


sign_list_Eunique_dir <- list(unique_erato_5th_strict, unique_erato_D1_strict, unique_erato_D2_strict)
sign_list_Munique_dir <- list(unique_melp_5th_strict, unique_melp_D1_strict, unique_melp_D2_strict)

#unique peaks
plot(NULL, xlim=c(1600,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(1600,0,by=-100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

ptest <- c()
for(e in 1:21){
  
  countpeaks_E <- c()
  countpeaks_M <- c()
  for(p in 1:length(sign_list_Eunique_dir)){
    
    signs_E <- subset(sign_list_Eunique_dir[[p]], sign_list_Eunique_dir[[p]]$start > chromTable$start_pan[e] & sign_list_Eunique_dir[[p]]$end < chromTable$end_pan[e])
    signs_M <- subset(sign_list_Munique_dir[[p]], sign_list_Munique_dir[[p]]$start > chromTable$start_pan[e] & sign_list_Munique_dir[[p]]$end < chromTable$end_pan[e])
    
    
    countpeaks_E <- c(countpeaks_E, nrow(signs_E))
    countpeaks_M <- c(countpeaks_M, nrow(signs_M))
  }
  countpeaks <- as.table(cbind(countpeaks_E,countpeaks_M))
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(1600,0), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  pval <- c()
  for(k in 1:3) {
    pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
    pval <- c(pval, pstat$p.value)
  }
  ptest <- rbind(ptest, pval)
}

plot(NULL, xlim=c(1600,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "unique")
####


plot(NULL, xlim = c(-2000000,max(chromTable$end_pan-chromTable$start_pan)+1000000), ylim = c(0,1), axes=F, ylab = '', xlab = '', xaxs="i", yaxs="i")
axis(1,seq(0,max(chromTable$end_pan-chromTable$start_pan),by=3000000),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)


par(mar=c(0,0,0,1))

plot(NULL, xlim = c(-2000000,max(chromTable$end_pan-chromTable$start_pan)+1000000), ylim = c(0,21), axes=F, ylab = '', xlab = '', xaxs="i", yaxs="i", xaxt='none', yaxt='none')
# axis(seq(1,30,1))

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
  text(-500000, e-0.75, 'e', cex=1)
  text(-500000, e-0.25, 'm', cex=1)
  
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

plot(NULL, xlim=c(0,10), ylim = c(0,1), axes=FALSE, ann=FALSE)
legend(5, 0.5, c("5th","D1","D2"), fill=c('hotpink2', 'orange', 'dodgerblue'), border=NA, horiz=TRUE, xjust = 0.5, yjust = 0.5, box.col = NA)










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


#### write tables

load("sign_E_5th_D1D2_p05fc1_5thup.rda")
load("sign_M_5th_D1D2_p05fc1_5thup.rda")

load("sign_E_D1_5thD2_p05fc1_D1up.rda")
load("sign_M_D1_5thD2_p05fc1_D1up.rda")

load("sign_E_D2_5thD1_p05fc1_D2up.rda")
load("sign_M_D2_5thD1_p05fc1_D2up.rda")

load("sign_E_5th_D1D2_p05fc1_5thdown.rda")
load("sign_M_5th_D1D2_p05fc1_5thdown.rda")

load("sign_E_D1_5thD2_p05fc1_D1down.rda")
load("sign_M_D1_5thD2_p05fc1_D1down.rda")

load("sign_E_D2_5thD1_p05fc1_D2down.rda")
load("sign_M_D2_5thD1_p05fc1_D2down.rda")

head(sign_E_5th_D1D2_p05fc1)


sign_list_E_dir <- list(sign_E_5th_D1D2_5thup, sign_E_D1_5thD2_D1up, sign_E_D2_5thD1_D2up, sign_E_5th_D1D2_5thdown, sign_E_D1_5thD2_D1down, sign_E_D2_5thD1_D2down)
sign_list_M_dir <- list(sign_M_5th_D1D2_5thup, sign_M_D1_5thD2_D1up, sign_M_D2_5thD1_D2up, sign_M_5th_D1D2_5thdown, sign_M_D1_5thD2_D1down, sign_M_D2_5thD1_D2down)

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


write.table(sign_list_E_df_all_dir[[1]][,c("scafE","startE","endE")], file = 'development_5th_p05FC1_erato_5thup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[2]][,c("scafE","startE","endE")], file = 'development_D1_p05FC1_erato_D1up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[3]][,c("scafE","startE","endE")], file = 'development_D2_p05FC1_erato_D2up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[4]][,c("scafE","startE","endE")], file = 'development_5th_p05FC1_erato_5thdown.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[5]][,c("scafE","startE","endE")], file = 'development_D1_p05FC1_erato_D1down.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[6]][,c("scafE","startE","endE")], file = 'development_D2_p05FC1_erato_D2down.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("scafM","startM","endM")], file = 'development_5th_p05FC1_melp_5thup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[2]][,c("scafM","startM","endM")], file = 'development_D1_p05FC1_melp_D1up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[3]][,c("scafM","startM","endM")], file = 'development_D2_p05FC1_melp_D2up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[4]][,c("scafM","startM","endM")], file = 'development_5th_p05FC1_melp_5thdown.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[5]][,c("scafM","startM","endM")], file = 'development_D1_p05FC1_melp_D1down.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[6]][,c("scafM","startM","endM")], file = 'development_D2_p05FC1_melp_D2down.bed', col.names = F, row.names = F, quote = F)


write.table(subset(sign_list_E_df_all_dir[[1]], sign_list_E_df_all_dir[[1]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_5th_p05FC1_eratoUnique_5thup.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_all_dir[[2]], sign_list_E_df_all_dir[[2]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_D1_p05FC1_eratoUnique_D1up.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_all_dir[[3]], sign_list_E_df_all_dir[[3]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_D2_p05FC1_eratoUnique_D2up.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_all_dir[[4]], sign_list_E_df_all_dir[[4]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_5th_p05FC1_eratoUnique_5thdown.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_all_dir[[5]], sign_list_E_df_all_dir[[5]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_D1_p05FC1_eratoUnique_D1down.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_E_df_all_dir[[6]], sign_list_E_df_all_dir[[6]]$scafM == "NA")[,c("scafE","startE","endE")], file = 'development_D2_p05FC1_eratoUnique_D2down.bed', col.names = F, row.names = F, quote = F)

write.table(subset(sign_list_M_df_all_dir[[1]], sign_list_M_df_all_dir[[1]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_5th_p05FC1_melpUnique_5thup.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_all_dir[[2]], sign_list_M_df_all_dir[[2]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_D1_p05FC1_melpUnique_D1up.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_all_dir[[3]], sign_list_M_df_all_dir[[3]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_D2_p05FC1_melpUnique_D2up.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_all_dir[[4]], sign_list_M_df_all_dir[[4]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_5th_p05FC1_melpUnique_5thdown.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_all_dir[[5]], sign_list_M_df_all_dir[[5]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_D1_p05FC1_melpUnique_D1down.bed', col.names = F, row.names = F, quote = F)
write.table(subset(sign_list_M_df_all_dir[[6]], sign_list_M_df_all_dir[[6]]$scafE == "NA")[,c("scafM","startM","endM")], file = 'development_D2_p05FC1_melpUnique_D2down.bed', col.names = F, row.names = F, quote = F)


development_5thup_shared <- merge(sign_list_E_df_all_dir[[1]],sign_list_M_df_all_dir[[1]], by=c('start','end'))
development_D1up_shared <- merge(sign_list_E_df_all_dir[[2]],sign_list_M_df_all_dir[[2]], by=c('start','end'))
development_D2up_shared <- merge(sign_list_E_df_all_dir[[3]],sign_list_M_df_all_dir[[3]], by=c('start','end'))

development_5thdown_shared <- merge(sign_list_E_df_all_dir[[4]],sign_list_M_df_all_dir[[4]], by=c('start','end'))
development_D1down_shared <- merge(sign_list_E_df_all_dir[[5]],sign_list_M_df_all_dir[[5]], by=c('start','end'))
development_D2down_shared <- merge(sign_list_E_df_all_dir[[6]],sign_list_M_df_all_dir[[6]], by=c('start','end'))


write.table(development_5thup_shared[,c("scafE.x","startE.x","endE.x")], file = "development_5th_p05FC1_shared_5thup.bed", col.names = F, row.names = F, quote = F)
write.table(development_D1up_shared[,c("scafE.x","startE.x","endE.x")], file = "development_D1_p05FC1_shared_D1up.bed", col.names = F, row.names = F, quote = F)
write.table(development_D2up_shared[,c("scafE.x","startE.x","endE.x")], file = "development_D2_p05FC1_shared_D2up.bed", col.names = F, row.names = F, quote = F)

write.table(development_5thdown_shared[,c("scafE.x","startE.x","endE.x")], file = "development_5th_p05FC1_shared_5thdown.bed", col.names = F, row.names = F, quote = F)
write.table(development_D1down_shared[,c("scafE.x","startE.x","endE.x")], file = "development_D1_p05FC1_shared_D1down.bed", col.names = F, row.names = F, quote = F)
write.table(development_D2down_shared[,c("scafE.x","startE.x","endE.x")], file = "development_D2_p05FC1_shared_D2down.bed", col.names = F, row.names = F, quote = F)


write.table(development_5thup_shared[,c("scafM.x","startM.x","endM.x")], file = "development_5th_p05FC1_shared_5thup_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(development_D1up_shared[,c("scafM.x","startM.x","endM.x")], file = "development_D1_p05FC1_shared_D1up_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(development_D2up_shared[,c("scafM.x","startM.x","endM.x")], file = "development_D2_p05FC1_shared_D2up_melpCoords.bed", col.names = F, row.names = F, quote = F)

write.table(development_5thdown_shared[,c("scafM.x","startM.x","endM.x")], file = "development_5th_p05FC1_shared_5thdown_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(development_D1down_shared[,c("scafM.x","startM.x","endM.x")], file = "development_D1_p05FC1_shared_D1down_melpCoords.bed", col.names = F, row.names = F, quote = F)
write.table(development_D2down_shared[,c("scafM.x","startM.x","endM.x")], file = "development_D2_p05FC1_shared_D2down_melpCoords.bed", col.names = F, row.names = F, quote = F)



###############################
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

load('sign_E_5th_D1D2_p05fc1_ALL.rda')
load('sign_M_5th_D1D2_p05fc1_ALL.rda')

load('sign_E_D1_5thD2_p05fc1_ALL.rda')
load('sign_M_D1_5thD2_p05fc1_ALL.rda')

load('sign_E_D2_5thD1_p05fc1_ALL.rda')
load('sign_M_D2_5thD1_p05fc1_ALL.rda')

sign_E_5th_D1D2_5thup <- subset(res_E_5th_D1D2, res_E_5th_D1D2$padj < 0.05 & (res_E_5th_D1D2$log2FoldChange) < -1)
sign_M_5th_D1D2_5thup <- subset(res_M_5th_D1D2, res_M_5th_D1D2$padj < 0.05 & (res_M_5th_D1D2$log2FoldChange) < -1)

sign_E_5th_D1D2_5thdown <- subset(res_E_5th_D1D2, res_E_5th_D1D2$padj < 0.05 & (res_E_5th_D1D2$log2FoldChange) > 1)
sign_M_5th_D1D2_5thdown <- subset(res_M_5th_D1D2, res_M_5th_D1D2$padj < 0.05 & (res_M_5th_D1D2$log2FoldChange) > 1)

sign_E_D1_5thD2_D1up <- subset(res_E_D1_5thD2, res_E_D1_5thD2$padj < 0.05 & (res_E_D1_5thD2$log2FoldChange) > 1)
sign_M_D1_5thD2_D1up <- subset(res_M_D1_5thD2, res_M_D1_5thD2$padj < 0.05 & (res_M_D1_5thD2$log2FoldChange) > 1)

sign_E_D1_5thD2_D1down <- subset(res_E_D1_5thD2, res_E_D1_5thD2$padj < 0.05 & (res_E_D1_5thD2$log2FoldChange) < -1)
sign_M_D1_5thD2_D1down <- subset(res_M_D1_5thD2, res_M_D1_5thD2$padj < 0.05 & (res_M_D1_5thD2$log2FoldChange) < -1)

sign_E_D2_5thD1_D2up <- subset(res_E_D2_5thD1, res_E_D2_5thD1$padj < 0.05 & (res_E_D2_5thD1$log2FoldChange) > 1)
sign_M_D2_5thD1_D2up <- subset(res_M_D2_5thD1, res_M_D2_5thD1$padj < 0.05 & (res_M_D2_5thD1$log2FoldChange) > 1)

sign_E_D2_5thD1_D2down <- subset(res_E_D2_5thD1, res_E_D2_5thD1$padj < 0.05 & (res_E_D2_5thD1$log2FoldChange) < -1)
sign_M_D2_5thD1_D2down <- subset(res_M_D2_5thD1, res_M_D2_5thD1$padj < 0.05 & (res_M_D2_5thD1$log2FoldChange) < -1)


# 5th up

merged_EM_5thup <- merge(as.data.frame(res_E_5th_D1D2), as.data.frame(sign_M_5th_D1D2_5thup), by = 'row.names')
merged_ME_5thup <- merge(as.data.frame(res_M_5th_D1D2), as.data.frame(sign_E_5th_D1D2_5thup), by = 'row.names')



merged_shared_5thup <- merge(as.data.frame(sign_E_5th_D1D2_5thup), as.data.frame(sign_M_5th_D1D2_5thup), by = 'row.names')

merged_EM_5thup$log2FoldChange.x <- merged_EM_5thup$log2FoldChange.x*-1
merged_EM_5thup$log2FoldChange.y <- merged_EM_5thup$log2FoldChange.y*-1
merged_ME_5thup$log2FoldChange.x <- merged_ME_5thup$log2FoldChange.x*-1
merged_ME_5thup$log2FoldChange.y <- merged_ME_5thup$log2FoldChange.y*-1
merged_shared_5thup$log2FoldChange.x <- merged_shared_5thup$log2FoldChange.x*-1
merged_shared_5thup$log2FoldChange.y <- merged_shared_5thup$log2FoldChange.y*-1

merged_EM_5thup_unique <- merged_EM_5thup[!merged_EM_5thup$Row.names %in% merged_shared_5thup$Row.names, ]
merged_ME_5thup_unique <- merged_ME_5thup[!merged_ME_5thup$Row.names %in% merged_shared_5thup$Row.names, ]

par(mfrow=c(3,2))

plot(merged_EM_5thup_unique$log2FoldChange.x, merged_EM_5thup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "5th instar up")
par(new=T)
plot(merged_ME_5thup_unique$log2FoldChange.x, merged_ME_5thup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_5thup$log2FoldChange.x, merged_shared_5thup$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_5thup_unique$log2FoldChange.y ~ merged_EM_5thup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_5thup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_5thup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_5thup_unique$log2FoldChange.y ~ merged_ME_5thup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_5thup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_5thup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_5thup$log2FoldChange.y ~ merged_shared_5thup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_5thup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_5thup$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                      as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                      "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_5thup_unique$log2FoldChange.y ~ merged_EM_5thup_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_5thup_unique$log2FoldChange.y ~ merged_ME_5thup_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_5thup$log2FoldChange.y ~ merged_shared_5thup$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_5th_D1D2_5thup$Row.names <- rownames(sign_E_5th_D1D2_5thup)
sign_E_5th_D1D2_5thup <- separate(data = as.data.frame(sign_E_5th_D1D2_5thup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_5th_D1D2_5thup[,c(2:5,7,8)] <- sapply(sign_E_5th_D1D2_5thup[,c(2:5,7,8)],as.numeric) 
sign_E_5th_D1D2_5thup$Genome <- 'pan'

write.table(sign_E_5th_D1D2_5thup[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_D1D2_5thup.bed', col.names = F, row.names = F, quote = F)

sign_M_5th_D1D2_5thup$Row.names <- rownames(sign_M_5th_D1D2_5thup)
sign_M_5th_D1D2_5thup <- separate(data = as.data.frame(sign_M_5th_D1D2_5thup), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_D1D2_5thup[,c(2:5,7,8)] <- sapply(sign_M_5th_D1D2_5thup[,c(2:5,7,8)],as.numeric) 
sign_M_5th_D1D2_5thup$Genome <- 'pan'

sign_M_5th_D1D2_5thup_write <- sign_M_5th_D1D2_5thup[,c(18,10,11)]
sign_M_5th_D1D2_5thup_write <- subset(sign_M_5th_D1D2_5thup_write, as.numeric(sign_M_5th_D1D2_5thup_write$end2)-as.numeric(sign_M_5th_D1D2_5thup_write$start2) < 10000)
write.table(sign_M_5th_D1D2_5thup_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_D1D2_5thup.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_5thup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_D1D2_5thup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_D1D2_5thup_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_D1D2_5thup_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_D1D2_5thup_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_5th_D1D2_5thup)
nrow(sign_M_5th_D1D2_5thup)
nrow(merged_shared_5thup)

#
nrow(sign_E_5th_D1D2_5thup)-nrow(merged_ME_5thup[grep("NA", merged_ME_5thup$Row.names), ])
nrow(sign_M_5th_D1D2_5thup)-nrow(merged_EM_5thup[grep("NA", merged_EM_5thup$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_5th_D1D2_5thup_bed <- sign_E_5th_D1D2_5thup[,c(18,8,9)]
colnames(sign_E_5th_D1D2_5thup_bed) <- c('chrom','start','end')
sign_E_5th_D1D2_5thup_bed$start <- as.numeric(sign_E_5th_D1D2_5thup_bed$start)
sign_E_5th_D1D2_5thup_bed$end <- as.numeric(sign_E_5th_D1D2_5thup_bed$end)

sign_M_5th_D1D2_5thup_bed <- sign_M_5th_D1D2_5thup_write
colnames(sign_M_5th_D1D2_5thup_bed) <- c('chrom','start','end')
sign_M_5th_D1D2_5thup_bed$start <- as.numeric(sign_M_5th_D1D2_5thup_bed$start)
sign_M_5th_D1D2_5thup_bed$end <- as.numeric(sign_M_5th_D1D2_5thup_bed$end)

merged_shared_5thup_bed <- merged_shared_5thup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_5thup_bed) <- c('chrom','start','end')
merged_shared_5thup_bed$start <- as.numeric(merged_shared_5thup_bed$start)
merged_shared_5thup_bed$end <- as.numeric(merged_shared_5thup_bed$end)

sign_E_5th_D1D2_5thup_bed_match <- bed_intersect(sign_E_5th_D1D2_5thup_bed,conserved_E)
sign_M_5th_D1D2_5thup_bed_match <- bed_intersect(sign_M_5th_D1D2_5thup_bed,conserved_M)
merged_shared_5thup_bed_match <- bed_intersect(merged_shared_5thup_bed,conserved_S)

nrow(sign_E_5th_D1D2_5thup_bed_match)
nrow(sign_M_5th_D1D2_5thup_bed_match)
nrow(merged_shared_5thup_bed_match)

sign_E_5th_D1D2_5thup_cons <- merge(sign_E_5th_D1D2_5thup, sign_E_5th_D1D2_5thup_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_5th_D1D2_5thup_cons <- merge(sign_M_5th_D1D2_5thup, sign_M_5th_D1D2_5thup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_5thup_cons <- merge(merged_ME_5thup, sign_E_5th_D1D2_5thup_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_5thup_cons <- merge(merged_EM_5thup, sign_M_5th_D1D2_5thup_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_5th_D1D2_5thup_cons)-nrow(merged_ME_5thup_cons[grep("NA", merged_ME_5thup_cons$Row.names), ])
nrow(sign_M_5th_D1D2_5thup_cons)-nrow(merged_EM_5thup_cons[grep("NA", merged_EM_5thup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_5th_D1D2_5thup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_5th_D1D2_5thup_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_5thup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_5th_D1D2_5thup_bed_match)
nrow(sign_M_5th_D1D2_5thup_bed_match)
nrow(merged_shared_5thup_bed_match)

nrow(sign_E_5th_D1D2_5thup)-nrow(merged_ME_5thup[grep("NA", merged_ME_5thup$Row.names), ])
nrow(sign_M_5th_D1D2_5thup)-nrow(merged_EM_5thup[grep("NA", merged_EM_5thup$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))



#########################


# 5th down

merged_EM_5thdown <- merge(as.data.frame(res_E_5th_D1D2), as.data.frame(sign_M_5th_D1D2_5thdown), by = 'row.names')
merged_ME_5thdown <- merge(as.data.frame(res_M_5th_D1D2), as.data.frame(sign_E_5th_D1D2_5thdown), by = 'row.names')

merged_shared_5thdown <- merge(as.data.frame(sign_E_5th_D1D2_5thdown), as.data.frame(sign_M_5th_D1D2_5thdown), by = 'row.names')

merged_EM_5thdown_unique <- merged_EM_5thdown[!merged_EM_5thdown$Row.names %in% merged_shared_5thdown$Row.names, ]
merged_ME_5thdown_unique <- merged_ME_5thdown[!merged_ME_5thdown$Row.names %in% merged_shared_5thdown$Row.names, ]

plot(merged_EM_5thdown_unique$log2FoldChange.x, merged_EM_5thdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "5th instar down")
par(new=T)
plot(merged_ME_5thdown_unique$log2FoldChange.x, merged_ME_5thdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_5thdown$log2FoldChange.x, merged_shared_5thdown$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_5thdown_unique$log2FoldChange.y ~ merged_EM_5thdown_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_5thdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_5thdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_5thdown_unique$log2FoldChange.y ~ merged_ME_5thdown_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_5thdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_5thdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_5thdown$log2FoldChange.y ~ merged_shared_5thdown$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_5thdown$log2FoldChange.x)), x2 = max(na.omit(merged_shared_5thdown$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                        as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                        "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

######################## % completely missing

sign_E_5th_D1D2_5thdown$Row.names <- rownames(sign_E_5th_D1D2_5thdown)
sign_E_5th_D1D2_5thdown <- separate(data = as.data.frame(sign_E_5th_D1D2_5thdown), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_5th_D1D2_5thdown[,c(2:5,7,8)] <- sapply(sign_E_5th_D1D2_5thdown[,c(2:5,7,8)],as.numeric) 
sign_E_5th_D1D2_5thdown$Genome <- 'pan'

write.table(sign_E_5th_D1D2_5thdown[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_D1D2_5thdown.bed', col.names = F, row.names = F, quote = F)

sign_M_5th_D1D2_5thdown$Row.names <- rownames(sign_M_5th_D1D2_5thdown)
sign_M_5th_D1D2_5thdown <- separate(data = as.data.frame(sign_M_5th_D1D2_5thdown), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_5th_D1D2_5thdown[,c(2:5,7,8)] <- sapply(sign_M_5th_D1D2_5thdown[,c(2:5,7,8)],as.numeric) 
sign_M_5th_D1D2_5thdown$Genome <- 'pan'

sign_M_5th_D1D2_5thdown_write <- sign_M_5th_D1D2_5thdown[,c(18,10,11)]
sign_M_5th_D1D2_5thdown_write <- subset(sign_M_5th_D1D2_5thdown_write, as.numeric(sign_M_5th_D1D2_5thdown_write$end2)-as.numeric(sign_M_5th_D1D2_5thdown_write$start2) < 10000)
write.table(sign_M_5th_D1D2_5thdown_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_D1D2_5thdown.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_5thdown[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_D1D2_5thdown.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_5th_D1D2_5thdown_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_5th_D1D2_5thdown_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_5th_D1D2_5thdown_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_5th_D1D2_5thdown)
nrow(sign_M_5th_D1D2_5thdown)
nrow(merged_shared_5thdown)

#
nrow(sign_E_5th_D1D2_5thdown)-nrow(merged_ME_5thdown[grep("NA", merged_ME_5thdown$Row.names), ])
nrow(sign_M_5th_D1D2_5thdown)-nrow(merged_EM_5thdown[grep("NA", merged_EM_5thdown$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_5th_D1D2_5thdown_bed <- sign_E_5th_D1D2_5thdown[,c(18,8,9)]
colnames(sign_E_5th_D1D2_5thdown_bed) <- c('chrom','start','end')
sign_E_5th_D1D2_5thdown_bed$start <- as.numeric(sign_E_5th_D1D2_5thdown_bed$start)
sign_E_5th_D1D2_5thdown_bed$end <- as.numeric(sign_E_5th_D1D2_5thdown_bed$end)

sign_M_5th_D1D2_5thdown_bed <- sign_M_5th_D1D2_5thdown_write
colnames(sign_M_5th_D1D2_5thdown_bed) <- c('chrom','start','end')
sign_M_5th_D1D2_5thdown_bed$start <- as.numeric(sign_M_5th_D1D2_5thdown_bed$start)
sign_M_5th_D1D2_5thdown_bed$end <- as.numeric(sign_M_5th_D1D2_5thdown_bed$end)

merged_shared_5thdown_bed <- merged_shared_5thdown[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_5thdown_bed) <- c('chrom','start','end')
merged_shared_5thdown_bed$start <- as.numeric(merged_shared_5thdown_bed$start)
merged_shared_5thdown_bed$end <- as.numeric(merged_shared_5thdown_bed$end)

sign_E_5th_D1D2_5thdown_bed_match <- bed_intersect(sign_E_5th_D1D2_5thdown_bed,conserved_E)
sign_M_5th_D1D2_5thdown_bed_match <- bed_intersect(sign_M_5th_D1D2_5thdown_bed,conserved_M)
merged_shared_5thdown_bed_match <- bed_intersect(merged_shared_5thdown_bed,conserved_S)


sign_E_5th_D1D2_5thdown_cons <- merge(sign_E_5th_D1D2_5thdown, sign_E_5th_D1D2_5thdown_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_5th_D1D2_5thdown_cons <- merge(sign_M_5th_D1D2_5thdown, sign_M_5th_D1D2_5thdown_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_5thdown_cons <- merge(merged_ME_5thdown, sign_E_5th_D1D2_5thdown_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_5thdown_cons <- merge(merged_EM_5thdown, sign_M_5th_D1D2_5thdown_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_5th_D1D2_5thdown_cons)-nrow(merged_ME_5thdown_cons[grep("NA", merged_ME_5thdown_cons$Row.names), ])
nrow(sign_M_5th_D1D2_5thdown_cons)-nrow(merged_EM_5thdown_cons[grep("NA", merged_EM_5thdown_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_5th_D1D2_5thdown_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_5th_D1D2_5thdown_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_5thdown_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_5th_D1D2_5thdown_bed_match)
nrow(sign_M_5th_D1D2_5thdown_bed_match)
nrow(merged_shared_5thdown_bed_match)

nrow(sign_E_5th_D1D2_5thdown_cons)-nrow(merged_ME_5thdown_cons[grep("NA", merged_ME_5thdown_cons$Row.names), ])
nrow(sign_M_5th_D1D2_5thdown_cons)-nrow(merged_EM_5thdown_cons[grep("NA", merged_EM_5thdown_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))


#########################

# D1 up

merged_EM_D1up <- merge(as.data.frame(res_E_D1_5thD2), as.data.frame(sign_M_D1_5thD2_D1up), by = 'row.names')
merged_ME_D1up <- merge(as.data.frame(res_M_D1_5thD2), as.data.frame(sign_E_D1_5thD2_D1up), by = 'row.names')



merged_shared_D1up <- merge(as.data.frame(sign_E_D1_5thD2_D1up), as.data.frame(sign_M_D1_5thD2_D1up), by = 'row.names')


merged_EM_D1up_unique <- merged_EM_D1up[!merged_EM_D1up$Row.names %in% merged_shared_D1up$Row.names, ]
merged_ME_D1up_unique <- merged_ME_D1up[!merged_ME_D1up$Row.names %in% merged_shared_D1up$Row.names, ]


plot(merged_EM_D1up_unique$log2FoldChange.x, merged_EM_D1up_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h pupa up")
par(new=T)
plot(merged_ME_D1up_unique$log2FoldChange.x, merged_ME_D1up_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D1up$log2FoldChange.x, merged_shared_D1up$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_D1up_unique$log2FoldChange.y ~ merged_EM_D1up_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D1up_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D1up_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D1up_unique$log2FoldChange.y ~ merged_ME_D1up_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D1up_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D1up_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D1up$log2FoldChange.y ~ merged_shared_D1up$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D1up$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D1up$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                        as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                        "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

######################## % completely missing

sign_E_D1_5thD2_D1up$Row.names <- rownames(sign_E_D1_5thD2_D1up)
sign_E_D1_5thD2_D1up <- separate(data = as.data.frame(sign_E_D1_5thD2_D1up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_5thD2_D1up[,c(2:5,7,8)] <- sapply(sign_E_D1_5thD2_D1up[,c(2:5,7,8)],as.numeric) 
sign_E_D1_5thD2_D1up$Genome <- 'pan'

write.table(sign_E_D1_5thD2_D1up[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_5thD2_D1up.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_5thD2_D1up$Row.names <- rownames(sign_M_D1_5thD2_D1up)
sign_M_D1_5thD2_D1up <- separate(data = as.data.frame(sign_M_D1_5thD2_D1up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_5thD2_D1up[,c(2:5,7,8)] <- sapply(sign_M_D1_5thD2_D1up[,c(2:5,7,8)],as.numeric) 
sign_M_D1_5thD2_D1up$Genome <- 'pan'

sign_M_D1_5thD2_D1up_write <- sign_M_D1_5thD2_D1up[,c(18,10,11)]
sign_M_D1_5thD2_D1up_write <- subset(sign_M_D1_5thD2_D1up_write, as.numeric(sign_M_D1_5thD2_D1up_write$end2)-as.numeric(sign_M_D1_5thD2_D1up_write$start2) < 10000)
write.table(sign_M_D1_5thD2_D1up_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_5thD2_D1up.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_D1up[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_5thD2_D1up.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_5thD2_D1up_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_5thD2_D1up_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_5thD2_D1up_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_5thD2_D1up)
nrow(sign_M_D1_5thD2_D1up)
nrow(merged_shared_D1up)

#
nrow(sign_E_D1_5thD2_D1up)-nrow(merged_ME_D1up[grep("NA", merged_ME_D1up$Row.names), ])
nrow(sign_M_D1_5thD2_D1up)-nrow(merged_EM_D1up[grep("NA", merged_EM_D1up$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))


######################## Calculate for conserved set
sign_E_D1_5thD2_D1up_bed <- sign_E_D1_5thD2_D1up[,c(18,8,9)]
colnames(sign_E_D1_5thD2_D1up_bed) <- c('chrom','start','end')
sign_E_D1_5thD2_D1up_bed$start <- as.numeric(sign_E_D1_5thD2_D1up_bed$start)
sign_E_D1_5thD2_D1up_bed$end <- as.numeric(sign_E_D1_5thD2_D1up_bed$end)

sign_M_D1_5thD2_D1up_bed <- sign_M_D1_5thD2_D1up_write
colnames(sign_M_D1_5thD2_D1up_bed) <- c('chrom','start','end')
sign_M_D1_5thD2_D1up_bed$start <- as.numeric(sign_M_D1_5thD2_D1up_bed$start)
sign_M_D1_5thD2_D1up_bed$end <- as.numeric(sign_M_D1_5thD2_D1up_bed$end)

merged_shared_D1up_bed <- merged_shared_D1up[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D1up_bed) <- c('chrom','start','end')
merged_shared_D1up_bed$start <- as.numeric(merged_shared_D1up_bed$start)
merged_shared_D1up_bed$end <- as.numeric(merged_shared_D1up_bed$end)

sign_E_D1_5thD2_D1up_bed_match <- bed_intersect(sign_E_D1_5thD2_D1up_bed,conserved_E)
sign_M_D1_5thD2_D1up_bed_match <- bed_intersect(sign_M_D1_5thD2_D1up_bed,conserved_M)
merged_shared_D1up_bed_match <- bed_intersect(merged_shared_D1up_bed,conserved_S)


sign_E_D1_5thD2_D1up_cons <- merge(sign_E_D1_5thD2_D1up, sign_E_D1_5thD2_D1up_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_5thD2_D1up_cons <- merge(sign_M_D1_5thD2_D1up, sign_M_D1_5thD2_D1up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D1up_cons <- merge(merged_ME_D1up, sign_E_D1_5thD2_D1up_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D1up_cons <- merge(merged_EM_D1up, sign_M_D1_5thD2_D1up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_D1_5thD2_D1up_cons)-nrow(merged_ME_D1up_cons[grep("NA", merged_ME_D1up_cons$Row.names), ])
nrow(sign_M_D1_5thD2_D1up_cons)-nrow(merged_EM_D1up_cons[grep("NA", merged_EM_D1up_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_5thD2_D1up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_5thD2_D1up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D1up_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_5thD2_D1up_bed_match)
nrow(sign_M_D1_5thD2_D1up_bed_match)
nrow(merged_shared_D1up_bed_match)

nrow(sign_E_D1_5thD2_D1up_cons)-nrow(merged_ME_D1up_cons[grep("NA", merged_ME_D1up_cons$Row.names), ])
nrow(sign_M_D1_5thD2_D1up_cons)-nrow(merged_EM_D1up_cons[grep("NA", merged_EM_D1up_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))
#########################


# D1 down

merged_EM_D1down <- merge(as.data.frame(res_E_D1_5thD2), as.data.frame(sign_M_D1_5thD2_D1down), by = 'row.names')
merged_ME_D1down <- merge(as.data.frame(res_M_D1_5thD2), as.data.frame(sign_E_D1_5thD2_D1down), by = 'row.names')

merged_shared_D1down <- merge(as.data.frame(sign_E_D1_5thD2_D1down), as.data.frame(sign_M_D1_5thD2_D1down), by = 'row.names')

merged_EM_D1down$log2FoldChange.x <- merged_EM_D1down$log2FoldChange.x*-1
merged_EM_D1down$log2FoldChange.y <- merged_EM_D1down$log2FoldChange.y*-1
merged_ME_D1down$log2FoldChange.x <- merged_ME_D1down$log2FoldChange.x*-1
merged_ME_D1down$log2FoldChange.y <- merged_ME_D1down$log2FoldChange.y*-1
merged_shared_D1down$log2FoldChange.x <- merged_shared_D1down$log2FoldChange.x*-1
merged_shared_D1down$log2FoldChange.y <- merged_shared_D1down$log2FoldChange.y*-1

merged_EM_D1down_unique <- merged_EM_D1down[!merged_EM_D1down$Row.names %in% merged_shared_D1down$Row.names, ]
merged_ME_D1down_unique <- merged_ME_D1down[!merged_ME_D1down$Row.names %in% merged_shared_D1down$Row.names, ]



plot(merged_EM_D1down_unique$log2FoldChange.x, merged_EM_D1down_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h pupa down")
par(new=T)
plot(merged_ME_D1down_unique$log2FoldChange.x, merged_ME_D1down_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D1down$log2FoldChange.x, merged_shared_D1down$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_D1down_unique$log2FoldChange.y ~ merged_EM_D1down_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D1down_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D1down_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D1down_unique$log2FoldChange.y ~ merged_ME_D1down_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D1down_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D1down_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D1down$log2FoldChange.y ~ merged_shared_D1down$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D1down$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D1down$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                        as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                        "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

######################## % completely missing

sign_E_D1_5thD2_D1down$Row.names <- rownames(sign_E_D1_5thD2_D1down)
sign_E_D1_5thD2_D1down <- separate(data = as.data.frame(sign_E_D1_5thD2_D1down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_5thD2_D1down[,c(2:5,7,8)] <- sapply(sign_E_D1_5thD2_D1down[,c(2:5,7,8)],as.numeric) 
sign_E_D1_5thD2_D1down$Genome <- 'pan'

write.table(sign_E_D1_5thD2_D1down[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_5thD2_D1down.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_5thD2_D1down$Row.names <- rownames(sign_M_D1_5thD2_D1down)
sign_M_D1_5thD2_D1down <- separate(data = as.data.frame(sign_M_D1_5thD2_D1down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_5thD2_D1down[,c(2:5,7,8)] <- sapply(sign_M_D1_5thD2_D1down[,c(2:5,7,8)],as.numeric) 
sign_M_D1_5thD2_D1down$Genome <- 'pan'

sign_M_D1_5thD2_D1down_write <- sign_M_D1_5thD2_D1down[,c(18,10,11)]
sign_M_D1_5thD2_D1down_write <- subset(sign_M_D1_5thD2_D1down_write, as.numeric(sign_M_D1_5thD2_D1down_write$end2)-as.numeric(sign_M_D1_5thD2_D1down_write$start2) < 10000)
write.table(sign_M_D1_5thD2_D1down_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_5thD2_D1down.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_D1down[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_5thD2_D1down.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_5thD2_D1down_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_5thD2_D1down_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D1_5thD2_D1down_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_5thD2_D1down)
nrow(sign_M_D1_5thD2_D1down)
nrow(merged_shared_D1down)

#
nrow(sign_E_D1_5thD2_D1down)-nrow(merged_ME_D1down[grep("NA", merged_ME_D1down$Row.names), ])
nrow(sign_M_D1_5thD2_D1down)-nrow(merged_EM_D1down[grep("NA", merged_EM_D1down$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_5thD2_D1down_bed <- sign_E_D1_5thD2_D1down[,c(18,8,9)]
colnames(sign_E_D1_5thD2_D1down_bed) <- c('chrom','start','end')
sign_E_D1_5thD2_D1down_bed$start <- as.numeric(sign_E_D1_5thD2_D1down_bed$start)
sign_E_D1_5thD2_D1down_bed$end <- as.numeric(sign_E_D1_5thD2_D1down_bed$end)

sign_M_D1_5thD2_D1down_bed <- sign_M_D1_5thD2_D1down_write
colnames(sign_M_D1_5thD2_D1down_bed) <- c('chrom','start','end')
sign_M_D1_5thD2_D1down_bed$start <- as.numeric(sign_M_D1_5thD2_D1down_bed$start)
sign_M_D1_5thD2_D1down_bed$end <- as.numeric(sign_M_D1_5thD2_D1down_bed$end)

merged_shared_D1down_bed <- merged_shared_D1down[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D1down_bed) <- c('chrom','start','end')
merged_shared_D1down_bed$start <- as.numeric(merged_shared_D1down_bed$start)
merged_shared_D1down_bed$end <- as.numeric(merged_shared_D1down_bed$end)

sign_E_D1_5thD2_D1down_bed_match <- bed_intersect(sign_E_D1_5thD2_D1down_bed,conserved_E)
sign_M_D1_5thD2_D1down_bed_match <- bed_intersect(sign_M_D1_5thD2_D1down_bed,conserved_M)
merged_shared_D1down_bed_match <- bed_intersect(merged_shared_D1down_bed,conserved_S)


sign_E_D1_5thD2_D1down_cons <- merge(sign_E_D1_5thD2_D1down, sign_E_D1_5thD2_D1down_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_5thD2_D1down_cons <- merge(sign_M_D1_5thD2_D1down, sign_M_D1_5thD2_D1down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D1down_cons <- merge(merged_ME_D1down, sign_E_D1_5thD2_D1down_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D1down_cons <- merge(merged_EM_D1down, sign_M_D1_5thD2_D1down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_D1_5thD2_D1down_cons)-nrow(merged_ME_D1down_cons[grep("NA", merged_ME_D1down_cons$Row.names), ])
nrow(sign_M_D1_5thD2_D1down_cons)-nrow(merged_EM_D1down_cons[grep("NA", merged_EM_D1down_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_5thD2_D1down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_5thD2_D1down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D1down_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_5thD2_D1down_bed_match)
nrow(sign_M_D1_5thD2_D1down_bed_match)
nrow(merged_shared_D1down_bed_match)

nrow(sign_E_D1_5thD2_D1down_cons)-nrow(merged_ME_D1down_cons[grep("NA", merged_ME_D1down_cons$Row.names), ])
nrow(sign_M_D1_5thD2_D1down_cons)-nrow(merged_EM_D1down_cons[grep("NA", merged_EM_D1down_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))
#########################
#########################

#########################

# D2 up

merged_EM_D2up <- merge(as.data.frame(res_E_D2_5thD1), as.data.frame(sign_M_D2_5thD1_D2up), by = 'row.names')
merged_ME_D2up <- merge(as.data.frame(res_M_D2_5thD1), as.data.frame(sign_E_D2_5thD1_D2up), by = 'row.names')



merged_shared_D2up <- merge(as.data.frame(sign_E_D2_5thD1_D2up), as.data.frame(sign_M_D2_5thD1_D2up), by = 'row.names')


merged_EM_D2up_unique <- merged_EM_D2up[!merged_EM_D2up$Row.names %in% merged_shared_D2up$Row.names, ]
merged_ME_D2up_unique <- merged_ME_D2up[!merged_ME_D2up$Row.names %in% merged_shared_D2up$Row.names, ]


plot(merged_EM_D2up_unique$log2FoldChange.x, merged_EM_D2up_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "60h pupa up")
par(new=T)
plot(merged_ME_D2up_unique$log2FoldChange.x, merged_ME_D2up_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D2up$log2FoldChange.x, merged_shared_D2up$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_D2up_unique$log2FoldChange.y ~ merged_EM_D2up_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D2up_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D2up_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D2up_unique$log2FoldChange.y ~ merged_ME_D2up_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D2up_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D2up_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D2up$log2FoldChange.y ~ merged_shared_D2up$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D2up$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D2up$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                        as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                        "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

######################## % completely missing

sign_E_D2_5thD1_D2up$Row.names <- rownames(sign_E_D2_5thD1_D2up)
sign_E_D2_5thD1_D2up <- separate(data = as.data.frame(sign_E_D2_5thD1_D2up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_5thD1_D2up[,c(2:5,7,8)] <- sapply(sign_E_D2_5thD1_D2up[,c(2:5,7,8)],as.numeric) 
sign_E_D2_5thD1_D2up$Genome <- 'pan'

write.table(sign_E_D2_5thD1_D2up[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_5thD1_D2up.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_5thD1_D2up$Row.names <- rownames(sign_M_D2_5thD1_D2up)
sign_M_D2_5thD1_D2up <- separate(data = as.data.frame(sign_M_D2_5thD1_D2up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_5thD1_D2up[,c(2:5,7,8)] <- sapply(sign_M_D2_5thD1_D2up[,c(2:5,7,8)],as.numeric) 
sign_M_D2_5thD1_D2up$Genome <- 'pan'

sign_M_D2_5thD1_D2up_write <- sign_M_D2_5thD1_D2up[,c(18,10,11)]
sign_M_D2_5thD1_D2up_write <- subset(sign_M_D2_5thD1_D2up_write, as.numeric(sign_M_D2_5thD1_D2up_write$end2)-as.numeric(sign_M_D2_5thD1_D2up_write$start2) < 10000)
write.table(sign_M_D2_5thD1_D2up_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_5thD1_D2up.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_D2up[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_5thD1_D2up.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_5thD1_D2up_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_5thD1_D2up_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_5thD1_D2up_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D2_5thD1_D2up)
nrow(sign_M_D2_5thD1_D2up)
nrow(merged_shared_D2up)

#
nrow(sign_E_D2_5thD1_D2up)-nrow(merged_ME_D2up[grep("NA", merged_ME_D2up$Row.names), ])
nrow(sign_M_D2_5thD1_D2up)-nrow(merged_EM_D2up[grep("NA", merged_EM_D2up$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))


######################## Calculate for conserved set
sign_E_D2_5thD1_D2up_bed <- sign_E_D2_5thD1_D2up[,c(18,8,9)]
colnames(sign_E_D2_5thD1_D2up_bed) <- c('chrom','start','end')
sign_E_D2_5thD1_D2up_bed$start <- as.numeric(sign_E_D2_5thD1_D2up_bed$start)
sign_E_D2_5thD1_D2up_bed$end <- as.numeric(sign_E_D2_5thD1_D2up_bed$end)

sign_M_D2_5thD1_D2up_bed <- sign_M_D2_5thD1_D2up_write
colnames(sign_M_D2_5thD1_D2up_bed) <- c('chrom','start','end')
sign_M_D2_5thD1_D2up_bed$start <- as.numeric(sign_M_D2_5thD1_D2up_bed$start)
sign_M_D2_5thD1_D2up_bed$end <- as.numeric(sign_M_D2_5thD1_D2up_bed$end)

merged_shared_D2up_bed <- merged_shared_D2up[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D2up_bed) <- c('chrom','start','end')
merged_shared_D2up_bed$start <- as.numeric(merged_shared_D2up_bed$start)
merged_shared_D2up_bed$end <- as.numeric(merged_shared_D2up_bed$end)

sign_E_D2_5thD1_D2up_bed_match <- bed_intersect(sign_E_D2_5thD1_D2up_bed,conserved_E)
sign_M_D2_5thD1_D2up_bed_match <- bed_intersect(sign_M_D2_5thD1_D2up_bed,conserved_M)
merged_shared_D2up_bed_match <- bed_intersect(merged_shared_D2up_bed,conserved_S)


sign_E_D2_5thD1_D2up_cons <- merge(sign_E_D2_5thD1_D2up, sign_E_D2_5thD1_D2up_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D2_5thD1_D2up_cons <- merge(sign_M_D2_5thD1_D2up, sign_M_D2_5thD1_D2up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D2up_cons <- merge(merged_ME_D2up, sign_E_D2_5thD1_D2up_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D2up_cons <- merge(merged_EM_D2up, sign_M_D2_5thD1_D2up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_D2_5thD1_D2up_cons)-nrow(merged_ME_D2up_cons[grep("NA", merged_ME_D2up_cons$Row.names), ])
nrow(sign_M_D2_5thD1_D2up_cons)-nrow(merged_EM_D2up_cons[grep("NA", merged_EM_D2up_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D2_5thD1_D2up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D2_5thD1_D2up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D2up_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D2_5thD1_D2up_bed_match)
nrow(sign_M_D2_5thD1_D2up_bed_match)
nrow(merged_shared_D2up_bed_match)

nrow(sign_E_D2_5thD1_D2up_cons)-nrow(merged_ME_D2up_cons[grep("NA", merged_ME_D2up_cons$Row.names), ])
nrow(sign_M_D2_5thD1_D2up_cons)-nrow(merged_EM_D2up_cons[grep("NA", merged_EM_D2up_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################

# D2 down

merged_EM_D2down <- merge(as.data.frame(res_E_D2_5thD1), as.data.frame(sign_M_D2_5thD1_D2down), by = 'row.names')
merged_ME_D2down <- merge(as.data.frame(res_M_D2_5thD1), as.data.frame(sign_E_D2_5thD1_D2down), by = 'row.names')

merged_shared_D2down <- merge(as.data.frame(sign_E_D2_5thD1_D2down), as.data.frame(sign_M_D2_5thD1_D2down), by = 'row.names')

merged_EM_D2down$log2FoldChange.x <- merged_EM_D2down$log2FoldChange.x*-1
merged_EM_D2down$log2FoldChange.y <- merged_EM_D2down$log2FoldChange.y*-1
merged_ME_D2down$log2FoldChange.x <- merged_ME_D2down$log2FoldChange.x*-1
merged_ME_D2down$log2FoldChange.y <- merged_ME_D2down$log2FoldChange.y*-1
merged_shared_D2down$log2FoldChange.x <- merged_shared_D2down$log2FoldChange.x*-1
merged_shared_D2down$log2FoldChange.y <- merged_shared_D2down$log2FoldChange.y*-1

merged_EM_D2down_unique <- merged_EM_D2down[!merged_EM_D2down$Row.names %in% merged_shared_D2down$Row.names, ]
merged_ME_D2down_unique <- merged_ME_D2down[!merged_ME_D2down$Row.names %in% merged_shared_D2down$Row.names, ]


plot(merged_EM_D2down_unique$log2FoldChange.x, merged_EM_D2down_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('blue', alpha.f = 0.1), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "60h pupa down")
par(new=T)
plot(merged_ME_D2down_unique$log2FoldChange.x, merged_ME_D2down_unique$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('orange', alpha.f = 0.1), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_D2down$log2FoldChange.x, merged_shared_D2down$log2FoldChange.y, pch = 19, xlim = c(-1,7), ylim = c(0,8), col = adjustcolor('black', alpha.f = 0.1), xlab = '', ylab = '')


ablineclip(lm(merged_EM_D2down_unique$log2FoldChange.y ~ merged_EM_D2down_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_D2down_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_D2down_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_D2down_unique$log2FoldChange.y ~ merged_ME_D2down_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_D2down_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_D2down_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_D2down$log2FoldChange.y ~ merged_shared_D2down$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_D2down$log2FoldChange.x)), x2 = max(na.omit(merged_shared_D2down$log2FoldChange.x)))

legend(-1, 8, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                        as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                        "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

######################## % completely missing

sign_E_D2_5thD1_D2down$Row.names <- rownames(sign_E_D2_5thD1_D2down)
sign_E_D2_5thD1_D2down <- separate(data = as.data.frame(sign_E_D2_5thD1_D2down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_5thD1_D2down[,c(2:5,7,8)] <- sapply(sign_E_D2_5thD1_D2down[,c(2:5,7,8)],as.numeric) 
sign_E_D2_5thD1_D2down$Genome <- 'pan'

write.table(sign_E_D2_5thD1_D2down[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_5thD1_D2down.bed', col.names = F, row.names = F, quote = F)

sign_M_D2_5thD1_D2down$Row.names <- rownames(sign_M_D2_5thD1_D2down)
sign_M_D2_5thD1_D2down <- separate(data = as.data.frame(sign_M_D2_5thD1_D2down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D2_5thD1_D2down[,c(2:5,7,8)] <- sapply(sign_M_D2_5thD1_D2down[,c(2:5,7,8)],as.numeric) 
sign_M_D2_5thD1_D2down$Genome <- 'pan'

sign_M_D2_5thD1_D2down_write <- sign_M_D2_5thD1_D2down[,c(18,10,11)]
sign_M_D2_5thD1_D2down_write <- subset(sign_M_D2_5thD1_D2down_write, as.numeric(sign_M_D2_5thD1_D2down_write$end2)-as.numeric(sign_M_D2_5thD1_D2down_write$start2) < 10000)
write.table(sign_M_D2_5thD1_D2down_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_5thD1_D2down.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_D2down[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_5thD1_D2down.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D2_5thD1_D2down_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D2_5thD1_D2down_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_D2_5thD1_D2down_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D2_5thD1_D2down)
nrow(sign_M_D2_5thD1_D2down)
nrow(merged_shared_D2down)

#
nrow(sign_E_D2_5thD1_D2down)-nrow(merged_ME_D2down[grep("NA", merged_ME_D2down$Row.names), ])
nrow(sign_M_D2_5thD1_D2down)-nrow(merged_EM_D2down[grep("NA", merged_EM_D2down$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D2_5thD1_D2down_bed <- sign_E_D2_5thD1_D2down[,c(18,8,9)]
colnames(sign_E_D2_5thD1_D2down_bed) <- c('chrom','start','end')
sign_E_D2_5thD1_D2down_bed$start <- as.numeric(sign_E_D2_5thD1_D2down_bed$start)
sign_E_D2_5thD1_D2down_bed$end <- as.numeric(sign_E_D2_5thD1_D2down_bed$end)

sign_M_D2_5thD1_D2down_bed <- sign_M_D2_5thD1_D2down_write
colnames(sign_M_D2_5thD1_D2down_bed) <- c('chrom','start','end')
sign_M_D2_5thD1_D2down_bed$start <- as.numeric(sign_M_D2_5thD1_D2down_bed$start)
sign_M_D2_5thD1_D2down_bed$end <- as.numeric(sign_M_D2_5thD1_D2down_bed$end)

merged_shared_D2down_bed <- merged_shared_D2down[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_D2down_bed) <- c('chrom','start','end')
merged_shared_D2down_bed$start <- as.numeric(merged_shared_D2down_bed$start)
merged_shared_D2down_bed$end <- as.numeric(merged_shared_D2down_bed$end)

sign_E_D2_5thD1_D2down_bed_match <- bed_intersect(sign_E_D2_5thD1_D2down_bed,conserved_E)
sign_M_D2_5thD1_D2down_bed_match <- bed_intersect(sign_M_D2_5thD1_D2down_bed,conserved_M)
merged_shared_D2down_bed_match <- bed_intersect(merged_shared_D2down_bed,conserved_S)


sign_E_D2_5thD1_D2down_cons <- merge(sign_E_D2_5thD1_D2down, sign_E_D2_5thD1_D2down_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D2_5thD1_D2down_cons <- merge(sign_M_D2_5thD1_D2down, sign_M_D2_5thD1_D2down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_D2down_cons <- merge(merged_ME_D2down, sign_E_D2_5thD1_D2down_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_D2down_cons <- merge(merged_EM_D2down, sign_M_D2_5thD1_D2down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(sign_E_D2_5thD1_D2down_cons)-nrow(merged_ME_D2down_cons[grep("NA", merged_ME_D2down_cons$Row.names), ])
nrow(sign_M_D2_5thD1_D2down_cons)-nrow(merged_EM_D2down_cons[grep("NA", merged_EM_D2down_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D2_5thD1_D2down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D2_5thD1_D2down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_D2down_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D2_5thD1_D2down_bed_match)
nrow(sign_M_D2_5thD1_D2down_bed_match)
nrow(merged_shared_D2down_bed_match)

nrow(sign_E_D2_5thD1_D2down_cons)-nrow(merged_ME_D2down_cons[grep("NA", merged_ME_D2down_cons$Row.names), ])
nrow(sign_M_D2_5thD1_D2down_cons)-nrow(merged_EM_D2down_cons[grep("NA", merged_EM_D2down_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

#########################