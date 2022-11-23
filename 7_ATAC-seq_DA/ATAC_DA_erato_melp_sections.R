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
match_EM <- read.table("MACS2/match_melp_to_erato.0.5.txt", h=F)
match_ME <- read.table("MACS2/match_erato_to_melp.0.5.txt", h=F)

colnames(match_EM) <- c('genome1', 'pan_start_E', 'pan_end_E', 'genome2', 'pan_start_M', 'pan_end_M', 'length')
colnames(match_ME) <- c('genome1', 'pan_start_M', 'pan_end_M', 'genome2', 'pan_start_E', 'pan_end_E', 'length')

match_ME_A <- subset(match_ME, match_ME$genome2 != '.')
match_EM_A <- subset(match_EM, match_EM$genome2 != '.')

nrow(counts_erato_pan)
nrow(counts_melp_pan)

nrow(match_ME)
nrow(match_EM)

nrow(match_ME_A)
nrow(match_EM_A)

# merge intersect files with peak counts
counts_erato_pan_match <- merge(counts_erato_pan, match_EM, by.x = c('pan_start','pan_end'), by.y =c('pan_start_E', 'pan_end_E'), all.x = TRUE, all.y = TRUE)
# counts_melp_pan_match <- merge(counts_melp_pan, match_EM, by.x = c('pan_start','pan_end'), by.y =c('pan_start_M', 'pan_end_M'))
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

nrow(counts_erato[is.na(counts_erato$BR10_Demophoon_Brain),])

ctsF[is.na(ctsF)] <- 0

barplot(colSums(ctsF2), las=2)

###
# Gradient analysis
###

# D1 FW erato
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))
sampleInfo_subE$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] == 'FP'){sampleInfo_subE$tissueNumeric[e] <- 1}
  if(sampleInfo_subE$tissue[e] == 'FM'){sampleInfo_subE$tissueNumeric[e] <- 2}
  if(sampleInfo_subE$tissue[e] == 'FD'){sampleInfo_subE$tissueNumeric[e] <- 3}
}

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissueNumeric)

atacDDS_E_grad <- DESeq(ddsE)
res_E_grad <- results(atacDDS_E_grad)


save(res_E_grad, file='res_E_grad.rda')


sign_E_D1_grad_pos <- subset(res_E_grad, res_E_grad$padj < 0.05 & (res_E_grad$log2FoldChange) > 1)
sign_E_D1_grad_neg <- subset(res_E_grad, res_E_grad$padj < 0.05 & (res_E_grad$log2FoldChange) < -1)

sign_E_D1_grad_pos_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_grad_pos)), col = 'rownames(sign_E_D1_grad_pos)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_grad_neg_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_grad_neg)), col = 'rownames(sign_E_D1_grad_neg)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(sign_E_D1_grad_pos_splits[,c("genome","start","end")], file = "D1_grad_pos_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_E_D1_grad_neg_splits[,c("genome","start","end")], file = "D1_grad_neg_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')

# D2 FW erato
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))
sampleInfo_subE$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] == 'FP'){sampleInfo_subE$tissueNumeric[e] <- 1}
  if(sampleInfo_subE$tissue[e] == 'FM'){sampleInfo_subE$tissueNumeric[e] <- 2}
  if(sampleInfo_subE$tissue[e] == 'FD'){sampleInfo_subE$tissueNumeric[e] <- 3}
}

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissueNumeric)

atacDDS_E_grad_D2 <- DESeq(ddsE)
res_E_grad_D2 <- results(atacDDS_E_grad_D2)

sign_E_D2_grad_pos <- subset(res_E_grad_D2, res_E_grad_D2$padj < 0.05 & (res_E_grad_D2$log2FoldChange) > 1)
sign_E_D2_grad_neg <- subset(res_E_grad_D2, res_E_grad_D2$padj < 0.05 & (res_E_grad_D2$log2FoldChange) < -1)

sign_E_D2_grad_pos_splits <- separate(data =  as.data.frame(rownames(sign_E_D2_grad_pos)), col = 'rownames(sign_E_D2_grad_pos)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D2_grad_neg_splits <- separate(data =  as.data.frame(rownames(sign_E_D2_grad_neg)), col = 'rownames(sign_E_D2_grad_neg)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(sign_E_D2_grad_pos_splits[,c("genome","start","end")], file = "D2_grad_pos_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_E_D2_grad_neg_splits[,c("genome","start","end")], file = "D2_grad_neg_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')


sign_E_D1_grad_pos$erato_gene <- rownames(sign_E_D1_grad_pos)
sign_E_D1_grad_neg$erato_gene <- rownames(sign_E_D1_grad_neg)
sign_E_D2_grad_pos$erato_gene <- rownames(sign_E_D2_grad_pos)
sign_E_D2_grad_neg$erato_gene <- rownames(sign_E_D2_grad_neg)

merge(as.data.frame(sign_E_D1_grad_pos), as.data.frame(sign_E_D2_grad_pos), by = 'erato_gene')
merge(as.data.frame(sign_E_D1_grad_neg), as.data.frame(sign_E_D2_grad_neg), by = 'erato_gene')
# The results table will then contain the regression coeffcient in the table log2FoldChange. That the column is called "log2FoldChange" is then, of course, a bit a misnomer; it is the slope of the regression line, i.e., change in expression divided by change in the numeric predictor. (The expression is measured in units of log2 normalized counts, the predictor in whatever unit you use; in your case log10 of the dose.).

# D1 FW melpomene
sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))
sampleInfo_subM$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] == 'FP'){sampleInfo_subM$tissueNumeric[e] <- 1}
  if(sampleInfo_subM$tissue[e] == 'FM'){sampleInfo_subM$tissueNumeric[e] <- 2}
  if(sampleInfo_subM$tissue[e] == 'FD'){sampleInfo_subM$tissueNumeric[e] <- 3}
}

ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissueNumeric)

atacDDS_M_grad <- DESeq(ddsM)
res_M_grad <- results(atacDDS_M_grad)

save(res_M_grad, file='res_M_grad.rda')


sign_M_D1_grad_pos <- subset(res_M_grad, res_M_grad$padj < 0.05 & (res_M_grad$log2FoldChange) > 1)
sign_M_D1_grad_neg <- subset(res_M_grad, res_M_grad$padj < 0.05 & (res_M_grad$log2FoldChange) < -1)

sign_M_D1_grad_pos_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_grad_pos)), col = 'rownames(sign_M_D1_grad_pos)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_grad_neg_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_grad_neg)), col = 'rownames(sign_M_D1_grad_neg)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(sign_M_D1_grad_pos_splits[,c("genome","start","end")], file = "D1_grad_pos_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_M_D1_grad_neg_splits[,c("genome","start","end")], file = "D1_grad_neg_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')


####

# extract normalized counts of sign peaks


# erato D1 Pos
test_counts <- counts(atacDDS_E_grad[rownames(sign_E_D1_grad_pos),], normalized = TRUE)
# test_counts <- test_counts[!rownames(test_counts)%in%rownames(sign_E_D1_FP_rest_down),]
# test_counts <- test_counts[!rownames(test_counts)%in%rownames(sign_E_D1_FD_rest_up),]
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_D1_pos_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_D1_pos_era <- aggregate(test_counts_m_D1_pos_era, list(test_counts_m_D1_pos_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_D1_pos_era <- scale(test_counts_m_m_D1_pos_era[,c(1:(ncol(test_counts_m_m_D1_pos_era)-1))])
test_counts_m_m_s_D1_pos_era <- cbind(as.data.frame(test_counts_m_m_s_D1_pos_era), test_counts_m_m_D1_pos_era$tissueNumeric)


# erato D1 Neg
test_counts <- counts(atacDDS_E_grad[rownames(sign_E_D1_grad_neg),], normalized = TRUE)
erato_grad_D1_neg_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_D1_neg_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_D1_neg_era <- aggregate(test_counts_m_D1_neg_era, list(test_counts_m_D1_neg_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_D1_neg_era <- scale(test_counts_m_m_D1_neg_era[,c(1:(ncol(test_counts_m_m_D1_neg_era)-1))])
test_counts_m_m_s_D1_neg_era <- cbind(as.data.frame(test_counts_m_m_s_D1_neg_era), test_counts_m_m_D1_neg_era$tissueNumeric)


# melp D1 Pos
test_counts <- counts(atacDDS_M_grad[rownames(sign_M_D1_grad_pos),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_D1_pos_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_D1_pos_melp <- aggregate(test_counts_m_D1_pos_melp, list(test_counts_m_D1_pos_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_D1_pos_melp <- scale(test_counts_m_m_D1_pos_melp[,c(1:(ncol(test_counts_m_m_D1_pos_melp)-1))])
test_counts_m_m_s_D1_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_D1_pos_melp), test_counts_m_m_D1_pos_melp$tissueNumeric)


# melp D1 Neg
test_counts <- counts(atacDDS_M_grad[rownames(sign_M_D1_grad_neg),], normalized = TRUE)
melp_grad_D1_neg_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_D1_neg_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_D1_neg_melp <- aggregate(test_counts_m_D1_neg_melp, list(test_counts_m_D1_neg_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_D1_neg_melp <- scale(test_counts_m_m_D1_neg_melp[,c(1:(ncol(test_counts_m_m_D1_neg_melp)-1))])
test_counts_m_m_s_D1_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_D1_neg_melp), test_counts_m_m_D1_neg_melp$tissueNumeric)

## shared

test_counts_m_m_s_D1_pos_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_D1_pos_era)[c(3:ncol(test_counts_m_m_s_D1_pos_era)-1)]  ), col = 'colnames(test_counts_m_m_s_D1_pos_era)[c(3:ncol(test_counts_m_m_s_D1_pos_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_D1_neg_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_D1_neg_era)[c(3:ncol(test_counts_m_m_s_D1_neg_era)-1)]  ), col = 'colnames(test_counts_m_m_s_D1_neg_era)[c(3:ncol(test_counts_m_m_s_D1_neg_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_D1_pos_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_D1_pos_melp)[c(3:ncol(test_counts_m_m_s_D1_pos_melp)-1)]), col = 'colnames(test_counts_m_m_s_D1_pos_melp)[c(3:ncol(test_counts_m_m_s_D1_pos_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_D1_neg_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_D1_neg_melp)[c(3:ncol(test_counts_m_m_s_D1_neg_melp)-1)]), col = 'colnames(test_counts_m_m_s_D1_neg_melp)[c(3:ncol(test_counts_m_m_s_D1_neg_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

test_counts_m_m_s_D1_pos_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_D1_pos_era)[c(3:ncol(test_counts_m_m_s_D1_pos_era)-1)])  ,test_counts_m_m_s_D1_pos_era_splits )
test_counts_m_m_s_D1_neg_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_D1_neg_era)[c(3:ncol(test_counts_m_m_s_D1_neg_era)-1)])  ,test_counts_m_m_s_D1_neg_era_splits )
test_counts_m_m_s_D1_pos_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_D1_pos_melp)[c(3:ncol(test_counts_m_m_s_D1_pos_melp)-1)]), test_counts_m_m_s_D1_pos_melp_splits)
test_counts_m_m_s_D1_neg_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_D1_neg_melp)[c(3:ncol(test_counts_m_m_s_D1_neg_melp)-1)]), test_counts_m_m_s_D1_neg_melp_splits)

test_counts_m_m_s_D1_pos_era_split  <- subset(test_counts_m_m_s_D1_pos_era_split , test_counts_m_m_s_D1_pos_era_split$scafM  != "NA")
test_counts_m_m_s_D1_neg_era_split  <- subset(test_counts_m_m_s_D1_neg_era_split , test_counts_m_m_s_D1_neg_era_split$scafM  != "NA")
test_counts_m_m_s_D1_pos_melp_split <- subset(test_counts_m_m_s_D1_pos_melp_split, test_counts_m_m_s_D1_pos_melp_split$scafE != "NA")
test_counts_m_m_s_D1_neg_melp_split <- subset(test_counts_m_m_s_D1_neg_melp_split, test_counts_m_m_s_D1_neg_melp_split$scafE != "NA")

test_counts_m_m_s_D1_pos_era_sharedn <- merge(test_counts_m_m_s_D1_pos_era_split, test_counts_m_m_s_D1_pos_melp_split, by = c('start','end'))[,3]
test_counts_m_m_s_D1_neg_era_sharedn <- merge(test_counts_m_m_s_D1_neg_era_split, test_counts_m_m_s_D1_neg_melp_split, by = c('start','end'))[,3]

test_counts_m_m_s_D1_pos_era_shared <-  as.data.frame(test_counts_m_m_s_D1_pos_era[, as.character(test_counts_m_m_s_D1_pos_era_sharedn)])
test_counts_m_m_s_D1_neg_era_shared <-  as.data.frame(test_counts_m_m_s_D1_neg_era[, as.character(test_counts_m_m_s_D1_neg_era_sharedn)])
test_counts_m_m_s_D1_pos_melp_shared <- as.data.frame(test_counts_m_m_s_D1_pos_melp[, as.character(test_counts_m_m_s_D1_pos_era_sharedn)])
test_counts_m_m_s_D1_neg_melp_shared <- as.data.frame(test_counts_m_m_s_D1_neg_melp[, as.character(test_counts_m_m_s_D1_neg_era_sharedn)])

write.table(test_counts_m_m_s_D1_pos_era_splits[,c("scafE","startE","endE")], file = "D1_grad_pos_era_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_D1_neg_era_splits[,c("scafE","startE","endE")], file = "D1_grad_neg_era_06122021.bed", col.names = F, row.names = F, quote = F)

write.table(test_counts_m_m_s_D1_pos_melp_splits[,c("scafM","startM","endM")], file = "D1_grad_pos_melp_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_D1_neg_melp_splits[,c("scafM","startM","endM")], file = "D1_grad_neg_melp_06122021.bed", col.names = F, row.names = F, quote = F)

test_counts_m_m_s_D1_pos_era_splits_split <- separate(data =  as.data.frame((test_counts_m_m_s_D1_pos_era_sharedn)), col = '(test_counts_m_m_s_D1_pos_era_sharedn)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_D1_neg_era_splits_split <- separate(data =  as.data.frame((test_counts_m_m_s_D1_neg_era_sharedn)), col = '(test_counts_m_m_s_D1_neg_era_sharedn)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(test_counts_m_m_s_D1_pos_era_splits_split[,c("scafE","startE","endE")], file = "D1_grad_pos_era_shared_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_D1_neg_era_splits_split[,c("scafE","startE","endE")], file = "D1_grad_neg_era_shared_06122021.bed", col.names = F, row.names = F, quote = F)


##
par(mfcol=c(3,2))

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_pos_era)-1)){
  lines(test_counts_m_m_s_D1_pos_era[,ncol(test_counts_m_m_s_D1_pos_era)],  test_counts_m_m_s_D1_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_pos_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_D1_pos_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_neg_era)-1)){
  lines(test_counts_m_m_s_D1_neg_era[,ncol(test_counts_m_m_s_D1_neg_era)],  test_counts_m_m_s_D1_neg_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_neg_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_D1_neg_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,type='n',axes=FALSE,ann=FALSE)


plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_pos_melp)-1)){
  lines(test_counts_m_m_s_D1_pos_melp[,ncol(test_counts_m_m_s_D1_pos_melp)],  test_counts_m_m_s_D1_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_pos_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_D1_pos_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_neg_melp)-1)){
  lines(test_counts_m_m_s_D1_neg_melp[,ncol(test_counts_m_m_s_D1_neg_melp)],  test_counts_m_m_s_D1_neg_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_neg_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_D1_neg_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,type='n',axes=FALSE,ann=FALSE)
####







sign_list_E_dir <- list(sign_E_D1_grad_pos, sign_E_D1_grad_neg)
sign_list_M_dir <- list(sign_M_D1_grad_pos, sign_M_D1_grad_neg)


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
  sign_list_M_df_dir[[i]] <- separate(data = sign_list_M_df_dir[[i]], col = 'rownames(sign_list_M_dir[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_M_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_M_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_M_df_all_dir[[i]] <- cbind(sign_list_M_dir[[i]], sign_list_M_df_dir[[i]])
}

sections_D1_grad_pos_shared_E <- merge(sign_list_E_df_all_dir[[1]],sign_list_M_df_all_dir[[1]], by=c('start','end'))
sections_D1_grad_neg_shared_E <- merge(sign_list_E_df_all_dir[[2]],sign_list_M_df_all_dir[[2]], by=c('start','end'))

write.table(sign_list_E_df_all_dir[[1]][,c("scafE","startE","endE")], file = 'sections_erato_grad_up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[2]][,c("scafE","startE","endE")], file = 'sections_erato_grad_down.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("scafM","startM","endM")], file = 'sections_melp_grad_up.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[2]][,c("scafM","startM","endM")], file = 'sections_melp_grad_down.bed', col.names = F, row.names = F, quote = F)

write.table(sections_D1_grad_pos_shared_E[,c("scafE.x","startE.x","endE.x")], file = 'sections_shared_erato_grad_up.bed', col.names = F, row.names = F, quote = F)
write.table(sections_D1_grad_neg_shared_E[,c("scafE.x","startE.x","endE.x")], file = 'sections_shared_erato_grad_down.bed', col.names = F, row.names = F, quote = F)






######################
# D1 FM vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FM'){
     sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E_FM_rest <- DESeq(ddsE)

res_E_FM_rest <- results(atacDDS_E_FM_rest)

save(res_E_FM_rest, file='res_E_FM_rest.rda')

sign_E_D1_FM_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D1_FM_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D1_FM_rest_down, file='sign_E_D1_FM_rest_down.rda')
save(sign_E_D1_FM_rest_up, file='sign_E_D1_FM_rest_up.rda')

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FM'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M_FM_rest <- DESeq(ddsM)

res_M_FM_rest <- results(atacDDS_M_FM_rest)

save(res_M_FM_rest, file='res_M_FM_rest.rda')

sign_M_D1_FM_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D1_FM_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D1_FM_rest_down, file='sign_M_D1_FM_rest_down.rda')
save(sign_M_D1_FM_rest_up, file='sign_M_D1_FM_rest_up.rda')

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$tissue == 'FM' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2A_FM <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))

ddsA_FM <- DESeqDataSetFromMatrix(countData = ctsF2A_FM,
                               colData = sampleInfo_subA,
                               design = ~species)

atacDDS_A_FM <- DESeq(ddsA_FM)
res_A_D1_FM <- results(atacDDS_A_FM)
sign_A_D1_FM <- subset(res_A_D1_FM, res_A_D1_FM$padj < 0.05 & abs(res_A_D1_FM$log2FoldChange > 1))

nrow(sign_A_D1_FM)
nrow(merge(as.data.frame(sign_A_D1_FM), as.data.frame(sign_E_D1_FM_rest_up), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FM), as.data.frame(sign_M_D1_FM_rest_up), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_D1_FM), as.data.frame(sign_E_D1_FM_rest_down), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FM), as.data.frame(sign_M_D1_FM_rest_down), by = 'row.names'))

######################
# D1 FP vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FP'){
    sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E_FP_rest <- DESeq(ddsE)

res_E_FP_rest <- results(atacDDS_E_FP_rest)

save(res_E_FP_rest, file='res_E_FP_rest.rda')

sign_E_D1_FP_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D1_FP_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D1_FP_rest_down, file='sign_E_D1_FP_rest_down.rda')
save(sign_E_D1_FP_rest_up, file='sign_E_D1_FP_rest_up.rda')


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FP'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M_FP_rest <- DESeq(ddsM)

res_M_FP_rest <- results(atacDDS_M_FP_rest)

save(res_M_FP_rest, file='res_M_FP_rest.rda')

sign_M_D1_FP_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D1_FP_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D1_FP_rest_down, file='sign_M_D1_FP_rest_down.rda')
save(sign_M_D1_FP_rest_up, file='sign_M_D1_FP_rest_up.rda')

sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$tissue == 'FP' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2A_FP <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))

ddsA_FP <- DESeqDataSetFromMatrix(countData = ctsF2A_FP,
                                  colData = sampleInfo_subA,
                                  design = ~species)

atacDDS_A_FP <- DESeq(ddsA_FP)
res_A_D1_FP <- results(atacDDS_A_FP)
sign_A_D1_FP <- subset(res_A_D1_FP, res_A_D1_FP$padj < 0.05 & abs(res_A_D1_FP$log2FoldChange > 1))

nrow(sign_A_D1_FP)
nrow(merge(as.data.frame(sign_A_D1_FP), as.data.frame(sign_E_D1_FP_rest_up), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FP), as.data.frame(sign_M_D1_FP_rest_up), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_D1_FP), as.data.frame(sign_E_D1_FP_rest_down), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FP), as.data.frame(sign_M_D1_FP_rest_down), by = 'row.names'))

######################
# D1 FD vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FD'){
    sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E_FD_rest <- DESeq(ddsE)

res_E_FD_rest <- results(atacDDS_E_FD_rest)

save(res_E_FD_rest, file='res_E_FD_rest.rda')

sign_E_D1_FD_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D1_FD_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D1_FD_rest_down, file='sign_E_D1_FD_rest_down.rda')
save(sign_E_D1_FD_rest_up, file='sign_E_D1_FD_rest_up.rda')


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FD'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M_FD_rest <- DESeq(ddsM)

res_M_FD_rest <- results(atacDDS_M_FD_rest)

save(res_M_FD_rest, file='res_M_FD_rest.rda')

sign_M_D1_FD_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D1_FD_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D1_FD_rest_down, file='sign_M_D1_FD_rest_down.rda')
save(sign_M_D1_FD_rest_up, file='sign_M_D1_FD_rest_up.rda')


sampleInfo_subA <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$tissue == 'FD' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2A_FD <- select(ctsF2A, c(as.character(sampleInfo_subA$id)))

ddsA_FD <- DESeqDataSetFromMatrix(countData = ctsF2A_FD,
                                  colData = sampleInfo_subA,
                                  design = ~species)

atacDDS_A_FD <- DESeq(ddsA_FD)
res_A_D1_FD <- results(atacDDS_A_FD)
sign_A_D1_FD <- subset(res_A_D1_FD, res_A_D1_FD$padj < 0.05 & abs(res_A_D1_FD$log2FoldChange > 1))

nrow(sign_A_D1_FD)
nrow(merge(as.data.frame(sign_A_D1_FD), as.data.frame(sign_E_D1_FD_rest_up), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FD), as.data.frame(sign_M_D1_FD_rest_up), by = 'row.names'))

nrow(merge(as.data.frame(sign_A_D1_FD), as.data.frame(sign_E_D1_FD_rest_down), by = 'row.names'))
nrow(merge(as.data.frame(sign_A_D1_FD), as.data.frame(sign_M_D1_FD_rest_down), by = 'row.names'))
######################
# D2 FM vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FM'){
    sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_D2_FM_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D2_FM_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D2_FM_rest_down, file='sign_E_D2_FM_rest_down.rda')
save(sign_E_D2_FM_rest_up, file='sign_E_D2_FM_rest_up.rda')

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FM'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M <- DESeq(ddsM)

res_M <- results(atacDDS_M)

sign_M_D2_FM_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D2_FM_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D2_FM_rest_down, file='sign_M_D2_FM_rest_down.rda')
save(sign_M_D2_FM_rest_up, file='sign_M_D2_FM_rest_up.rda')

######################
# D2 FP vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FP'){
    sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_D2_FP_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D2_FP_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D2_FP_rest_down, file='sign_E_D2_FP_rest_down.rda')
save(sign_E_D2_FP_rest_up, file='sign_E_D2_FP_rest_up.rda')


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FP'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M <- DESeq(ddsM)

res_M <- results(atacDDS_M)

sign_M_D2_FP_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D2_FP_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D2_FP_rest_down, file='sign_M_D2_FP_rest_down.rda')
save(sign_M_D2_FP_rest_up, file='sign_M_D2_FP_rest_up.rda')

######################
# D2 FD vs rest
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] != 'FD'){
    sampleInfo_subE$tissue[e] <- 'rest'
  }
}
ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~tissue)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_D2_FD_rest_down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_E_D2_FD_rest_up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)

save(sign_E_D2_FD_rest_down, file='sign_E_D2_FD_rest_down.rda')
save(sign_E_D2_FD_rest_up, file='sign_E_D2_FD_rest_up.rda')


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] != 'FD'){
    sampleInfo_subM$tissue[e] <- 'rest'
  }
}
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~tissue)
atacDDS_M <- DESeq(ddsM)

res_M <- results(atacDDS_M)

sign_M_D2_FD_rest_down <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)
sign_M_D2_FD_rest_up <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_M_D2_FD_rest_down, file='sign_M_D2_FD_rest_down.rda')
save(sign_M_D2_FD_rest_up, file='sign_M_D2_FD_rest_up.rda')

######################
# D1 vs D2 erato FM
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue == 'FM') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~stage)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_FM_D1D2_D1up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_E_FM_D1D2_D1down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)

save(sign_E_FM_D1D2_D1up, file='sign_E_FM_D1D2_D1up.rda')
save(sign_E_FM_D1D2_D1down, file='sign_E_FM_D1D2_D1down.rda')

######################
# D1 vs D2 erato FP
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue == 'FP') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~stage)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_FP_D1D2_D1up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_E_FP_D1D2_D1down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)

save(sign_E_FP_D1D2_D1up, file='sign_E_FP_D1D2_D1up.rda')
save(sign_E_FP_D1D2_D1down, file='sign_E_FP_D1D2_D1down.rda')

######################
# D1 vs D2 erato FD
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$tissue == 'FD') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~stage)
atacDDS_E <- DESeq(ddsE)

res_E <- results(atacDDS_E)

sign_E_FD_D1D2_D1up <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_E_FD_D1D2_D1down <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)

save(sign_E_FD_D1D2_D1up, file='sign_E_FD_D1D2_D1up.rda')
save(sign_E_FD_D1D2_D1down, file='sign_E_FD_D1D2_D1down.rda')

######################
# D1 FP vs FM  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FD'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FD'))



ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D1_FP_FM_FPup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_FP_FM_FPup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D1_FP_FM_FPdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_FP_FM_FPdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D1_FP_FM_FPup, file='sign_E_D1_FP_FM_p05fc1_FPup.rda')
save(sign_E_D1_FP_FM_FPdown, file='sign_E_D1_FP_FM_p05fc1_FPdown.rda')
save(sign_M_D1_FP_FM_FPup, file='sign_M_D1_FP_FM_p05fc1_FPup.rda')
save(sign_M_D1_FP_FM_FPdown, file='sign_M_D1_FP_FM_p05fc1_FPdown.rda')


######################
# D1 FP vs FD  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FM'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FM'))



ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D1_FP_FD_FPup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_FP_FD_FPup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D1_FP_FD_FPdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_FP_FD_FPdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D1_FP_FD_FPup, file='sign_E_D1_FP_FD_p05fc1_FPup.rda')
save(sign_E_D1_FP_FD_FPdown, file='sign_E_D1_FP_FD_p05fc1_FPdown.rda')
save(sign_M_D1_FP_FD_FPup, file='sign_M_D1_FP_FD_p05fc1_FPup.rda')
save(sign_M_D1_FP_FD_FPdown, file='sign_M_D1_FP_FD_p05fc1_FPdown.rda')


######################
# D1 FM vs FD  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FP'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FP'))



ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D1_FM_FD_FMup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_FM_FD_FMup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D1_FM_FD_FMdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_FM_FD_FMdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D1_FM_FD_FMup, file='sign_E_D1_FM_FD_p05fc1_FMup.rda')
save(sign_E_D1_FM_FD_FMdown, file='sign_E_D1_FM_FD_p05fc1_FMdown.rda')
save(sign_M_D1_FM_FD_FMup, file='sign_M_D1_FM_FD_p05fc1_FMup.rda')
save(sign_M_D1_FM_FD_FMdown, file='sign_M_D1_FM_FD_p05fc1_FMdown.rda')










######################
# D2 FP vs FM  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FD' ))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FD'))



ctsF2E <- select(ctsF2, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF2, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FP_FM_FPup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FP_FM_FPup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FP_FM_FPdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FP_FM_FPdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D2_FP_FM_FPup, file='sign_E_D2_FP_FM_p05fc1_FPup.rda')
save(sign_E_D2_FP_FM_FPdown, file='sign_E_D2_FP_FM_p05fc1_FPdown.rda')
save(sign_M_D2_FP_FM_FPup, file='sign_M_D2_FP_FM_p05fc1_FPup.rda')
save(sign_M_D2_FP_FM_FPdown, file='sign_M_D2_FP_FM_p05fc1_FPdown.rda')

atac_Rlog <- vst(atacDDS_M)
plotPCA(atac_Rlog, intgroup = c("tissue"), ntop = nrow(atac_Rlog))

######################
# D2 FP vs FD  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FM'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FM'))



ctsF2E <- select(ctsF2, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF2, c(as.character(sampleInfo_subM$id)))

barplot(colSums(ctsF2M), las=2)

##
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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FP_FD_FPup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FP_FD_FPup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FP_FD_FPdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FP_FD_FPdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D2_FP_FD_FPup, file='sign_E_D2_FP_FD_p05fc1_FPup.rda')
save(sign_E_D2_FP_FD_FPdown, file='sign_E_D2_FP_FD_p05fc1_FPdown.rda')
save(sign_M_D2_FP_FD_FPup, file='sign_M_D2_FP_FD_p05fc1_FPup.rda')
save(sign_M_D2_FP_FD_FPdown, file='sign_M_D2_FP_FD_p05fc1_FPdown.rda')


######################
# D2 FM vs FD  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FP'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP' & sampleInfo_sub$tissue != 'FP'))



ctsF2E <- select(ctsF2, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF2, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FM_FD_FMup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FM_FD_FMup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FM_FD_FMdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FM_FD_FMdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D2_FM_FD_FMup, file='sign_E_D2_FM_FD_p05fc1_FMup.rda')
save(sign_E_D2_FM_FD_FMdown, file='sign_E_D2_FM_FD_p05fc1_FMdown.rda')
save(sign_M_D2_FM_FD_FMup, file='sign_M_D2_FM_FD_p05fc1_FMup.rda')
save(sign_M_D2_FM_FD_FMdown, file='sign_M_D2_FM_FD_p05fc1_FMdown.rda')





######################
# D1 FW vs HW 
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' ))

sampleInfo_subE$wing <- 'forewing'
for(x in 1:nrow(sampleInfo_subE)){ 
  if(sampleInfo_subE$tissue[[x]] == 'HA' | sampleInfo_subE$tissue[[x]] == 'HP') {sampleInfo_subE$wing[[x]] <- 'hindwing'}}

sampleInfo_subM$wing <- 'forewing'
for(x in 1:nrow(sampleInfo_subM)){ 
  if(sampleInfo_subM$tissue[[x]] == 'HA' | sampleInfo_subM$tissue[[x]] == 'HP') {sampleInfo_subM$wing[[x]] <- 'hindwing'}}

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

###
# setup model
###

dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                              colData = sampleInfo_sub,
                              design = ~wing)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~wing)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~wing)


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D1_FWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D1_FWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

sign_E_D1_FWdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D1_FWdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

save(sign_E_D1_FWup, file='sign_E_D1_p05fc1_FWup.rda')
save(sign_E_D1_FWdown, file='sign_E_D1_p05fc1_FWdown.rda')
save(sign_M_D1_FWup, file='sign_M_D1_p05fc1_FWup.rda')
save(sign_M_D1_FWdown, file='sign_M_D1_p05fc1_FWdown.rda')


######################
# D2 FW vs HW 
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' ))

sampleInfo_subE$wing <- 'forewing'
for(x in 1:nrow(sampleInfo_subE)){ 
  if(sampleInfo_subE$tissue[[x]] == 'HA' | sampleInfo_subE$tissue[[x]] == 'HP') {sampleInfo_subE$wing[[x]] <- 'hindwing'}}

sampleInfo_subM$wing <- 'forewing'
for(x in 1:nrow(sampleInfo_subM)){ 
  if(sampleInfo_subM$tissue[[x]] == 'HA' | sampleInfo_subM$tissue[[x]] == 'HP') {sampleInfo_subM$wing[[x]] <- 'hindwing'}}

ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

###
# setup model
###

dds <- DESeqDataSetFromMatrix(countData = ctsF2,
                              colData = sampleInfo_sub,
                              design = ~wing)

ddsE <- DESeqDataSetFromMatrix(countData = ctsF2E,
                               colData = sampleInfo_subE,
                               design = ~wing)

ddsM <- DESeqDataSetFromMatrix(countData = ctsF2M,
                               colData = sampleInfo_subM,
                               design = ~wing)


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FWup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FWup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

sign_E_D2_FWdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FWdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

save(sign_E_D2_FWup, file='sign_E_D2_p05fc1_FWup.rda')
save(sign_E_D2_FWdown, file='sign_E_D2_p05fc1_FWdown.rda')
save(sign_M_D2_FWup, file='sign_M_D2_p05fc1_FWup.rda')
save(sign_M_D2_FWdown, file='sign_M_D2_p05fc1_FWdown.rda')





############################# Tissue vs all
######################
# D2 FP  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subE)){ if(sampleInfo_subE$tissue[e] != 'FP'){ sampleInfo_subE$tissue[e] <- 'other'}}
for(e in 1:nrow(sampleInfo_subM)){ if(sampleInfo_subM$tissue[e] != 'FP'){ sampleInfo_subM$tissue[e] <- 'other'}}

ctsF2E <- select(ctsF2, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF2, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FPup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FPup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FPdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FPdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D2_FPup, file='sign_E_D2_p05fc1_FPup.rda')
save(sign_E_D2_FPdown, file='sign_E_D2_p05fc1_FPdown.rda')
save(sign_M_D2_FPup, file='sign_M_D2_p05fc1_FPup.rda')
save(sign_M_D2_FPdown, file='sign_M_D2_p05fc1_FPdown.rda')

######################
# D2 FP  
sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY2') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))

sampleInfo_subE$tissue <- as.character(sampleInfo_subE$tissue)
sampleInfo_subM$tissue <- as.character(sampleInfo_subM$tissue)
for(e in 1:nrow(sampleInfo_subE)){ if(sampleInfo_subE$tissue[e] != 'FM'){ sampleInfo_subE$tissue[e] <- 'other'}}
for(e in 1:nrow(sampleInfo_subM)){ if(sampleInfo_subM$tissue[e] != 'FM'){ sampleInfo_subM$tissue[e] <- 'other'}}

ctsF2E <- select(ctsF2, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF2, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)
res_E <- results(atacDDS_E)
res_M <- results(atacDDS_M)

sign_E_D2_FMup <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) > 1)
sign_M_D2_FMup <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) > 1)

sign_E_D2_FMdown <- subset(res_E, res_E$padj < 0.05 & (res_E$log2FoldChange) < -1)
sign_M_D2_FMdown <- subset(res_M, res_M$padj < 0.05 & (res_M$log2FoldChange) < -1)

save(sign_E_D2_FMup, file='sign_E_D2_p05fc1_FMup.rda')
save(sign_E_D2_FMdown, file='sign_E_D2_p05fc1_FMdown.rda')
save(sign_M_D2_FMup, file='sign_M_D2_p05fc1_FMup.rda')
save(sign_M_D2_FMdown, file='sign_M_D2_p05fc1_FMdown.rda')

#############################

load('sign_E_D1_FP_FM_p05fc1_FPup.rda')
load('sign_E_D1_FP_FM_p05fc1_FPdown.rda')
load('sign_M_D1_FP_FM_p05fc1_FPup.rda')
load('sign_M_D1_FP_FM_p05fc1_FPdown.rda')

load('sign_E_D1_FP_FD_p05fc1_FPup.rda')
load('sign_E_D1_FP_FD_p05fc1_FPdown.rda')
load('sign_M_D1_FP_FD_p05fc1_FPup.rda')
load('sign_M_D1_FP_FD_p05fc1_FPdown.rda')

load('sign_E_D1_FM_FD_p05fc1_FMup.rda')
load('sign_E_D1_FM_FD_p05fc1_FMdown.rda')
load('sign_M_D1_FM_FD_p05fc1_FMup.rda')
load('sign_M_D1_FM_FD_p05fc1_FMdown.rda')

load('sign_E_D2_FP_FM_p05fc1_FPup.rda')
load('sign_E_D2_FP_FM_p05fc1_FPdown.rda')
load('sign_M_D2_FP_FM_p05fc1_FPup.rda')
load('sign_M_D2_FP_FM_p05fc1_FPdown.rda')

load('sign_E_D2_FP_FD_p05fc1_FPup.rda')
load('sign_E_D2_FP_FD_p05fc1_FPdown.rda')
load('sign_M_D2_FP_FD_p05fc1_FPup.rda')
load('sign_M_D2_FP_FD_p05fc1_FPdown.rda')

load('sign_E_D2_FM_FD_p05fc1_FMup.rda')
load('sign_E_D2_FM_FD_p05fc1_FMdown.rda')
load('sign_M_D2_FM_FD_p05fc1_FMup.rda')
load('sign_M_D2_FM_FD_p05fc1_FMdown.rda')


load('sign_E_D1_FM_rest_down.rda')
load('sign_E_D1_FM_rest_up.rda')
load('sign_M_D1_FM_rest_down.rda')
load('sign_M_D1_FM_rest_up.rda')

load('sign_E_D1_FP_rest_down.rda')
load('sign_E_D1_FP_rest_up.rda')
load('sign_M_D1_FP_rest_down.rda')
load('sign_M_D1_FP_rest_up.rda')

load('sign_E_D1_FD_rest_down.rda')
load('sign_E_D1_FD_rest_up.rda')
load('sign_M_D1_FD_rest_down.rda')
load('sign_M_D1_FD_rest_up.rda')


load('sign_E_D2_FM_rest_down.rda')
load('sign_E_D2_FM_rest_up.rda')
load('sign_M_D2_FM_rest_down.rda')
load('sign_M_D2_FM_rest_up.rda')

load('sign_E_D2_FP_rest_down.rda')
load('sign_E_D2_FP_rest_up.rda')
load('sign_M_D2_FP_rest_down.rda')
load('sign_M_D2_FP_rest_up.rda')

load('sign_E_D2_FD_rest_down.rda')
load('sign_E_D2_FD_rest_up.rda')
load('sign_M_D2_FD_rest_down.rda')
load('sign_M_D2_FD_rest_up.rda')


load('sign_E_FM_D1D2_D1up.rda')
load('sign_E_FM_D1D2_D1down.rda')

load('sign_E_FP_D1D2_D1up.rda')
load('sign_E_FP_D1D2_D1down.rda')

load('sign_E_FD_D1D2_D1up.rda')
load('sign_E_FD_D1D2_D1down.rda')

sign_E_D1_FM_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FM_rest_up)), col = 'rownames(sign_E_D1_FM_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FM_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FM_rest_down)), col = 'rownames(sign_E_D1_FM_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(sign_E_D1_FM_rest_up_splits[,c("genome","start","end")], file = "D1_FMup_pos_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_E_D1_FM_rest_down_splits[,c("genome","start","end")], file = "D1_FMdown_neg_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')

sign_M_D1_FM_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FM_rest_up)), col = 'rownames(sign_M_D1_FM_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FM_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FM_rest_down)), col = 'rownames(sign_M_D1_FM_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(sign_M_D1_FM_rest_up_splits[,c("genome","start","end")], file = "D1_FMup_pos_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_M_D1_FM_rest_down_splits[,c("genome","start","end")], file = "D1_FMdown_neg_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')


sign_E_D1_FP_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FP_rest_up)), col = 'rownames(sign_E_D1_FP_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_E_D1_FP_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FP_rest_down)), col = 'rownames(sign_E_D1_FP_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")

write.table(sign_E_D1_FP_rest_up_splits[,c("genome","start","end")], file = "D1_FPup_pos_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_E_D1_FP_rest_down_splits[,c("genome","start","end")], file = "D1_FPdown_neg_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')

sign_M_D1_FP_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FP_rest_up)), col = 'rownames(sign_M_D1_FP_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_M_D1_FP_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FP_rest_down)), col = 'rownames(sign_M_D1_FP_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")

write.table(sign_M_D1_FP_rest_up_splits[,c("genome","start","end")], file = "D1_FPup_pos_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_M_D1_FP_rest_down_splits[,c("genome","start","end")], file = "D1_FPdown_neg_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')


sign_E_D1_FD_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FD_rest_up)), col = 'rownames(sign_E_D1_FD_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFD", "startM", "endM"), sep = "_")
sign_E_D1_FD_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_E_D1_FD_rest_down)), col = 'rownames(sign_E_D1_FD_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFD", "startM", "endM"), sep = "_")

write.table(sign_E_D1_FD_rest_up_splits[,c("genome","start","end")], file = "D1_FDup_pos_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_E_D1_FD_rest_down_splits[,c("genome","start","end")], file = "D1_FDdown_neg_era_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')

sign_M_D1_FD_rest_up_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FD_rest_up)), col = 'rownames(sign_M_D1_FD_rest_up)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFD", "startM", "endM"), sep = "_")
sign_M_D1_FD_rest_down_splits <- separate(data =  as.data.frame(rownames(sign_M_D1_FD_rest_down)), col = 'rownames(sign_M_D1_FD_rest_down)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scaFD", "startM", "endM"), sep = "_")

write.table(sign_M_D1_FD_rest_up_splits[,c("genome","start","end")], file = "D1_FDup_pos_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sign_M_D1_FD_rest_down_splits[,c("genome","start","end")], file = "D1_FDdown_neg_melp_shared_PAN.bed", col.names = F, row.names = F, quote = F, sep = '\t')

####
# extract normalized counts of sign peaks

sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))
sampleInfo_subE$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_subE)){
  if(sampleInfo_subE$tissue[e] == 'FP'){sampleInfo_subE$tissueNumeric[e] <- 1}
  if(sampleInfo_subE$tissue[e] == 'FM'){sampleInfo_subE$tissueNumeric[e] <- 2}
  if(sampleInfo_subE$tissue[e] == 'FD'){sampleInfo_subE$tissueNumeric[e] <- 3}
}

sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))
sampleInfo_subM$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_subM)){
  if(sampleInfo_subM$tissue[e] == 'FP'){sampleInfo_subM$tissueNumeric[e] <- 1}
  if(sampleInfo_subM$tissue[e] == 'FM'){sampleInfo_subM$tissueNumeric[e] <- 2}
  if(sampleInfo_subM$tissue[e] == 'FD'){sampleInfo_subM$tissueNumeric[e] <- 3}
}



# erato FP up
test_counts <- counts(atacDDS_E_FP_rest[rownames(sign_E_D1_FP_rest_up),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FP_pos_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FP_pos_era <- aggregate(test_counts_m_FP_pos_era, list(test_counts_m_FP_pos_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FP_pos_era <- scale(test_counts_m_m_FP_pos_era[,c(1:(ncol(test_counts_m_m_FP_pos_era)-1))])
test_counts_m_m_s_FP_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FP_pos_era), test_counts_m_m_FP_pos_era$tissueNumeric)


# erato FP down
test_counts <- counts(atacDDS_E_FP_rest[rownames(sign_E_D1_FP_rest_down),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FP_neg_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FP_neg_era <- aggregate(test_counts_m_FP_neg_era, list(test_counts_m_FP_neg_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FP_neg_era <- scale(test_counts_m_m_FP_neg_era[,c(1:(ncol(test_counts_m_m_FP_neg_era)-1))])
test_counts_m_m_s_FP_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FP_neg_era), test_counts_m_m_FP_neg_era$tissueNumeric)


# erato FM up
test_counts <- counts(atacDDS_E_FM_rest[rownames(sign_E_D1_FM_rest_up),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FM_pos_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FM_pos_era <- aggregate(test_counts_m_FM_pos_era, list(test_counts_m_FM_pos_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FM_pos_era <- scale(test_counts_m_m_FM_pos_era[,c(1:(ncol(test_counts_m_m_FM_pos_era)-1))])
test_counts_m_m_s_FM_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FM_pos_era), test_counts_m_m_FM_pos_era$tissueNumeric)


# erato FM down
test_counts <- counts(atacDDS_E_FM_rest[rownames(sign_E_D1_FM_rest_down),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FM_neg_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FM_neg_era <- aggregate(test_counts_m_FM_neg_era, list(test_counts_m_FM_neg_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FM_neg_era <- scale(test_counts_m_m_FM_neg_era[,c(1:(ncol(test_counts_m_m_FM_neg_era)-1))])
test_counts_m_m_s_FM_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FM_neg_era), test_counts_m_m_FM_neg_era$tissueNumeric)



# erato FD up
test_counts <- counts(atacDDS_E_FD_rest[rownames(sign_E_D1_FD_rest_up),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FD_pos_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FD_pos_era <- aggregate(test_counts_m_FD_pos_era, list(test_counts_m_FD_pos_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FD_pos_era <- scale(test_counts_m_m_FD_pos_era[,c(1:(ncol(test_counts_m_m_FD_pos_era)-1))])
test_counts_m_m_s_FD_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FD_pos_era), test_counts_m_m_FD_pos_era$tissueNumeric)


# erato FD down
test_counts <- counts(atacDDS_E_FD_rest[rownames(sign_E_D1_FD_rest_down),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FD_neg_era <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FD_neg_era <- aggregate(test_counts_m_FD_neg_era, list(test_counts_m_FD_neg_era$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FD_neg_era <- scale(test_counts_m_m_FD_neg_era[,c(1:(ncol(test_counts_m_m_FD_neg_era)-1))])
test_counts_m_m_s_FD_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FD_neg_era), test_counts_m_m_FD_neg_era$tissueNumeric)




# melp FP up
test_counts <- counts(atacDDS_M_FP_rest[rownames(sign_M_D1_FP_rest_up),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FP_pos_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_FP_pos_melp <- aggregate(test_counts_m_FP_pos_melp, list(test_counts_m_FP_pos_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FP_pos_melp <- scale(test_counts_m_m_FP_pos_melp[,c(1:(ncol(test_counts_m_m_FP_pos_melp)-1))])
test_counts_m_m_s_FP_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FP_pos_melp), test_counts_m_m_FP_pos_melp$tissueNumeric)


# melp FP down
test_counts <- counts(atacDDS_M_FD_rest[rownames(sign_M_D1_FP_rest_down),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FP_neg_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_FP_neg_melp <- aggregate(test_counts_m_FP_neg_melp, list(test_counts_m_FP_neg_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FP_neg_melp <- scale(test_counts_m_m_FP_neg_melp[,c(1:(ncol(test_counts_m_m_FP_neg_melp)-1))])
test_counts_m_m_s_FP_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FP_neg_melp), test_counts_m_m_FP_neg_melp$tissueNumeric)


# melp FM up
test_counts <- counts(atacDDS_M_FM_rest[rownames(sign_M_D1_FM_rest_up),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FM_pos_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_FM_pos_melp <- aggregate(test_counts_m_FM_pos_melp, list(test_counts_m_FM_pos_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FM_pos_melp <- scale(test_counts_m_m_FM_pos_melp[,c(1:(ncol(test_counts_m_m_FM_pos_melp)-1))])
test_counts_m_m_s_FM_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FM_pos_melp), test_counts_m_m_FM_pos_melp$tissueNumeric)


# melp FM down
test_counts <- counts(atacDDS_M_FM_rest[rownames(sign_M_D1_FM_rest_down),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FM_neg_melp <- merge(test_countst, sampleInfo_subE, by = 'id')

test_counts_m_m_FM_neg_melp <- aggregate(test_counts_m_FM_neg_melp, list(test_counts_m_FM_neg_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FM_neg_melp <- scale(test_counts_m_m_FM_neg_melp[,c(1:(ncol(test_counts_m_m_FM_neg_melp)-1))])
test_counts_m_m_s_FM_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FM_neg_melp), test_counts_m_m_FM_neg_melp$tissueNumeric)



# melp FD up
test_counts <- counts(atacDDS_M_FD_rest[rownames(sign_M_D1_FD_rest_up),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FD_pos_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_FD_pos_melp <- aggregate(test_counts_m_FD_pos_melp, list(test_counts_m_FD_pos_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FD_pos_melp <- scale(test_counts_m_m_FD_pos_melp[,c(1:(ncol(test_counts_m_m_FD_pos_melp)-1))])
test_counts_m_m_s_FD_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FD_pos_melp), test_counts_m_m_FD_pos_melp$tissueNumeric)



# melp FD down
test_counts <- counts(atacDDS_M_FD_rest[rownames(sign_M_D1_FD_rest_down),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$id <- rownames(test_countst)

test_counts_m_FD_neg_melp <- merge(test_countst, sampleInfo_subM, by = 'id')

test_counts_m_m_FD_neg_melp <- aggregate(test_counts_m_FD_neg_melp, list(test_counts_m_FD_neg_melp$tissueNumeric), FUN=mean)[,c(2:(2+nrow(test_counts)),(2+nrow(test_counts))+5)]

test_counts_m_m_s_FD_neg_melp <- scale(test_counts_m_m_FD_neg_melp[,c(1:(ncol(test_counts_m_m_FD_neg_melp)-1))])
test_counts_m_m_s_FD_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FD_neg_melp), test_counts_m_m_FD_neg_melp$tissueNumeric)

## shared

test_counts_m_m_s_FP_pos_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FP_pos_era)[c(3:ncol(test_counts_m_m_s_FP_pos_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FP_pos_era)[c(3:ncol(test_counts_m_m_s_FP_pos_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FP_neg_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FP_neg_era)[c(3:ncol(test_counts_m_m_s_FP_neg_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FP_neg_era)[c(3:ncol(test_counts_m_m_s_FP_neg_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FP_pos_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FP_pos_melp)[c(3:ncol(test_counts_m_m_s_FP_pos_melp)-1)]), col = 'colnames(test_counts_m_m_s_FP_pos_melp)[c(3:ncol(test_counts_m_m_s_FP_pos_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FP_neg_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FP_neg_melp)[c(3:ncol(test_counts_m_m_s_FP_neg_melp)-1)]), col = 'colnames(test_counts_m_m_s_FP_neg_melp)[c(3:ncol(test_counts_m_m_s_FP_neg_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

test_counts_m_m_s_FP_pos_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FP_pos_era)[c(3:ncol(test_counts_m_m_s_FP_pos_era)-1)])  ,test_counts_m_m_s_FP_pos_era_splits )
test_counts_m_m_s_FP_neg_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FP_neg_era)[c(3:ncol(test_counts_m_m_s_FP_neg_era)-1)])  ,test_counts_m_m_s_FP_neg_era_splits )
test_counts_m_m_s_FP_pos_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FP_pos_melp)[c(3:ncol(test_counts_m_m_s_FP_pos_melp)-1)]), test_counts_m_m_s_FP_pos_melp_splits)
test_counts_m_m_s_FP_neg_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FP_neg_melp)[c(3:ncol(test_counts_m_m_s_FP_neg_melp)-1)]), test_counts_m_m_s_FP_neg_melp_splits)

test_counts_m_m_s_FP_pos_era_split  <- subset(test_counts_m_m_s_FP_pos_era_split , test_counts_m_m_s_FP_pos_era_split$scafM  != "NA")
test_counts_m_m_s_FP_neg_era_split  <- subset(test_counts_m_m_s_FP_neg_era_split , test_counts_m_m_s_FP_neg_era_split$scafM  != "NA")
test_counts_m_m_s_FP_pos_melp_split <- subset(test_counts_m_m_s_FP_pos_melp_split, test_counts_m_m_s_FP_pos_melp_split$scafE != "NA")
test_counts_m_m_s_FP_neg_melp_split <- subset(test_counts_m_m_s_FP_neg_melp_split, test_counts_m_m_s_FP_neg_melp_split$scafE != "NA")

test_counts_m_m_s_FP_pos_era_sharedn <- merge(test_counts_m_m_s_FP_pos_era_split, test_counts_m_m_s_FP_pos_melp_split, by = c('start','end'))[,3]
test_counts_m_m_s_FP_neg_era_sharedn <- merge(test_counts_m_m_s_FP_neg_era_split, test_counts_m_m_s_FP_neg_melp_split, by = c('start','end'))[,3]

test_counts_m_m_s_FP_pos_era_shared <-  as.data.frame(test_counts_m_m_s_FP_pos_era[, as.character(test_counts_m_m_s_FP_pos_era_sharedn)])
test_counts_m_m_s_FP_neg_era_shared <-  as.data.frame(test_counts_m_m_s_FP_neg_era[, as.character(test_counts_m_m_s_FP_neg_era_sharedn)])
test_counts_m_m_s_FP_pos_melp_shared <- as.data.frame(test_counts_m_m_s_FP_pos_melp[, as.character(test_counts_m_m_s_FP_pos_era_sharedn)])
test_counts_m_m_s_FP_neg_melp_shared <- as.data.frame(test_counts_m_m_s_FP_neg_melp[, as.character(test_counts_m_m_s_FP_neg_era_sharedn)])


test_counts_m_m_s_FM_pos_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FM_pos_era)[c(3:ncol(test_counts_m_m_s_FM_pos_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FM_pos_era)[c(3:ncol(test_counts_m_m_s_FM_pos_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FM_neg_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FM_neg_era)[c(3:ncol(test_counts_m_m_s_FM_neg_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FM_neg_era)[c(3:ncol(test_counts_m_m_s_FM_neg_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FM_pos_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FM_pos_melp)[c(3:ncol(test_counts_m_m_s_FM_pos_melp)-1)]), col = 'colnames(test_counts_m_m_s_FM_pos_melp)[c(3:ncol(test_counts_m_m_s_FM_pos_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FM_neg_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FM_neg_melp)[c(3:ncol(test_counts_m_m_s_FM_neg_melp)-1)]), col = 'colnames(test_counts_m_m_s_FM_neg_melp)[c(3:ncol(test_counts_m_m_s_FM_neg_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

test_counts_m_m_s_FM_pos_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FM_pos_era)[c(3:ncol(test_counts_m_m_s_FM_pos_era)-1)])  ,test_counts_m_m_s_FM_pos_era_splits )
test_counts_m_m_s_FM_neg_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FM_neg_era)[c(3:ncol(test_counts_m_m_s_FM_neg_era)-1)])  ,test_counts_m_m_s_FM_neg_era_splits )
test_counts_m_m_s_FM_pos_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FM_pos_melp)[c(3:ncol(test_counts_m_m_s_FM_pos_melp)-1)]), test_counts_m_m_s_FM_pos_melp_splits)
test_counts_m_m_s_FM_neg_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FM_neg_melp)[c(3:ncol(test_counts_m_m_s_FM_neg_melp)-1)]), test_counts_m_m_s_FM_neg_melp_splits)

test_counts_m_m_s_FM_pos_era_split  <- subset(test_counts_m_m_s_FM_pos_era_split , test_counts_m_m_s_FM_pos_era_split$scafM  != "NA")
test_counts_m_m_s_FM_neg_era_split  <- subset(test_counts_m_m_s_FM_neg_era_split , test_counts_m_m_s_FM_neg_era_split$scafM  != "NA")
test_counts_m_m_s_FM_pos_melp_split <- subset(test_counts_m_m_s_FM_pos_melp_split, test_counts_m_m_s_FM_pos_melp_split$scafE != "NA")
test_counts_m_m_s_FM_neg_melp_split <- subset(test_counts_m_m_s_FM_neg_melp_split, test_counts_m_m_s_FM_neg_melp_split$scafE != "NA")

test_counts_m_m_s_FM_pos_era_sharedn <- merge(test_counts_m_m_s_FM_pos_era_split, test_counts_m_m_s_FM_pos_melp_split, by = c('start','end'))[,3]
test_counts_m_m_s_FM_neg_era_sharedn <- merge(test_counts_m_m_s_FM_neg_era_split, test_counts_m_m_s_FM_neg_melp_split, by = c('start','end'))[,3]

test_counts_m_m_s_FM_pos_era_shared <-  as.data.frame(test_counts_m_m_s_FM_pos_era[, as.character(test_counts_m_m_s_FM_pos_era_sharedn)])
test_counts_m_m_s_FM_neg_era_shared <-  as.data.frame(test_counts_m_m_s_FM_neg_era[, as.character(test_counts_m_m_s_FM_neg_era_sharedn)])
test_counts_m_m_s_FM_pos_melp_shared <- as.data.frame(test_counts_m_m_s_FM_pos_melp[, as.character(test_counts_m_m_s_FM_pos_era_sharedn)])
test_counts_m_m_s_FM_neg_melp_shared <- as.data.frame(test_counts_m_m_s_FM_neg_melp[, as.character(test_counts_m_m_s_FM_neg_era_sharedn)])


test_counts_m_m_s_FD_pos_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FD_pos_era)[c(3:ncol(test_counts_m_m_s_FD_pos_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FD_pos_era)[c(3:ncol(test_counts_m_m_s_FD_pos_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FD_neg_era_splits <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FD_neg_era)[c(3:ncol(test_counts_m_m_s_FD_neg_era)-1)]  ), col = 'colnames(test_counts_m_m_s_FD_neg_era)[c(3:ncol(test_counts_m_m_s_FD_neg_era) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FD_pos_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FD_pos_melp)[c(3:ncol(test_counts_m_m_s_FD_pos_melp)-1)]), col = 'colnames(test_counts_m_m_s_FD_pos_melp)[c(3:ncol(test_counts_m_m_s_FD_pos_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
test_counts_m_m_s_FD_neg_melp_splits <- separate(data = as.data.frame(colnames(test_counts_m_m_s_FD_neg_melp)[c(3:ncol(test_counts_m_m_s_FD_neg_melp)-1)]), col = 'colnames(test_counts_m_m_s_FD_neg_melp)[c(3:ncol(test_counts_m_m_s_FD_neg_melp) - 1)]', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

test_counts_m_m_s_FD_pos_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FD_pos_era)[c(3:ncol(test_counts_m_m_s_FD_pos_era)-1)])  ,test_counts_m_m_s_FD_pos_era_splits )
test_counts_m_m_s_FD_neg_era_split  <- cbind(as.data.frame(colnames(test_counts_m_m_s_FD_neg_era)[c(3:ncol(test_counts_m_m_s_FD_neg_era)-1)])  ,test_counts_m_m_s_FD_neg_era_splits )
test_counts_m_m_s_FD_pos_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FD_pos_melp)[c(3:ncol(test_counts_m_m_s_FD_pos_melp)-1)]), test_counts_m_m_s_FD_pos_melp_splits)
test_counts_m_m_s_FD_neg_melp_split <- cbind(as.data.frame(colnames(test_counts_m_m_s_FD_neg_melp)[c(3:ncol(test_counts_m_m_s_FD_neg_melp)-1)]), test_counts_m_m_s_FD_neg_melp_splits)

test_counts_m_m_s_FD_pos_era_split  <- subset(test_counts_m_m_s_FD_pos_era_split , test_counts_m_m_s_FD_pos_era_split$scafM  != "NA")
test_counts_m_m_s_FD_neg_era_split  <- subset(test_counts_m_m_s_FD_neg_era_split , test_counts_m_m_s_FD_neg_era_split$scafM  != "NA")
test_counts_m_m_s_FD_pos_melp_split <- subset(test_counts_m_m_s_FD_pos_melp_split, test_counts_m_m_s_FD_pos_melp_split$scafE != "NA")
test_counts_m_m_s_FD_neg_melp_split <- subset(test_counts_m_m_s_FD_neg_melp_split, test_counts_m_m_s_FD_neg_melp_split$scafE != "NA")

test_counts_m_m_s_FD_pos_era_sharedn <- merge(test_counts_m_m_s_FD_pos_era_split, test_counts_m_m_s_FD_pos_melp_split, by = c('start','end'))[,3]
test_counts_m_m_s_FD_neg_era_sharedn <- merge(test_counts_m_m_s_FD_neg_era_split, test_counts_m_m_s_FD_neg_melp_split, by = c('start','end'))[,3]

test_counts_m_m_s_FD_pos_era_shared <-  as.data.frame(test_counts_m_m_s_FD_pos_era[, as.character(test_counts_m_m_s_FD_pos_era_sharedn)])
test_counts_m_m_s_FD_neg_era_shared <-  as.data.frame(test_counts_m_m_s_FD_neg_era[, as.character(test_counts_m_m_s_FD_neg_era_sharedn)])
test_counts_m_m_s_FD_pos_melp_shared <- as.data.frame(test_counts_m_m_s_FD_pos_melp[, as.character(test_counts_m_m_s_FD_pos_era_sharedn)])
test_counts_m_m_s_FD_neg_melp_shared <- as.data.frame(test_counts_m_m_s_FD_neg_melp[, as.character(test_counts_m_m_s_FD_neg_era_sharedn)])

test_counts_m_m_s_FP_pos_era_shared_split <- separate(data =  as.data.frame(colnames(test_counts_m_m_s_FP_pos_era_shared)), col = 'colnames(test_counts_m_m_s_FP_pos_era_shared)', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")

write.table(test_counts_m_m_s_FP_pos_era_splits[,c("scafE","startE","endE")], file = "FP_pos_era_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_FM_pos_era_splits[,c("scafE","startE","endE")], file = "FM_pos_era_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_FD_pos_era_splits[,c("scafE","startE","endE")], file = "FD_pos_era_06122021.bed", col.names = F, row.names = F, quote = F)

write.table(test_counts_m_m_s_FP_pos_melp_splits[,c("scafM","startM","endM")], file = "FP_pos_melp_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_FM_pos_melp_splits[,c("scafM","startM","endM")], file = "FM_pos_melp_06122021.bed", col.names = F, row.names = F, quote = F)
write.table(test_counts_m_m_s_FD_pos_melp_splits[,c("scafM","startM","endM")], file = "FD_pos_melp_06122021.bed", col.names = F, row.names = F, quote = F)

write.table(test_counts_m_m_s_FP_pos_era_shared_split[,c("scafE","startE","endE")], file = "FP_pos_era_shared_06122021.bed", col.names = F, row.names = F, quote = F)

##
par(mfcol=c(3,2))

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FP_pos_era)-1)){
  lines(test_counts_m_m_s_FP_pos_era[,ncol(test_counts_m_m_s_FP_pos_era)],  test_counts_m_m_s_FP_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_pos_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FP_pos_era_shared[,i], col = adjustcolor("green",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_era)-1)){
  lines(test_counts_m_m_s_FP_neg_era[,ncol(test_counts_m_m_s_FP_neg_era)],  test_counts_m_m_s_FP_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FP_neg_era_shared[,i], col = adjustcolor("green",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FM_pos_era)-1)){
  lines(test_counts_m_m_s_FM_pos_era[,ncol(test_counts_m_m_s_FM_pos_era)],  test_counts_m_m_s_FM_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_pos_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FM_pos_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_era)-1)){
  lines(test_counts_m_m_s_FM_neg_era[,ncol(test_counts_m_m_s_FM_neg_era)],  test_counts_m_m_s_FM_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FM_neg_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FD_pos_era)-1)){
  lines(test_counts_m_m_s_FD_pos_era[,ncol(test_counts_m_m_s_FD_pos_era)],  test_counts_m_m_s_FD_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_pos_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FD_pos_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_era)-1)){
  lines(test_counts_m_m_s_FD_neg_era[,ncol(test_counts_m_m_s_FD_neg_era)],  test_counts_m_m_s_FD_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_era_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FD_neg_era_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FP_pos_melp)-1)){
  lines(test_counts_m_m_s_FP_pos_melp[,ncol(test_counts_m_m_s_FP_pos_melp)],  test_counts_m_m_s_FP_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_pos_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FP_pos_melp_shared[,i], col = adjustcolor("green",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_melp)-1)){
  lines(test_counts_m_m_s_FP_neg_melp[,ncol(test_counts_m_m_s_FP_neg_melp)],  test_counts_m_m_s_FP_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FP_neg_melp_shared[,i], col = adjustcolor("green",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FM_pos_melp)-1)){
  lines(test_counts_m_m_s_FM_pos_melp[,ncol(test_counts_m_m_s_FM_pos_melp)],  test_counts_m_m_s_FM_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_pos_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FM_pos_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_melp)-1)){
  lines(test_counts_m_m_s_FM_neg_melp[,ncol(test_counts_m_m_s_FM_neg_melp)],  test_counts_m_m_s_FM_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FM_neg_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FD_pos_melp)-1)){
  lines(test_counts_m_m_s_FD_pos_melp[,ncol(test_counts_m_m_s_FD_pos_melp)],  test_counts_m_m_s_FD_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_pos_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FD_pos_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_melp)-1)){
  lines(test_counts_m_m_s_FD_neg_melp[,ncol(test_counts_m_m_s_FD_neg_melp)],  test_counts_m_m_s_FD_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_melp_shared))){
  lines(c(1,2,3),  test_counts_m_m_s_FD_neg_melp_shared[,i], col = adjustcolor("green",1), type = 'b')
}



####

sign_E_D1_FPup <- rbind(sign_E_D1_FP_FD_FPup, sign_E_D1_FP_FM_FPup)
sign_E_D1_FPup <- sign_E_D1_FPup[unique(c(rownames(sign_E_D1_FP_FD_FPup),rownames(sign_E_D1_FP_FM_FPup))),]

sign_E_D1_FPdown <- rbind(sign_E_D1_FP_FD_FPdown, sign_E_D1_FP_FM_FPdown)
sign_E_D1_FPdown <- sign_E_D1_FPdown[unique(c(rownames(sign_E_D1_FP_FD_FPdown),rownames(sign_E_D1_FP_FM_FPdown))),]

sign_E_D1_FPup_intersect <- merge(as.data.frame(sign_E_D1_FP_FD_FPup), as.data.frame(sign_E_D1_FP_FM_FPup), by = "row.names")
sign_E_D1_FPdown_intersect <- merge(as.data.frame(sign_E_D1_FP_FD_FPdown), as.data.frame(sign_E_D1_FP_FM_FPdown), by = "row.names")

nrow(sign_E_D1_FPup)
nrow(sign_E_D1_FPdown)
nrow(sign_E_D1_FPup_intersect)
nrow(sign_E_D1_FPdown_intersect)
nrow(sign_E_D1_FD_rest_up)
nrow(sign_E_D1_FD_rest_down)

sign_E_D1_FMup <- rbind(sign_E_D1_FP_FM_FPdown, sign_E_D1_FM_FD_FMup)
sign_E_D1_FMup <- sign_E_D1_FMup[unique(c(rownames(sign_E_D1_FP_FM_FPdown),rownames(sign_E_D1_FM_FD_FMup))),]

sign_E_D1_FMdown <- rbind(sign_E_D1_FP_FM_FPup, sign_E_D1_FM_FD_FMdown)
sign_E_D1_FMdown <- sign_E_D1_FMdown[unique(c(rownames(sign_E_D1_FP_FM_FPup),rownames(sign_E_D1_FM_FD_FMdown))),]

sign_E_D1_FMup_intersect <- merge(as.data.frame(sign_E_D1_FP_FM_FPdown), as.data.frame(sign_E_D1_FM_FD_FMup), by = "row.names")
sign_E_D1_FMdown_intersect <- merge(as.data.frame(sign_E_D1_FP_FM_FPup), as.data.frame(sign_E_D1_FM_FD_FMdown), by = "row.names")

nrow(sign_E_D1_FMup)
nrow(sign_E_D1_FMdown)
nrow(sign_E_D1_FMup_intersect)
nrow(sign_E_D1_FMdown_intersect)
nrow(sign_E_D1_FM_rest_up)
nrow(sign_E_D1_FM_rest_down)

sign_E_D1_FDup <- rbind(sign_E_D1_FP_FD_FPdown, sign_E_D1_FM_FD_FMdown)
sign_E_D1_FDup <- sign_E_D1_FDup[unique(c(rownames(sign_E_D1_FP_FD_FPdown),rownames(sign_E_D1_FM_FD_FMdown))),]

sign_E_D1_FDdown <- rbind(sign_E_D1_FP_FD_FPup, sign_E_D1_FM_FD_FMup)
sign_E_D1_FDdown <- sign_E_D1_FDdown[unique(c(rownames(sign_E_D1_FP_FD_FPup),rownames(sign_E_D1_FM_FD_FMup))),]

sign_E_D1_FDup_intersect <- merge(as.data.frame(sign_E_D1_FP_FD_FPdown), as.data.frame(sign_E_D1_FM_FD_FMdown), by = "row.names")
sign_E_D1_FDdown_intersect <- merge(as.data.frame(sign_E_D1_FP_FD_FPup), as.data.frame(sign_E_D1_FM_FD_FMup), by = "row.names")

nrow(sign_E_D1_FDup)
nrow(sign_E_D1_FDdown)
nrow(sign_E_D1_FDup_intersect)
nrow(sign_E_D1_FDdown_intersect)
nrow(sign_E_D1_FD_rest_up)
nrow(sign_E_D1_FD_rest_down)

sign_M_D1_FPup <- rbind(sign_M_D1_FP_FD_FPup, sign_M_D1_FP_FM_FPup)
sign_M_D1_FPup <- sign_M_D1_FPup[unique(c(rownames(sign_M_D1_FP_FD_FPup),rownames(sign_M_D1_FP_FM_FPup))),]

sign_M_D1_FPdown <- rbind(sign_M_D1_FP_FD_FPdown, sign_M_D1_FP_FM_FPdown)
sign_M_D1_FPdown <- sign_M_D1_FPdown[unique(c(rownames(sign_M_D1_FP_FD_FPdown),rownames(sign_M_D1_FP_FM_FPdown))),]

sign_M_D1_FPup_intersect <- merge(as.data.frame(sign_M_D1_FP_FD_FPup), as.data.frame(sign_M_D1_FP_FM_FPup), by = "row.names")
sign_M_D1_FPdown_intersect <- merge(as.data.frame(sign_M_D1_FP_FD_FPdown), as.data.frame(sign_M_D1_FP_FM_FPdown), by = "row.names")

nrow(sign_M_D1_FPup)
nrow(sign_M_D1_FPdown)
nrow(sign_M_D1_FPup_intersect)
nrow(sign_M_D1_FPdown_intersect)
nrow(sign_M_D1_FP_rest_up)
nrow(sign_M_D1_FP_rest_down)

sign_M_D1_FMup <- rbind(sign_M_D1_FP_FM_FPdown, sign_M_D1_FM_FD_FMup)
sign_M_D1_FMup <- sign_M_D1_FMup[unique(c(rownames(sign_M_D1_FP_FM_FPdown),rownames(sign_M_D1_FM_FD_FMup))),]

sign_M_D1_FMdown <- rbind(sign_M_D1_FP_FM_FPup, sign_M_D1_FM_FD_FMdown)
sign_M_D1_FMdown <- sign_M_D1_FMdown[unique(c(rownames(sign_M_D1_FP_FM_FPup),rownames(sign_M_D1_FM_FD_FMdown))),]

sign_M_D1_FMup_intersect <- merge(as.data.frame(sign_M_D1_FP_FM_FPdown), as.data.frame(sign_M_D1_FM_FD_FMup), by = "row.names")
sign_M_D1_FMdown_intersect <- merge(as.data.frame(sign_M_D1_FP_FM_FPup), as.data.frame(sign_M_D1_FM_FD_FMdown), by = "row.names")

nrow(sign_M_D1_FMup)
nrow(sign_M_D1_FMdown)
nrow(sign_M_D1_FMup_intersect)
nrow(sign_M_D1_FMdown_intersect)
nrow(sign_M_D1_FM_rest_up)
nrow(sign_M_D1_FM_rest_down)

sign_M_D1_FDup <- rbind(sign_M_D1_FP_FD_FPdown, sign_M_D1_FM_FD_FMdown)
sign_M_D1_FDup <- sign_M_D1_FDup[unique(c(rownames(sign_M_D1_FP_FD_FPdown),rownames(sign_M_D1_FM_FD_FMdown))),]

sign_M_D1_FDdown <- rbind(sign_M_D1_FP_FD_FPup, sign_M_D1_FM_FD_FMup)
sign_M_D1_FDdown <- sign_M_D1_FDdown[unique(c(rownames(sign_M_D1_FP_FD_FPup),rownames(sign_M_D1_FM_FD_FMup))),]

sign_M_D1_FDup_intersect <- merge(as.data.frame(sign_M_D1_FP_FD_FPdown), as.data.frame(sign_M_D1_FM_FD_FMdown), by = "row.names")
sign_M_D1_FDdown_intersect <- merge(as.data.frame(sign_M_D1_FP_FD_FPup), as.data.frame(sign_M_D1_FM_FD_FMup), by = "row.names")

nrow(sign_M_D1_FDup)
nrow(sign_M_D1_FDdown)
nrow(sign_M_D1_FDup_intersect)
nrow(sign_M_D1_FDdown_intersect)
nrow(sign_M_D1_FD_rest_up)
nrow(sign_M_D1_FD_rest_down)


sign_E_D2_FPup <- rbind(sign_E_D2_FP_FD_FPup, sign_E_D2_FP_FM_FPup)
sign_E_D2_FPup <- sign_E_D2_FPup[unique(c(rownames(sign_E_D2_FP_FD_FPup),rownames(sign_E_D2_FP_FM_FPup))),]

sign_E_D2_FPdown <- rbind(sign_E_D2_FP_FD_FPdown, sign_E_D2_FP_FM_FPdown)
sign_E_D2_FPdown <- sign_E_D2_FPdown[unique(c(rownames(sign_E_D2_FP_FD_FPdown),rownames(sign_E_D2_FP_FM_FPdown))),]

sign_E_D2_FPup_intersect <- merge(as.data.frame(sign_E_D2_FP_FD_FPup), as.data.frame(sign_E_D2_FP_FM_FPup), by = "row.names")
sign_E_D2_FPdown_intersect <- merge(as.data.frame(sign_E_D2_FP_FD_FPdown), as.data.frame(sign_E_D2_FP_FM_FPdown), by = "row.names")

nrow(sign_E_D2_FPup)
nrow(sign_E_D2_FPdown)
nrow(sign_E_D2_FPup_intersect)
nrow(sign_E_D2_FPdown_intersect)
nrow(sign_E_D2_FP_rest_up)
nrow(sign_E_D2_FP_rest_down)

sign_E_D2_FMup <- rbind(sign_E_D2_FP_FM_FPdown, sign_E_D2_FM_FD_FMup)
sign_E_D2_FMup <- sign_E_D2_FMup[unique(c(rownames(sign_E_D2_FP_FM_FPdown),rownames(sign_E_D2_FM_FD_FMup))),]

sign_E_D2_FMdown <- rbind(sign_E_D2_FP_FM_FPup, sign_E_D2_FM_FD_FMdown)
sign_E_D2_FMdown <- sign_E_D2_FMdown[unique(c(rownames(sign_E_D2_FP_FM_FPup),rownames(sign_E_D2_FM_FD_FMdown))),]

sign_E_D2_FMup_intersect <- merge(as.data.frame(sign_E_D2_FP_FM_FPdown), as.data.frame(sign_E_D2_FM_FD_FMup), by = "row.names")
sign_E_D2_FMdown_intersect <- merge(as.data.frame(sign_E_D2_FP_FM_FPup), as.data.frame(sign_E_D2_FM_FD_FMdown), by = "row.names")

nrow(sign_E_D2_FMup)
nrow(sign_E_D2_FMdown)
nrow(sign_E_D2_FMup_intersect)
nrow(sign_E_D2_FMdown_intersect)
nrow(sign_E_D2_FM_rest_up)
nrow(sign_E_D2_FM_rest_down)

sign_E_D2_FDup <- rbind(sign_E_D2_FP_FD_FPdown, sign_E_D2_FM_FD_FMdown)
sign_E_D2_FDup <- sign_E_D2_FDup[unique(c(rownames(sign_E_D2_FP_FD_FPdown),rownames(sign_E_D2_FM_FD_FMdown))),]

sign_E_D2_FDdown <- rbind(sign_E_D2_FP_FD_FPup, sign_E_D2_FM_FD_FMup)
sign_E_D2_FDdown <- sign_E_D2_FDdown[unique(c(rownames(sign_E_D2_FP_FD_FPup),rownames(sign_E_D2_FM_FD_FMup))),]

sign_E_D2_FDup_intersect <- merge(as.data.frame(sign_E_D2_FP_FD_FPdown), as.data.frame(sign_E_D2_FM_FD_FMdown), by = "row.names")
sign_E_D2_FDdown_intersect <- merge(as.data.frame(sign_E_D2_FP_FD_FPup), as.data.frame(sign_E_D2_FM_FD_FMup), by = "row.names")

nrow(sign_E_D2_FDup)
nrow(sign_E_D2_FDdown)
nrow(sign_E_D2_FDup_intersect)
nrow(sign_E_D2_FDdown_intersect)
nrow(sign_E_D2_FD_rest_up)
nrow(sign_E_D2_FD_rest_down)


sign_M_D2_FPup <- rbind(sign_M_D2_FP_FD_FPup, sign_M_D2_FP_FM_FPup)
sign_M_D2_FPup <- sign_M_D2_FPup[unique(c(rownames(sign_M_D2_FP_FD_FPup),rownames(sign_M_D2_FP_FM_FPup))),]

sign_M_D2_FPup_intersect <- merge(as.data.frame(sign_M_D2_FP_FD_FPup), as.data.frame(sign_M_D2_FP_FM_FPup), by = "row.names")
sign_M_D2_FPdown_intersect <- merge(as.data.frame(sign_M_D2_FP_FD_FPdown), as.data.frame(sign_M_D2_FP_FM_FPdown), by = "row.names")

nrow(sign_M_D2_FPup_intersect)
nrow(sign_M_D2_FPdown_intersect)
nrow(sign_M_D2_FP_rest_up)
nrow(sign_M_D2_FP_rest_down)

sign_M_D2_FMup <- rbind(sign_M_D2_FP_FM_FPdown, sign_M_D2_FM_FD_FMup)
sign_M_D2_FMup <- sign_M_D2_FMup[unique(c(rownames(sign_M_D2_FP_FM_FPdown),rownames(sign_M_D2_FM_FD_FMup))),]

sign_M_D2_FMup_intersect <- merge(as.data.frame(sign_M_D2_FP_FM_FPdown), as.data.frame(sign_M_D2_FM_FD_FMup), by = "row.names")
sign_M_D2_FMdown_intersect <- merge(as.data.frame(sign_M_D2_FP_FM_FPup), as.data.frame(sign_M_D2_FM_FD_FMdown), by = "row.names")

nrow(sign_M_D2_FMup_intersect)
nrow(sign_M_D2_FMdown_intersect)
nrow(sign_M_D2_FM_rest_up)
nrow(sign_M_D2_FM_rest_down)

sign_M_D2_FDup <- rbind(sign_M_D2_FP_FD_FPdown, sign_M_D2_FM_FD_FMdown)
sign_M_D2_FDup <- sign_M_D2_FDup[unique(c(rownames(sign_M_D2_FP_FD_FPdown),rownames(sign_M_D2_FM_FD_FMdown))),]

sign_M_D2_FDup_intersect <- merge(as.data.frame(sign_M_D2_FP_FD_FPdown), as.data.frame(sign_M_D2_FM_FD_FMdown), by = "row.names")
sign_M_D2_FDdown_intersect <- merge(as.data.frame(sign_M_D2_FP_FD_FPup), as.data.frame(sign_M_D2_FM_FD_FMup), by = "row.names")

nrow(sign_M_D2_FDup_intersect)
nrow(sign_M_D2_FDdown_intersect)
nrow(sign_M_D2_FD_rest_up)
nrow(sign_M_D2_FD_rest_down)



sign_E_D1D2_FPup_intersect <- merge(sign_E_D1_FPup_intersect,sign_E_D2_FPup_intersect, by = "Row.names")
sign_E_D1D2_FPdown_intersect <- merge(sign_E_D1_FPdown_intersect,sign_E_D2_FPdown_intersect, by = "Row.names")

sign_E_D1D2_FMup_intersect <- merge(sign_E_D1_FMup_intersect,sign_E_D2_FMup_intersect, by = "Row.names")
sign_E_D1D2_FMdown_intersect <- merge(sign_E_D1_FMdown_intersect,sign_E_D2_FMdown_intersect, by = "Row.names")

sign_E_D1D2_FDup_intersect <- merge(sign_E_D1_FDup_intersect,sign_E_D2_FDup_intersect, by = "Row.names")
sign_E_D1D2_FDdown_intersect <- merge(sign_E_D1_FDdown_intersect,sign_E_D2_FDdown_intersect, by = "Row.names")

nrow(sign_E_D1D2_FPup_intersect)
nrow(sign_E_D1D2_FPdown_intersect)
nrow(sign_E_D1D2_FMup_intersect)
nrow(sign_E_D1D2_FMdown_intersect)
nrow(sign_E_D1D2_FDup_intersect)
nrow(sign_E_D1D2_FDdown_intersect)

sign_E_D1D2_FPup_intersect <- merge(sign_E_D1_FPup_intersect,sign_E_D2_FPup_intersect, by = "Row.names")
sign_E_D1D2_FPdown_intersect <- merge(sign_E_D1_FPdown_intersect,sign_E_D2_FPdown_intersect, by = "Row.names")

sign_E_D1D2_FMup_intersect <- merge(sign_E_D1_FMup_intersect,sign_E_D2_FMup_intersect, by = "Row.names")
sign_E_D1D2_FMdown_intersect <- merge(sign_E_D1_FMdown_intersect,sign_E_D2_FMdown_intersect, by = "Row.names")

sign_E_D1D2_FDup_intersect <- merge(sign_E_D1_FDup_intersect,sign_E_D2_FDup_intersect, by = "Row.names")
sign_E_D1D2_FDdown_intersect <- merge(sign_E_D1_FDdown_intersect,sign_E_D2_FDdown_intersect, by = "Row.names")

nrow(sign_E_D1D2_FPup_intersect)
nrow(sign_E_D1D2_FPdown_intersect)
nrow(sign_E_D1D2_FMup_intersect)
nrow(sign_E_D1D2_FMdown_intersect)
nrow(sign_E_D1D2_FDup_intersect)
nrow(sign_E_D1D2_FDdown_intersect)

sign_E_D1D2_FPup <- merge(sign_E_D1_FPup, sign_E_D2_FPup, by = "row.names")
sign_E_D1D2_FPdown <- merge(sign_E_D1_FPdown, sign_E_D2_FPdown, by = "row.names")

sign_E_D1D2_FMup <- merge(sign_E_D1_FMup, sign_E_D2_FMup, by = "row.names")
sign_E_D1D2_FMdown <- merge(sign_E_D1_FMdown, sign_E_D2_FMdown, by = "row.names")

sign_E_D1D2_FDup <- merge(sign_E_D1_FDup, sign_E_D2_FDup, by = "row.names")
sign_E_D1D2_FDdown <- merge(sign_E_D1_FDdown, sign_E_D2_FDdown, by = "row.names")

nrow(sign_E_D1D2_FPup)
nrow(sign_E_D1D2_FPdown)
nrow(sign_E_D1D2_FMup)
nrow(sign_E_D1D2_FMdown)
nrow(sign_E_D1D2_FDup)
nrow(sign_E_D1D2_FDdown)


sign_E_D1D2_FP_rest_up <- merge(as.data.frame(sign_E_D1_FP_rest_up), as.data.frame(sign_E_D2_FP_rest_up), by = "row.names")
sign_E_D1D2_FP_rest_down <- merge(as.data.frame(sign_E_D1_FP_rest_down), as.data.frame(sign_E_D2_FP_rest_down), by = "row.names")

sign_E_D1D2_FM_rest_up <- merge(as.data.frame(sign_E_D1_FM_rest_up), as.data.frame(sign_E_D2_FM_rest_up), by = "row.names")
sign_E_D1D2_FM_rest_down <- merge(as.data.frame(sign_E_D1_FM_rest_down), as.data.frame(sign_E_D2_FM_rest_down), by = "row.names")

sign_E_D1D2_FD_rest_up <- merge(as.data.frame(sign_E_D1_FD_rest_up), as.data.frame(sign_E_D2_FD_rest_up), by = "row.names")
sign_E_D1D2_FD_rest_down <- merge(as.data.frame(sign_E_D1_FD_rest_down), as.data.frame(sign_E_D2_FD_rest_down), by = "row.names")

nrow(sign_E_D1D2_FP_rest_up)
nrow(sign_E_D1D2_FP_rest_down)
nrow(sign_E_D1D2_FM_rest_up)
nrow(sign_E_D1D2_FM_rest_down)
nrow(sign_E_D1D2_FD_rest_up)
nrow(sign_E_D1D2_FD_rest_down)




sign_list_E_dir <- list(sign_E_D1_FPup, sign_E_D1_FMup, sign_E_D1_FDup, sign_E_D1_FPdown, sign_E_D1_FMdown, sign_E_D1_FDdown, sign_E_D2_FPup, sign_E_D2_FMup, sign_E_D2_FDup)
sign_list_M_dir <- list(sign_M_D1_FPup, sign_M_D1_FMup, sign_M_D1_FDup, sign_M_D1_FPdown, sign_M_D1_FMdown, sign_M_D1_FDdown, sign_M_D2_FPup, sign_M_D2_FMup, sign_M_D2_FDup)

sign_list_E_dir <- list(sign_E_D1_FP_rest_up, sign_E_D1_FM_rest_up, sign_E_D1_FD_rest_up, sign_E_D1_FP_rest_down, sign_E_D1_FM_rest_down, sign_E_D1_FD_rest_down, sign_E_D2_FP_rest_up, sign_E_D2_FM_rest_up, sign_E_D2_FD_rest_up, sign_E_D2_FP_rest_down, sign_E_D2_FM_rest_down, sign_E_D2_FD_rest_down)
sign_list_M_dir <- list(sign_M_D1_FP_rest_up, sign_M_D1_FM_rest_up, sign_M_D1_FD_rest_up, sign_M_D1_FP_rest_down, sign_M_D1_FM_rest_down, sign_M_D1_FD_rest_down, sign_M_D2_FP_rest_up, sign_M_D2_FM_rest_up, sign_M_D2_FD_rest_up, sign_M_D2_FP_rest_down, sign_M_D2_FM_rest_down, sign_M_D2_FD_rest_down)



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

write.table(sign_list_E_df_all_dir[[1]][,c("scafE","startE","endE")], file = 'sections_D1_FPup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[2]][,c("scafE","startE","endE")], file = 'sections_D1_FMup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[3]][,c("scafE","startE","endE")], file = 'sections_D1_FDup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[4]][,c("scafE","startE","endE")], file = 'sections_D1_FPdown_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[5]][,c("scafE","startE","endE")], file = 'sections_D1_FMdown_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[6]][,c("scafE","startE","endE")], file = 'sections_D1_FDdown_rest_erato.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_E_df_all_dir[[7]][,c("scafE","startE","endE")], file = 'sections_D2_FPup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[8]][,c("scafE","startE","endE")], file = 'sections_D2_FMup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[9]][,c("scafE","startE","endE")], file = 'sections_D2_FDup_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[10]][,c("scafE","startE","endE")], file = 'sections_D2_FPdown_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[11]][,c("scafE","startE","endE")], file = 'sections_D2_FMdown_rest_erato.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[12]][,c("scafE","startE","endE")], file = 'sections_D2_FDdown_rest_erato.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("scafM","startM","endM")], file = 'sections_D1_FPup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[2]][,c("scafM","startM","endM")], file = 'sections_D1_FMup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[3]][,c("scafM","startM","endM")], file = 'sections_D1_FDup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[4]][,c("scafM","startM","endM")], file = 'sections_D1_FPdown_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[5]][,c("scafM","startM","endM")], file = 'sections_D1_FMdown_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[6]][,c("scafM","startM","endM")], file = 'sections_D1_FDdown_rest_melp.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[7]][,c("scafM","startM","endM")], file = 'sections_D2_FPup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[8]][,c("scafM","startM","endM")], file = 'sections_D2_FMup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[9]][,c("scafM","startM","endM")], file = 'sections_D2_FDup_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[10]][,c("scafM","startM","endM")], file = 'sections_D2_FPdown_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[11]][,c("scafM","startM","endM")], file = 'sections_D2_FMdown_rest_melp.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[12]][,c("scafM","startM","endM")], file = 'sections_D2_FDdown_rest_melp.bed', col.names = F, row.names = F, quote = F)

sections_D1_FPup_shared_E <- merge(sign_list_E_df_all_dir[[1]],sign_list_M_df_all_dir[[1]], by=c('start','end'))
sections_D1_FPdown_shared_E <- merge(sign_list_E_df_all_dir[[4]],sign_list_M_df_all_dir[[4]], by=c('start','end'))

nrow(sections_D1_FPup_shared_E)
nrow(sections_D1_FPdown_shared_E)

write.table(sections_D1_FPup_shared_E[,c("scafE.x","startE.x","endE.x")], file = "sections_D1_FPup_rest_shared.bed", col.names = F, row.names = F, quote = F)
write.table(sections_D1_FPdown_shared_E[,c("scafE.x","startE.x","endE.x")], file = "sections_D1_FPdown_rest_shared.bed", col.names = F, row.names = F, quote = F)

sections_D1_FMup_shared_E <- merge(sign_list_E_df_all_dir[[2]],sign_list_M_df_all_dir[[2]], by=c('start','end'))
sections_D1_FMdown_shared_E <- merge(sign_list_E_df_all_dir[[5]],sign_list_M_df_all_dir[[5]], by=c('start','end'))

nrow(sections_D1_FMup_shared_E)
nrow(sections_D1_FMdown_shared_E)

sections_D1_FDup_shared_E <- merge(sign_list_E_df_all_dir[[3]],sign_list_M_df_all_dir[[3]], by=c('start','end'))
sections_D1_FDdown_shared_E <- merge(sign_list_E_df_all_dir[[6]],sign_list_M_df_all_dir[[6]], by=c('start','end'))

nrow(sections_D1_FDup_shared_E)
nrow(sections_D1_FDdown_shared_E)

## for intersected
sign_list_E_dir <- list(sign_E_D1_FPup_intersect, sign_E_D1_FMup_intersect, sign_E_D1_FDup_intersect, sign_E_D1_FPdown_intersect, sign_E_D2_FPup_intersect, sign_E_D2_FMup_intersect, sign_E_D2_FDup_intersect)
sign_list_M_dir <- list(sign_M_D1_FPup_intersect, sign_M_D1_FMup_intersect, sign_M_D1_FDup_intersect, sign_M_D1_FPdown_intersect, sign_M_D2_FPup_intersect, sign_M_D2_FMup_intersect, sign_M_D2_FDup_intersect)


## for intersected
sign_list_E_df_dir <- list()
sign_list_E_df_all_dir <- list()
for(i in 1:length(sign_list_E_dir)){
  sign_list_E_df_dir[[i]] <- as.data.frame(sign_list_E_dir[[i]][,1])
  colnames(sign_list_E_df_dir[[i]]) <- 'rownames(sign_list_E_dir[[i]])'
  sign_list_E_df_dir[[i]] <- separate(data = sign_list_E_df_dir[[i]], col = 'rownames(sign_list_E_dir[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_E_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_E_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_E_df_all_dir[[i]] <- cbind(sign_list_E_dir[[i]], sign_list_E_df_dir[[i]])
}

sign_list_M_df_dir <- list()
sign_list_M_df_all_dir <- list()
for(i in 1:length(sign_list_M_dir)){
  sign_list_M_df_dir[[i]] <- as.data.frame(sign_list_M_dir[[i]][,1])
  colnames(sign_list_M_df_dir[[i]]) <- 'rownames(sign_list_M_dir[[i]])'
  sign_list_M_df_dir[[i]] <- separate(data = sign_list_M_df_dir[[i]], col = 'rownames(sign_list_M_dir[[i]])', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_M_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_M_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_M_df_all_dir[[i]] <- cbind(sign_list_M_dir[[i]], sign_list_M_df_dir[[i]])
} 

sign_list_E_dir <- list(sign_E_D1_FM_rest_up)

sign_list_E_df_dir <- list()
sign_list_E_df_all_dir <- list()
for(i in 1:length(sign_list_E_dir)){
  sign_list_E_df_dir[[i]] <- as.data.frame(rownames(sign_list_E_dir[[i]]))
  sign_list_E_df_dir[[i]] <- separate(data = sign_list_E_df_dir[[i]], col = 'rownames(sign_list_E_dir[[i]])', into = c("genome", "start", "end", "start2", "end2", "scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
  sign_list_E_df_dir[[i]][,c(2:5,7,8)] <- sapply(sign_list_E_df_dir[[i]][,c(2:5,7,8)],as.numeric) 
  sign_list_E_df_all_dir[[i]] <- cbind(sign_list_E_dir[[i]], sign_list_E_df_dir[[i]])
}
write.table(sign_list_E_df_all_dir[[1]][,c("scafE","startE","endE")], file = 'sections_D1_p05FC1_erato_FMupVSrest.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_E_df_all_dir[[1]][,c("scafE","startE","endE")], file = 'sections_D1_p05FC1_erato_FPup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[2]][,c("scafE","startE","endE")], file = 'sections_D1_p05FC1_erato_FMup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[3]][,c("scafE","startE","endE")], file = 'sections_D1_p05FC1_erato_FDup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[4]][,c("scafE","startE","endE")], file = 'sections_D2_p05FC1_erato_FPup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[5]][,c("scafE","startE","endE")], file = 'sections_D2_p05FC1_erato_FMup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[6]][,c("scafE","startE","endE")], file = 'sections_D2_p05FC1_erato_FDup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[7]][,c("scafE","startE","endE")], file = 'D1_p05FC1_erato_FWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[8]][,c("scafE","startE","endE")], file = 'D1_p05FC1_erato_HWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[9]][,c("scafE","startE","endE")], file = 'D2_p05FC1_erato_FWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[10]][,c("scafE","startE","endE")], file = 'D2_p05FC1_erato_HWup.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("scafM","startM","endM")], file = 'sections_D1_p05FC1_melp_FPup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[2]][,c("scafM","startM","endM")], file = 'sections_D1_p05FC1_melp_FMup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[3]][,c("scafM","startM","endM")], file = 'sections_D1_p05FC1_melp_FDup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[4]][,c("scafM","startM","endM")], file = 'sections_D2_p05FC1_melp_FPup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[5]][,c("scafM","startM","endM")], file = 'sections_D2_p05FC1_melp_FMup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[6]][,c("scafM","startM","endM")], file = 'sections_D2_p05FC1_melp_FDup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[7]][,c("scafM","startM","endM")], file = 'D2_p05FC1_melp_FWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[8]][,c("scafM","startM","endM")], file = 'D2_p05FC1_melp_HWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[9]][,c("scafM","startM","endM")], file = 'D2_p05FC1_melp_FWup.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[10]][,c("scafM","startM","endM")], file = 'D2_p05FC1_melp_HWup.bed', col.names = F, row.names = F, quote = F)


sections_D1_FPup_shared_E <- merge(sign_list_E_df_all_dir[[1]],sign_list_M_df_all_dir[[1]], by=c('start','end'))
sections_D1_FMup_shared_E <- merge(sign_list_E_df_all_dir[[2]],sign_list_M_df_all_dir[[2]], by=c('start','end'))
sections_D1_FDup_shared_E <- merge(sign_list_E_df_all_dir[[3]],sign_list_M_df_all_dir[[3]], by=c('start','end'))

sections_D1_FPdown_shared_E <- merge(sign_list_E_df_all_dir[[4]],sign_list_M_df_all_dir[[4]], by=c('start','end'))

nrow(sections_D1_FPup_shared_E)
nrow(sections_D1_FMup_shared_E)
nrow(sections_D1_FDup_shared_E)
nrow(sections_D1_FPdown_shared_E)

write.table(sections_D1_FPup_shared_E[,c("scafE.x","startE.x","endE.x")], file = "sections_D1_p05FC1_erato_FPup_shared.bed", col.names = F, row.names = F, quote = F)
write.table(sections_D1_FMup_shared_E[,c("scafE.x","startE.x","endE.x")], file = "sections_D1_p05FC1_erato_FMup_shared.bed", col.names = F, row.names = F, quote = F)
write.table(sections_D1_FDup_shared_E[,c("scafE.x","startE.x","endE.x")], file = "sections_D1_p05FC1_erato_FDup_shared.bed", col.names = F, row.names = F, quote = F)


sections_D1_FPup_shared_M <- merge(sign_list_E_df_all_dir[[1]],sign_list_M_df_all_dir[[1]], by=c('start','end'))
sections_D1_FMup_shared_M <- merge(sign_list_E_df_all_dir[[2]],sign_list_M_df_all_dir[[2]], by=c('start','end'))
sections_D1_FDup_shared_M <- merge(sign_list_E_df_all_dir[[3]],sign_list_M_df_all_dir[[3]], by=c('start','end'))

write.table(sections_D1_FPup_shared_M[,c("scafM.x","startM.x","endM.x")], file = "sections_D1_p05FC1_melp_FPup_shared.bed", col.names = F, row.names = F, quote = F)
write.table(sections_D1_FMup_shared_M[,c("scafM.x","startM.x","endM.x")], file = "sections_D1_p05FC1_melp_FMup_shared.bed", col.names = F, row.names = F, quote = F)
write.table(sections_D1_FDup_shared_M[,c("scafM.x","startM.x","endM.x")], file = "sections_D1_p05FC1_melp_FDup_shared.bed", col.names = F, row.names = F, quote = F)




write.table(sign_list_E_df_all_dir[[1]][,c("genome","start","end")], file = 'sections_D1_p05FC1_erato_FPup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[2]][,c("genome","start","end")], file = 'sections_D1_p05FC1_erato_FMup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[3]][,c("genome","start","end")], file = 'sections_D1_p05FC1_erato_FDup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[4]][,c("genome","start","end")], file = 'sections_D2_p05FC1_erato_FPup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[5]][,c("genome","start","end")], file = 'sections_D2_p05FC1_erato_FMup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[6]][,c("genome","start","end")], file = 'sections_D2_p05FC1_erato_FDup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[7]][,c("genome","start","end")], file = 'D1_p05FC1_erato_FWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[8]][,c("genome","start","end")], file = 'D1_p05FC1_erato_HWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[9]][,c("genome","start","end")], file = 'D2_p05FC1_erato_FWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_E_df_all_dir[[10]][,c("genome","start","end")], file = 'D2_p05FC1_erato_HWup_pan.bed', col.names = F, row.names = F, quote = F)

write.table(sign_list_M_df_all_dir[[1]][,c("genome","start","end")], file = 'sections_D1_p05FC1_melp_FPup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[2]][,c("genome","start","end")], file = 'sections_D1_p05FC1_melp_FMup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[3]][,c("genome","start","end")], file = 'sections_D1_p05FC1_melp_FDup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[4]][,c("genome","start","end")], file = 'sections_D2_p05FC1_melp_FPup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[5]][,c("genome","start","end")], file = 'sections_D2_p05FC1_melp_FMup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[6]][,c("genome","start","end")], file = 'sections_D2_p05FC1_melp_FDup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[7]][,c("genome","start","end")], file = 'D1_p05FC1_melp_FWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[8]][,c("genome","start","end")], file = 'D1_p05FC1_melp_HWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[9]][,c("genome","start","end")], file = 'D2_p05FC1_melp_FWup_pan.bed', col.names = F, row.names = F, quote = F)
write.table(sign_list_M_df_all_dir[[10]][,c("genome","start","end")], file = 'D2_p05FC1_melp_HWup_pan.bed', col.names = F, row.names = F, quote = F)





subset(sign_list_E_df_all_dir[[3]], sign_list_E_df_all_dir[[3]]$start == 381431014)











































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

sign_list_E <- list(sign_E_5th_D1D2_gradup, sign_E_D1_5thD2_D1up, sign_E_D2_5thD1_D2up)
sign_list_M <- list(sign_M_5th_D1D2_gradup, sign_M_D1_5thD2_D1up, sign_M_D2_5thD1_D2up)

colList <- c('#ff0000b5', '#009e73ff', '#0072b2ff')

sign_list_E_df <- list()
sign_list_E_df_all <- list()
for(i in 1:length(sign_list_E)){
  sign_list_E_df[[i]] <- as.data.frame(rownames(sign_list_E[[i]]))
  sign_list_E_df[[i]] <- separate(data = sign_list_E_df[[i]], col = 'rownames(sign_list_E[[i]])', into = c("genome", "start", "end", "start2", "end2"), sep = "_")
  sign_list_E_df[[i]][,c(2:5)] <- sapply(sign_list_E_df[[i]][,c(2:5)],as.numeric) 
  sign_list_E_df_all[[i]] <- cbind(sign_list_E[[i]], sign_list_E_df[[i]])
}

sign_list_M_df <- list()
sign_list_M_df_all <- list()
for(i in 1:length(sign_list_M)){
  sign_list_M_df[[i]] <- as.data.frame(rownames(sign_list_M[[i]]))
  sign_list_M_df[[i]] <- separate(data = sign_list_M_df[[i]], col = 'rownames(sign_list_M[[i]])', into = c("genome", "start", "end", "start2", "end2"), sep = "_")
  sign_list_M_df[[i]][,c(2:5)] <- sapply(sign_list_M_df[[i]][,c(2:5)],as.numeric) 
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

#shared peaks
plot(NULL, xlim=c(0,1000), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(0,1000,by=100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

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
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(0,1000), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  pval <- c()
  for(k in 1:3) {
    pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
    pval <- c(pval, pstat$p.value)
  }
  ptest <- rbind(ptest, pval)
}

plot(NULL, xlim=c(0,1000), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "shared")

#unique peaks
plot(NULL, xlim=c(1000,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1,seq(1000,0,by=-100),line=-1,col="black",col.ticks="black",col.axis="black", cex.axis=1, pos = 1)

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
  barplot(countpeaks, col=c(colList[c(1:3)]), horiz = TRUE, xlim=c(1000,0), axes=F, yaxt='none', xaxt='none', border=NA, xaxs="r")
  
  pval <- c()
  for(k in 1:3) {
    pstat <- chisq.test(cbind(countpeaks_M[k],countpeaks_E[k]))
    pval <- c(pval, pstat$p.value)
  }
  ptest <- rbind(ptest, pval)
}

plot(NULL, xlim=c(1000,0), ylim = c(0,1), axes=FALSE, ann=FALSE)
text(50,0.5, "unique")



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
legend(5, 0.5, c("5th","D1","D2"), fill=c('hotpink2', 'cornflowerblue', 'green'), border=NA, horiz=TRUE, xjust = 0.5, yjust = 0.5, box.col = NA)










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








sampleInfo_subE <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'erato' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))


sampleInfo_subM <- subset(sampleInfo_sub, ((sampleInfo_sub$stage == 'pupalDAY1') &
                                             sampleInfo_sub$species == 'melp' & sampleInfo_sub$tissue != 'HA' & sampleInfo_sub$tissue != 'HP'))



ctsF2E <- select(ctsF, c(as.character(sampleInfo_subE$id)))
ctsF2M <- select(ctsF, c(as.character(sampleInfo_subM$id)))

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


###
# run model
###
dds
atacDDS_E <- DESeq(ddsE)
atacDDS_M <- DESeq(ddsM)

atac_Rlog <- vst(atacDDS_M)
plotPCA(atac_Rlog, intgroup = c("tissue"), ntop = nrow(atac_Rlog)) + geom_text(aes(label=colnames(atac_Rlog)),vjust=2)


################################

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


load('res_E_grad.rda')
load('res_M_grad.rda')

load('res_E_FD_rest.rda')
load('res_M_FD_rest.rda')

load('res_E_FM_rest.rda')
load('res_M_FM_rest.rda')

load('res_E_FP_rest.rda')
load('res_M_FP_rest.rda')


sign_E_D1_grad_pos <- subset(res_E_grad, res_E_grad$padj < 0.05 & (res_E_grad$log2FoldChange) > 1)
sign_E_D1_grad_neg <- subset(res_E_grad, res_E_grad$padj < 0.05 & (res_E_grad$log2FoldChange) < -1)

sign_M_D1_grad_pos <- subset(res_M_grad, res_M_grad$padj < 0.05 & (res_M_grad$log2FoldChange) > 1)
sign_M_D1_grad_neg <- subset(res_M_grad, res_M_grad$padj < 0.05 & (res_M_grad$log2FoldChange) < -1)

sign_E_D1_FD_rest_up <- subset(res_E_FD_rest, res_E_FD_rest$padj < 0.05 & (res_E_FD_rest$log2FoldChange) < -1)
sign_E_D1_FD_rest_down <- subset(res_E_FD_rest, res_E_FD_rest$padj < 0.05 & (res_E_FD_rest$log2FoldChange) > 1)

sign_M_D1_FD_rest_up <- subset(res_M_FD_rest, res_M_FD_rest$padj < 0.05 & (res_M_FD_rest$log2FoldChange) < -1)
sign_M_D1_FD_rest_down <- subset(res_M_FD_rest, res_M_FD_rest$padj < 0.05 & (res_M_FD_rest$log2FoldChange) > 1)

sign_E_D1_FM_rest_up <- subset(res_E_FM_rest, res_E_FM_rest$padj < 0.05 & (res_E_FM_rest$log2FoldChange) < -1)
sign_E_D1_FM_rest_down <- subset(res_E_FM_rest, res_E_FM_rest$padj < 0.05 & (res_E_FM_rest$log2FoldChange) > 1)

sign_M_D1_FM_rest_up <- subset(res_M_FM_rest, res_M_FM_rest$padj < 0.05 & (res_M_FM_rest$log2FoldChange) < -1)
sign_M_D1_FM_rest_down <- subset(res_M_FM_rest, res_M_FM_rest$padj < 0.05 & (res_M_FM_rest$log2FoldChange) > 1)

sign_E_D1_FP_rest_up <- subset(res_E_FP_rest, res_E_FP_rest$padj < 0.05 & (res_E_FP_rest$log2FoldChange) < -1)
sign_E_D1_FP_rest_down <- subset(res_E_FP_rest, res_E_FP_rest$padj < 0.05 & (res_E_FP_rest$log2FoldChange) > 1)

sign_M_D1_FP_rest_up <- subset(res_M_FP_rest, res_M_FP_rest$padj < 0.05 & (res_M_FP_rest$log2FoldChange) < -1)
sign_M_D1_FP_rest_down <- subset(res_M_FP_rest, res_M_FP_rest$padj < 0.05 & (res_M_FP_rest$log2FoldChange) > 1)



# grad+

merged_EM_gradup <- merge(as.data.frame(res_E_grad), as.data.frame(sign_M_D1_grad_pos), by = 'row.names')
merged_ME_gradup <- merge(as.data.frame(res_M_grad), as.data.frame(sign_E_D1_grad_pos), by = 'row.names')



merged_shared_gradup <- merge(as.data.frame(sign_E_D1_grad_pos), as.data.frame(sign_M_D1_grad_pos), by = 'row.names')

# merged_EM_gradup$log2FoldChange.x <- merged_EM_gradup$log2FoldChange.x*-1
# merged_EM_gradup$log2FoldChange.y <- merged_EM_gradup$log2FoldChange.y*-1
# merged_ME_gradup$log2FoldChange.x <- merged_ME_gradup$log2FoldChange.x*-1
# merged_ME_gradup$log2FoldChange.y <- merged_ME_gradup$log2FoldChange.y*-1
# merged_shared_gradup$log2FoldChange.x <- merged_shared_gradup$log2FoldChange.x*-1
# merged_shared_gradup$log2FoldChange.y <- merged_shared_gradup$log2FoldChange.y*-1

merged_EM_gradup_unique <- merged_EM_gradup[!merged_EM_gradup$Row.names %in% merged_shared_gradup$Row.names, ]
merged_ME_gradup_unique <- merged_ME_gradup[!merged_ME_gradup$Row.names %in% merged_shared_gradup$Row.names, ]

par(mfrow=c(4,2))

plot(merged_EM_gradup_unique$log2FoldChange.x, merged_EM_gradup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h grad+")
par(new=T)
plot(merged_ME_gradup_unique$log2FoldChange.x, merged_ME_gradup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_gradup$log2FoldChange.x, merged_shared_gradup$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_gradup_unique$log2FoldChange.y ~ merged_EM_gradup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_gradup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_gradup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_gradup_unique$log2FoldChange.y ~ merged_ME_gradup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_gradup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_gradup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_gradup$log2FoldChange.y ~ merged_shared_gradup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_gradup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_gradup$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_gradup_unique$log2FoldChange.y ~ merged_EM_gradup_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_gradup_unique$log2FoldChange.y ~ merged_ME_gradup_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_gradup$log2FoldChange.y ~ merged_shared_gradup$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_grad_pos$Row.names <- rownames(sign_E_D1_grad_pos)
sign_E_D1_grad_pos <- separate(data = as.data.frame(sign_E_D1_grad_pos), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_grad_pos[,c(2:5,7,8)] <- sapply(sign_E_D1_grad_pos[,c(2:5,7,8)],as.numeric) 
sign_E_D1_grad_pos$Genome <- 'pan'

write.table(sign_E_D1_grad_pos[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_grad_pos.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_grad_pos$Row.names <- rownames(sign_M_D1_grad_pos)
sign_M_D1_grad_pos <- separate(data = as.data.frame(sign_M_D1_grad_pos), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_grad_pos[,c(2:5,7,8)] <- sapply(sign_M_D1_grad_pos[,c(2:5,7,8)],as.numeric) 
sign_M_D1_grad_pos$Genome <- 'pan'

sign_M_D1_grad_pos_write <- sign_M_D1_grad_pos[,c(18,10,11)]
sign_M_D1_grad_pos_write <- subset(sign_M_D1_grad_pos_write, as.numeric(sign_M_D1_grad_pos_write$end2)-as.numeric(sign_M_D1_grad_pos_write$start2) < 10000)
write.table(sign_M_D1_grad_pos_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_grad_pos.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_gradup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_gradup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_grad_pos_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_grad_pos_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_gradup_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_grad_pos)
nrow(sign_M_D1_grad_pos)
nrow(merged_shared_gradup)

#
nrow(sign_E_D1_grad_pos)-nrow(merged_ME_gradup[grep("NA", merged_ME_gradup$Row.names), ])
nrow(sign_M_D1_grad_pos)-nrow(merged_EM_gradup[grep("NA", merged_EM_gradup$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))


######################## Calculate for conserved set
sign_E_D1_grad_pos_bed <- sign_E_D1_grad_pos[,c(18,8,9)]
colnames(sign_E_D1_grad_pos_bed) <- c('chrom','start','end')
sign_E_D1_grad_pos_bed$start <- as.numeric(sign_E_D1_grad_pos_bed$start)
sign_E_D1_grad_pos_bed$end <- as.numeric(sign_E_D1_grad_pos_bed$end)

sign_M_D1_grad_pos_bed <- sign_M_D1_grad_pos_write
colnames(sign_M_D1_grad_pos_bed) <- c('chrom','start','end')
sign_M_D1_grad_pos_bed$start <- as.numeric(sign_M_D1_grad_pos_bed$start)
sign_M_D1_grad_pos_bed$end <- as.numeric(sign_M_D1_grad_pos_bed$end)

merged_shared_gradup_bed <- merged_shared_gradup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_gradup_bed) <- c('chrom','start','end')
merged_shared_gradup_bed$start <- as.numeric(merged_shared_gradup_bed$start)
merged_shared_gradup_bed$end <- as.numeric(merged_shared_gradup_bed$end)

sign_E_D1_grad_pos_bed_match <- bed_intersect(sign_E_D1_grad_pos_bed,conserved_E)
sign_M_D1_grad_pos_bed_match <- bed_intersect(sign_M_D1_grad_pos_bed,conserved_M)
merged_shared_gradup_bed_match <- bed_intersect(merged_shared_gradup_bed,conserved_S)


sign_E_D1_grad_pos_cons <- merge(sign_E_D1_grad_pos, sign_E_D1_grad_pos_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_grad_pos_cons <- merge(sign_M_D1_grad_pos, sign_M_D1_grad_pos_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_gradup_cons <- merge(merged_ME_gradup, sign_E_D1_grad_pos_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_gradup_cons <- merge(merged_EM_gradup, sign_M_D1_grad_pos_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D2HWup_cons[grep("NA", merged_ME_D2HWup_cons$Row.names), ])
nrow(merged_EM_D2HWup_cons[grep("NA", merged_EM_D2HWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_grad_pos_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_grad_pos_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_gradup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_grad_pos_bed_match)
nrow(sign_M_D1_grad_pos_bed_match)
nrow(merged_shared_gradup_bed_match)

nrow(merged_ME_gradup_cons[grep("NA", merged_ME_gradup_cons$Row.names), ])
nrow(merged_EM_gradup_cons[grep("NA", merged_EM_gradup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# grad-

merged_EM_graddown <- merge(as.data.frame(res_E_grad), as.data.frame(sign_M_D1_grad_neg), by = 'row.names')
merged_ME_graddown <- merge(as.data.frame(res_M_grad), as.data.frame(sign_E_D1_grad_neg), by = 'row.names')



merged_shared_graddown <- merge(as.data.frame(sign_E_D1_grad_neg), as.data.frame(sign_M_D1_grad_neg), by = 'row.names')

merged_EM_graddown$log2FoldChange.x <- merged_EM_graddown$log2FoldChange.x*-1
merged_EM_graddown$log2FoldChange.y <- merged_EM_graddown$log2FoldChange.y*-1
merged_ME_graddown$log2FoldChange.x <- merged_ME_graddown$log2FoldChange.x*-1
merged_ME_graddown$log2FoldChange.y <- merged_ME_graddown$log2FoldChange.y*-1
merged_shared_graddown$log2FoldChange.x <- merged_shared_graddown$log2FoldChange.x*-1
merged_shared_graddown$log2FoldChange.y <- merged_shared_graddown$log2FoldChange.y*-1

merged_EM_graddown_unique <- merged_EM_graddown[!merged_EM_graddown$Row.names %in% merged_shared_graddown$Row.names, ]
merged_ME_graddown_unique <- merged_ME_graddown[!merged_ME_graddown$Row.names %in% merged_shared_graddown$Row.names, ]

plot(merged_EM_graddown_unique$log2FoldChange.x, merged_EM_graddown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h grad-")
par(new=T)
plot(merged_ME_graddown_unique$log2FoldChange.x, merged_ME_graddown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_graddown$log2FoldChange.x, merged_shared_graddown$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_graddown_unique$log2FoldChange.y ~ merged_EM_graddown_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_graddown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_graddown_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_graddown_unique$log2FoldChange.y ~ merged_ME_graddown_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_graddown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_graddown_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_graddown$log2FoldChange.y ~ merged_shared_graddown$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_graddown$log2FoldChange.x)), x2 = max(na.omit(merged_shared_graddown$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_graddown_unique$log2FoldChange.y ~ merged_EM_graddown_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_graddown_unique$log2FoldChange.y ~ merged_ME_graddown_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_graddown$log2FoldChange.y ~ merged_shared_graddown$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_grad_neg$Row.names <- rownames(sign_E_D1_grad_neg)
sign_E_D1_grad_neg <- separate(data = as.data.frame(sign_E_D1_grad_neg), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_grad_neg[,c(2:5,7,8)] <- sapply(sign_E_D1_grad_neg[,c(2:5,7,8)],as.numeric) 
sign_E_D1_grad_neg$Genome <- 'pan'

write.table(sign_E_D1_grad_neg[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_grad_neg.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_grad_neg$Row.names <- rownames(sign_M_D1_grad_neg)
sign_M_D1_grad_neg <- separate(data = as.data.frame(sign_M_D1_grad_neg), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_grad_neg[,c(2:5,7,8)] <- sapply(sign_M_D1_grad_neg[,c(2:5,7,8)],as.numeric) 
sign_M_D1_grad_neg$Genome <- 'pan'

sign_M_D1_grad_neg_write <- sign_M_D1_grad_neg[,c(18,10,11)]
sign_M_D1_grad_neg_write <- subset(sign_M_D1_grad_neg_write, as.numeric(sign_M_D1_grad_neg_write$end2)-as.numeric(sign_M_D1_grad_neg_write$start2) < 10000)
write.table(sign_M_D1_grad_neg_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_grad_neg.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_graddown[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_graddown.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_grad_neg_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_grad_neg_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_graddown_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_grad_neg)
nrow(sign_M_D1_grad_neg)
nrow(merged_shared_graddown)

#
nrow(sign_E_D1_grad_neg)-nrow(merged_ME_graddown[grep("NA", merged_ME_graddown$Row.names), ])
nrow(sign_M_D1_grad_neg)-nrow(merged_EM_graddown[grep("NA", merged_EM_graddown$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_grad_neg_bed <- sign_E_D1_grad_neg[,c(18,8,9)]
colnames(sign_E_D1_grad_neg_bed) <- c('chrom','start','end')
sign_E_D1_grad_neg_bed$start <- as.numeric(sign_E_D1_grad_neg_bed$start)
sign_E_D1_grad_neg_bed$end <- as.numeric(sign_E_D1_grad_neg_bed$end)

sign_M_D1_grad_neg_bed <- sign_M_D1_grad_neg_write
colnames(sign_M_D1_grad_neg_bed) <- c('chrom','start','end')
sign_M_D1_grad_neg_bed$start <- as.numeric(sign_M_D1_grad_neg_bed$start)
sign_M_D1_grad_neg_bed$end <- as.numeric(sign_M_D1_grad_neg_bed$end)

merged_shared_graddown_bed <- merged_shared_graddown[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_graddown_bed) <- c('chrom','start','end')
merged_shared_graddown_bed$start <- as.numeric(merged_shared_graddown_bed$start)
merged_shared_graddown_bed$end <- as.numeric(merged_shared_graddown_bed$end)

sign_E_D1_grad_neg_bed_match <- bed_intersect(sign_E_D1_grad_neg_bed,conserved_E)
sign_M_D1_grad_neg_bed_match <- bed_intersect(sign_M_D1_grad_neg_bed,conserved_M)
merged_shared_graddown_bed_match <- bed_intersect(merged_shared_graddown_bed,conserved_S)


sign_E_D1_grad_neg_cons <- merge(sign_E_D1_grad_neg, sign_E_D1_grad_neg_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_grad_neg_cons <- merge(sign_M_D1_grad_neg, sign_M_D1_grad_neg_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_graddown_cons <- merge(merged_ME_graddown, sign_E_D1_grad_neg_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_graddown_cons <- merge(merged_EM_graddown, sign_M_D1_grad_neg_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_D2HWup_cons[grep("NA", merged_ME_D2HWup_cons$Row.names), ])
nrow(merged_EM_D2HWup_cons[grep("NA", merged_EM_D2HWup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_grad_neg_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_grad_neg_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_graddown_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_grad_neg_bed_match)
nrow(sign_M_D1_grad_neg_bed_match)
nrow(merged_shared_graddown_bed_match)

nrow(merged_ME_graddown_cons[grep("NA", merged_ME_graddown_cons$Row.names), ])
nrow(merged_EM_graddown_cons[grep("NA", merged_EM_graddown_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# FD up 

merged_EM_FDup <- merge(as.data.frame(res_E_FD_rest), as.data.frame(sign_M_D1_FD_rest_up), by = 'row.names')
merged_ME_FDup <- merge(as.data.frame(res_M_FD_rest), as.data.frame(sign_E_D1_FD_rest_up), by = 'row.names')



merged_shared_FDup <- merge(as.data.frame(sign_E_D1_FD_rest_up), as.data.frame(sign_M_D1_FD_rest_up), by = 'row.names')

merged_EM_FDup$log2FoldChange.x <- merged_EM_FDup$log2FoldChange.x*-1
merged_EM_FDup$log2FoldChange.y <- merged_EM_FDup$log2FoldChange.y*-1
merged_ME_FDup$log2FoldChange.x <- merged_ME_FDup$log2FoldChange.x*-1
merged_ME_FDup$log2FoldChange.y <- merged_ME_FDup$log2FoldChange.y*-1
merged_shared_FDup$log2FoldChange.x <- merged_shared_FDup$log2FoldChange.x*-1
merged_shared_FDup$log2FoldChange.y <- merged_shared_FDup$log2FoldChange.y*-1

merged_EM_FDup_unique <- merged_EM_FDup[!merged_EM_FDup$Row.names %in% merged_shared_FDup$Row.names, ]
merged_ME_FDup_unique <- merged_ME_FDup[!merged_ME_FDup$Row.names %in% merged_shared_FDup$Row.names, ]

plot(merged_EM_FDup_unique$log2FoldChange.x, merged_EM_FDup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FD up")
par(new=T)
plot(merged_ME_FDup_unique$log2FoldChange.x, merged_ME_FDup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FDup$log2FoldChange.x, merged_shared_FDup$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FDup_unique$log2FoldChange.y ~ merged_EM_FDup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FDup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FDup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FDup_unique$log2FoldChange.y ~ merged_ME_FDup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FDup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FDup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FDup$log2FoldChange.y ~ merged_shared_FDup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FDup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FDup$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FDup_unique$log2FoldChange.y ~ merged_EM_FDup_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FDup_unique$log2FoldChange.y ~ merged_ME_FDup_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FDup$log2FoldChange.y ~ merged_shared_FDup$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FD_rest_up$Row.names <- rownames(sign_E_D1_FD_rest_up)
sign_E_D1_FD_rest_up <- separate(data = as.data.frame(sign_E_D1_FD_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FD_rest_up[,c(2:5,7,8)] <- sapply(sign_E_D1_FD_rest_up[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FD_rest_up$Genome <- 'pan'

write.table(sign_E_D1_FD_rest_up[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FD_rest_up.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FD_rest_up$Row.names <- rownames(sign_M_D1_FD_rest_up)
sign_M_D1_FD_rest_up <- separate(data = as.data.frame(sign_M_D1_FD_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FD_rest_up[,c(2:5,7,8)] <- sapply(sign_M_D1_FD_rest_up[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FD_rest_up$Genome <- 'pan'

sign_M_D1_FD_rest_up_write <- sign_M_D1_FD_rest_up[,c(18,10,11)]
sign_M_D1_FD_rest_up_write <- subset(sign_M_D1_FD_rest_up_write, as.numeric(sign_M_D1_FD_rest_up_write$end2)-as.numeric(sign_M_D1_FD_rest_up_write$start2) < 10000)
write.table(sign_M_D1_FD_rest_up_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FD_rest_up.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FDup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FDup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FD_rest_up_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FD_rest_up_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FDup_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FD_rest_up)
nrow(sign_M_D1_FD_rest_up)
nrow(merged_shared_FDup)

#
nrow(sign_E_D1_FD_rest_up)-nrow(merged_ME_FDup[grep("NA", merged_ME_FDup$Row.names), ])
nrow(sign_M_D1_FD_rest_up)-nrow(merged_EM_FDup[grep("NA", merged_EM_FDup$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FD_rest_up_bed <- sign_E_D1_FD_rest_up[,c(18,8,9)]
colnames(sign_E_D1_FD_rest_up_bed) <- c('chrom','start','end')
sign_E_D1_FD_rest_up_bed$start <- as.numeric(sign_E_D1_FD_rest_up_bed$start)
sign_E_D1_FD_rest_up_bed$end <- as.numeric(sign_E_D1_FD_rest_up_bed$end)

sign_M_D1_FD_rest_up_bed <- sign_M_D1_FD_rest_up_write
colnames(sign_M_D1_FD_rest_up_bed) <- c('chrom','start','end')
sign_M_D1_FD_rest_up_bed$start <- as.numeric(sign_M_D1_FD_rest_up_bed$start)
sign_M_D1_FD_rest_up_bed$end <- as.numeric(sign_M_D1_FD_rest_up_bed$end)

merged_shared_shared_FDup_bed <- merged_shared_shared_FDup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_shared_FDup_bed) <- c('chrom','start','end')
merged_shared_shared_FDup_bed$start <- as.numeric(merged_shared_shared_FDup_bed$start)
merged_shared_shared_FDup_bed$end <- as.numeric(merged_shared_shared_FDup_bed$end)

sign_E_D1_FD_rest_up_bed_match <- bed_intersect(sign_E_D1_FD_rest_up_bed,conserved_E)
sign_M_D1_FD_rest_up_bed_match <- bed_intersect(sign_M_D1_FD_rest_up_bed,conserved_M)
merged_shared_shared_FDup_bed_match <- bed_intersect(merged_shared_shared_FDup_bed,conserved_S)


sign_E_D1_FD_rest_up_cons <- merge(sign_E_D1_FD_rest_up, sign_E_D1_FD_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FD_rest_up_cons <- merge(sign_M_D1_FD_rest_up, sign_M_D1_FD_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FDup_cons <- merge(merged_ME_FDup, sign_E_D1_FD_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FDup_cons <- merge(merged_EM_FDup, sign_M_D1_FD_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FDup_cons[grep("NA", merged_ME_FDup_cons$Row.names), ])
nrow(merged_EM_FDup_cons[grep("NA", merged_EM_FDup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FD_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FD_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_shared_FDup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FD_rest_up_bed_match)
nrow(sign_M_D1_FD_rest_up_bed_match)
nrow(merged_shared_shared_FDup_bed_match)

nrow(merged_ME_FDup_cons[grep("NA", merged_ME_FDup_cons$Row.names), ])
nrow(merged_EM_FDup_cons[grep("NA", merged_EM_FDup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# FD down 

merged_EM_FDdown <- merge(as.data.frame(res_E_FD_rest), as.data.frame(sign_M_D1_FD_rest_down), by = 'row.names')
merged_ME_FDdown <- merge(as.data.frame(res_M_FD_rest), as.data.frame(sign_E_D1_FD_rest_down), by = 'row.names')



merged_shared_FDdown <- merge(as.data.frame(sign_E_D1_FD_rest_down), as.data.frame(sign_M_D1_FD_rest_down), by = 'row.names')

# merged_EM_FDdown$log2FoldChange.x <- merged_EM_FDdown$log2FoldChange.x*-1
# merged_EM_FDdown$log2FoldChange.y <- merged_EM_FDdown$log2FoldChange.y*-1
# merged_ME_FDdown$log2FoldChange.x <- merged_ME_FDdown$log2FoldChange.x*-1
# merged_ME_FDdown$log2FoldChange.y <- merged_ME_FDdown$log2FoldChange.y*-1
# merged_shared_FDdown$log2FoldChange.x <- merged_shared_FDdown$log2FoldChange.x*-1
# merged_shared_FDdown$log2FoldChange.y <- merged_shared_FDdown$log2FoldChange.y*-1

merged_EM_FDdown_unique <- merged_EM_FDdown[!merged_EM_FDdown$Row.names %in% merged_shared_FDdown$Row.names, ]
merged_ME_FDdown_unique <- merged_ME_FDdown[!merged_ME_FDdown$Row.names %in% merged_shared_FDdown$Row.names, ]

plot(merged_EM_FDdown_unique$log2FoldChange.x, merged_EM_FDdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FD down")
par(new=T)
plot(merged_ME_FDdown_unique$log2FoldChange.x, merged_ME_FDdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FDdown$log2FoldChange.x, merged_shared_FDdown$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FDdown_unique$log2FoldChange.y ~ merged_EM_FDdown_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FDdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FDdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FDdown_unique$log2FoldChange.y ~ merged_ME_FDdown_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FDdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FDdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FDdown$log2FoldChange.y ~ merged_shared_FDdown$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FDdown$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FDdown$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FDdown_unique$log2FoldChange.y ~ merged_EM_FDdown_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FDdown_unique$log2FoldChange.y ~ merged_ME_FDdown_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FDdown$log2FoldChange.y ~ merged_shared_FDdown$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FD_rest_down$Row.names <- rownames(sign_E_D1_FD_rest_down)
sign_E_D1_FD_rest_down <- separate(data = as.data.frame(sign_E_D1_FD_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FD_rest_down[,c(2:5,7,8)] <- sapply(sign_E_D1_FD_rest_down[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FD_rest_down$Genome <- 'pan'

write.table(sign_E_D1_FD_rest_down[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FD_rest_down.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FD_rest_down$Row.names <- rownames(sign_M_D1_FD_rest_down)
sign_M_D1_FD_rest_down <- separate(data = as.data.frame(sign_M_D1_FD_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FD_rest_down[,c(2:5,7,8)] <- sapply(sign_M_D1_FD_rest_down[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FD_rest_down$Genome <- 'pan'

sign_M_D1_FD_rest_down_write <- sign_M_D1_FD_rest_down[,c(18,10,11)]
sign_M_D1_FD_rest_down_write <- subset(sign_M_D1_FD_rest_down_write, as.numeric(sign_M_D1_FD_rest_down_write$end2)-as.numeric(sign_M_D1_FD_rest_down_write$start2) < 10000)
write.table(sign_M_D1_FD_rest_down_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FD_rest_down.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FDdown[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FDdown.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FD_rest_down_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FD_rest_down_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FDdown_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FD_rest_down)
nrow(sign_M_D1_FD_rest_down)
nrow(merged_shared_FDdown)

#
nrow(sign_E_D1_FD_rest_down)-nrow(merged_ME_FDdown[grep("NA", merged_ME_FDdown$Row.names), ])
nrow(sign_M_D1_FD_rest_down)-nrow(merged_EM_FDdown[grep("NA", merged_EM_FDdown$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FD_rest_down_bed <- sign_E_D1_FD_rest_down[,c(18,8,9)]
colnames(sign_E_D1_FD_rest_down_bed) <- c('chrom','start','end')
sign_E_D1_FD_rest_down_bed$start <- as.numeric(sign_E_D1_FD_rest_down_bed$start)
sign_E_D1_FD_rest_down_bed$end <- as.numeric(sign_E_D1_FD_rest_down_bed$end)

sign_M_D1_FD_rest_down_bed <- sign_M_D1_FD_rest_down_write
colnames(sign_M_D1_FD_rest_down_bed) <- c('chrom','start','end')
sign_M_D1_FD_rest_down_bed$start <- as.numeric(sign_M_D1_FD_rest_down_bed$start)
sign_M_D1_FD_rest_down_bed$end <- as.numeric(sign_M_D1_FD_rest_down_bed$end)

merged_shared_shared_FDdown_bed <- merged_shared_shared_FDdown[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_shared_FDdown_bed) <- c('chrom','start','end')
merged_shared_shared_FDdown_bed$start <- as.numeric(merged_shared_shared_FDdown_bed$start)
merged_shared_shared_FDdown_bed$end <- as.numeric(merged_shared_shared_FDdown_bed$end)

sign_E_D1_FD_rest_down_bed_match <- bed_intersect(sign_E_D1_FD_rest_down_bed,conserved_E)
sign_M_D1_FD_rest_down_bed_match <- bed_intersect(sign_M_D1_FD_rest_down_bed,conserved_M)
merged_shared_shared_FDdown_bed_match <- bed_intersect(merged_shared_shared_FDdown_bed,conserved_S)


sign_E_D1_FD_rest_down_cons <- merge(sign_E_D1_FD_rest_down, sign_E_D1_FD_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FD_rest_down_cons <- merge(sign_M_D1_FD_rest_down, sign_M_D1_FD_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FDdown_cons <- merge(merged_ME_FDdown, sign_E_D1_FD_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FDdown_cons <- merge(merged_EM_FDdown, sign_M_D1_FD_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FDdown_cons[grep("NA", merged_ME_FDdown_cons$Row.names), ])
nrow(merged_EM_FDdown_cons[grep("NA", merged_EM_FDdown_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FD_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FD_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_shared_FDdown_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FD_rest_down_bed_match)
nrow(sign_M_D1_FD_rest_down_bed_match)
nrow(merged_shared_shared_FDdown_bed_match)

nrow(merged_ME_FDdown_cons[grep("NA", merged_ME_FDdown_cons$Row.names), ])
nrow(merged_EM_FDdown_cons[grep("NA", merged_EM_FDdown_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# FM up 

merged_EM_FMup <- merge(as.data.frame(res_E_FM_rest), as.data.frame(sign_M_D1_FM_rest_up), by = 'row.names')
merged_ME_FMup <- merge(as.data.frame(res_M_FM_rest), as.data.frame(sign_E_D1_FM_rest_up), by = 'row.names')



merged_shared_FMup <- merge(as.data.frame(sign_E_D1_FM_rest_up), as.data.frame(sign_M_D1_FM_rest_up), by = 'row.names')

merged_EM_FMup$log2FoldChange.x <- merged_EM_FMup$log2FoldChange.x*-1
merged_EM_FMup$log2FoldChange.y <- merged_EM_FMup$log2FoldChange.y*-1
merged_ME_FMup$log2FoldChange.x <- merged_ME_FMup$log2FoldChange.x*-1
merged_ME_FMup$log2FoldChange.y <- merged_ME_FMup$log2FoldChange.y*-1
merged_shared_FMup$log2FoldChange.x <- merged_shared_FMup$log2FoldChange.x*-1
merged_shared_FMup$log2FoldChange.y <- merged_shared_FMup$log2FoldChange.y*-1

merged_EM_FMup_unique <- merged_EM_FMup[!merged_EM_FMup$Row.names %in% merged_shared_FMup$Row.names, ]
merged_ME_FMup_unique <- merged_ME_FMup[!merged_ME_FMup$Row.names %in% merged_shared_FMup$Row.names, ]

plot(merged_EM_FMup_unique$log2FoldChange.x, merged_EM_FMup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FM up")
par(new=T)
plot(merged_ME_FMup_unique$log2FoldChange.x, merged_ME_FMup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FMup$log2FoldChange.x, merged_shared_FMup$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FMup_unique$log2FoldChange.y ~ merged_EM_FMup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FMup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FMup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FMup_unique$log2FoldChange.y ~ merged_ME_FMup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FMup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FMup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FMup$log2FoldChange.y ~ merged_shared_FMup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FMup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FMup$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FMup_unique$log2FoldChange.y ~ merged_EM_FMup_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FMup_unique$log2FoldChange.y ~ merged_ME_FMup_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FMup$log2FoldChange.y ~ merged_shared_FMup$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FM_rest_up$Row.names <- rownames(sign_E_D1_FM_rest_up)
sign_E_D1_FM_rest_up <- separate(data = as.data.frame(sign_E_D1_FM_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FM_rest_up[,c(2:5,7,8)] <- sapply(sign_E_D1_FM_rest_up[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FM_rest_up$Genome <- 'pan'

write.table(sign_E_D1_FM_rest_up[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FM_rest_up.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FM_rest_up$Row.names <- rownames(sign_M_D1_FM_rest_up)
sign_M_D1_FM_rest_up <- separate(data = as.data.frame(sign_M_D1_FM_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FM_rest_up[,c(2:5,7,8)] <- sapply(sign_M_D1_FM_rest_up[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FM_rest_up$Genome <- 'pan'

sign_M_D1_FM_rest_up_write <- sign_M_D1_FM_rest_up[,c(18,10,11)]
sign_M_D1_FM_rest_up_write <- subset(sign_M_D1_FM_rest_up_write, as.numeric(sign_M_D1_FM_rest_up_write$end2)-as.numeric(sign_M_D1_FM_rest_up_write$start2) < 10000)
write.table(sign_M_D1_FM_rest_up_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FM_rest_up.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FMup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FMup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FM_rest_up_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FM_rest_up_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FMup_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FM_rest_up)
nrow(sign_M_D1_FM_rest_up)
nrow(merged_shared_FMup)

#
nrow(sign_E_D1_FM_rest_up)-nrow(merged_ME_FMup[grep("NA", merged_ME_FMup$Row.names), ])
nrow(sign_M_D1_FM_rest_up)-nrow(merged_EM_FMup[grep("NA", merged_EM_FMup$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FM_rest_up_bed <- sign_E_D1_FM_rest_up[,c(18,8,9)]
colnames(sign_E_D1_FM_rest_up_bed) <- c('chrom','start','end')
sign_E_D1_FM_rest_up_bed$start <- as.numeric(sign_E_D1_FM_rest_up_bed$start)
sign_E_D1_FM_rest_up_bed$end <- as.numeric(sign_E_D1_FM_rest_up_bed$end)

sign_M_D1_FM_rest_up_bed <- sign_M_D1_FM_rest_up_write
colnames(sign_M_D1_FM_rest_up_bed) <- c('chrom','start','end')
sign_M_D1_FM_rest_up_bed$start <- as.numeric(sign_M_D1_FM_rest_up_bed$start)
sign_M_D1_FM_rest_up_bed$end <- as.numeric(sign_M_D1_FM_rest_up_bed$end)

merged_shared_shared_FMup_bed <- merged_shared_shared_FMup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_shared_FMup_bed) <- c('chrom','start','end')
merged_shared_shared_FMup_bed$start <- as.numeric(merged_shared_shared_FMup_bed$start)
merged_shared_shared_FMup_bed$end <- as.numeric(merged_shared_shared_FMup_bed$end)

sign_E_D1_FM_rest_up_bed_match <- bed_intersect(sign_E_D1_FM_rest_up_bed,conserved_E)
sign_M_D1_FM_rest_up_bed_match <- bed_intersect(sign_M_D1_FM_rest_up_bed,conserved_M)
merged_shared_shared_FMup_bed_match <- bed_intersect(merged_shared_shared_FMup_bed,conserved_S)


sign_E_D1_FM_rest_up_cons <- merge(sign_E_D1_FM_rest_up, sign_E_D1_FM_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FM_rest_up_cons <- merge(sign_M_D1_FM_rest_up, sign_M_D1_FM_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FMup_cons <- merge(merged_ME_FMup, sign_E_D1_FM_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FMup_cons <- merge(merged_EM_FMup, sign_M_D1_FM_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FMup_cons[grep("NA", merged_ME_FMup_cons$Row.names), ])
nrow(merged_EM_FMup_cons[grep("NA", merged_EM_FMup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FM_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FM_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_shared_FMup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FM_rest_up_bed_match)
nrow(sign_M_D1_FM_rest_up_bed_match)
nrow(merged_shared_shared_FMup_bed_match)

nrow(merged_ME_FMup_cons[grep("NA", merged_ME_FMup_cons$Row.names), ])
nrow(merged_EM_FMup_cons[grep("NA", merged_EM_FMup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))


# FM down 

merged_EM_FMdown <- merge(as.data.frame(res_E_FM_rest), as.data.frame(sign_M_D1_FM_rest_down), by = 'row.names')
merged_ME_FMdown <- merge(as.data.frame(res_M_FM_rest), as.data.frame(sign_E_D1_FM_rest_down), by = 'row.names')



merged_shared_FMdown <- merge(as.data.frame(sign_E_D1_FM_rest_down), as.data.frame(sign_M_D1_FM_rest_down), by = 'row.names')

# merged_EM_FMdown$log2FoldChange.x <- merged_EM_FMdown$log2FoldChange.x*-1
# merged_EM_FMdown$log2FoldChange.y <- merged_EM_FMdown$log2FoldChange.y*-1
# merged_ME_FMdown$log2FoldChange.x <- merged_ME_FMdown$log2FoldChange.x*-1
# merged_ME_FMdown$log2FoldChange.y <- merged_ME_FMdown$log2FoldChange.y*-1
# merged_shared_FMdown$log2FoldChange.x <- merged_shared_FMdown$log2FoldChange.x*-1
# merged_shared_FMdown$log2FoldChange.y <- merged_shared_FMdown$log2FoldChange.y*-1

merged_EM_FMdown_unique <- merged_EM_FMdown[!merged_EM_FMdown$Row.names %in% merged_shared_FMdown$Row.names, ]
merged_ME_FMdown_unique <- merged_ME_FMdown[!merged_ME_FMdown$Row.names %in% merged_shared_FMdown$Row.names, ]

plot(merged_EM_FMdown_unique$log2FoldChange.x, merged_EM_FMdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.8), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FM down")
par(new=T)
plot(merged_ME_FMdown_unique$log2FoldChange.x, merged_ME_FMdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.8), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FMdown$log2FoldChange.x, merged_shared_FMdown$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.8), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FMdown_unique$log2FoldChange.y ~ merged_EM_FMdown_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FMdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FMdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FMdown_unique$log2FoldChange.y ~ merged_ME_FMdown_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FMdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FMdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FMdown$log2FoldChange.y ~ merged_shared_FMdown$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FMdown$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FMdown$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FMdown_unique$log2FoldChange.y ~ merged_EM_FMdown_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FMdown_unique$log2FoldChange.y ~ merged_ME_FMdown_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FMdown$log2FoldChange.y ~ merged_shared_FMdown$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FM_rest_down$Row.names <- rownames(sign_E_D1_FM_rest_down)
sign_E_D1_FM_rest_down <- separate(data = as.data.frame(sign_E_D1_FM_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_E_D1_FM_rest_down[,c(2:5,7,8)] <- sapply(sign_E_D1_FM_rest_down[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FM_rest_down$Genome <- 'pan'

write.table(sign_E_D1_FM_rest_down[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FM_rest_down.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FM_rest_down$Row.names <- rownames(sign_M_D1_FM_rest_down)
sign_M_D1_FM_rest_down <- separate(data = as.data.frame(sign_M_D1_FM_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scafM", "startM", "endM"), sep = "_")
sign_M_D1_FM_rest_down[,c(2:5,7,8)] <- sapply(sign_M_D1_FM_rest_down[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FM_rest_down$Genome <- 'pan'

sign_M_D1_FM_rest_down_write <- sign_M_D1_FM_rest_down[,c(18,10,11)]
sign_M_D1_FM_rest_down_write <- subset(sign_M_D1_FM_rest_down_write, as.numeric(sign_M_D1_FM_rest_down_write$end2)-as.numeric(sign_M_D1_FM_rest_down_write$start2) < 10000)
write.table(sign_M_D1_FM_rest_down_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FM_rest_down.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FMdown[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FMdown.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FM_rest_down_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FM_rest_down_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FMdown_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FM_rest_down)
nrow(sign_M_D1_FM_rest_down)
nrow(merged_shared_FMdown)

#
nrow(sign_E_D1_FM_rest_down)-nrow(merged_ME_FMdown[grep("NA", merged_ME_FMdown$Row.names), ])
nrow(sign_M_D1_FM_rest_down)-nrow(merged_EM_FMdown[grep("NA", merged_EM_FMdown$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FM_rest_down_bed <- sign_E_D1_FM_rest_down[,c(18,8,9)]
colnames(sign_E_D1_FM_rest_down_bed) <- c('chrom','start','end')
sign_E_D1_FM_rest_down_bed$start <- as.numeric(sign_E_D1_FM_rest_down_bed$start)
sign_E_D1_FM_rest_down_bed$end <- as.numeric(sign_E_D1_FM_rest_down_bed$end)

sign_M_D1_FM_rest_down_bed <- sign_M_D1_FM_rest_down_write
colnames(sign_M_D1_FM_rest_down_bed) <- c('chrom','start','end')
sign_M_D1_FM_rest_down_bed$start <- as.numeric(sign_M_D1_FM_rest_down_bed$start)
sign_M_D1_FM_rest_down_bed$end <- as.numeric(sign_M_D1_FM_rest_down_bed$end)

merged_shared_shared_FMdown_bed <- merged_shared_shared_FMdown[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_shared_FMdown_bed) <- c('chrom','start','end')
merged_shared_shared_FMdown_bed$start <- as.numeric(merged_shared_shared_FMdown_bed$start)
merged_shared_shared_FMdown_bed$end <- as.numeric(merged_shared_shared_FMdown_bed$end)

sign_E_D1_FM_rest_down_bed_match <- bed_intersect(sign_E_D1_FM_rest_down_bed,conserved_E)
sign_M_D1_FM_rest_down_bed_match <- bed_intersect(sign_M_D1_FM_rest_down_bed,conserved_M)
merged_shared_shared_FMdown_bed_match <- bed_intersect(merged_shared_shared_FMdown_bed,conserved_S)


sign_E_D1_FM_rest_down_cons <- merge(sign_E_D1_FM_rest_down, sign_E_D1_FM_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FM_rest_down_cons <- merge(sign_M_D1_FM_rest_down, sign_M_D1_FM_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FMdown_cons <- merge(merged_ME_FMdown, sign_E_D1_FM_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FMdown_cons <- merge(merged_EM_FMdown, sign_M_D1_FM_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FMdown_cons[grep("NA", merged_ME_FMdown_cons$Row.names), ])
nrow(merged_EM_FMdown_cons[grep("NA", merged_EM_FMdown_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FM_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FM_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_shared_FMdown_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FM_rest_down_bed_match)
nrow(sign_M_D1_FM_rest_down_bed_match)
nrow(merged_shared_shared_FMdown_bed_match)

nrow(merged_ME_FMdown_cons[grep("NA", merged_ME_FMdown_cons$Row.names), ])
nrow(merged_EM_FMdown_cons[grep("NA", merged_EM_FMdown_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# FP up 

merged_EM_FPup <- merge(as.data.frame(res_E_FP_rest), as.data.frame(sign_M_D1_FP_rest_up), by = 'row.names')
merged_ME_FPup <- merge(as.data.frame(res_M_FP_rest), as.data.frame(sign_E_D1_FP_rest_up), by = 'row.names')



merged_shared_FPup <- merge(as.data.frame(sign_E_D1_FP_rest_up), as.data.frame(sign_M_D1_FP_rest_up), by = 'row.names')

merged_EM_FPup$log2FoldChange.x <- merged_EM_FPup$log2FoldChange.x*-1
merged_EM_FPup$log2FoldChange.y <- merged_EM_FPup$log2FoldChange.y*-1
merged_ME_FPup$log2FoldChange.x <- merged_ME_FPup$log2FoldChange.x*-1
merged_ME_FPup$log2FoldChange.y <- merged_ME_FPup$log2FoldChange.y*-1
merged_shared_FPup$log2FoldChange.x <- merged_shared_FPup$log2FoldChange.x*-1
merged_shared_FPup$log2FoldChange.y <- merged_shared_FPup$log2FoldChange.y*-1

merged_EM_FPup_unique <- merged_EM_FPup[!merged_EM_FPup$Row.names %in% merged_shared_FPup$Row.names, ]
merged_ME_FPup_unique <- merged_ME_FPup[!merged_ME_FPup$Row.names %in% merged_shared_FPup$Row.names, ]

plot(merged_EM_FPup_unique$log2FoldChange.x, merged_EM_FPup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.5), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FP up")
par(new=T)
plot(merged_ME_FPup_unique$log2FoldChange.x, merged_ME_FPup_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.5), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FPup$log2FoldChange.x, merged_shared_FPup$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.5), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FPup_unique$log2FoldChange.y ~ merged_EM_FPup_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FPup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FPup_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FPup_unique$log2FoldChange.y ~ merged_ME_FPup_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FPup_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FPup_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FPup$log2FoldChange.y ~ merged_shared_FPup$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FPup$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FPup$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FPup_unique$log2FoldChange.y ~ merged_EM_FPup_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FPup_unique$log2FoldChange.y ~ merged_ME_FPup_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FPup$log2FoldChange.y ~ merged_shared_FPup$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FP_rest_up$Row.names <- rownames(sign_E_D1_FP_rest_up)
sign_E_D1_FP_rest_up <- separate(data = as.data.frame(sign_E_D1_FP_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_E_D1_FP_rest_up[,c(2:5,7,8)] <- sapply(sign_E_D1_FP_rest_up[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FP_rest_up$Genome <- 'pan'

write.table(sign_E_D1_FP_rest_up[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FP_rest_up.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FP_rest_up$Row.names <- rownames(sign_M_D1_FP_rest_up)
sign_M_D1_FP_rest_up <- separate(data = as.data.frame(sign_M_D1_FP_rest_up), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_M_D1_FP_rest_up[,c(2:5,7,8)] <- sapply(sign_M_D1_FP_rest_up[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FP_rest_up$Genome <- 'pan'

sign_M_D1_FP_rest_up_write <- sign_M_D1_FP_rest_up[,c(18,10,11)]
sign_M_D1_FP_rest_up_write <- subset(sign_M_D1_FP_rest_up_write, as.numeric(sign_M_D1_FP_rest_up_write$end2)-as.numeric(sign_M_D1_FP_rest_up_write$start2) < 10000)
write.table(sign_M_D1_FP_rest_up_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FP_rest_up.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FPup[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FPup.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FP_rest_up_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FP_rest_up_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FPup_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FP_rest_up)
nrow(sign_M_D1_FP_rest_up)
nrow(merged_shared_FPup)

#
nrow(sign_E_D1_FP_rest_up)-nrow(merged_ME_FPup[grep("NA", merged_ME_FPup$Row.names), ])
nrow(sign_M_D1_FP_rest_up)-nrow(merged_EM_FPup[grep("NA", merged_EM_FPup$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FP_rest_up_bed <- sign_E_D1_FP_rest_up[,c(18,8,9)]
colnames(sign_E_D1_FP_rest_up_bed) <- c('chrom','start','end')
sign_E_D1_FP_rest_up_bed$start <- as.numeric(sign_E_D1_FP_rest_up_bed$start)
sign_E_D1_FP_rest_up_bed$end <- as.numeric(sign_E_D1_FP_rest_up_bed$end)

sign_M_D1_FP_rest_up_bed <- sign_M_D1_FP_rest_up_write
colnames(sign_M_D1_FP_rest_up_bed) <- c('chrom','start','end')
sign_M_D1_FP_rest_up_bed$start <- as.numeric(sign_M_D1_FP_rest_up_bed$start)
sign_M_D1_FP_rest_up_bed$end <- as.numeric(sign_M_D1_FP_rest_up_bed$end)

merged_shared_FPup_bed <- merged_shared_FPup[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_FPup_bed) <- c('chrom','start','end')
merged_shared_FPup_bed$start <- as.numeric(merged_shared_FPup_bed$start)
merged_shared_FPup_bed$end <- as.numeric(merged_shared_FPup_bed$end)

sign_E_D1_FP_rest_up_bed_match <- bed_intersect(sign_E_D1_FP_rest_up_bed,conserved_E)
sign_M_D1_FP_rest_up_bed_match <- bed_intersect(sign_M_D1_FP_rest_up_bed,conserved_M)
merged_shared_FPup_bed_match <- bed_intersect(merged_shared_FPup_bed,conserved_S)


sign_E_D1_FP_rest_up_cons <- merge(sign_E_D1_FP_rest_up, sign_E_D1_FP_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FP_rest_up_cons <- merge(sign_M_D1_FP_rest_up, sign_M_D1_FP_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FPup_cons <- merge(merged_ME_FPup, sign_E_D1_FP_rest_up_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FPup_cons <- merge(merged_EM_FPup, sign_M_D1_FP_rest_up_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FPup_cons[grep("NA", merged_ME_FPup_cons$Row.names), ])
nrow(merged_EM_FPup_cons[grep("NA", merged_EM_FPup_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FP_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FP_rest_up_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_FPup_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FP_rest_up_bed_match)
nrow(sign_M_D1_FP_rest_up_bed_match)
nrow(merged_shared_FPup_bed_match)

nrow(merged_ME_FPup_cons[grep("NA", merged_ME_FPup_cons$Row.names), ])
nrow(merged_EM_FPup_cons[grep("NA", merged_EM_FPup_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))

# FP down 

merged_EM_FPdown <- merge(as.data.frame(res_E_FP_rest), as.data.frame(sign_M_D1_FP_rest_down), by = 'row.names')
merged_ME_FPdown <- merge(as.data.frame(res_M_FP_rest), as.data.frame(sign_E_D1_FP_rest_down), by = 'row.names')



merged_shared_FPdown <- merge(as.data.frame(sign_E_D1_FP_rest_down), as.data.frame(sign_M_D1_FP_rest_down), by = 'row.names')

# merged_EM_FPdown$log2FoldChange.x <- merged_EM_FPdown$log2FoldChange.x*-1
# merged_EM_FPdown$log2FoldChange.y <- merged_EM_FPdown$log2FoldChange.y*-1
# merged_ME_FPdown$log2FoldChange.x <- merged_ME_FPdown$log2FoldChange.x*-1
# merged_ME_FPdown$log2FoldChange.y <- merged_ME_FPdown$log2FoldChange.y*-1
# merged_shared_FPdown$log2FoldChange.x <- merged_shared_FPdown$log2FoldChange.x*-1
# merged_shared_FPdown$log2FoldChange.y <- merged_shared_FPdown$log2FoldChange.y*-1

merged_EM_FPdown_unique <- merged_EM_FPdown[!merged_EM_FPdown$Row.names %in% merged_shared_FPdown$Row.names, ]
merged_ME_FPdown_unique <- merged_ME_FPdown[!merged_ME_FPdown$Row.names %in% merged_shared_FPdown$Row.names, ]

plot(merged_EM_FPdown_unique$log2FoldChange.x, merged_EM_FPdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('blue', alpha.f = 0.5), xlab = 'Foldchange species 2', ylab = ('Foldchange species 1 (significant)'), main = "36h FP down")
par(new=T)
plot(merged_ME_FPdown_unique$log2FoldChange.x, merged_ME_FPdown_unique$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('orange', alpha.f = 0.5), xlab = '', ylab = '')
par(new=T)
plot(merged_shared_FPdown$log2FoldChange.x, merged_shared_FPdown$log2FoldChange.y, pch = 19, xlim = c(-1,2), ylim = c(0,3), col = adjustcolor('black', alpha.f = 0.5), xlab = '', ylab = '')


ablineclip(lm(merged_EM_FPdown_unique$log2FoldChange.y ~ merged_EM_FPdown_unique$log2FoldChange.x), col = 'blue', x1 = min(na.omit(merged_EM_FPdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_EM_FPdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_ME_FPdown_unique$log2FoldChange.y ~ merged_ME_FPdown_unique$log2FoldChange.x), col = 'orange', x1 = min(na.omit(merged_ME_FPdown_unique$log2FoldChange.x)), x2 = max(na.omit(merged_ME_FPdown_unique$log2FoldChange.x)))
ablineclip(lm(merged_shared_FPdown$log2FoldChange.y ~ merged_shared_FPdown$log2FoldChange.x), col = 'black', x1 = min(na.omit(merged_shared_FPdown$log2FoldChange.x)), x2 = max(na.omit(merged_shared_FPdown$log2FoldChange.x)))

legend(-1, 3, legend=c(as.expression(bquote(paste('significant in ', italic('H. erato')))), 
                       as.expression(bquote(paste('significant in ', italic('H. melpomene')))),
                       "significant in both (shared)"),
       col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1, bty = "n")

# legend(3.5, 8, legend=c(paste('significant in ', 'H. erato', ' R2= ', round(summary(lm(merged_EM_FPdown_unique$log2FoldChange.y ~ merged_EM_FPdown_unique$log2FoldChange.x))$r.squared, 3)), 
#                         paste('significant in ', 'H. melpomene', ' R2= ', round(summary(lm(merged_ME_FPdown_unique$log2FoldChange.y ~ merged_ME_FPdown_unique$log2FoldChange.x))$r.squared, 3)),
#                         paste("significant in both (shared) R2= ", round(summary(lm(merged_shared_FPdown$log2FoldChange.y ~ merged_shared_FPdown$log2FoldChange.x))$r.squared, 3))),
#        col=c("orange", "blue", "black"), pch= c(19,19,19), cex=1)
######################## % completely missing

sign_E_D1_FP_rest_down$Row.names <- rownames(sign_E_D1_FP_rest_down)
sign_E_D1_FP_rest_down <- separate(data = as.data.frame(sign_E_D1_FP_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_E_D1_FP_rest_down[,c(2:5,7,8)] <- sapply(sign_E_D1_FP_rest_down[,c(2:5,7,8)],as.numeric) 
sign_E_D1_FP_rest_down$Genome <- 'pan'

write.table(sign_E_D1_FP_rest_down[,c(18,8,9)], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FP_rest_down.bed', col.names = F, row.names = F, quote = F)

sign_M_D1_FP_rest_down$Row.names <- rownames(sign_M_D1_FP_rest_down)
sign_M_D1_FP_rest_down <- separate(data = as.data.frame(sign_M_D1_FP_rest_down), col = 'Row.names', into = c("genome", "start", "end", "start2", "end2","scafE", "startE", "endE", "scaFP", "startM", "endM"), sep = "_")
sign_M_D1_FP_rest_down[,c(2:5,7,8)] <- sapply(sign_M_D1_FP_rest_down[,c(2:5,7,8)],as.numeric) 
sign_M_D1_FP_rest_down$Genome <- 'pan'

sign_M_D1_FP_rest_down_write <- sign_M_D1_FP_rest_down[,c(18,10,11)]
sign_M_D1_FP_rest_down_write <- subset(sign_M_D1_FP_rest_down_write, as.numeric(sign_M_D1_FP_rest_down_write$end2)-as.numeric(sign_M_D1_FP_rest_down_write$start2) < 10000)
write.table(sign_M_D1_FP_rest_down_write, '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FP_rest_down.bed', col.names = F, row.names = F, quote = F)

write.table(merged_shared_FPdown[,c('Genome.x','start.x','end.x')], '9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FPdown.bed', col.names = F, row.names = F, quote = F)

IDY_unique_E <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_E_D1_FP_rest_down_PAN_IDY.txt', h=T)
IDY_unique_M <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_M_D1_FP_rest_down_PAN_IDY.txt', h=T)
IDY_unique_S <- read.table('9_ATAC-seq_peak_conservation/IDY_shared_unique/sign_SHARED_FPdown_PAN_IDY.txt', h=T)

# Total peaks
nrow(sign_E_D1_FP_rest_down)
nrow(sign_M_D1_FP_rest_down)
nrow(merged_shared_FPdown)

#
nrow(sign_E_D1_FP_rest_down)-nrow(merged_ME_FPdown[grep("NA", merged_ME_FPdown$Row.names), ])
nrow(sign_M_D1_FP_rest_down)-nrow(merged_EM_FPdown[grep("NA", merged_EM_FPdown$Row.names), ])

nrow(subset(IDY_unique_E,IDY_unique_E$IDY == 0))
nrow(subset(IDY_unique_E,IDY_unique_E$IDY <= 0.5))

nrow(subset(IDY_unique_M,IDY_unique_M$IDY == 0))
nrow(subset(IDY_unique_M,IDY_unique_M$IDY <= 0.5))

nrow(subset(IDY_unique_S,IDY_unique_S$IDY == 0))
nrow(subset(IDY_unique_S,IDY_unique_S$IDY <= 0.5))

######################## Calculate for conserved set
sign_E_D1_FP_rest_down_bed <- sign_E_D1_FP_rest_down[,c(18,8,9)]
colnames(sign_E_D1_FP_rest_down_bed) <- c('chrom','start','end')
sign_E_D1_FP_rest_down_bed$start <- as.numeric(sign_E_D1_FP_rest_down_bed$start)
sign_E_D1_FP_rest_down_bed$end <- as.numeric(sign_E_D1_FP_rest_down_bed$end)

sign_M_D1_FP_rest_down_bed <- sign_M_D1_FP_rest_down_write
colnames(sign_M_D1_FP_rest_down_bed) <- c('chrom','start','end')
sign_M_D1_FP_rest_down_bed$start <- as.numeric(sign_M_D1_FP_rest_down_bed$start)
sign_M_D1_FP_rest_down_bed$end <- as.numeric(sign_M_D1_FP_rest_down_bed$end)

merged_shared_FPdown_bed <- merged_shared_FPdown[,c('Genome.x','start.x','end.x')]
colnames(merged_shared_FPdown_bed) <- c('chrom','start','end')
merged_shared_FPdown_bed$start <- as.numeric(merged_shared_FPdown_bed$start)
merged_shared_FPdown_bed$end <- as.numeric(merged_shared_FPdown_bed$end)

sign_E_D1_FP_rest_down_bed_match <- bed_intersect(sign_E_D1_FP_rest_down_bed,conserved_E)
sign_M_D1_FP_rest_down_bed_match <- bed_intersect(sign_M_D1_FP_rest_down_bed,conserved_M)
merged_shared_FPdown_bed_match <- bed_intersect(merged_shared_FPdown_bed,conserved_S)


sign_E_D1_FP_rest_down_cons <- merge(sign_E_D1_FP_rest_down, sign_E_D1_FP_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
sign_M_D1_FP_rest_down_cons <- merge(sign_M_D1_FP_rest_down, sign_M_D1_FP_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

merged_ME_FPdown_cons <- merge(merged_ME_FPdown, sign_E_D1_FP_rest_down_bed_match,  by.x = c('start'), by.y = c('start.x'))
merged_EM_FPdown_cons <- merge(merged_EM_FPdown, sign_M_D1_FP_rest_down_bed_match,  by.x = c('start2'), by.y = c('start.x'))

nrow(merged_ME_FPdown_cons[grep("NA", merged_ME_FPdown_cons$Row.names), ])
nrow(merged_EM_FPdown_cons[grep("NA", merged_EM_FPdown_cons$Row.names), ])


IDY_unique_E_cons <- merge(IDY_unique_E, sign_E_D1_FP_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_M_cons <- merge(IDY_unique_M, sign_M_D1_FP_rest_down_bed_match, by.x = c('start'), by.y = c('start.x'))
IDY_unique_S_cons <- merge(IDY_unique_S, merged_shared_FPdown_bed_match, by.x = c('start'), by.y = c('start.x'))

nrow(sign_E_D1_FP_rest_down_bed_match)
nrow(sign_M_D1_FP_rest_down_bed_match)
nrow(merged_shared_FPdown_bed_match)

nrow(merged_ME_FPdown_cons[grep("NA", merged_ME_FPdown_cons$Row.names), ])
nrow(merged_EM_FPdown_cons[grep("NA", merged_EM_FPdown_cons$Row.names), ])

nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY == 0))
nrow(subset(IDY_unique_E_cons,IDY_unique_E_cons$IDY <= 0.5))

nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY == 0))
nrow(subset(IDY_unique_M_cons,IDY_unique_M_cons$IDY <= 0.5))

nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY == 0))
nrow(subset(IDY_unique_S_cons,IDY_unique_S_cons$IDY <= 0.5))
