library(dplyr)
library(DESeq2)

### H. erato

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                                      (sampleInfo$day == 'larval'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(1,1,1,1))
barplot(colSums(counts_erato2), las=2)

dds_5th <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~wing)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == '1'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(1,1,1,1))
barplot(colSums(counts_erato2), las=2)

dds_D1 <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                  colData = sampleInfo_sub,
                                  design = ~wing)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == '2'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(1,1,1,1))
barplot(colSums(counts_erato2), las=2)

dds_D2 <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

###
# run model
###

rnaDDS_5th <- DESeq(dds_5th)
rnaDDS_D1 <- DESeq(dds_D1)
rnaDDS_D2 <- DESeq(dds_D2)

res <- results(rnaDDS)
res

write.table(res, file = "FWHW/FWHW_diff_erato.txt", quote = F, sep='\t')

subset(res, res$padj < 0.05 & abs(res$log2FoldChange) >= 1)
sign_upHW <- subset(res, res$padj < 0.1 & (res$log2FoldChange) > 1)
sign_upFW <- subset(res, res$padj < 0.1 & (res$log2FoldChange) < -1)


res["evm.model.Herato1202.713",]

####################################################################################

TFs <- read.table("FWHW/TF_IDs_FWHW.txt", h=T)

enrichment <- read.table("FWHW/motif_enrichment_summary.txt", h=T)

enrich_sub_FWup <- subset(enrichment, enrichment$comp == "FWHW_5th_p05FC1_erato_FWup")
enrich_sub_HWup <- subset(enrichment, enrichment$comp == "FWHW_5th_p05FC1_erato_HWup")

enrich_sub_FWup_agg <- aggregate(enrich_sub_FWup, by = list(enrich_sub_FWup$TF), FUN = mean)
enrich_sub_HWup_agg <- aggregate(enrich_sub_HWup, by = list(enrich_sub_HWup$TF), FUN = mean)

colnames(enrich_sub_FWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')
colnames(enrich_sub_HWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')

TF_merg_FWup <- merge(TFs, enrich_sub_FWup_agg, by = "TF", all.x = F)
TF_merg_HWup <- merge(TFs, enrich_sub_HWup_agg, by = "TF", all.x = F)


res_5th <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1 <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2 <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))

# res2 <- results(rnaDDS, contrast=c('day','larval','2'))

res_5th_sub <- res_5th[row.names(res_5th) %in% TFs$eratoID, ]
res_D1_sub <- res_D1[row.names(res_D1) %in% TFs$eratoID, ]
res_D2_sub <- res_D2[row.names(res_D2) %in% TFs$eratoID, ]

res_5th_sub$eratoID <- rownames(res_5th_sub)
res_D1_sub$eratoID <- rownames(res_D1_sub)
res_D2_sub$eratoID <- rownames(res_D2_sub)


res_merg_FWup <- merge(as.data.frame(res_5th_sub), TF_merg_FWup, by='eratoID')
res_merg_HWup <- merge(as.data.frame(res_5th_sub), TF_merg_HWup, by='eratoID')

res_merg_FWup[is.na(res_merg_FWup)] <- 0
res_merg_HWup[is.na(res_merg_HWup)] <- 0

plot(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, ylim=c(-8,2), xlim=c(-100,0), col='black', pch=19)
text(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, res_merg_FWup$TF, cex=0.6, pos=4, col="black")
par(new=T)
plot(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, ylim=c(-8,2), xlim=c(-100,0), col='red', pch=19)
text(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, res_merg_HWup$TF, cex=0.6, pos=4, col="red")


###
enrich_sub_FWup <- subset(enrichment, enrichment$comp == "FWHW_D1_p05FC1_erato_FWup")
enrich_sub_HWup <- subset(enrichment, enrichment$comp == "FWHW_D1_p05FC1_erato_HWup")

enrich_sub_FWup_agg <- aggregate(enrich_sub_FWup, by = list(enrich_sub_FWup$TF), FUN = mean)
enrich_sub_HWup_agg <- aggregate(enrich_sub_HWup, by = list(enrich_sub_HWup$TF), FUN = mean)

colnames(enrich_sub_FWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')
colnames(enrich_sub_HWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')

TF_merg_FWup <- merge(TFs, enrich_sub_FWup_agg, by = "TF", all.x = F)
TF_merg_HWup <- merge(TFs, enrich_sub_HWup_agg, by = "TF", all.x = F)


res_5th <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1 <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2 <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))

# res2 <- results(rnaDDS, contrast=c('day','larval','2'))

res_5th_sub <- res_5th[row.names(res_5th) %in% TFs$eratoID, ]
res_D1_sub <- res_D1[row.names(res_D1) %in% TFs$eratoID, ]
res_D2_sub <- res_D2[row.names(res_D2) %in% TFs$eratoID, ]

res_5th_sub$eratoID <- rownames(res_5th_sub)
res_D1_sub$eratoID <- rownames(res_D1_sub)
res_D2_sub$eratoID <- rownames(res_D2_sub)


res_merg_FWup <- merge(as.data.frame(res_D1_sub), TF_merg_FWup, by='eratoID')
res_merg_HWup <- merge(as.data.frame(res_D1_sub), TF_merg_HWup, by='eratoID')

res_merg_FWup[is.na(res_merg_FWup)] <- 0
res_merg_HWup[is.na(res_merg_HWup)] <- 0

plot(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, ylim=c(-8,8), xlim=c(-100,0), col='black', pch=19)
text(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, res_merg_FWup$TF, cex=0.6, pos=4, col="black")
par(new=T)
plot(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, ylim=c(-8,8), xlim=c(-100,0), col='red', pch=19)
text(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, res_merg_HWup$TF, cex=0.6, pos=4, col="red")
###



###
enrich_sub_FWup <- subset(enrichment, enrichment$comp == "FWHW_D2_p05FC1_erato_FWup")
enrich_sub_HWup <- subset(enrichment, enrichment$comp == "FWHW_D2_p05FC1_erato_HWup")

enrich_sub_FWup_agg <- aggregate(enrich_sub_FWup, by = list(enrich_sub_FWup$TF), FUN = mean)
enrich_sub_HWup_agg <- aggregate(enrich_sub_HWup, by = list(enrich_sub_HWup$TF), FUN = mean)

colnames(enrich_sub_FWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')
colnames(enrich_sub_HWup_agg) <- c( 'TF', 'comp', 'motif', 'evalue', 'flybase', 'tomtom', 'X')

TF_merg_FWup <- merge(TFs, enrich_sub_FWup_agg, by = "TF", all.x = F)
TF_merg_HWup <- merge(TFs, enrich_sub_HWup_agg, by = "TF", all.x = F)


res_5th <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1 <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2 <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))

# res2 <- results(rnaDDS, contrast=c('day','larval','2'))

res_5th_sub <- res_5th[row.names(res_5th) %in% TFs$eratoID, ]
res_D1_sub <- res_D1[row.names(res_D1) %in% TFs$eratoID, ]
res_D2_sub <- res_D2[row.names(res_D2) %in% TFs$eratoID, ]

res_5th_sub$eratoID <- rownames(res_5th_sub)
res_D1_sub$eratoID <- rownames(res_D1_sub)
res_D2_sub$eratoID <- rownames(res_D2_sub)


res_merg_FWup <- merge(as.data.frame(res_D2_sub), TF_merg_FWup, by='eratoID')
res_merg_HWup <- merge(as.data.frame(res_D2_sub), TF_merg_HWup, by='eratoID')

res_merg_FWup[is.na(res_merg_FWup)] <- 0
res_merg_HWup[is.na(res_merg_HWup)] <- 0

plot(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, ylim=c(-3,3), xlim=c(-100,0), col='black', pch=19)
text(log(res_merg_FWup$evalue), res_merg_FWup$log2FoldChange, res_merg_FWup$TF, cex=0.6, pos=4, col="black")
par(new=T)
plot(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, ylim=c(-3,3), xlim=c(-100,0), col='red', pch=19)
text(log(res_merg_HWup$evalue), res_merg_HWup$log2FoldChange, res_merg_HWup$TF, cex=0.6, pos=4, col="red")
###












res_5th_sub$day <- 0
res_D1_sub$day <- 1
res_D2_sub$day <- 2

res_merg1 <- rbind(as.data.frame(res_5th_sub), as.data.frame(res_D1_sub), as.data.frame(res_D2_sub))

res_merg2 <- merge(res_merg1, TFs, by='eratoID')

res_merg2[is.na(res_merg2)] <- 0

par(mfrow=c(30,1), mar=c(0.1,0,0.1,0))
for(e in 1:nrow(TFs)){
  subdata <- subset(res_merg2, res_merg2$eratoID == TFs$eratoID[e])
  
  sub0 <- subset(subdata, subdata$day == 0)
  sub1 <- subset(subdata, subdata$day == 1)
  sub2 <- subset(subdata, subdata$day == 2)
  # if(sub0$log2FoldChange == 0 & sub1$log2FoldChange == 0 & sub2$log2FoldChange == 0){next}
  plot(NULL, axes=FALSE, ann=FALSE, xlim=c(-1,3), ylim=c(-10,10))
  rect(-1,-10,10,0, col='gray65', border=NA)
  rect(-1,0,3,10, col='gray95', border=NA)
  par(new=T)
  subdata <- subdata[order(subdata$day),]
  plot(subdata$day, subdata$log2FoldChange, axes=FALSE, ann=FALSE, xlim=c(-1,3), ylim=c(-10,10), type='l', col=' red', lwd=2)
  text(-0.5,0, as.character(subdata$TF[1]))
}

write.table(res_merg2, file="FWHW/erato_DF_TFs.txt", quote = F, sep='\t')



######
# Volcano


par(mfrow=c(1,6), mar=c(2,2,2,0))
plot(res_5th$log2FoldChange, -log10(res_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.05))
res_merg2_5th <- subset(res_merg2, res_merg2$day == 0)
par(new=T)
# plot(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1))
text(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), labels=res_merg2_5th$TF, col='red')


plot(res_D1$log2FoldChange, -log10(res_D1$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.05))
res_merg2_5th <- subset(res_merg2, res_merg2$day == 1)
par(new=T)
# plot(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1))
text(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), labels=res_merg2_5th$TF, col='red')


plot(res_D2$log2FoldChange, -log10(res_D2$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.05))
res_merg2_5th <- subset(res_merg2, res_merg2$day == 2)
par(new=T)
# plot(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1))
text(res_merg2_5th$log2FoldChange, -log10(res_merg2_5th$pvalue), labels=res_merg2_5th$TF, col='red')



### H. melpomene

counts_melp <- read.table('counts_nonStranded/counts_table_melp_nonStranded.txt', h=T)

head(counts_melp)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == 'larval'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_melp2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_melp2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_melp2,
                              colData = sampleInfo_sub,
                              design = ~wing)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/FWHW_diff_melp.txt", quote = F, sep='\t')

##
counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == 'larval'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~wing)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/FWHW_diff_erato.txt", quote = F, sep='\t')

subset(res, res$padj < 0.05 & abs(res$log2FoldChange) >= 1)
sign_upHW <- subset(res, res$padj < 0.05 & (res$log2FoldChange) > 1)
sign_upFW <- subset(res, res$padj < 0.05 & (res$log2FoldChange) < -1)


res["evm.model.Hmelp1202.713",]




######################################################################################

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato'))

sampleInfo_sub$comp <- 'pupa'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == 'larval'){sampleInfo_sub$comp[e] <- '5thinstar' }}

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/5th_D1D2_diff_erato.txt", quote = F, sep='\t')


######################################################################################

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato'))

sampleInfo_sub$comp <- '5thD2'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == '1'){sampleInfo_sub$comp[e] <- 'D1' }
  }

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/D1_5thD2_diff_erato.txt", quote = F, sep='\t')


######################################################################################

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato'))

sampleInfo_sub$comp <- '5thD1'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == '2'){sampleInfo_sub$comp[e] <- 'D2' }
}

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/D2_5thD1_diff_erato.txt", quote = F, sep='\t')



######################################################################################

counts_erato <- read.table('counts_nonStranded/counts_table_melp_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene'))

sampleInfo_sub$comp <- 'pupa'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == 'larval'){sampleInfo_sub$comp[e] <- '5thinstar' }}

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/5th_D1D2_diff_melp.txt", quote = F, sep='\t')


######################################################################################

counts_erato <- read.table('counts_nonStranded/counts_table_melp_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene'))

sampleInfo_sub$comp <- '5thD2'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == '1'){sampleInfo_sub$comp[e] <- 'D1' }
}

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/D1_5thD2_diff_melp.txt", quote = F, sep='\t')


######################################################################################

counts_erato <- read.table('counts_table_melp_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene'))

sampleInfo_sub$comp <- '5thD1'
for(e in 1:nrow(sampleInfo_sub)){
  if(sampleInfo_sub$day[e] == '2'){sampleInfo_sub$comp[e] <- 'D2' }
}

sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

par(mar=c(12,6,4,4))
barplot(colSums(counts_erato2), las=2)

dds <- DESeqDataSetFromMatrix(countData = counts_erato2,
                              colData = sampleInfo_sub,
                              design = ~comp)

###
# run model
###
dds
rnaDDS <- DESeq(dds)
res <- results(rnaDDS)
res

write.table(res, file = "FWHW/D2_5thD1_diff_melp.txt", quote = F, sep='\t')

