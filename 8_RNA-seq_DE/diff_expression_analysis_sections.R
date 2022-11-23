library(dplyr)
library(DESeq2)

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)
counts_melp <- read.table('counts_nonStranded/counts_table_melp_nonStranded.txt', h=T)


sampleInfo <- read.table('counts_nonStranded/sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

sampleInfo$proximal <- 'rest'
sampleInfo$medial <- 'rest'
sampleInfo$distal <- 'rest'
for(e in 1:nrow(sampleInfo)){if(sampleInfo$compartment[e] == 'proximal'){sampleInfo$proximal[e] <- 'proximal'}}
for(e in 1:nrow(sampleInfo)){if(sampleInfo$compartment[e] == 'medial'){sampleInfo$medial[e] <- 'medial'}}
for(e in 1:nrow(sampleInfo)){if(sampleInfo$compartment[e] == 'distal'){sampleInfo$distal[e] <- 'distal'}}

sampleInfo_erato_D1 <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '1') & (sampleInfo$wing == 'forewing'))
sampleInfo_melp_D1 <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '1') & (sampleInfo$wing == 'forewing'))

sampleInfo_erato_D2 <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '2') & (sampleInfo$wing == 'forewing'))
sampleInfo_melp_D2 <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '2') & (sampleInfo$wing == 'forewing'))


sampleInfo_erato_D1$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_erato_D1)){
  if(sampleInfo_erato_D1$compartment[e] == 'proximal'){sampleInfo_erato_D1$tissueNumeric[e] <- 1}
  if(sampleInfo_erato_D1$compartment[e] == 'medial'){sampleInfo_erato_D1$tissueNumeric[e] <- 2}
  if(sampleInfo_erato_D1$compartment[e] == 'distal'){sampleInfo_erato_D1$tissueNumeric[e] <- 3}
}

sampleInfo_melp_D1$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_melp_D1)){
  if(sampleInfo_melp_D1$compartment[e] == 'proximal'){sampleInfo_melp_D1$tissueNumeric[e] <- 1}
  if(sampleInfo_melp_D1$compartment[e] == 'medial'){sampleInfo_melp_D1$tissueNumeric[e] <- 2}
  if(sampleInfo_melp_D1$compartment[e] == 'distal'){sampleInfo_melp_D1$tissueNumeric[e] <- 3}
}

sampleInfo_erato_D2$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_erato_D2)){
  if(sampleInfo_erato_D2$compartment[e] == 'proximal'){sampleInfo_erato_D2$tissueNumeric[e] <- 1}
  if(sampleInfo_erato_D2$compartment[e] == 'medial'){sampleInfo_erato_D2$tissueNumeric[e] <- 2}
  if(sampleInfo_erato_D2$compartment[e] == 'distal'){sampleInfo_erato_D2$tissueNumeric[e] <- 3}
}

sampleInfo_melp_D2$tissueNumeric <- NA
for(e in 1:nrow(sampleInfo_melp_D2)){
  if(sampleInfo_melp_D2$compartment[e] == 'proximal'){sampleInfo_melp_D2$tissueNumeric[e] <- 1}
  if(sampleInfo_melp_D2$compartment[e] == 'medial'){sampleInfo_melp_D2$tissueNumeric[e] <- 2}
  if(sampleInfo_melp_D2$compartment[e] == 'distal'){sampleInfo_melp_D2$tissueNumeric[e] <- 3}
}

counts_eratoF <- counts_erato[,-1]
counts_melpF <- counts_melp[,-1]

rownames(counts_eratoF) <- counts_erato[,1]
rownames(counts_melpF) <- counts_melp[,1]

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_erato_D1$sample_id)))
counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_melp_D1$sample_id)))

counts_erato2_D2 <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2$sample_id)))
counts_melp2_D2 <- select(counts_melpF, c(as.character(sampleInfo_melp_D2$sample_id)))

dds_D1_FD_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~distal)
dds_D1_FM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~medial)
dds_D1_FP_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~proximal)

dds_D1_FD_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D1, design = ~distal)
dds_D1_FM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D1, design = ~medial)
dds_D1_FP_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D1, design = ~proximal)

dds_D1_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~tissueNumeric)
dds_D1_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D1, design = ~tissueNumeric)

dds_D2_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_D2, colData = sampleInfo_erato_D2, design = ~tissueNumeric)
dds_D2_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_D2, colData = sampleInfo_melp_D2, design = ~tissueNumeric)


rnaDDS_E_grad <- DESeq(dds_D1_erato)
rnaDDS_M_grad <- DESeq(dds_D1_melp)

res_rna_E_grad <- results(rnaDDS_E_grad)
res_rna_M_grad <- results(rnaDDS_M_grad)

sign_E_rna_D1_grad_pos <- subset(res_rna_E_grad, res_rna_E_grad$padj < 0.05 & (res_rna_E_grad$log2FoldChange) > 1)
sign_E_rna_D1_grad_neg <- subset(res_rna_E_grad, res_rna_E_grad$padj < 0.05 & (res_rna_E_grad$log2FoldChange) < -1)

sign_M_rna_D1_grad_pos <- subset(res_rna_M_grad, res_rna_M_grad$padj < 0.05 & (res_rna_M_grad$log2FoldChange) > 1)
sign_M_rna_D1_grad_neg <- subset(res_rna_M_grad, res_rna_M_grad$padj < 0.05 & (res_rna_M_grad$log2FoldChange) < -1)

resultsNames(rnaDDS_E_grad)

as.data.frame(DESeq2::results(rnaDDS_E_grad, format = "DataFrame", name = "tissueNumeric"))

rnaDDS_E_grad_D2 <- DESeq(dds_D2_erato)
rnaDDS_M_grad_D2 <- DESeq(dds_D2_melp)

res_rna_E_grad_D2 <- results(rnaDDS_E_grad_D2)
res_rna_M_grad_D2 <- results(rnaDDS_M_grad_D2)

sign_E_rna_D2_grad_pos <- subset(res_rna_E_grad_D2, res_rna_E_grad_D2$padj < 0.05 & (res_rna_E_grad_D2$log2FoldChange) > 1)
sign_E_rna_D2_grad_neg <- subset(res_rna_E_grad_D2, res_rna_E_grad_D2$padj < 0.05 & (res_rna_E_grad_D2$log2FoldChange) < -1)

sign_M_rna_D2_grad_pos <- subset(res_rna_M_grad_D2, res_rna_M_grad_D2$padj < 0.05 & (res_rna_M_grad_D2$log2FoldChange) > 1)
sign_M_rna_D2_grad_neg <- subset(res_rna_M_grad_D2, res_rna_M_grad_D2$padj < 0.05 & (res_rna_M_grad_D2$log2FoldChange) < -1)


# extract normalized counts of sign genes


# erato D1 Pos
test_counts <- counts(rnaDDS_E_grad[rownames(sign_E_rna_D1_grad_pos),], normalized = TRUE)
erato_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_D1_pos_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_D1_pos_era <- aggregate(test_counts_m_D1_pos_era, list(test_counts_m_D1_pos_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_D1_pos_era <- scale(test_counts_m_m_D1_pos_era[,c(1:(ncol(test_counts_m_m_D1_pos_era)-1))])
test_counts_m_m_s_D1_pos_era <- cbind(as.data.frame(test_counts_m_m_s_D1_pos_era), test_counts_m_m_D1_pos_era$tissueNumeric)


# erato D1 Neg
test_counts <- counts(rnaDDS_E_grad[rownames(sign_E_rna_D1_grad_neg),], normalized = TRUE)
erato_grad_D1_neg_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_D1_neg_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_D1_neg_era <- aggregate(test_counts_m_D1_neg_era, list(test_counts_m_D1_neg_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_D1_neg_era <- scale(test_counts_m_m_D1_neg_era[,c(1:(ncol(test_counts_m_m_D1_neg_era)-1))])
test_counts_m_m_s_D1_neg_era <- cbind(as.data.frame(test_counts_m_m_s_D1_neg_era), test_counts_m_m_D1_neg_era$tissueNumeric)


# melp D1 Pos
test_counts <- counts(rnaDDS_M_grad[rownames(sign_M_rna_D1_grad_pos),], normalized = TRUE)
melp_grad_D1_pos_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_D1_pos_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_D1_pos_melp <- aggregate(test_counts_m_D1_pos_melp, list(test_counts_m_D1_pos_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_D1_pos_melp <- scale(test_counts_m_m_D1_pos_melp[,c(1:(ncol(test_counts_m_m_D1_pos_melp)-1))])
test_counts_m_m_s_D1_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_D1_pos_melp), test_counts_m_m_D1_pos_melp$tissueNumeric)




# melp D1 Neg
test_counts <- counts(rnaDDS_M_grad[rownames(sign_M_rna_D1_grad_neg),], normalized = TRUE)
melp_grad_D1_neg_genes <- rownames(test_counts)
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_D1_neg_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_D1_neg_melp <- aggregate(test_counts_m_D1_neg_melp, list(test_counts_m_D1_neg_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_D1_neg_melp <- scale(test_counts_m_m_D1_neg_melp[,c(1:(ncol(test_counts_m_m_D1_neg_melp)-1))])
test_counts_m_m_s_D1_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_D1_neg_melp), test_counts_m_m_D1_neg_melp$tissueNumeric)

## identidfy shared genes
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'era_gene')
head(homologs)

D1_pos_era <- as.data.frame(colnames(test_counts_m_m_s_D1_pos_era))
D1_pos_melp <- as.data.frame(colnames(test_counts_m_m_s_D1_pos_melp))

D1_neg_era <- as.data.frame(colnames(test_counts_m_m_s_D1_neg_era))
D1_neg_melp <- as.data.frame(colnames(test_counts_m_m_s_D1_neg_melp))

colnames(D1_pos_era) <- 'era_gene'
colnames(D1_pos_melp) <- 'melp_gene'

colnames(D1_neg_era) <- 'era_gene'
colnames(D1_neg_melp) <- 'melp_gene'

D1_pos_era_m <- merge(D1_pos_era, homologs, by = 'era_gene')
D1_pos_melp_m <- merge(D1_pos_melp, homologs, by = 'melp_gene')

D1_neg_era_m <- merge(D1_neg_era, homologs, by = 'era_gene')
D1_neg_melp_m <- merge(D1_neg_melp, homologs, by = 'melp_gene')

D1_pos_shared <- merge(D1_pos_era_m, D1_pos_melp_m, by = 'era_gene')
D1_neg_shared <- merge(D1_neg_era_m, D1_neg_melp_m, by = 'era_gene')
##

par(mfcol=c(3,2))

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_pos_era)-1)){
  lines(test_counts_m_m_s_D1_pos_era[,ncol(test_counts_m_m_s_D1_pos_era)],  test_counts_m_m_s_D1_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
##
plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_neg_era)-1)){
  lines(test_counts_m_m_s_D1_neg_era[,ncol(test_counts_m_m_s_D1_neg_era)],  test_counts_m_m_s_D1_neg_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_neg_era[,D1_neg_shared$era_gene,])-1)){
  lines(c(1,2,3),  test_counts_m_m_s_D1_neg_era[,D1_neg_shared$era_gene,][,i], col = adjustcolor("green",1), type = 'b')
}
##
plot(0,type='n',axes=FALSE,ann=FALSE)

##
plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_pos_melp)-1)){
  lines(test_counts_m_m_s_D1_pos_melp[,ncol(test_counts_m_m_s_D1_pos_melp)],  test_counts_m_m_s_D1_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}

##
plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_D1_neg_melp)-1)){
  lines(test_counts_m_m_s_D1_neg_melp[,ncol(test_counts_m_m_s_D1_neg_melp)],  test_counts_m_m_s_D1_neg_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_D1_neg_melp[,D1_neg_shared$melp_gene.y,])-1)){
  lines(c(1,2,3),  test_counts_m_m_s_D1_neg_melp[,D1_neg_shared$melp_gene.y,][,i], col = adjustcolor("green",1), type = 'b')
}
##
plot(0,type='n',axes=FALSE,ann=FALSE)

#### TFs

TF_list <- read.table('sections/TF_enriched_melp_era.txt', h=F)
head(TF_list)
colnames(TF_list) <- c('TF', 'melp', 'era')

TF_grad_era <- as.data.frame(res_rna_E_grad)[as.character(TF_list$era),]
TF_grad_melp <- as.data.frame(res_rna_M_grad)[as.character(TF_list$melp),]

write.table(TF_grad_era, file = 'sections/TF_grad_era_08122021.txt', quote = F, sep = '\t')
write.table(TF_grad_melp, file = 'sections/TF_grad_melp_08122021.txt', quote = F, sep = '\t')
####

as.data.frame(rownames(sign_E_rna_D1_grad_pos))
as.data.frame(rownames(sign_E_rna_D1_grad_neg))

as.data.frame(rownames(sign_M_rna_D1_grad_pos))
as.data.frame(rownames(sign_M_rna_D1_grad_neg))

homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')
head(homologs)

sign_E_rna_D1_grad_pos$erato_gene <- rownames(sign_E_rna_D1_grad_pos)
sign_E_rna_D1_grad_neg$erato_gene <- rownames(sign_E_rna_D1_grad_neg)
sign_M_rna_D1_grad_pos$melp_gene <- rownames(sign_M_rna_D1_grad_pos)
sign_M_rna_D1_grad_neg$melp_gene <- rownames(sign_M_rna_D1_grad_neg)

sign_E_rna_D2_grad_pos$erato_gene <- rownames(sign_E_rna_D2_grad_pos)
sign_E_rna_D2_grad_neg$erato_gene <- rownames(sign_E_rna_D2_grad_neg)

merge(as.data.frame(sign_E_rna_D1_grad_pos), as.data.frame(sign_E_rna_D2_grad_pos), by = 'erato_gene')
merge(as.data.frame(sign_E_rna_D1_grad_neg), as.data.frame(sign_E_rna_D2_grad_neg), by = 'erato_gene')

sign_E_rna_D1_grad_pos_homol <- merge(as.data.frame(sign_E_rna_D1_grad_pos), homologs, by = 'erato_gene')
sign_E_rna_D1_grad_neg_homol <- merge(as.data.frame(sign_E_rna_D1_grad_neg), homologs, by = 'erato_gene')
sign_M_rna_D1_grad_pos_homol <- merge(as.data.frame(sign_M_rna_D1_grad_pos), homologs, by = 'melp_gene')
sign_M_rna_D1_grad_neg_homol <- merge(as.data.frame(sign_M_rna_D1_grad_neg), homologs, by = 'melp_gene')


nrow(merge(sign_E_rna_D1_grad_pos_homol, sign_M_rna_D1_grad_pos_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_rna_D1_grad_neg_homol, sign_M_rna_D1_grad_neg_homol, by = 'melp_gene', all = F))



###
par(mfrow=c(3,2))
counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_erato_D1$sample_id)))
dds_D1_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~compartment)
rnaDDS_D1_erato <- DESeq(dds_D1_erato)

wntA_counts <- counts(rnaDDS_D1_erato['evm.model.Herato1001.161',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "WntA D1 erato", ylim = c(0,500), ylab = "RNA read count", xlab = "wing section")


counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_melp_D1$sample_id)))
dds_D1_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D1, design = ~compartment)
rnaDDS_D1_melp <- DESeq(dds_D1_melp)

wntA_counts <- counts(rnaDDS_D1_melp['HMEL018100g1.t1',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "WntA D1 melp", ylim = c(0,500), ylab = "RNA read count", xlab = "wing section")

##

wntA_counts <- counts(rnaDDS_D1_erato['evm.model.Herato1505.85',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "cortex D1 erato", ylim = c(0,150), ylab = "RNA read count", xlab = "wing section")



wntA_counts <- counts(rnaDDS_D1_melp['HMEL000025-RG',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "cortex D1 melp", ylim = c(0,150), ylab = "RNA read count", xlab = "wing section")

##


wntA_counts <- counts(rnaDDS_D1_erato['evm.model.Herato1801.64',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "optix D1 erato", ylim = c(0,10), ylab = "RNA read count", xlab = "wing section")



wntA_counts <- counts(rnaDDS_D1_melp['HMEL001028-RA',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D1, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "optix D1 melp", ylim = c(0,10), ylab = "RNA read count", xlab = "wing section")
###


###
par(mfrow=c(3,2))
counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2$sample_id)))
dds_D2_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D2, design = ~compartment)
rnaDDS_D2_erato <- DESeq(dds_D2_erato)

wntA_counts <- counts(rnaDDS_D2_erato['evm.model.Herato1001.161',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "wntA D2 erato", ylim = c(0,500), ylab = "RNA read count", xlab = "wing section")


counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_melp_D2$sample_id)))
dds_D2_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D2, design = ~compartment)
rnaDDS_D2_melp <- DESeq(dds_D2_melp)

wntA_counts <- counts(rnaDDS_D2_melp['HMEL018100g1.t1',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "wntA D2 melp", ylim = c(0,500), ylab = "RNA read count", xlab = "wing section")

##

wntA_counts <- counts(rnaDDS_D2_erato['evm.model.Herato1505.85',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "cortex D2 erato", ylim = c(0,150), ylab = "RNA read count", xlab = "wing section")



wntA_counts <- counts(rnaDDS_D2_melp['HMEL000025-RG',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "cortex D2 melp", ylim = c(0,150), ylab = "RNA read count", xlab = "wing section")

##


wntA_counts <- counts(rnaDDS_D2_erato['evm.model.Herato1801.64',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_erato_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~as.character(wntA_counts_merged$compartment), main = "optix D2 erato", ylim = c(0,10), ylab = "RNA read count", xlab = "wing section")



wntA_counts <- counts(rnaDDS_D2_melp['HMEL001028-RA',], normalized = TRUE)

wntA_counts <- as.data.frame(t(wntA_counts))
wntA_counts$sample_id <- rownames(wntA_counts)

wntA_counts_merged <- merge(wntA_counts, sampleInfo_melp_D2, by = 'sample_id')

boxplot(wntA_counts_merged[,2] ~ as.character(wntA_counts_merged$compartment), main = "optix D2 melp", ylim = c(0,10), ylab = "RNA read count", xlab = "wing section")
###







rnaDDS_D1_FD_erato <- DESeq(dds_D1_FD_erato)
rnaDDS_D1_FM_erato <- DESeq(dds_D1_FM_erato)
rnaDDS_D1_FP_erato <- DESeq(dds_D1_FP_erato)

rnaDDS_D1_FD_melp <- DESeq(dds_D1_FD_melp)
rnaDDS_D1_FM_melp <- DESeq(dds_D1_FM_melp)
rnaDDS_D1_FP_melp <- DESeq(dds_D1_FP_melp)

res_D1_FD_erato <- results(rnaDDS_D1_FD_erato, contrast=c('distal','distal','rest'))
res_D1_FM_erato <- results(rnaDDS_D1_FM_erato, contrast=c('medial','medial','rest'))
res_D1_FP_erato <- results(rnaDDS_D1_FP_erato, contrast=c('proximal','proximal','rest'))

res_D1_FD_melp <- results(rnaDDS_D1_FD_melp, contrast=c('distal','distal','rest'))
res_D1_FM_melp <- results(rnaDDS_D1_FM_melp, contrast=c('medial','medial','rest'))
res_D1_FP_melp <- results(rnaDDS_D1_FP_melp, contrast=c('proximal','proximal','rest'))

sign_D1_FD_erato <- subset(res_D1_FD_erato, res_D1_FD_erato$padj <= 0.05)
sign_D1_FM_erato <- subset(res_D1_FM_erato, res_D1_FM_erato$padj <= 0.05)
sign_D1_FP_erato <- subset(res_D1_FP_erato, res_D1_FP_erato$padj <= 0.05)

sign_D1_FD_melp <- subset(res_D1_FD_melp, res_D1_FD_melp$padj <= 0.05)
sign_D1_FM_melp <- subset(res_D1_FM_melp, res_D1_FM_melp$padj <= 0.05)
sign_D1_FP_melp <- subset(res_D1_FP_melp, res_D1_FP_melp$padj <= 0.05)

sign_D1_FDup_erato <- subset(res_D1_FD_erato, res_D1_FD_erato$padj <= 0.05 & res_D1_FD_erato$log2FoldChange > 1)
sign_D1_FMup_erato <- subset(res_D1_FM_erato, res_D1_FM_erato$padj <= 0.05 & res_D1_FM_erato$log2FoldChange > 1)
sign_D1_FPup_erato <- subset(res_D1_FP_erato, res_D1_FP_erato$padj <= 0.05 & res_D1_FP_erato$log2FoldChange > 1)

nrow(sign_D1_FDup_erato)
nrow(sign_D1_FMup_erato)
nrow(sign_D1_FPup_erato)

sign_D1_FDup_melp <- subset(res_D1_FD_melp, res_D1_FD_melp$padj <= 0.05 & res_D1_FD_melp$log2FoldChange > 1)
sign_D1_FMup_melp <- subset(res_D1_FM_melp, res_D1_FM_melp$padj <= 0.05 & res_D1_FM_melp$log2FoldChange > 1)
sign_D1_FPup_melp <- subset(res_D1_FP_melp, res_D1_FP_melp$padj <= 0.05 & res_D1_FP_melp$log2FoldChange > 1)

nrow(sign_D1_FDup_melp)
nrow(sign_D1_FMup_melp)
nrow(sign_D1_FPup_melp)

sign_D1_FDdown_erato <- subset(res_D1_FD_erato, res_D1_FD_erato$padj <= 0.05 & res_D1_FD_erato$log2FoldChange < -1)
sign_D1_FMdown_erato <- subset(res_D1_FM_erato, res_D1_FM_erato$padj <= 0.05 & res_D1_FM_erato$log2FoldChange < -1)
sign_D1_FPdown_erato <- subset(res_D1_FP_erato, res_D1_FP_erato$padj <= 0.05 & res_D1_FP_erato$log2FoldChange < -1)

nrow(sign_D1_FDdown_erato)
nrow(sign_D1_FMdown_erato)
nrow(sign_D1_FPdown_erato)

sign_D1_FDdown_melp <- subset(res_D1_FD_melp, res_D1_FD_melp$padj <= 0.05 & res_D1_FD_melp$log2FoldChange < -1)
sign_D1_FMdown_melp <- subset(res_D1_FM_melp, res_D1_FM_melp$padj <= 0.05 & res_D1_FM_melp$log2FoldChange < -1)
sign_D1_FPdown_melp <- subset(res_D1_FP_melp, res_D1_FP_melp$padj <= 0.05 & res_D1_FP_melp$log2FoldChange < -1)

nrow(sign_D1_FDdown_melp)
nrow(sign_D1_FMdown_melp)
nrow(sign_D1_FPdown_melp)



# erato D1 FP
test_counts <- counts(rnaDDS_D1_FP_erato[rownames(sign_D1_FPup_erato),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
erato_grad_D1_pos_genes 
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FP_pos_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FP_pos_era <- aggregate(test_counts_m_FP_pos_era, list(test_counts_m_FP_pos_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FP_pos_era <- scale(test_counts_m_m_FP_pos_era[,c(1:(ncol(test_counts_m_m_FP_pos_era)-1))])
test_counts_m_m_s_FP_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FP_pos_era), test_counts_m_m_FP_pos_era$tissueNumeric)


# erato D1 FP
test_counts <- counts(rnaDDS_D1_FP_erato[rownames(sign_D1_FPdown_erato),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FP_neg_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FP_neg_era <- aggregate(test_counts_m_FP_neg_era, list(test_counts_m_FP_neg_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FP_neg_era <- scale(test_counts_m_m_FP_neg_era[,c(1:(ncol(test_counts_m_m_FP_neg_era)-1))])
test_counts_m_m_s_FP_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FP_neg_era), test_counts_m_m_FP_neg_era$tissueNumeric)


# erato D1 FM
test_counts <- counts(rnaDDS_D1_FM_erato[rownames(sign_D1_FMup_erato),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FM_pos_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FM_pos_era <- aggregate(test_counts_m_FM_pos_era, list(test_counts_m_FM_pos_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FM_pos_era <- scale(test_counts_m_m_FM_pos_era[,c(1:(ncol(test_counts_m_m_FM_pos_era)-1))])
test_counts_m_m_s_FM_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FM_pos_era), test_counts_m_m_FM_pos_era$tissueNumeric)


# erato D1 FM
test_counts <- counts(rnaDDS_D1_FM_erato[rownames(sign_D1_FMdown_erato),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FM_neg_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FM_neg_era <- aggregate(test_counts_m_FM_neg_era, list(test_counts_m_FM_neg_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FM_neg_era <- scale(test_counts_m_m_FM_neg_era[,c(1:(ncol(test_counts_m_m_FM_neg_era)-1))])
test_counts_m_m_s_FM_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FM_neg_era), test_counts_m_m_FM_neg_era$tissueNumeric)


# erato D1 FD
test_counts <- counts(rnaDDS_D1_FD_erato[rownames(sign_D1_FDup_erato),], normalized = TRUE)
test_counts <- as.data.frame(test_counts)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FD_pos_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FD_pos_era <- aggregate(test_counts_m_FD_pos_era, list(test_counts_m_FD_pos_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FD_pos_era <- scale(test_counts_m_m_FD_pos_era[,c(1:(ncol(test_counts_m_m_FD_pos_era)-1))])
test_counts_m_m_s_FD_pos_era <- cbind(as.data.frame(test_counts_m_m_s_FD_pos_era), test_counts_m_m_FD_pos_era$tissueNumeric)


# erato D1 FD
test_counts <- counts(rnaDDS_D1_FD_erato[rownames(sign_D1_FDdown_erato),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%erato_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FD_neg_era <- merge(test_countst, sampleInfo_erato_D1, by = 'sample_id')

test_counts_m_m_FD_neg_era <- aggregate(test_counts_m_FD_neg_era, list(test_counts_m_FD_neg_era$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FD_neg_era <- scale(test_counts_m_m_FD_neg_era[,c(1:(ncol(test_counts_m_m_FD_neg_era)-1))])
test_counts_m_m_s_FD_neg_era <- cbind(as.data.frame(test_counts_m_m_s_FD_neg_era), test_counts_m_m_FD_neg_era$tissueNumeric)


# melp D1 FP
test_counts <- counts(rnaDDS_D1_FP_melp[rownames(sign_D1_FPup_melp),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FP_pos_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FP_pos_melp <- aggregate(test_counts_m_FP_pos_melp, list(test_counts_m_FP_pos_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FP_pos_melp <- scale(test_counts_m_m_FP_pos_melp[,c(1:(ncol(test_counts_m_m_FP_pos_melp)-1))])
test_counts_m_m_s_FP_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FP_pos_melp), test_counts_m_m_FP_pos_melp$tissueNumeric)


# melp D1 FP
test_counts <- counts(rnaDDS_D1_FP_melp[rownames(sign_D1_FPdown_melp),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FP_neg_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FP_neg_melp <- aggregate(test_counts_m_FP_neg_melp, list(test_counts_m_FP_neg_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FP_neg_melp <- scale(test_counts_m_m_FP_neg_melp[,c(1:(ncol(test_counts_m_m_FP_neg_melp)-1))])
test_counts_m_m_s_FP_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FP_neg_melp), test_counts_m_m_FP_neg_melp$tissueNumeric)





# melp D1 FM
test_counts <- counts(rnaDDS_D1_FM_melp[rownames(sign_D1_FMup_melp),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FM_pos_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FM_pos_melp <- aggregate(test_counts_m_FM_pos_melp, list(test_counts_m_FM_pos_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FM_pos_melp <- scale(test_counts_m_m_FM_pos_melp[,c(1:(ncol(test_counts_m_m_FM_pos_melp)-1))])
test_counts_m_m_s_FM_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FM_pos_melp), test_counts_m_m_FM_pos_melp$tissueNumeric)



# melp D1 FM
test_counts <- counts(rnaDDS_D1_FM_melp[rownames(sign_D1_FMdown_melp),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FM_neg_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FM_neg_melp <- aggregate(test_counts_m_FM_neg_melp, list(test_counts_m_FM_neg_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FM_neg_melp <- scale(test_counts_m_m_FM_neg_melp[,c(1:(ncol(test_counts_m_m_FM_neg_melp)-1))])
test_counts_m_m_s_FM_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FM_neg_melp), test_counts_m_m_FM_neg_melp$tissueNumeric)


# melp D1 FD
test_counts <- counts(rnaDDS_D1_FD_melp[rownames(sign_D1_FDup_melp),], normalized = TRUE)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FD_pos_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FD_pos_melp <- aggregate(test_counts_m_FD_pos_melp, list(test_counts_m_FD_pos_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FD_pos_melp <- scale(test_counts_m_m_FD_pos_melp[,c(1:(ncol(test_counts_m_m_FD_pos_melp)-1))])
test_counts_m_m_s_FD_pos_melp <- cbind(as.data.frame(test_counts_m_m_s_FD_pos_melp), test_counts_m_m_FD_pos_melp$tissueNumeric)


# melp D1 FD
test_counts <- counts(rnaDDS_D1_FD_melp[rownames(sign_D1_FDdown_melp),], normalized = TRUE)
test_counts <- as.data.frame(test_counts)
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_pos_genes,]
test_counts <- test_counts[!rownames(test_counts)%in%melp_grad_D1_neg_genes,]
test_countst <- as.data.frame(t(test_counts))
test_countst$sample_id <- rownames(test_countst)

test_counts_m_FD_neg_melp <- merge(test_countst, sampleInfo_melp_D1, by = 'sample_id')

test_counts_m_m_FD_neg_melp <- aggregate(test_counts_m_FD_neg_melp, list(test_counts_m_FD_neg_melp$tissueNumeric), FUN=mean)[,c(3:(2+nrow(test_counts)),(2+nrow(test_counts))+10)]

test_counts_m_m_s_FD_neg_melp <- scale(test_counts_m_m_FD_neg_melp[,c(1:(ncol(test_counts_m_m_FD_neg_melp)-1))])
test_counts_m_m_s_FD_neg_melp <- cbind(as.data.frame(test_counts_m_m_s_FD_neg_melp), test_counts_m_m_FD_neg_melp$tissueNumeric)


## identidfy shared genes
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'era_gene')
head(homologs)

FP_pos_era  <- as.data.frame(colnames(test_counts_m_m_s_FP_pos_era))
FP_pos_melp <- as.data.frame(colnames(test_counts_m_m_s_FP_pos_melp))

FP_neg_era  <- as.data.frame(colnames(test_counts_m_m_s_FP_neg_era))
FP_neg_melp <- as.data.frame(colnames(test_counts_m_m_s_FP_neg_melp))

FM_pos_era  <- as.data.frame(colnames(test_counts_m_m_s_FM_pos_era))
FM_pos_melp <- as.data.frame(colnames(test_counts_m_m_s_FM_pos_melp))

FM_neg_era  <- as.data.frame(colnames(test_counts_m_m_s_FM_neg_era))
FM_neg_melp <- as.data.frame(colnames(test_counts_m_m_s_FM_neg_melp))

FD_pos_era  <- as.data.frame(colnames(test_counts_m_m_s_FD_pos_era))
FD_pos_melp <- as.data.frame(colnames(test_counts_m_m_s_FD_pos_melp))

FD_neg_era  <- as.data.frame(colnames(test_counts_m_m_s_FD_neg_era))
FD_neg_melp <- as.data.frame(colnames(test_counts_m_m_s_FD_neg_melp))

colnames(FP_pos_era) <- 'era_gene'
colnames(FP_pos_melp) <- 'melp_gene'

colnames(FP_neg_era) <- 'era_gene'
colnames(FP_neg_melp) <- 'melp_gene'

colnames(FM_pos_era) <- 'era_gene'
colnames(FM_pos_melp) <- 'melp_gene'

colnames(FM_neg_era) <- 'era_gene'
colnames(FM_neg_melp) <- 'melp_gene'

colnames(FD_pos_era) <- 'era_gene'
colnames(FD_pos_melp) <- 'melp_gene'

colnames(FD_neg_era) <- 'era_gene'
colnames(FD_neg_melp) <- 'melp_gene'

FP_pos_era_m  <- merge(FP_pos_era, homologs, by = 'era_gene')
FP_pos_melp_m <- merge(FP_pos_melp, homologs, by = 'melp_gene')

FP_neg_era_m  <- merge(FP_neg_era, homologs, by = 'era_gene')
FP_neg_melp_m <- merge(FP_neg_melp, homologs, by = 'melp_gene')

FP_pos_shared <- merge(FP_pos_era_m, FP_pos_melp_m, by = 'era_gene')
FP_neg_shared <- merge(FP_neg_era_m, FP_neg_melp_m, by = 'era_gene')

FM_pos_era_m  <- merge(FM_pos_era, homologs, by = 'era_gene')
FM_pos_melp_m <- merge(FM_pos_melp, homologs, by = 'melp_gene')

FM_neg_era_m  <- merge(FM_neg_era, homologs, by = 'era_gene')
FM_neg_melp_m <- merge(FM_neg_melp, homologs, by = 'melp_gene')

FM_pos_shared <- merge(FM_pos_era_m, FM_pos_melp_m, by = 'era_gene')
FM_neg_shared <- merge(FM_neg_era_m, FM_neg_melp_m, by = 'era_gene')

FD_pos_era_m  <- merge(FD_pos_era, homologs, by = 'era_gene')
FD_pos_melp_m <- merge(FD_pos_melp, homologs, by = 'melp_gene')

FD_neg_era_m  <- merge(FD_neg_era, homologs, by = 'era_gene')
FD_neg_melp_m <- merge(FD_neg_melp, homologs, by = 'melp_gene')

FD_pos_shared <- merge(FD_pos_era_m, FD_pos_melp_m, by = 'era_gene')
FD_neg_shared <- merge(FD_neg_era_m, FD_neg_melp_m, by = 'era_gene')
###

par(mfcol=c(3,2))

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FP_pos_era)-1)){
  lines(test_counts_m_m_s_FP_pos_era[,ncol(test_counts_m_m_s_FP_pos_era)],  test_counts_m_m_s_FP_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_pos_era)-1)){
  lines(c(1,2,3),  test_counts_m_m_s_FP_pos_era[,FP_pos_shared$era_gene][,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_era)-1)){
  lines(test_counts_m_m_s_FP_neg_era[,ncol(test_counts_m_m_s_FP_neg_era)],  test_counts_m_m_s_FP_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FM_pos_era)-1)){
  lines(test_counts_m_m_s_FM_pos_era[,ncol(test_counts_m_m_s_FM_pos_era)],  test_counts_m_m_s_FM_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_era)-1)){
  lines(test_counts_m_m_s_FM_neg_era[,ncol(test_counts_m_m_s_FM_neg_era)],  test_counts_m_m_s_FM_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FD_pos_era)-1)){
  lines(test_counts_m_m_s_FD_pos_era[,ncol(test_counts_m_m_s_FD_pos_era)],  test_counts_m_m_s_FD_pos_era[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_era)-1)){
  lines(test_counts_m_m_s_FD_neg_era[,ncol(test_counts_m_m_s_FD_neg_era)],  test_counts_m_m_s_FD_neg_era[,i], col = adjustcolor("darkred",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FP_pos_melp)-1)){
  lines(test_counts_m_m_s_FP_pos_melp[,ncol(test_counts_m_m_s_FP_pos_melp)],  test_counts_m_m_s_FP_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_pos_melp)-1)){
  lines(c(1,2,3),  test_counts_m_m_s_FP_pos_melp[,FP_pos_shared$melp_gene.y][,i], col = adjustcolor("green",1), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FP_neg_melp)-1)){
  lines(test_counts_m_m_s_FP_neg_melp[,ncol(test_counts_m_m_s_FP_neg_melp)],  test_counts_m_m_s_FP_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FM_pos_melp)-1)){
  lines(test_counts_m_m_s_FM_pos_melp[,ncol(test_counts_m_m_s_FM_pos_melp)],  test_counts_m_m_s_FM_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FM_neg_melp)-1)){
  lines(test_counts_m_m_s_FM_neg_melp[,ncol(test_counts_m_m_s_FM_neg_melp)],  test_counts_m_m_s_FM_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}

plot(0,0,xlim = c(1,3), ylim = c(-1.5,1.5),type = "n", ylab= "", xlab = "")

for (i in 1:(ncol(test_counts_m_m_s_FD_pos_melp)-1)){
  lines(test_counts_m_m_s_FD_pos_melp[,ncol(test_counts_m_m_s_FD_pos_melp)],  test_counts_m_m_s_FD_pos_melp[,i], col = adjustcolor("black",0.3), type = 'b')
}
for (i in 1:(ncol(test_counts_m_m_s_FD_neg_melp)-1)){
  lines(test_counts_m_m_s_FD_neg_melp[,ncol(test_counts_m_m_s_FD_neg_melp)],  test_counts_m_m_s_FD_neg_melp[,i], col = adjustcolor("darkred",0.3), type = 'b')
}


#### TFs

TF_list <- read.table('sections/TF_enriched_melp_era.txt', h=F)
head(TF_list)
colnames(TF_list) <- c('TF', 'melp', 'era')

TF_FP_era <- as.data.frame(res_D1_FP_erato)[as.character(TF_list$era),]
TF_FM_era <- as.data.frame(res_D1_FM_erato)[as.character(TF_list$era),]
TF_FD_era <- as.data.frame(res_D1_FD_erato)[as.character(TF_list$era),]

TF_FP_melp <- as.data.frame(res_D1_FP_melp)[as.character(TF_list$melp),]
TF_FM_melp <- as.data.frame(res_D1_FM_melp)[as.character(TF_list$melp),]
TF_FD_melp <- as.data.frame(res_D1_FD_melp)[as.character(TF_list$melp),]

write.table(TF_FP_era, file = 'sections/TF_FP_era_08122021.txt', quote = F, sep = '\t')
write.table(TF_FM_era, file = 'sections/TF_FM_era_08122021.txt', quote = F, sep = '\t')
write.table(TF_FD_era, file = 'sections/TF_FD_era_08122021.txt', quote = F, sep = '\t')

write.table(TF_FP_melp, file = 'sections/TF_FP_melp_08122021.txt', quote = F, sep = '\t')
write.table(TF_FM_melp, file = 'sections/TF_FM_melp_08122021.txt', quote = F, sep = '\t')
write.table(TF_FD_melp, file = 'sections/TF_FD_melp_08122021.txt', quote = F, sep = '\t')
###

par(mfrow=c(2,3))

plot(res_D1_FP_erato$log2FoldChange, -log10(res_D1_FP_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FP_erato$log2FoldChange, -log10(sign_D1_FP_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')

plot(res_D1_FM_erato$log2FoldChange, -log10(res_D1_FM_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FM_erato$log2FoldChange, -log10(sign_D1_FM_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')

plot(res_D1_FD_erato$log2FoldChange, -log10(res_D1_FD_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FD_erato$log2FoldChange, -log10(sign_D1_FD_erato$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')


plot(res_D1_FP_melp$log2FoldChange, -log10(res_D1_FP_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FP_melp$log2FoldChange, -log10(sign_D1_FP_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')

plot(res_D1_FM_melp$log2FoldChange, -log10(res_D1_FM_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FM_melp$log2FoldChange, -log10(sign_D1_FM_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')

plot(res_D1_FD_melp$log2FoldChange, -log10(res_D1_FD_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F, xlab='',ylab='')
par(new=T)
plot(sign_D1_FD_melp$log2FoldChange, -log10(sign_D1_FD_melp$pvalue), pch=19, xlim=c(-25,25), ylim=c(0,15), col=adjustcolor('black', alpha=1), labels = F, xlab='',ylab='')



res_D1_FP_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]
res_D1_FM_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]
res_D1_FD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]


dds_D1_comp_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D1, design = ~compartment)
rnaDDS_D1_comp_erato <- DESeq(dds_D1_comp_erato)
res_D1_PvsD_erato <- results(rnaDDS_D1_comp_erato, contrast=c('compartment','proximal','distal'))
res_D1_MvsD_erato <- results(rnaDDS_D1_comp_erato, contrast=c('compartment','medial','distal'))
res_D1_MvsP_erato <- results(rnaDDS_D1_comp_erato, contrast=c('compartment','medial','proximal'))

res_D1_PvsD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]
res_D1_MvsD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]
res_D1_MvsP_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57"),]

############ shared DE expressed with shared peaks
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')
head(homologs)

sign_D1_FDup_erato$erato_gene <- rownames(sign_D1_FDup_erato)
sign_D1_FMup_erato$erato_gene <- rownames(sign_D1_FMup_erato)
sign_D1_FPup_erato$erato_gene <- rownames(sign_D1_FPup_erato)

sign_D1_FDup_erato_homol <- merge(as.data.frame(sign_D1_FDup_erato), homologs, by = 'erato_gene')
sign_D1_FMup_erato_homol <- merge(as.data.frame(sign_D1_FMup_erato), homologs, by = 'erato_gene')
sign_D1_FPup_erato_homol <- merge(as.data.frame(sign_D1_FPup_erato), homologs, by = 'erato_gene')

sign_D1_FDup_melp$melp_gene <- rownames(sign_D1_FDup_melp)
sign_D1_FMup_melp$melp_gene <- rownames(sign_D1_FMup_melp)
sign_D1_FPup_melp$melp_gene <- rownames(sign_D1_FPup_melp)

sign_D1_FDup_melp_homol <- merge(as.data.frame(sign_D1_FDup_melp), homologs, by = 'melp_gene')
sign_D1_FMup_melp_homol <- merge(as.data.frame(sign_D1_FMup_melp), homologs, by = 'melp_gene')
sign_D1_FPup_melp_homol <- merge(as.data.frame(sign_D1_FPup_melp), homologs, by = 'melp_gene')

nrow(merge(sign_D1_FDup_erato_homol, sign_D1_FDup_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D1_FMup_erato_homol, sign_D1_FMup_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D1_FPup_erato_homol, sign_D1_FPup_melp_homol, by = 'melp_gene', all = F))


sign_D1_FDdown_erato$erato_gene <- rownames(sign_D1_FDdown_erato)
sign_D1_FMdown_erato$erato_gene <- rownames(sign_D1_FMdown_erato)
sign_D1_FPdown_erato$erato_gene <- rownames(sign_D1_FPdown_erato)

sign_D1_FDdown_erato_homol <- merge(as.data.frame(sign_D1_FDdown_erato), homologs, by = 'erato_gene')
sign_D1_FMdown_erato_homol <- merge(as.data.frame(sign_D1_FMdown_erato), homologs, by = 'erato_gene')
sign_D1_FPdown_erato_homol <- merge(as.data.frame(sign_D1_FPdown_erato), homologs, by = 'erato_gene')

sign_D1_FDdown_melp$melp_gene <- rownames(sign_D1_FDdown_melp)
sign_D1_FMdown_melp$melp_gene <- rownames(sign_D1_FMdown_melp)
sign_D1_FPdown_melp$melp_gene <- rownames(sign_D1_FPdown_melp)

sign_D1_FDdown_melp_homol <- merge(as.data.frame(sign_D1_FDdown_melp), homologs, by = 'melp_gene')
sign_D1_FMdown_melp_homol <- merge(as.data.frame(sign_D1_FMdown_melp), homologs, by = 'melp_gene')
sign_D1_FPdown_melp_homol <- merge(as.data.frame(sign_D1_FPdown_melp), homologs, by = 'melp_gene')

nrow(merge(sign_D1_FDdown_erato_homol, sign_D1_FDdown_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D1_FMdown_erato_homol, sign_D1_FMdown_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D1_FPdown_erato_homol, sign_D1_FPdown_melp_homol, by = 'melp_gene', all = F))






### D2

sampleInfo_erato_D2 <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '2') & (sampleInfo$wing == 'forewing'))
sampleInfo_melp_D2 <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '2') & (sampleInfo$wing == 'forewing'))

counts_eratoF <- counts_erato[,-1]
counts_melpF <- counts_melp[,-1]

rownames(counts_eratoF) <- counts_erato[,1]
rownames(counts_melpF) <- counts_melp[,1]

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2$sample_id)))
counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_melp_D2$sample_id)))

dds_D2_FD_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D2, design = ~distal)
dds_D2_FM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D2, design = ~medial)
dds_D2_FP_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D2, design = ~proximal)

dds_D2_FD_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D2, design = ~distal)
dds_D2_FM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D2, design = ~medial)
dds_D2_FP_melp <- DESeqDataSetFromMatrix(countData = counts_melp2, colData = sampleInfo_melp_D2, design = ~proximal)

rnaDDS_D2_FD_erato <- DESeq(dds_D2_FD_erato)
rnaDDS_D2_FM_erato <- DESeq(dds_D2_FM_erato)
rnaDDS_D2_FP_erato <- DESeq(dds_D2_FP_erato)

rnaDDS_D2_FD_melp <- DESeq(dds_D2_FD_melp)
rnaDDS_D2_FM_melp <- DESeq(dds_D2_FM_melp)
rnaDDS_D2_FP_melp <- DESeq(dds_D2_FP_melp)

res_D2_FD_erato <- results(rnaDDS_D2_FD_erato, contrast=c('distal','distal','rest'))
res_D2_FM_erato <- results(rnaDDS_D2_FM_erato, contrast=c('medial','medial','rest'))
res_D2_FP_erato <- results(rnaDDS_D2_FP_erato, contrast=c('proximal','proximal','rest'))

res_D2_FD_melp <- results(rnaDDS_D2_FD_melp, contrast=c('distal','distal','rest'))
res_D2_FM_melp <- results(rnaDDS_D2_FM_melp, contrast=c('medial','medial','rest'))
res_D2_FP_melp <- results(rnaDDS_D2_FP_melp, contrast=c('proximal','proximal','rest'))

sign_D2_FD_erato <- subset(res_D2_FD_erato, res_D2_FD_erato$padj <= 0.05)
sign_D2_FM_erato <- subset(res_D2_FM_erato, res_D2_FM_erato$padj <= 0.05)
sign_D2_FP_erato <- subset(res_D2_FP_erato, res_D2_FP_erato$padj <= 0.05)

sign_D2_FD_melp <- subset(res_D2_FD_melp, res_D2_FD_melp$padj <= 0.05)
sign_D2_FM_melp <- subset(res_D2_FM_melp, res_D2_FM_melp$padj <= 0.05)
sign_D2_FP_melp <- subset(res_D2_FP_melp, res_D2_FP_melp$padj <= 0.05)

sign_D2_FDup_erato <- subset(res_D2_FD_erato, res_D2_FD_erato$padj <= 0.05 & res_D2_FD_erato$log2FoldChange > 1)
sign_D2_FMup_erato <- subset(res_D2_FM_erato, res_D2_FM_erato$padj <= 0.05 & res_D2_FM_erato$log2FoldChange > 1)
sign_D2_FPup_erato <- subset(res_D2_FP_erato, res_D2_FP_erato$padj <= 0.05 & res_D2_FP_erato$log2FoldChange > 1)

nrow(sign_D2_FDup_erato)
nrow(sign_D2_FMup_erato)
nrow(sign_D2_FPup_erato)

sign_D2_FDup_melp <- subset(res_D2_FD_melp, res_D2_FD_melp$padj <= 0.05 & res_D2_FD_melp$log2FoldChange > 1)
sign_D2_FMup_melp <- subset(res_D2_FM_melp, res_D2_FM_melp$padj <= 0.05 & res_D2_FM_melp$log2FoldChange > 1)
sign_D2_FPup_melp <- subset(res_D2_FP_melp, res_D2_FP_melp$padj <= 0.05 & res_D2_FP_melp$log2FoldChange > 1)

nrow(sign_D2_FDup_melp)
nrow(sign_D2_FMup_melp)
nrow(sign_D2_FPup_melp)

sign_D2_FDdown_erato <- subset(res_D2_FD_erato, res_D2_FD_erato$padj <= 0.05 & res_D2_FD_erato$log2FoldChange < -1)
sign_D2_FMdown_erato <- subset(res_D2_FM_erato, res_D2_FM_erato$padj <= 0.05 & res_D2_FM_erato$log2FoldChange < -1)
sign_D2_FPdown_erato <- subset(res_D2_FP_erato, res_D2_FP_erato$padj <= 0.05 & res_D2_FP_erato$log2FoldChange < -1)

nrow(sign_D2_FDdown_erato)
nrow(sign_D2_FMdown_erato)
nrow(sign_D2_FPdown_erato)

sign_D2_FDdown_melp <- subset(res_D2_FD_melp, res_D2_FD_melp$padj <= 0.05 & res_D2_FD_melp$log2FoldChange < -1)
sign_D2_FMdown_melp <- subset(res_D2_FM_melp, res_D2_FM_melp$padj <= 0.05 & res_D2_FM_melp$log2FoldChange < -1)
sign_D2_FPdown_melp <- subset(res_D2_FP_melp, res_D2_FP_melp$padj <= 0.05 & res_D2_FP_melp$log2FoldChange < -1)

nrow(sign_D2_FDdown_melp)
nrow(sign_D2_FMdown_melp)
nrow(sign_D2_FPdown_melp)

res_D2_FP_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234"),]
res_D2_FM_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234"),]
res_D2_FD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234"),]


dds_D2_comp_erato <- DESeqDataSetFromMatrix(countData = counts_erato2, colData = sampleInfo_erato_D2, design = ~compartment)
rnaDDS_D2_comp_erato <- DESeq(dds_D2_comp_erato)
res_D2_PvsD_erato <- results(rnaDDS_D2_comp_erato, contrast=c('compartment','proximal','distal'))
res_D2_MvsD_erato <- results(rnaDDS_D2_comp_erato, contrast=c('compartment','medial','distal'))
res_D2_MvsP_erato <- results(rnaDDS_D2_comp_erato, contrast=c('compartment','medial','proximal'))

res_D2_PvsD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234","evm.model.Herato0701.757","evm.model.Herato1701.314","evm.model.Herato1904.162","evm.model.Herato1904.21"),]
res_D2_MvsD_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234","evm.model.Herato0701.757","evm.model.Herato1701.314","evm.model.Herato1904.162","evm.model.Herato1904.21"),]
res_D2_MvsP_erato[c("evm.model.Herato1904.247","evm.model.Herato0901.271","evm.model.Herato0101.514","evm.model.Herato2101.308","evm.model.Herato1901.57","evm.model.Herato0503.234","evm.model.Herato0701.757","evm.model.Herato1701.314","evm.model.Herato1904.162","evm.model.Herato1904.21"),]




############ shared DE expressed with shared peaks
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')
head(homologs)

sign_D2_FDup_erato$erato_gene <- rownames(sign_D2_FDup_erato)
sign_D2_FMup_erato$erato_gene <- rownames(sign_D2_FMup_erato)
sign_D2_FPup_erato$erato_gene <- rownames(sign_D2_FPup_erato)

sign_D2_FDup_erato_homol <- merge(as.data.frame(sign_D2_FDup_erato), homologs, by = 'erato_gene')
sign_D2_FMup_erato_homol <- merge(as.data.frame(sign_D2_FMup_erato), homologs, by = 'erato_gene')
sign_D2_FPup_erato_homol <- merge(as.data.frame(sign_D2_FPup_erato), homologs, by = 'erato_gene')

sign_D2_FDup_melp$melp_gene <- rownames(sign_D2_FDup_melp)
sign_D2_FMup_melp$melp_gene <- rownames(sign_D2_FMup_melp)
sign_D2_FPup_melp$melp_gene <- rownames(sign_D2_FPup_melp)

sign_D2_FDup_melp_homol <- merge(as.data.frame(sign_D2_FDup_melp), homologs, by = 'melp_gene')
sign_D2_FMup_melp_homol <- merge(as.data.frame(sign_D2_FMup_melp), homologs, by = 'melp_gene')
sign_D2_FPup_melp_homol <- merge(as.data.frame(sign_D2_FPup_melp), homologs, by = 'melp_gene')

nrow(merge(sign_D2_FDup_erato_homol, sign_D2_FDup_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D2_FMup_erato_homol, sign_D2_FMup_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D2_FPup_erato_homol, sign_D2_FPup_melp_homol, by = 'melp_gene', all = F))


sign_D2_FDdown_erato$erato_gene <- rownames(sign_D2_FDdown_erato)
sign_D2_FMdown_erato$erato_gene <- rownames(sign_D2_FMdown_erato)
sign_D2_FPdown_erato$erato_gene <- rownames(sign_D2_FPdown_erato)

sign_D2_FDdown_erato_homol <- merge(as.data.frame(sign_D2_FDdown_erato), homologs, by = 'erato_gene')
sign_D2_FMdown_erato_homol <- merge(as.data.frame(sign_D2_FMdown_erato), homologs, by = 'erato_gene')
sign_D2_FPdown_erato_homol <- merge(as.data.frame(sign_D2_FPdown_erato), homologs, by = 'erato_gene')

sign_D2_FDdown_melp$melp_gene <- rownames(sign_D2_FDdown_melp)
sign_D2_FMdown_melp$melp_gene <- rownames(sign_D2_FMdown_melp)
sign_D2_FPdown_melp$melp_gene <- rownames(sign_D2_FPdown_melp)

sign_D2_FDdown_melp_homol <- merge(as.data.frame(sign_D2_FDdown_melp), homologs, by = 'melp_gene')
sign_D2_FMdown_melp_homol <- merge(as.data.frame(sign_D2_FMdown_melp), homologs, by = 'melp_gene')
sign_D2_FPdown_melp_homol <- merge(as.data.frame(sign_D2_FPdown_melp), homologs, by = 'melp_gene')

nrow(merge(sign_D2_FDdown_erato_homol, sign_D2_FDdown_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D2_FMdown_erato_homol, sign_D2_FMdown_melp_homol, by = 'melp_gene', all = F))
nrow(merge(sign_D2_FPdown_erato_homol, sign_D2_FPdown_melp_homol, by = 'melp_gene', all = F))







##### seperate model
sampleInfo_erato_D1_FPFM <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'medial')))
sampleInfo_melp_D1_FPFM <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'medial')))

sampleInfo_erato_D1_FPFD <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'distal')))
sampleInfo_melp_D1_FPFD <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'distal')))

sampleInfo_erato_D1_FDFM <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'distal') | (sampleInfo$compartment == 'medial')))
sampleInfo_melp_D1_FDFM <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '1') & ((sampleInfo$compartment == 'distal') | (sampleInfo$compartment == 'medial')))

counts_erato2_FPFM <- select(counts_eratoF, c(as.character(sampleInfo_erato_D1_FPFM$sample_id)))
counts_melp2_FPFM <- select(counts_melpF, c(as.character(sampleInfo_melp_D1_FPFM$sample_id)))

counts_erato2_FPFD <- select(counts_eratoF, c(as.character(sampleInfo_erato_D1_FPFD$sample_id)))
counts_melp2_FPFD <- select(counts_melpF, c(as.character(sampleInfo_melp_D1_FPFD$sample_id)))

counts_erato2_FDFM <- select(counts_eratoF, c(as.character(sampleInfo_erato_D1_FDFM$sample_id)))
counts_melp2_FDFM <- select(counts_melpF, c(as.character(sampleInfo_melp_D1_FDFM$sample_id)))

dds_D1_FPFM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FPFM, colData = sampleInfo_erato_D1_FPFM, design = ~compartment)
dds_D1_FPFM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FPFM, colData = sampleInfo_melp_D1_FPFM, design = ~compartment)

dds_D1_FPFD_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FPFD, colData = sampleInfo_erato_D1_FPFD, design = ~compartment)
dds_D1_FPFD_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FPFD, colData = sampleInfo_melp_D1_FPFD, design = ~compartment)

dds_D1_FDFM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FDFM, colData = sampleInfo_erato_D1_FDFM, design = ~compartment)
dds_D1_FDFM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FDFM, colData = sampleInfo_melp_D1_FDFM, design = ~compartment)


rnaDDS_D1_FPFM_erato <- DESeq(dds_D1_FPFM_erato)
rnaDDS_D1_FPFM_melp <- DESeq(dds_D1_FPFM_melp)

rnaDDS_D1_FPFD_erato <- DESeq(dds_D1_FPFD_erato)
rnaDDS_D1_FPFD_melp <- DESeq(dds_D1_FPFD_melp)

rnaDDS_D1_FDFM_erato <- DESeq(dds_D1_FDFM_erato)
rnaDDS_D1_FDFM_melp <- DESeq(dds_D1_FDFM_melp)

res_D1_FPFM_erato <- results(rnaDDS_D1_FPFM_erato, contrast=c('compartment','proximal','medial'))
res_D1_FPFM_melp <- results(rnaDDS_D1_FPFM_melp, contrast=c('compartment','proximal','medial'))

res_D1_FPFD_erato <- results(rnaDDS_D1_FPFD_erato, contrast=c('compartment','proximal','distal'))
res_D1_FPFD_melp <- results(rnaDDS_D1_FPFD_melp, contrast=c('compartment','proximal','distal'))

res_D1_FDFM_erato <- results(rnaDDS_D1_FDFM_erato, contrast=c('compartment','distal','medial'))
res_D1_FDFM_melp <- results(rnaDDS_D1_FDFM_melp, contrast=c('compartment','distal','medial'))

sign_D1_FPFM_FPup_erato <- subset(res_D1_FPFM_erato, res_D1_FPFM_erato$padj <= 0.05 & res_D1_FPFM_erato$log2FoldChange > 1)
sign_D1_FPFM_FPup_melp  <- subset(res_D1_FPFM_melp, res_D1_FPFM_melp$padj <= 0.05 & res_D1_FPFM_melp$log2FoldChange > 1)

sign_D1_FPFD_FPup_erato <- subset(res_D1_FPFD_erato, res_D1_FPFD_erato$padj <= 0.05 & res_D1_FPFD_erato$log2FoldChange > 1)
sign_D1_FPFD_FPup_melp  <- subset(res_D1_FPFD_melp, res_D1_FPFD_melp$padj <= 0.05 & res_D1_FPFD_melp$log2FoldChange > 1)

sign_D1_FDFM_FDup_erato <- subset(res_D1_FDFM_erato, res_D1_FDFM_erato$padj <= 0.05 & res_D1_FDFM_erato$log2FoldChange > 1)
sign_D1_FDFM_FDup_melp  <- subset(res_D1_FDFM_melp, res_D1_FDFM_melp$padj <= 0.05 & res_D1_FDFM_melp$log2FoldChange > 1)


sign_D1_FPFM_FMup_erato <- subset(res_D1_FPFM_erato, res_D1_FPFM_erato$padj <= 0.05 & res_D1_FPFM_erato$log2FoldChange < -1)
sign_D1_FPFM_FMup_melp  <- subset(res_D1_FPFM_melp, res_D1_FPFM_melp$padj <= 0.05 & res_D1_FPFM_melp$log2FoldChange < -1)

sign_D1_FPFD_FDup_erato <- subset(res_D1_FPFD_erato, res_D1_FPFD_erato$padj <= 0.05 & res_D1_FPFD_erato$log2FoldChange < -1)
sign_D1_FPFD_FDup_melp  <- subset(res_D1_FPFD_melp, res_D1_FPFD_melp$padj <= 0.05 & res_D1_FPFD_melp$log2FoldChange < -1)

sign_D1_FDFM_FMup_erato <- subset(res_D1_FDFM_erato, res_D1_FDFM_erato$padj <= 0.05 & res_D1_FDFM_erato$log2FoldChange < -1)
sign_D1_FDFM_FMup_melp  <- subset(res_D1_FDFM_melp, res_D1_FDFM_melp$padj <= 0.05 & res_D1_FDFM_melp$log2FoldChange < -1)



sign_E_D1_FPup <- rbind(sign_D1_FPFM_FPup_erato, sign_D1_FPFD_FPup_erato)
sign_E_D1_FPup <- sign_E_D1_FPup[unique(c(rownames(sign_D1_FPFM_FPup_erato),rownames(sign_D1_FPFD_FPup_erato))),]

sign_E_D1_FPdown <- rbind(sign_D1_FPFM_FMup_erato, sign_D1_FPFD_FDup_erato)
sign_E_D1_FPdown <- sign_E_D1_FPdown[unique(c(rownames(sign_D1_FPFM_FMup_erato),rownames(sign_D1_FPFD_FDup_erato))),]

sign_E_D1_FDup_intersect <- merge(as.data.frame(sign_D1_FPFM_FPup_erato), as.data.frame(sign_D1_FPFD_FPup_erato), by = "row.names")
sign_E_D1_FDdown_intersect <- merge(as.data.frame(sign_D1_FPFM_FMup_erato), as.data.frame(sign_D1_FPFD_FDup_erato), by = "row.names")

nrow(sign_E_D1_FPup)
nrow(sign_E_D1_FPdown)
nrow(sign_E_D1_FDup_intersect)
nrow(sign_E_D1_FDdown_intersect)

sign_M_D1_FPup <- rbind(sign_D1_FPFM_FPup_melp, sign_D1_FPFD_FPup_melp)
sign_M_D1_FPup <- sign_M_D1_FPup[unique(c(rownames(sign_D1_FPFM_FPup_melp),rownames(sign_D1_FPFD_FPup_melp))),]

sign_M_D1_FPdown <- rbind(sign_D1_FPFM_FMup_melp, sign_D1_FPFD_FDup_melp)
sign_M_D1_FPdown <- sign_M_D1_FPdown[unique(c(rownames(sign_D1_FPFM_FMup_melp),rownames(sign_D1_FPFD_FDup_melp))),]

sign_M_D1_FDup_intersect <- merge(as.data.frame(sign_D1_FPFM_FPup_melp), as.data.frame(sign_D1_FPFD_FPup_melp), by = "row.names")
sign_M_D1_FDdown_intersect <- merge(as.data.frame(sign_D1_FPFM_FMup_melp), as.data.frame(sign_D1_FPFD_FDup_melp), by = "row.names")

nrow(sign_M_D1_FPup)
nrow(sign_M_D1_FPdown)
nrow(sign_M_D1_FDup_intersect)
nrow(sign_M_D1_FDdown_intersect)


sign_E_D1_FMup <- rbind(sign_D1_FPFM_FMup_erato, sign_D1_FDFM_FMup_erato)
# sign_E_D1_FMup <- sign_E_D1_FPup[unique(c(rownames(sign_D1_FPFM_FMup_erato),rownames(sign_D1_FDFM_FMup_erato))),]

sign_E_D1_FMdown <- rbind(sign_D1_FPFM_FPup_erato, sign_D1_FDFM_FDup_erato)
# sign_E_D1_FMdown <- sign_E_D1_FPdown[unique(c(rownames(sign_D1_FPFM_FPup_erato),rownames(sign_D1_FDFM_FDup_erato))),]

sign_E_D1_FMup_intersect <- merge(as.data.frame(sign_D1_FPFM_FMup_erato), as.data.frame(sign_D1_FDFM_FMup_erato), by = "row.names")
sign_E_D1_FMdown_intersect <- merge(as.data.frame(sign_D1_FPFM_FPup_erato), as.data.frame(sign_D1_FDFM_FDup_erato), by = "row.names")

nrow(sign_E_D1_FMup)
nrow(sign_E_D1_FMdown)
nrow(sign_E_D1_FMup_intersect)
nrow(sign_E_D1_FMdown_intersect)


sign_M_D1_FMup <- rbind(sign_D1_FPFM_FMup_melp, sign_D1_FDFM_FMup_melp)
# sign_M_D1_FMup <- sign_E_D1_FPup[unique(c(rownames(sign_D1_FPFM_FMup_melp),rownames(sign_D1_FDFM_FMup_melp))),]

sign_M_D1_FMdown <- rbind(sign_D1_FPFM_FPup_melp, sign_D1_FDFM_FDup_melp)
# sign_M_D1_FMdown <- sign_E_D1_FPdown[unique(c(rownames(sign_D1_FPFM_FPup_melp),rownames(sign_D1_FDFM_FDup_melp))),]

sign_M_D1_FMup_intersect <- merge(as.data.frame(sign_D1_FPFM_FMup_melp), as.data.frame(sign_D1_FDFM_FMup_melp), by = "row.names")
sign_M_D1_FMdown_intersect <- merge(as.data.frame(sign_D1_FPFM_FPup_melp), as.data.frame(sign_D1_FDFM_FDup_melp), by = "row.names")

nrow(sign_M_D1_FMup)
nrow(sign_M_D1_FMdown)
nrow(sign_M_D1_FMup_intersect)
nrow(sign_M_D1_FMdown_intersect)


sign_E_D1_FDup <- rbind(sign_D1_FPFD_FDup_erato, sign_D1_FDFM_FDup_erato)
# sign_E_D1_FDup <- sign_E_D1_FPup[unique(c(rownames(sign_D1_FPFD_FDup_erato),rownames(sign_D1_FDFM_FDup_erato))),]

sign_E_D1_FDdown <- rbind(sign_D1_FPFD_FPup_erato, sign_D1_FDFM_FMup_erato)
# sign_E_D1_FDdown <- sign_E_D1_FPdown[unique(c(rownames(sign_D1_FPFD_FPup_erato),rownames(sign_D1_FDFM_FMup_erato))),]

sign_E_D1_FDup_intersect <- merge(as.data.frame(sign_D1_FPFD_FDup_erato), as.data.frame(sign_D1_FDFM_FDup_erato), by = "row.names")
sign_E_D1_FDdown_intersect <- merge(as.data.frame(sign_D1_FPFD_FPup_erato), as.data.frame(sign_D1_FDFM_FMup_erato), by = "row.names")

nrow(sign_E_D1_FDup)
nrow(sign_E_D1_FDdown)
nrow(sign_E_D1_FDup_intersect)
nrow(sign_E_D1_FDdown_intersect)



sign_M_D1_FDup <- rbind(sign_D1_FPFD_FDup_melp, sign_D1_FDFM_FDup_melp)
# sign_M_D1_FDup <- sign_E_D1_FPup[unique(c(rownames(sign_D1_FPFD_FDup_melp),rownames(sign_D1_FDFM_FDup_melp))),]

sign_M_D1_FDdown <- rbind(sign_D1_FPFD_FPup_melp, sign_D1_FDFM_FMup_melp)
# sign_M_D1_FDdown <- sign_E_D1_FPdown[unique(c(rownames(sign_D1_FPFD_FPup_melp),rownames(sign_D1_FDFM_FMup_melp))),]

sign_M_D1_FDup_intersect <- merge(as.data.frame(sign_D1_FPFD_FDup_melp), as.data.frame(sign_D1_FDFM_FDup_melp), by = "row.names")
sign_M_D1_FDdown_intersect <- merge(as.data.frame(sign_D1_FPFD_FPup_melp), as.data.frame(sign_D1_FDFM_FMup_melp), by = "row.names")

nrow(sign_M_D1_FDup)
nrow(sign_M_D1_FDdown)
nrow(sign_M_D1_FDup_intersect)
nrow(sign_M_D1_FDdown_intersect)



############ shared DE expressed with shared peaks
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')
head(homologs)

sign_E_D1_FDup$erato_gene <- rownames(sign_E_D1_FDup)
sign_E_D1_FMup$erato_gene <- rownames(sign_E_D1_FMup)
sign_E_D1_FPup$erato_gene <- rownames(sign_E_D1_FPup)

sign_E_D1_FDup_homol <- merge(as.data.frame(sign_E_D1_FDup), homologs, by = 'erato_gene')
sign_E_D1_FMup_homol <- merge(as.data.frame(sign_E_D1_FMup), homologs, by = 'erato_gene')
sign_E_D1_FPup_homol <- merge(as.data.frame(sign_E_D1_FPup), homologs, by = 'erato_gene')

sign_M_D1_FDup$melp_gene <- rownames(sign_M_D1_FDup)
sign_M_D1_FMup$melp_gene <- rownames(sign_M_D1_FMup)
sign_M_D1_FPup$melp_gene <- rownames(sign_M_D1_FPup)

sign_M_D1_FDup_homol <- merge(as.data.frame(sign_M_D1_FDup), homologs, by = 'melp_gene')
sign_M_D1_FMup_homol <- merge(as.data.frame(sign_M_D1_FMup), homologs, by = 'melp_gene')
sign_M_D1_FPup_homol <- merge(as.data.frame(sign_M_D1_FPup), homologs, by = 'melp_gene')

nrow(merge(sign_E_D1_FDup_homol, sign_M_D1_FDup_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D1_FMup_homol, sign_M_D1_FMup_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D1_FPup_homol, sign_M_D1_FPup_homol, by = 'melp_gene', all = F))


sign_E_D1_FDdown$erato_gene <- rownames(sign_E_D1_FDdown)
sign_E_D1_FMdown$erato_gene <- rownames(sign_E_D1_FMdown)
sign_E_D1_FPdown$erato_gene <- rownames(sign_E_D1_FPdown)

sign_E_D1_FDdown_homol <- merge(as.data.frame(sign_E_D1_FDdown), homologs, by = 'erato_gene')
sign_E_D1_FMdown_homol <- merge(as.data.frame(sign_E_D1_FMdown), homologs, by = 'erato_gene')
sign_E_D1_FPdown_homol <- merge(as.data.frame(sign_E_D1_FPdown), homologs, by = 'erato_gene')

sign_M_D1_FDdown$melp_gene <- rownames(sign_M_D1_FDdown)
sign_M_D1_FMdown$melp_gene <- rownames(sign_M_D1_FMdown)
sign_M_D1_FPdown$melp_gene <- rownames(sign_M_D1_FPdown)

sign_M_D1_FDdown_homol <- merge(as.data.frame(sign_M_D1_FDdown), homologs, by = 'melp_gene')
sign_M_D1_FMdown_homol <- merge(as.data.frame(sign_M_D1_FMdown), homologs, by = 'melp_gene')
sign_M_D1_FPdown_homol <- merge(as.data.frame(sign_M_D1_FPdown), homologs, by = 'melp_gene')

nrow(merge(sign_E_D1_FDdown_homol, sign_M_D1_FDdown_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D1_FMdown_homol, sign_M_D1_FMdown_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D1_FPdown_homol, sign_M_D1_FPdown_homol, by = 'melp_gene', all = F))



##### seperate model
sampleInfo_erato_D2_FPFM <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'medial')))
sampleInfo_melp_D2_FPFM <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'medial')))

sampleInfo_erato_D2_FPFD <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'distal')))
sampleInfo_melp_D2_FPFD <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'proximal') | (sampleInfo$compartment == 'distal')))

sampleInfo_erato_D2_FDFM <- subset(sampleInfo,  (sampleInfo$species == 'erato') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'distal') | (sampleInfo$compartment == 'medial')))
sampleInfo_melp_D2_FDFM <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & (sampleInfo$day == '2') & ((sampleInfo$compartment == 'distal') | (sampleInfo$compartment == 'medial')))

counts_erato2_FPFM <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2_FPFM$sample_id)))
counts_melp2_FPFM <- select(counts_melpF, c(as.character(sampleInfo_melp_D2_FPFM$sample_id)))

counts_erato2_FPFD <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2_FPFD$sample_id)))
counts_melp2_FPFD <- select(counts_melpF, c(as.character(sampleInfo_melp_D2_FPFD$sample_id)))

counts_erato2_FDFM <- select(counts_eratoF, c(as.character(sampleInfo_erato_D2_FDFM$sample_id)))
counts_melp2_FDFM <- select(counts_melpF, c(as.character(sampleInfo_melp_D2_FDFM$sample_id)))

dds_D2_FPFM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FPFM, colData = sampleInfo_erato_D2_FPFM, design = ~compartment)
dds_D2_FPFM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FPFM, colData = sampleInfo_melp_D2_FPFM, design = ~compartment)

dds_D2_FPFD_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FPFD, colData = sampleInfo_erato_D2_FPFD, design = ~compartment)
dds_D2_FPFD_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FPFD, colData = sampleInfo_melp_D2_FPFD, design = ~compartment)

dds_D2_FDFM_erato <- DESeqDataSetFromMatrix(countData = counts_erato2_FDFM, colData = sampleInfo_erato_D2_FDFM, design = ~compartment)
dds_D2_FDFM_melp <- DESeqDataSetFromMatrix(countData = counts_melp2_FDFM, colData = sampleInfo_melp_D2_FDFM, design = ~compartment)


rnaDDS_D2_FPFM_erato <- DESeq(dds_D2_FPFM_erato)
rnaDDS_D2_FPFM_melp <- DESeq(dds_D2_FPFM_melp)

rnaDDS_D2_FPFD_erato <- DESeq(dds_D2_FPFD_erato)
rnaDDS_D2_FPFD_melp <- DESeq(dds_D2_FPFD_melp)

rnaDDS_D2_FDFM_erato <- DESeq(dds_D2_FDFM_erato)
rnaDDS_D2_FDFM_melp <- DESeq(dds_D2_FDFM_melp)

res_D2_FPFM_erato <- results(rnaDDS_D2_FPFM_erato, contrast=c('compartment','proximal','medial'))
res_D2_FPFM_melp <- results(rnaDDS_D2_FPFM_melp, contrast=c('compartment','proximal','medial'))

res_D2_FPFD_erato <- results(rnaDDS_D2_FPFD_erato, contrast=c('compartment','proximal','distal'))
res_D2_FPFD_melp <- results(rnaDDS_D2_FPFD_melp, contrast=c('compartment','proximal','distal'))

res_D2_FDFM_erato <- results(rnaDDS_D2_FDFM_erato, contrast=c('compartment','distal','medial'))
res_D2_FDFM_melp <- results(rnaDDS_D2_FDFM_melp, contrast=c('compartment','distal','medial'))

sign_D2_FPFM_FPup_erato <- subset(res_D2_FPFM_erato, res_D2_FPFM_erato$padj <= 0.05 & res_D2_FPFM_erato$log2FoldChange > 1)
sign_D2_FPFM_FPup_melp  <- subset(res_D2_FPFM_melp, res_D2_FPFM_melp$padj <= 0.05 & res_D2_FPFM_melp$log2FoldChange > 1)

sign_D2_FPFD_FPup_erato <- subset(res_D2_FPFD_erato, res_D2_FPFD_erato$padj <= 0.05 & res_D2_FPFD_erato$log2FoldChange > 1)
sign_D2_FPFD_FPup_melp  <- subset(res_D2_FPFD_melp, res_D2_FPFD_melp$padj <= 0.05 & res_D2_FPFD_melp$log2FoldChange > 1)

sign_D2_FDFM_FDup_erato <- subset(res_D2_FDFM_erato, res_D2_FDFM_erato$padj <= 0.05 & res_D2_FDFM_erato$log2FoldChange > 1)
sign_D2_FDFM_FDup_melp  <- subset(res_D2_FDFM_melp, res_D2_FDFM_melp$padj <= 0.05 & res_D2_FDFM_melp$log2FoldChange > 1)


sign_D2_FPFM_FMup_erato <- subset(res_D2_FPFM_erato, res_D2_FPFM_erato$padj <= 0.05 & res_D2_FPFM_erato$log2FoldChange < -1)
sign_D2_FPFM_FMup_melp  <- subset(res_D2_FPFM_melp, res_D2_FPFM_melp$padj <= 0.05 & res_D2_FPFM_melp$log2FoldChange < -1)

sign_D2_FPFD_FDup_erato <- subset(res_D2_FPFD_erato, res_D2_FPFD_erato$padj <= 0.05 & res_D2_FPFD_erato$log2FoldChange < -1)
sign_D2_FPFD_FDup_melp  <- subset(res_D2_FPFD_melp, res_D2_FPFD_melp$padj <= 0.05 & res_D2_FPFD_melp$log2FoldChange < -1)

sign_D2_FDFM_FMup_erato <- subset(res_D2_FDFM_erato, res_D2_FDFM_erato$padj <= 0.05 & res_D2_FDFM_erato$log2FoldChange < -1)
sign_D2_FDFM_FMup_melp  <- subset(res_D2_FDFM_melp, res_D2_FDFM_melp$padj <= 0.05 & res_D2_FDFM_melp$log2FoldChange < -1)



sign_E_D2_FPup <- rbind(sign_D2_FPFM_FPup_erato, sign_D2_FPFD_FPup_erato)
sign_E_D2_FPup <- sign_E_D2_FPup[unique(c(rownames(sign_D2_FPFM_FPup_erato),rownames(sign_D2_FPFD_FPup_erato))),]

sign_E_D2_FPdown <- rbind(sign_D2_FPFM_FMup_erato, sign_D2_FPFD_FDup_erato)
sign_E_D2_FPdown <- sign_E_D2_FPdown[unique(c(rownames(sign_D2_FPFM_FMup_erato),rownames(sign_D2_FPFD_FDup_erato))),]

sign_E_D2_FDup_intersect <- merge(as.data.frame(sign_D2_FPFM_FPup_erato), as.data.frame(sign_D2_FPFD_FPup_erato), by = "row.names")
sign_E_D2_FDdown_intersect <- merge(as.data.frame(sign_D2_FPFM_FMup_erato), as.data.frame(sign_D2_FPFD_FDup_erato), by = "row.names")

nrow(sign_E_D2_FPup)
nrow(sign_E_D2_FPdown)
nrow(sign_E_D2_FDup_intersect)
nrow(sign_E_D2_FDdown_intersect)

sign_M_D2_FPup <- rbind(sign_D2_FPFM_FPup_melp, sign_D2_FPFD_FPup_melp)
sign_M_D2_FPup <- sign_E_D2_FPup[unique(c(rownames(sign_D2_FPFM_FPup_melp),rownames(sign_D2_FPFD_FPup_melp))),]

sign_M_D2_FPdown <- rbind(sign_D2_FPFM_FMup_melp, sign_D2_FPFD_FDup_melp)
# sign_M_D2_FPdown <- sign_E_D2_FPdown[unique(c(rownames(sign_D2_FPFM_FMup_melp),rownames(sign_D2_FPFD_FDup_melp))),]

sign_M_D2_FDup_intersect <- merge(as.data.frame(sign_D2_FPFM_FPup_melp), as.data.frame(sign_D2_FPFD_FPup_melp), by = "row.names")
sign_M_D2_FDdown_intersect <- merge(as.data.frame(sign_D2_FPFM_FMup_melp), as.data.frame(sign_D2_FPFD_FDup_melp), by = "row.names")

nrow(sign_M_D2_FPup)
nrow(sign_M_D2_FPdown)
nrow(sign_M_D2_FDup_intersect)
nrow(sign_M_D2_FDdown_intersect)


sign_E_D2_FMup <- rbind(sign_D2_FPFM_FMup_erato, sign_D2_FDFM_FMup_erato)
sign_E_D2_FMup <- sign_E_D2_FPup[unique(c(rownames(sign_D2_FPFM_FMup_erato),rownames(sign_D2_FDFM_FMup_erato))),]

sign_E_D2_FMdown <- rbind(sign_D2_FPFM_FPup_erato, sign_D2_FDFM_FDup_erato)
# sign_E_D2_FMdown <- sign_E_D2_FPdown[unique(c(rownames(sign_D2_FPFM_FPup_erato),rownames(sign_D2_FDFM_FDup_erato))),]

sign_E_D2_FMup_intersect <- merge(as.data.frame(sign_D2_FPFM_FMup_erato), as.data.frame(sign_D2_FDFM_FMup_erato), by = "row.names")
sign_E_D2_FMdown_intersect <- merge(as.data.frame(sign_D2_FPFM_FPup_erato), as.data.frame(sign_D2_FDFM_FDup_erato), by = "row.names")

nrow(sign_E_D2_FMup)
nrow(sign_E_D2_FMdown)
nrow(sign_E_D2_FMup_intersect)
nrow(sign_E_D2_FMdown_intersect)


sign_M_D2_FMup <- rbind(sign_D2_FPFM_FMup_melp, sign_D2_FDFM_FMup_melp)
# sign_M_D2_FMup <- sign_E_D2_FPup[unique(c(rownames(sign_D2_FPFM_FMup_melp),rownames(sign_D2_FDFM_FMup_melp))),]

sign_M_D2_FMdown <- rbind(sign_D2_FPFM_FPup_melp, sign_D2_FDFM_FDup_melp)
# sign_M_D2_FMdown <- sign_E_D2_FPdown[unique(c(rownames(sign_D2_FPFM_FPup_melp),rownames(sign_D2_FDFM_FDup_melp))),]

sign_M_D2_FMup_intersect <- merge(as.data.frame(sign_D2_FPFM_FMup_melp), as.data.frame(sign_D2_FDFM_FMup_melp), by = "row.names")
sign_M_D2_FMdown_intersect <- merge(as.data.frame(sign_D2_FPFM_FPup_melp), as.data.frame(sign_D2_FDFM_FDup_melp), by = "row.names")

nrow(sign_M_D2_FMup)
nrow(sign_M_D2_FMdown)
nrow(sign_M_D2_FMup_intersect)
nrow(sign_M_D2_FMdown_intersect)


sign_E_D2_FDup <- rbind(sign_D2_FPFD_FDup_erato, sign_D2_FDFM_FDup_erato)
# sign_E_D2_FDup <- sign_E_D2_FPup[unique(c(rownames(sign_D2_FPFD_FDup_erato),rownames(sign_D2_FDFM_FDup_erato))),]

sign_E_D2_FDdown <- rbind(sign_D2_FPFD_FPup_erato, sign_D2_FDFM_FMup_erato)
sign_E_D2_FDdown <- sign_E_D2_FPdown[unique(c(rownames(sign_D2_FPFD_FPup_erato),rownames(sign_D2_FDFM_FMup_erato))),]

sign_E_D2_FDup_intersect <- merge(as.data.frame(sign_D2_FPFD_FDup_erato), as.data.frame(sign_D2_FDFM_FDup_erato), by = "row.names")
sign_E_D2_FDdown_intersect <- merge(as.data.frame(sign_D2_FPFD_FPup_erato), as.data.frame(sign_D2_FDFM_FMup_erato), by = "row.names")

nrow(sign_E_D2_FDup)
nrow(sign_E_D2_FDdown)
nrow(sign_E_D2_FDup_intersect)
nrow(sign_E_D2_FDdown_intersect)



sign_M_D2_FDup <- rbind(sign_D2_FPFD_FDup_melp, sign_D2_FDFM_FDup_melp)
sign_M_D2_FDup <- sign_M_D2_FPup[unique(c(rownames(sign_D2_FPFD_FDup_melp),rownames(sign_D2_FDFM_FDup_melp))),]

sign_M_D2_FDdown <- rbind(sign_D2_FPFD_FPup_melp, sign_D2_FDFM_FMup_melp)
sign_M_D2_FDdown <- sign_M_D2_FPdown[unique(c(rownames(sign_D2_FPFD_FPup_melp),rownames(sign_D2_FDFM_FMup_melp))),]

sign_M_D2_FDup_intersect <- merge(as.data.frame(sign_D2_FPFD_FDup_melp), as.data.frame(sign_D2_FDFM_FDup_melp), by = "row.names")
sign_M_D2_FDdown_intersect <- merge(as.data.frame(sign_D2_FPFD_FPup_melp), as.data.frame(sign_D2_FDFM_FMup_melp), by = "row.names")

nrow(sign_M_D2_FDup)
nrow(sign_M_D2_FDdown)
nrow(sign_M_D2_FDup_intersect)
nrow(sign_M_D2_FDdown_intersect)





############ shared DE expressed with shared peaks
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')
head(homologs)

sign_E_D2_FDup$erato_gene <- rownames(sign_E_D2_FDup)
sign_E_D2_FMup$erato_gene <- rownames(sign_E_D2_FMup)
sign_E_D2_FPup$erato_gene <- rownames(sign_E_D2_FPup)

sign_E_D2_FDup_homol <- merge(as.data.frame(sign_E_D2_FDup), homologs, by = 'erato_gene')
sign_E_D2_FMup_homol <- merge(as.data.frame(sign_E_D2_FMup), homologs, by = 'erato_gene')
sign_E_D2_FPup_homol <- merge(as.data.frame(sign_E_D2_FPup), homologs, by = 'erato_gene')

sign_M_D2_FDup$melp_gene <- rownames(sign_M_D2_FDup)
sign_M_D2_FMup$melp_gene <- rownames(sign_M_D2_FMup)
sign_M_D2_FPup$melp_gene <- rownames(sign_M_D2_FPup)

sign_M_D2_FDup_homol <- merge(as.data.frame(sign_M_D2_FDup), homologs, by = 'melp_gene')
sign_M_D2_FMup_homol <- merge(as.data.frame(sign_M_D2_FMup), homologs, by = 'melp_gene')
sign_M_D2_FPup_homol <- merge(as.data.frame(sign_M_D2_FPup), homologs, by = 'melp_gene')

nrow(merge(sign_E_D2_FDup_homol, sign_M_D2_FDup_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D2_FMup_homol, sign_M_D2_FMup_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D2_FPup_homol, sign_M_D2_FPup_homol, by = 'melp_gene', all = F))


sign_E_D2_FDdown$erato_gene <- rownames(sign_E_D2_FDdown)
sign_E_D2_FMdown$erato_gene <- rownames(sign_E_D2_FMdown)
sign_E_D2_FPdown$erato_gene <- rownames(sign_E_D2_FPdown)

sign_E_D2_FDdown_homol <- merge(as.data.frame(sign_E_D2_FDdown), homologs, by = 'erato_gene')
sign_E_D2_FMdown_homol <- merge(as.data.frame(sign_E_D2_FMdown), homologs, by = 'erato_gene')
sign_E_D2_FPdown_homol <- merge(as.data.frame(sign_E_D2_FPdown), homologs, by = 'erato_gene')

sign_M_D2_FDdown$melp_gene <- rownames(sign_M_D2_FDdown)
sign_M_D2_FMdown$melp_gene <- rownames(sign_M_D2_FMdown)
sign_M_D2_FPdown$melp_gene <- rownames(sign_M_D2_FPdown)

sign_M_D2_FDdown_homol <- merge(as.data.frame(sign_M_D2_FDdown), homologs, by = 'melp_gene')
sign_M_D2_FMdown_homol <- merge(as.data.frame(sign_M_D2_FMdown), homologs, by = 'melp_gene')
sign_M_D2_FPdown_homol <- merge(as.data.frame(sign_M_D2_FPdown), homologs, by = 'melp_gene')

nrow(merge(sign_E_D2_FDdown_homol, sign_M_D2_FDdown_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D2_FMdown_homol, sign_M_D2_FMdown_homol, by = 'melp_gene', all = F))
nrow(merge(sign_E_D2_FPdown_homol, sign_M_D2_FPdown_homol, by = 'melp_gene', all = F))






sign_E_D1D2_FPup <- merge(sign_E_D1_FPup, sign_E_D2_FPup, by = "row.names")
sign_E_D1D2_FMup <- merge(as.data.frame(sign_E_D1_FMup), as.data.frame(na.omit(sign_E_D2_FMup)), by = "row.names")
sign_E_D1D2_FDup <- merge(as.data.frame(sign_E_D1_FDup), as.data.frame(na.omit(sign_E_D2_FDup)), by = "row.names")

sign_M_D1D2_FPup <- merge(sign_M_D1_FPup, sign_M_D2_FPup, by = "row.names") #0
sign_M_D1D2_FMup <- merge(as.data.frame(sign_M_D1_FMup), as.data.frame(sign_M_D2_FMup), by = "row.names")
sign_M_D1D2_FDup <- merge(sign_M_D1_FDup, sign_M_D2_FDup, by = "row.names") #0

nrow(sign_E_D1D2_FPup)
nrow(sign_E_D1D2_FMup)
nrow(sign_E_D1D2_FDup)

nrow(sign_M_D1D2_FPup)
nrow(sign_M_D1D2_FMup)
nrow(sign_M_D1D2_FDup)


sign_E_D1D2_FPdown <- merge(sign_E_D1_FPdown, sign_E_D2_FPdown, by = "row.names")
sign_E_D1D2_FMdown <- merge(as.data.frame(sign_E_D1_FMdown), as.data.frame(na.omit(sign_E_D2_FMdown)), by = "row.names")
sign_E_D1D2_FDdown <- merge(as.data.frame(sign_E_D1_FDdown), as.data.frame(na.omit(sign_E_D2_FDdown)), by = "row.names")

sign_M_D1D2_FPdown <- merge(sign_M_D1_FPdown, sign_M_D2_FPdown, by = "row.names") #0
sign_M_D1D2_FMdown <- merge(as.data.frame(sign_M_D1_FMdown), as.data.frame(sign_M_D2_FMdown), by = "row.names")
sign_M_D1D2_FDdown <- merge(sign_M_D1_FDdown, sign_M_D2_FDdown, by = "row.names") #0

nrow(sign_E_D1D2_FPdown)
nrow(sign_E_D1D2_FMdown)
nrow(sign_E_D1D2_FDdown)

nrow(sign_M_D1D2_FPdown)
nrow(sign_M_D1D2_FMdown)
nrow(sign_M_D1D2_FDdown)

sign_E_D1D2_FPup_intersect <- merge(sign_E_D1_FDup_intersect, sign_E_D2_FDup_intersect, by = "row.names")
nrow(sign_E_D1D2_FPup_intersect)


sign_D1D2_FPup_erato <- merge(as.data.frame(sign_D1_FPup_erato), as.data.frame(sign_D2_FPup_erato), by = "row.names")
sign_D1D2_FMup_erato <- merge(as.data.frame(sign_D1_FMup_erato), as.data.frame(sign_D2_FMup_erato), by = "row.names")
sign_D1D2_FDup_erato <- merge(as.data.frame(sign_D1_FDup_erato), as.data.frame(sign_D2_FDup_erato), by = "row.names")

nrow(sign_D1D2_FPup_erato)
nrow(sign_D1D2_FMup_erato)
nrow(sign_D1D2_FDup_erato)

sign_D1D2_FPup_melp <- merge(as.data.frame(sign_D1_FPup_melp), as.data.frame(sign_D2_FPup_melp), by = "row.names")
sign_D1D2_FMup_melp <- merge(as.data.frame(sign_D1_FMup_melp), as.data.frame(sign_D2_FMup_melp), by = "row.names")
sign_D1D2_FDup_melp <- merge(as.data.frame(sign_D1_FDup_melp), as.data.frame(sign_D2_FDup_melp), by = "row.names")

nrow(sign_D1D2_FPup_melp)
nrow(sign_D1D2_FMup_melp)
nrow(sign_D1D2_FDup_melp)


sign_D1D2_FPdown_erato <- merge(as.data.frame(sign_D1_FPdown_erato), as.data.frame(sign_D2_FPdown_erato), by = "row.names")
sign_D1D2_FMdown_erato <- merge(as.data.frame(sign_D1_FMdown_erato), as.data.frame(sign_D2_FMdown_erato), by = "row.names")
sign_D1D2_FDdown_erato <- merge(as.data.frame(sign_D1_FDdown_erato), as.data.frame(sign_D2_FDdown_erato), by = "row.names")

nrow(sign_D1D2_FPdown_erato)
nrow(sign_D1D2_FMdown_erato)
nrow(sign_D1D2_FDdown_erato)

sign_D1D2_FPdown_melp <- merge(as.data.frame(sign_D1_FPdown_melp), as.data.frame(sign_D2_FPdown_melp), by = "row.names")
sign_D1D2_FMdown_melp <- merge(as.data.frame(sign_D1_FMdown_melp), as.data.frame(sign_D2_FMdown_melp), by = "row.names")
sign_D1D2_FDdown_melp <- merge(as.data.frame(sign_D1_FDdown_melp), as.data.frame(sign_D2_FDdown_melp), by = "row.names")

nrow(sign_D1D2_FPdown_melp)
nrow(sign_D1D2_FMdown_melp)
nrow(sign_D1D2_FDdown_melp)

