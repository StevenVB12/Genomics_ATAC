library(dplyr)
library(DESeq2)

counts_erato <- read.table('counts_nonStranded/counts_table_erato_nonStranded.txt', h=T)

head(counts_erato)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

## 5th
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == 'larval'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))

dds_5th <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                  colData = sampleInfo_sub,
                                  design = ~wing)

## D1
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == '1'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

dds_D1 <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

## D2
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'erato') & 
                           (sampleInfo$day == '2'))
sampleInfo_sub

counts_eratoF <- counts_erato[,-1]
rownames(counts_eratoF) <- counts_erato[,1]
head(counts_eratoF)

counts_erato2 <- select(counts_eratoF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_erato2)

dds_D2 <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

###
# run model
###

rnaDDS_5th <- DESeq(dds_5th)
rnaDDS_D1 <- DESeq(dds_D1)
rnaDDS_D2 <- DESeq(dds_D2)

res_5th_erato <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1_erato <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2_erato <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))

sign_5th_FWup <- subset(res_5th_erato, res_5th_erato$padj <= 0.05 & res_5th_erato$log2FoldChange > 0)
sign_5th_HWup <- subset(res_5th_erato, res_5th_erato$padj <= 0.05 & res_5th_erato$log2FoldChange < 0)

sign_D1_FWup <- subset(res_D1_erato, res_D1_erato$padj <= 0.05 & res_D1_erato$log2FoldChange > 0)
sign_D1_HWup <- subset(res_D1_erato, res_D1_erato$padj <= 0.05 & res_D1_erato$log2FoldChange < 0)

sign_D2_FWup <- subset(res_D2_erato, res_D2_erato$padj <= 0.05 & res_D2_erato$log2FoldChange > 0)
sign_D2_HWup <- subset(res_D2_erato, res_D2_erato$padj <= 0.05 & res_D2_erato$log2FoldChange < 0)

ATAC_FWup_5th <- read.table('FWHW/FWHW_5th_p05FC1_erato_FWup_genes100kb.txt.filtered')
ATAC_HWup_5th <- read.table('FWHW/FWHW_5th_p05FC1_erato_HWup_genes100kb.txt.filtered')
ATAC_FWHW_5th <- rbind(ATAC_FWup_5th, ATAC_HWup_5th)

res_5th_ATAC_FWup <- res_5th_erato[row.names(res_5th_erato) %in% ATAC_FWup_5th$V2, ]
res_5th_ATAC_HWup <- res_5th_erato[row.names(res_5th_erato) %in% ATAC_HWup_5th$V2, ]
res_5th_ATAC_FWHW <- res_5th_erato[row.names(res_5th_erato) %in% ATAC_FWHW_5th$V2, ]

res_5th_ATAC_FWup_DE_erato <- subset(res_5th_ATAC_FWup, res_5th_ATAC_FWup$padj <=0.05)
res_5th_ATAC_HWup_DE_erato <- subset(res_5th_ATAC_HWup, res_5th_ATAC_HWup$padj <=0.05)

res_5th_ATAC_FWup_DE_erato <- sign_5th_FWup[row.names(sign_5th_FWup) %in% ATAC_FWHW_5th$V2, ]
res_5th_ATAC_HWup_DE_erato <- sign_5th_HWup[row.names(sign_5th_HWup) %in% ATAC_FWHW_5th$V2, ]

ATAC_FWup_D1 <- read.table('FWHW/FWHW_D1_p05FC1_erato_FWup_genes100kb.txt.filtered')
ATAC_HWup_D1 <- read.table('FWHW/FWHW_D1_p05FC1_erato_HWup_genes100kb.txt.filtered')
ATAC_FWHW_D1 <- rbind(ATAC_FWup_D1, ATAC_HWup_D1)

res_D1_ATAC_FWup <- res_D1_erato[row.names(res_D1_erato) %in% ATAC_FWup_D1$V2, ]
res_D1_ATAC_HWup <- res_D1_erato[row.names(res_D1_erato) %in% ATAC_HWup_D1$V2, ]

res_D1_ATAC_FWup_DE_erato <- subset(res_D1_ATAC_FWup, res_D1_ATAC_FWup$padj <=0.05)
res_D1_ATAC_HWup_DE_erato <- subset(res_D1_ATAC_HWup, res_D1_ATAC_HWup$padj <=0.05)

res_D1_ATAC_FWup_DE_erato <- sign_D1_FWup[row.names(sign_D1_FWup) %in% ATAC_FWHW_D1$V2, ]
res_D1_ATAC_HWup_DE_erato <- sign_D1_HWup[row.names(sign_D1_HWup) %in% ATAC_FWHW_D1$V2, ]

ATAC_FWup_D2 <- read.table('FWHW/FWHW_D2_p05FC1_erato_FWup_genes100kb.txt.filtered')
ATAC_HWup_D2 <- read.table('FWHW/FWHW_D2_p05FC1_erato_HWup_genes100kb.txt.filtered')
ATAC_FWHW_D2 <- rbind(ATAC_FWup_D2, ATAC_HWup_D2)

res_D2_ATAC_FWup <- res_D2_erato[row.names(res_D2_erato) %in% ATAC_FWup_D2$V2, ]
res_D2_ATAC_HWup <- res_D2_erato[row.names(res_D2_erato) %in% ATAC_HWup_D2$V2, ]

res_D2_ATAC_FWup_DE_erato <- subset(res_D2_ATAC_FWup, res_D2_ATAC_FWup$padj <=0.05)
res_D2_ATAC_HWup_DE_erato <- subset(res_D2_ATAC_HWup, res_D2_ATAC_HWup$padj <=0.05)

res_D2_ATAC_FWup_DE_erato <- sign_D2_FWup[row.names(sign_D2_FWup) %in% ATAC_FWHW_D2$V2, ]
res_D2_ATAC_HWup_DE_erato <- sign_D2_HWup[row.names(sign_D2_HWup) %in% ATAC_FWHW_D2$V2, ]

###########
counts_melp <- read.table('counts_nonStranded/counts_table_melp_nonStranded.txt', h=T)

head(counts_melp)

sampleInfo <- read.table('sample_info_Joe.txt', h = T)
sampleInfo <- unique(sampleInfo)

head(sampleInfo)

## 5th
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == 'larval'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))

dds_5th <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                  colData = sampleInfo_sub,
                                  design = ~wing)

## D1
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == '1'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_melp2)

dds_D1 <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

## D2
sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == '2'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_melp2)

dds_D2 <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

###
# run model
###

rnaDDS_5th <- DESeq(dds_5th)
rnaDDS_D1 <- DESeq(dds_D1)
rnaDDS_D2 <- DESeq(dds_D2)

res_5th_melp <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1_melp <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2_melp <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))

sign_5th_FWup <- subset(res_5th_melp, res_5th_melp$padj <= 0.05 & res_5th_melp$log2FoldChange > 0)
sign_5th_HWup <- subset(res_5th_melp, res_5th_melp$padj <= 0.05 & res_5th_melp$log2FoldChange < 0)

sign_D1_FWup <- subset(res_D1_melp, res_D1_melp$padj <= 0.05 & res_D1_melp$log2FoldChange > 0)
sign_D1_HWup <- subset(res_D1_melp, res_D1_melp$padj <= 0.05 & res_D1_melp$log2FoldChange < 0)

sign_D2_FWup <- subset(res_D2_melp, res_D2_melp$padj <= 0.05 & res_D2_melp$log2FoldChange > 0)
sign_D2_HWup <- subset(res_D2_melp, res_D2_melp$padj <= 0.05 & res_D2_melp$log2FoldChange < 0)


ATAC_FWup_5th <- read.table('FWHW/FWHW_5th_p05FC1_melp_FWup_genes100kb.txt.filtered')
ATAC_HWup_5th <- read.table('FWHW/FWHW_5th_p05FC1_melp_HWup_genes100kb.txt.filtered')
ATAC_FWHW_5th <- rbind(ATAC_FWup_5th, ATAC_HWup_5th)

res_5th_ATAC_FWup <- res_5th_melp[row.names(res_5th_melp) %in% ATAC_FWup_5th$V2, ]
res_5th_ATAC_HWup <- res_5th_melp[row.names(res_5th_melp) %in% ATAC_HWup_5th$V2, ]

res_5th_ATAC_FWup_DE_melp <- subset(res_5th_ATAC_FWup, res_5th_ATAC_FWup$padj <=0.05)
res_5th_ATAC_HWup_DE_melp <- subset(res_5th_ATAC_HWup, res_5th_ATAC_HWup$padj <=0.05)

res_5th_ATAC_FWup_DE_melp <- sign_5th_FWup[row.names(sign_5th_FWup) %in% ATAC_FWHW_5th$V2, ]
res_5th_ATAC_HWup_DE_melp <- sign_5th_HWup[row.names(sign_5th_HWup) %in% ATAC_FWHW_5th$V2, ]

ATAC_FWup_D1 <- read.table('FWHW/FWHW_D1_p05FC1_melp_FWup_genes100kb.txt.filtered')
ATAC_HWup_D1 <- read.table('FWHW/FWHW_D1_p05FC1_melp_HWup_genes100kb.txt.filtered')
ATAC_FWHW_D1 <- rbind(ATAC_FWup_D1, ATAC_HWup_D1)

res_D1_ATAC_FWup <- res_D1_melp[row.names(res_D1_melp) %in% ATAC_FWup_D1$V2, ]
res_D1_ATAC_HWup <- res_D1_melp[row.names(res_D1_melp) %in% ATAC_HWup_D1$V2, ]

res_D1_ATAC_FWup_DE_melp <- subset(res_D1_ATAC_FWup, res_D1_ATAC_FWup$padj <=0.05)
res_D1_ATAC_HWup_DE_melp <- subset(res_D1_ATAC_HWup, res_D1_ATAC_HWup$padj <=0.05)

res_D1_ATAC_FWup_DE_melp <- sign_D1_FWup[row.names(sign_D1_FWup) %in% ATAC_FWHW_D1$V2, ]
res_D2_ATAC_HWup_DE_melp <- sign_D1_HWup[row.names(sign_D1_HWup) %in% ATAC_FWHW_D1$V2, ]

ATAC_FWup_D2 <- read.table('FWHW/FWHW_D2_p05FC1_melp_FWup_genes100kb.txt.filtered')
ATAC_HWup_D2 <- read.table('FWHW/FWHW_D2_p05FC1_melp_HWup_genes100kb.txt.filtered')
ATAC_FWHW_D2 <- rbind(ATAC_FWup_D2, ATAC_HWup_D2)

res_D2_ATAC_FWup <- res_D2_melp[row.names(res_D2_melp) %in% ATAC_FWup_D2$V2, ]
res_D2_ATAC_HWup <- res_D2_melp[row.names(res_D2_melp) %in% ATAC_HWup_D2$V2, ]

res_D2_ATAC_FWup_DE_melp <- subset(res_D2_ATAC_FWup, res_D2_ATAC_FWup$padj <=0.05)
res_D2_ATAC_HWup_DE_melp <- subset(res_D2_ATAC_HWup, res_D2_ATAC_HWup$padj <=0.05)

res_D2_ATAC_FWup_DE_melp <- sign_D2_FWup[row.names(sign_D2_FWup) %in% ATAC_FWHW_D2$V2, ]
res_D2_ATAC_HWup_DE_melp <- sign_D2_HWup[row.names(sign_D2_HWup) %in% ATAC_FWHW_D2$V2, ]


############ shared DE expressed with shared peaks
homologs <- read.table('homologs_melp_erato.txt')
colnames(homologs) <- c('melp_gene', 'erato_gene')

res_5th_ATAC_FWup_DE_erato$erato_gene <- rownames(res_5th_ATAC_FWup_DE_erato)
res_5th_ATAC_HWup_DE_erato$erato_gene <- rownames(res_5th_ATAC_HWup_DE_erato)
res_D1_ATAC_FWup_DE_erato$erato_gene <- rownames(res_D1_ATAC_FWup_DE_erato)
res_D1_ATAC_HWup_DE_erato$erato_gene <- rownames(res_D1_ATAC_HWup_DE_erato)
res_D2_ATAC_FWup_DE_erato$erato_gene <- rownames(res_D2_ATAC_FWup_DE_erato)
res_D2_ATAC_HWup_DE_erato$erato_gene <- rownames(res_D2_ATAC_HWup_DE_erato)

res_5th_ATAC_FWup_DE_melp$melp_gene <- rownames(res_5th_ATAC_FWup_DE_melp)
res_5th_ATAC_HWup_DE_melp$melp_gene <- rownames(res_5th_ATAC_HWup_DE_melp)
res_D1_ATAC_FWup_DE_melp$melp_gene <- rownames(res_D1_ATAC_FWup_DE_melp)
res_D1_ATAC_HWup_DE_melp$melp_gene <- rownames(res_D1_ATAC_HWup_DE_melp)
res_D2_ATAC_FWup_DE_melp$melp_gene <- rownames(res_D2_ATAC_FWup_DE_melp)
res_D2_ATAC_HWup_DE_melp$melp_gene <- rownames(res_D2_ATAC_HWup_DE_melp)

res_5th_ATAC_FWup_DE_erato_homol <- merge(as.data.frame(res_5th_ATAC_FWup_DE_erato), homologs, by = 'erato_gene')
res_5th_ATAC_HWup_DE_erato_homol <- merge(as.data.frame(res_5th_ATAC_HWup_DE_erato), homologs, by = 'erato_gene')
res_D1_ATAC_FWup_DE_erato_homol <- merge(as.data.frame(res_D1_ATAC_FWup_DE_erato), homologs, by = 'erato_gene')
res_D1_ATAC_HWup_DE_erato_homol <- merge(as.data.frame(res_D1_ATAC_HWup_DE_erato), homologs, by = 'erato_gene')
res_D2_ATAC_FWup_DE_erato_homol <- merge(as.data.frame(res_D2_ATAC_FWup_DE_erato), homologs, by = 'erato_gene')
res_D2_ATAC_HWup_DE_erato_homol <- merge(as.data.frame(res_D2_ATAC_HWup_DE_erato), homologs, by = 'erato_gene')

res_5th_ATAC_FWup_DE_melp_homol <- merge(as.data.frame(res_5th_ATAC_FWup_DE_melp), homologs, by = 'melp_gene')
res_5th_ATAC_HWup_DE_melp_homol <- merge(as.data.frame(res_5th_ATAC_HWup_DE_melp), homologs, by = 'melp_gene')
res_D1_ATAC_FWup_DE_melp_homol <- merge(as.data.frame(res_D1_ATAC_FWup_DE_melp), homologs, by = 'melp_gene')
res_D1_ATAC_HWup_DE_melp_homol <- merge(as.data.frame(res_D1_ATAC_HWup_DE_melp), homologs, by = 'melp_gene')
res_D2_ATAC_FWup_DE_melp_homol <- merge(as.data.frame(res_D2_ATAC_FWup_DE_melp), homologs, by = 'melp_gene')
res_D2_ATAC_HWup_DE_melp_homol <- merge(as.data.frame(res_D2_ATAC_HWup_DE_melp), homologs, by = 'melp_gene')

nrow(merge(res_5th_ATAC_FWup_DE_erato_homol, res_5th_ATAC_FWup_DE_melp_homol, by = 'melp_gene', all = F))
nrow(merge(res_D1_ATAC_FWup_DE_erato_homol, res_D1_ATAC_FWup_DE_melp_homol, by = 'melp_gene', all = F))
nrow(merge(res_D2_ATAC_FWup_DE_erato_homol, res_D2_ATAC_FWup_DE_melp_homol, by = 'melp_gene', all = F))

nrow(merge(res_5th_ATAC_HWup_DE_erato_homol, res_5th_ATAC_HWup_DE_melp_homol, by = 'melp_gene', all = F))
nrow(merge(res_D1_ATAC_HWup_DE_erato_homol, res_D1_ATAC_HWup_DE_melp_homol, by = 'melp_gene', all = F))
nrow(merge(res_D2_ATAC_HWup_DE_erato_homol, res_D2_ATAC_HWup_DE_melp_homol, by = 'melp_gene', all = F))



#### shared DE expressed

res_5th_erato$erato_gene <- rownames(res_5th_erato)
res_D1_erato$erato_gene <- rownames(res_D1_erato)
res_D2_erato$erato_gene <- rownames(res_D2_erato)

res_5th_melp$melp_gene <- rownames(res_5th_melp)
res_D1_melp$melp_gene <- rownames(res_D1_melp)
res_D2_melp$melp_gene <- rownames(res_D2_melp)

res_5th_erato_DE_FWup <- subset(res_5th_erato, res_5th_erato$padj <= 0.05 & res_5th_erato$log2FoldChange > 0)
res_D1_erato_DE_FWup <- subset(res_D1_erato, res_D1_erato$padj <= 0.05 & res_D1_erato$log2FoldChange > 0)
res_D2_erato_DE_FWup <- subset(res_D2_erato, res_D2_erato$padj <= 0.05 & res_D2_erato$log2FoldChange > 0)

res_5th_melp_DE_FWup <- subset(res_5th_melp, res_5th_melp$padj <= 0.05 & res_5th_melp$log2FoldChange > 0)
res_D1_melp_DE_FWup <- subset(res_D1_melp, res_D1_melp$padj <= 0.05 & res_D1_melp$log2FoldChange > 0)
res_D2_melp_DE_FWup <- subset(res_D2_melp, res_D2_melp$padj <= 0.05 & res_D2_melp$log2FoldChange > 0)


res_5th_erato_DE_HWup <- subset(res_5th_erato, res_5th_erato$padj <= 0.05 & res_5th_erato$log2FoldChange < 0)
res_D1_erato_DE_HWup <- subset(res_D1_erato, res_D1_erato$padj <= 0.05 & res_D1_erato$log2FoldChange < 0)
res_D2_erato_DE_HWup <- subset(res_D2_erato, res_D2_erato$padj <= 0.05 & res_D2_erato$log2FoldChange < 0)

res_5th_melp_DE_HWup <- subset(res_5th_melp, res_5th_melp$padj <= 0.05 & res_5th_melp$log2FoldChange < 0)
res_D1_melp_DE_HWup <- subset(res_D1_melp, res_D1_melp$padj <= 0.05 & res_D1_melp$log2FoldChange < 0)
res_D2_melp_DE_HWup <- subset(res_D2_melp, res_D2_melp$padj <= 0.05 & res_D2_melp$log2FoldChange < 0)


res_5th_DE_erato_homol_FWup <- merge(as.data.frame(res_5th_erato_DE_FWup), homologs, by = 'erato_gene')
res_D1_DE_erato_homol_FWup <- merge(as.data.frame(res_D1_erato_DE_FWup), homologs, by = 'erato_gene')
res_D2_DE_erato_homol_FWup <- merge(as.data.frame(res_D2_erato_DE_FWup), homologs, by = 'erato_gene')

res_5th_DE_melp_homol_FWup <- merge(as.data.frame(res_5th_melp_DE_FWup), homologs, by = 'melp_gene')
res_D1_DE_melp_homol_FWup <- merge(as.data.frame(res_D1_melp_DE_FWup), homologs, by = 'melp_gene')
res_D2_DE_melp_homol_FWup <- merge(as.data.frame(res_D2_melp_DE_FWup), homologs, by = 'melp_gene')

nrow(merge(res_5th_DE_erato_homol_FWup, res_5th_DE_melp_homol_FWup, by = 'melp_gene', all = F))
nrow(merge(res_D1_DE_erato_homol_FWup, res_D1_DE_melp_homol_FWup, by = 'melp_gene', all = F))
nrow(merge(res_D2_DE_erato_homol_FWup, res_D2_DE_melp_homol_FWup, by = 'melp_gene', all = F))


res_5th_DE_erato_homol_HWup <- merge(as.data.frame(res_5th_erato_DE_HWup), homologs, by = 'erato_gene')
res_D1_DE_erato_homol_HWup <- merge(as.data.frame(res_D1_erato_DE_HWup), homologs, by = 'erato_gene')
res_D2_DE_erato_homol_HWup <- merge(as.data.frame(res_D2_erato_DE_HWup), homologs, by = 'erato_gene')

res_5th_DE_melp_homol_HWup <- merge(as.data.frame(res_5th_melp_DE_HWup), homologs, by = 'melp_gene')
res_D1_DE_melp_homol_HWup <- merge(as.data.frame(res_D1_melp_DE_HWup), homologs, by = 'melp_gene')
res_D2_DE_melp_homol_HWup <- merge(as.data.frame(res_D2_melp_DE_HWup), homologs, by = 'melp_gene')

nrow(merge(res_5th_DE_erato_homol_HWup, res_5th_DE_melp_homol_HWup, by = 'melp_gene', all = F))
nrow(merge(res_D1_DE_erato_homol_HWup, res_D1_DE_melp_homol_HWup, by = 'melp_gene', all = F))
nrow(merge(res_D2_DE_erato_homol_HWup, res_D2_DE_melp_homol_HWup, by = 'melp_gene', all = F))
