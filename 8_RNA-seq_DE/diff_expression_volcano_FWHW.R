library(dplyr)
library(DESeq2)

counts_erato <- read.table('counts_table_erato_nonStranded.txt', h=T)

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

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_erato2), las=2)

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

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_erato2), las=2)

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

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_erato2), las=2)

dds_D2 <- DESeqDataSetFromMatrix(countData = counts_erato2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

###
# run model
###

rnaDDS_5th <- DESeq(dds_5th)
rnaDDS_D1 <- DESeq(dds_D1)
rnaDDS_D2 <- DESeq(dds_D2)


res_5th <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1 <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2 <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))


DE_5th <- read.table('FWHW/FW_HW_genes_near_DE_peaks_5th.txt')

res_5th_ATAC <- res_5th[row.names(res_5th) %in% DE_5th$V1, ]

subset(res_5th_ATAC, res_5th_ATAC$pvalue <=0.05)

# boxplot(res_5th_ATAC$log2FoldChange)




DE_D1 <- read.table('FWHW/FW_HW_genes_near_DE_peaks_D1.txt')

res_D1_ATAC <- res_D1[row.names(res_D1) %in% DE_D1$V1, ]

(subset(res_D1_ATAC, res_D1_ATAC$pvalue <=0.05))

# boxplot(res_5th_ATAC$log2FoldChange)

res_5th_padj <- subset(res_5th, res_5th$padj < 0.05)
res_D1_padj <- subset(res_D1, res_D1$padj < 0.05)
res_D2_padj <- subset(res_D2, res_D2$padj < 0.05)

TFs <- read.table("FWHW/TF_IDs_FWHW.txt", h=T)
Inter <- read.table("FWHW/Interesting_genes_near_ATACpeaks_FWHW.txt", h=T, sep = '\t')

res_5th_ATACmotif <- res_5th[row.names(res_5th) %in% TFs$eratoID, ]
res_D1_ATACmotif <- res_D1[row.names(res_D1) %in% TFs$eratoID, ]
res_D2_ATACmotif <- res_D2[row.names(res_D2) %in% TFs$eratoID, ]

res_5th_Inter <- res_5th[row.names(res_5th) %in% Inter$era, ]
res_D1_Inter <- res_D1[row.names(res_D1) %in% Inter$era, ]
res_D2_Inter <- res_D2[row.names(res_D2) %in% Inter$era, ]

par(mfrow=c(1,6), mar=c(2,2,2,0))
plot(res_5th$log2FoldChange, -log10(res_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_5th_padj$log2FoldChange, -log10(res_5th_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_5th_ATACmotif$log2FoldChange, -log10(res_5th_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_5th_Inter$log2FoldChange, -log10(res_5th_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)

plot(res_D1$log2FoldChange, -log10(res_D1$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_D1_padj$log2FoldChange, -log10(res_D1_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_D1_ATACmotif$log2FoldChange, -log10(res_D1_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_D1_Inter$log2FoldChange, -log10(res_D1_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)

plot(res_D2$log2FoldChange, -log10(res_D2$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_D2_padj$log2FoldChange, -log10(res_D2_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_D2_ATACmotif$log2FoldChange, -log10(res_D2_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_D2_Inter$log2FoldChange, -log10(res_D2_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)




#######################################






counts_melp <- read.table('counts_table_melp_nonStranded.txt', h=T)

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

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_melp2), las=2)

dds_5th <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                  colData = sampleInfo_sub,
                                  design = ~wing)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == '1'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_melp2)

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_melp2), las=2)

dds_D1 <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

sampleInfo_sub <- subset(sampleInfo,  (sampleInfo$species == 'melpomene') & 
                           (sampleInfo$day == '2'))
sampleInfo_sub

counts_melpF <- counts_melp[,-1]
rownames(counts_melpF) <- counts_melp[,1]
head(counts_melpF)

counts_melp2 <- select(counts_melpF, c(as.character(sampleInfo_sub$sample_id)))
head(counts_melp2)

# par(mar=c(1,1,1,1))
# barplot(colSums(counts_melp2), las=2)

dds_D2 <- DESeqDataSetFromMatrix(countData = counts_melp2,
                                 colData = sampleInfo_sub,
                                 design = ~wing)

###
# run model
###

rnaDDS_5th <- DESeq(dds_5th)
rnaDDS_D1 <- DESeq(dds_D1)
rnaDDS_D2 <- DESeq(dds_D2)


res_5th <- results(rnaDDS_5th, contrast=c('wing','forewing','hindwing'))
res_D1 <- results(rnaDDS_D1, contrast=c('wing','forewing','hindwing'))
res_D2 <- results(rnaDDS_D2, contrast=c('wing','forewing','hindwing'))


DE_5th <- read.table('FWHW/FW_HW_genes_near_DE_peaks_5th_melp.txt')

res_5th_ATAC <- res_5th[row.names(res_5th) %in% DE_5th$V1, ]

subset(res_5th_ATAC, res_5th_ATAC$pvalue <=0.05)

# boxplot(res_5th_ATAC$log2FoldChange)




DE_D1 <- read.table('FWHW/FW_HW_genes_near_DE_peaks_D1_melp.txt')

res_D1_ATAC <- res_D1[row.names(res_D1) %in% DE_D1$V1, ]

(subset(res_D1_ATAC, res_D1_ATAC$pvalue <=0.05))

# boxplot(res_5th_ATAC$log2FoldChange)

res_5th_padj <- subset(res_5th, res_5th$padj < 0.05)
res_D1_padj <- subset(res_D1, res_D1$padj < 0.05)
res_D2_padj <- subset(res_D2, res_D2$padj < 0.05)

TFs <- read.table("FWHW/TF_IDs_FWHW.txt", h=T)
Inter <- read.table("FWHW/Interesting_genes_near_ATACpeaks_FWHW.txt", h=T, sep = '\t')

res_5th_ATACmotif <- res_5th[row.names(res_5th) %in% TFs$melpID, ]
res_D1_ATACmotif <- res_D1[row.names(res_D1) %in% TFs$melpID, ]
res_D2_ATACmotif <- res_D2[row.names(res_D2) %in% TFs$melpID, ]

res_5th_Inter <- res_5th[row.names(res_5th) %in% Inter$melp, ]
res_D1_Inter <- res_D1[row.names(res_D1) %in% Inter$melp, ]
res_D2_Inter <- res_D2[row.names(res_D2) %in% Inter$melp, ]

# par(mfrow=c(1,6), mar=c(2,2,2,0))
plot(res_5th$log2FoldChange, -log10(res_5th$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_5th_padj$log2FoldChange, -log10(res_5th_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_5th_ATACmotif$log2FoldChange, -log10(res_5th_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_5th_Inter$log2FoldChange, -log10(res_5th_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)

plot(res_D1$log2FoldChange, -log10(res_D1$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_D1_padj$log2FoldChange, -log10(res_D1_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_D1_ATACmotif$log2FoldChange, -log10(res_D1_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_D1_Inter$log2FoldChange, -log10(res_D1_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)

plot(res_D2$log2FoldChange, -log10(res_D2$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.04), labels = F)
par(new=T)
plot(res_D2_padj$log2FoldChange, -log10(res_D2_padj$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('black', alpha=0.5), axes = F)
par(new=T)
plot(res_D2_ATACmotif$log2FoldChange, -log10(res_D2_ATACmotif$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('red', alpha=1), axes=F)
# par(new=T)
# plot(res_D2_Inter$log2FoldChange, -log10(res_D2_Inter$pvalue), pch=19, xlim=c(-10,10), ylim=c(0,15), col=adjustcolor('green', alpha=1), axes=F)

