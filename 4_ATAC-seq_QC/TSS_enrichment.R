# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ATACseqQC")

library(GenomicFeatures)
txs <- makeTxDbFromGFF(file = 'G:\\My Drive\\ATAC-seq-share\\gff\\heliconius_erato_demophoon_v1_core_32_85_1.gff')
txs <- transcripts(txs)
summary(x)

library(ATACseqQC)

bamfile1 <- 'G:\\Shared drives\\Heliconius_genomics_UPRRP\\BAM_ATACseq_trimmomatic\\E4-FW_Herato.trim.filtered.sorted.nd.bam'
bamfile.labels1 <- gsub(".trim.filtered.sorted.nd.bam", "", basename(bamfile1))

baf1 <- readBamFile(bamfile1, asMates=TRUE, bigFile=TRUE)
# shift cut site
gal1 <- shiftGAlignmentsList(baf1)


# Transcription Start Site (TSS) Enrichment Score
tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore

plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
