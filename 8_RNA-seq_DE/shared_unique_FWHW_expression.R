

diff_erato <- read.table('FWHW/FWHW_diff_erato.txt')
diff_melp <- read.table('FWHW/FWHW_diff_melp.txt')
head(diff_melp,n=30)

# ## clean melp gene set
# melpToMask <- read.table('melptomask.names', sep='\t')
# head(melpToMask, n=30)
# nrow(melpToMask)
# melpToMask_genes <- subset(melpToMask, melpToMask$V3 == 'mRNA')
# melpToMask_genes2 <- separate(data = melpToMask_genes, col = 'V9', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
# nrow(melpToMask_genes2)
# 
# melpToMask_genes2$annot <- gsub("ID=", "", melpToMask_genes2$annot)
# head(melpToMask_genes2, n=30)

# diff_melp$rown <- rownames(diff_melp)
# results1 = setdiff(diff_melp$rown, as.character(melpToMask_genes2$annot))
# diff_melp <- diff_melp[results1,]
#

FWup_erato <- read.table('FWHW/FWHW_5th_p05FC1_eratoUnique_FWup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
HWup_erato <- read.table('FWHW/FWHW_5th_p05FC1_eratoUnique_HWup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

FWup_melp <- read.table('FWHW/FWHW_5th_p05FC1_melpUnique_FWup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
HWup_melp <- read.table('FWHW/FWHW_5th_p05FC1_melpUnique_HWup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

FWup_erato_shared <- read.table('FWHW/FWHW_5th_p05FC1_shared_FWup_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
HWup_erato_shared <- read.table('FWHW/FWHW_5th_p05FC1_shared_HWup_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

FWup_melp_shared <- read.table('FWHW/FWHW_5th_p05FC1_shared_FWup_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
HWup_melp_shared <- read.table('FWHW/FWHW_5th_p05FC1_shared_HWup_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


FWup_erato2 <- separate(data = FWup_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
FWup_erato2$annot <- gsub('ID=', '', FWup_erato2$annot)
FWup_erato2$annot <- gsub('-E\\d+', '', FWup_erato2$annot)
FWup_erato_FC <- diff_erato[FWup_erato2$annot,]
FWup_erato_FC$group <- 'openFW'
FWup_erato_FC$species <- 'erato'
FWup_erato_FC$cat <- 'unique'

HWup_erato2 <- separate(data = HWup_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
HWup_erato2$annot <- gsub('ID=', '', HWup_erato2$annot)
HWup_erato2$annot <- gsub('-E\\d+', '', HWup_erato2$annot)
HWup_erato_FC <- diff_erato[HWup_erato2$annot,]
HWup_erato_FC$group <- 'openHW'
HWup_erato_FC$species <- 'erato'
HWup_erato_FC$cat <- 'unique'

FWup_melp2 <- separate(data = FWup_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
FWup_melp2$annot <- gsub('ID=', '', FWup_melp2$annot)
FWup_melp2$annot <- gsub('-E\\d+', '', FWup_melp2$annot)
FWup_melp_FC <- diff_melp[FWup_melp2$annot,]
FWup_melp_FC$group <- 'openFW'
FWup_melp_FC$species <- 'melp'
FWup_melp_FC$cat <- 'unique'

HWup_melp2 <- separate(data = HWup_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
HWup_melp2$annot <- gsub('ID=', '', HWup_melp2$annot)
HWup_melp2$annot <- gsub('-E\\d+', '', HWup_melp2$annot)
HWup_melp_FC <- diff_melp[HWup_melp2$annot,]
HWup_melp_FC$group <- 'openHW'
HWup_melp_FC$species <- 'melp'
HWup_melp_FC$cat <- 'unique'


FWup_erato_shared2 <- separate(data = FWup_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
FWup_erato_shared2$annot <- gsub('ID=', '', FWup_erato_shared2$annot)
FWup_erato_shared2$annot <- gsub('-E\\d+', '', FWup_erato_shared2$annot)
FWup_erato_shared_FC <- diff_erato[FWup_erato_shared2$annot,]
FWup_erato_shared_FC$group <- 'openFW'
FWup_erato_shared_FC$species <- 'erato'
FWup_erato_shared_FC$cat <- 'shared'

HWup_erato_shared2 <- separate(data = HWup_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
HWup_erato_shared2$annot <- gsub('ID=', '', HWup_erato_shared2$annot)
HWup_erato_shared2$annot <- gsub('-E\\d+', '', HWup_erato_shared2$annot)
HWup_erato_shared_FC <- diff_erato[HWup_erato_shared2$annot,]
HWup_erato_shared_FC$group <- 'openHW'
HWup_erato_shared_FC$species <- 'erato'
HWup_erato_shared_FC$cat <- 'shared'

FWup_melp_shared2 <- separate(data = FWup_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
FWup_melp_shared2$annot <- gsub('ID=', '', FWup_melp_shared2$annot)
FWup_melp_shared2$annot <- gsub('-E\\d+', '', FWup_melp_shared2$annot)
FWup_melp_shared_FC <- diff_melp[FWup_melp_shared2$annot,]
FWup_melp_shared_FC$group <- 'openFW'
FWup_melp_shared_FC$species <- 'melp'
FWup_melp_shared_FC$cat <- 'shared'

HWup_melp_shared2 <- separate(data = HWup_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
HWup_melp_shared2$annot <- gsub('ID=', '', HWup_melp_shared2$annot)
HWup_melp_shared2$annot <- gsub('-E\\d+', '', HWup_melp_shared2$annot)
HWup_melp_shared_FC <- diff_melp[HWup_melp_shared2$annot,]
HWup_melp_shared_FC$group <- 'openHW'
HWup_melp_shared_FC$species <- 'melp'
HWup_melp_shared_FC$cat <- 'shared'

FC_all <- rbind(FWup_erato_FC, HWup_erato_FC, FWup_melp_FC, HWup_melp_FC,
                FWup_erato_shared_FC, HWup_erato_shared_FC, FWup_melp_shared_FC, HWup_melp_shared_FC)
FC_all[is.na(FC_all)] <- 0


library(ggplot2)
library(tidyverse)
library(viridis)

ggplot(FC_all, aes(x=group, y=log2FoldChange*-1, group=interaction(cat, group), fill=group)) +
  geom_boxplot(outlier.size=0) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.5, width=0.1) +
  theme(
    plot.title = element_text(size=11)
  ) +facet_grid(~cat*species)




