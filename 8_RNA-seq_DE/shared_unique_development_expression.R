library(dplyr)

diff_erato_5th_D1D2 <- read.table('development/5th_D1D2_diff_erato.txt')
diff_erato_D1_5thD2 <- read.table('development/D1_5thD2_diff_erato.txt')
diff_erato_D2_5thD1 <- read.table('development/D2_5thD1_diff_erato.txt')

diff_melp_5th_D1D2 <- read.table('development/5th_D1D2_diff_melp.txt')
diff_melp_D1_5thD2 <- read.table('development/D1_5thD2_diff_melp.txt')
diff_melp_D2_5thD1 <- read.table('development/D2_5thD1_diff_melp.txt')

##
# melpToMask <- read.table('melptomask.names', sep='\t')
# melpToMask_genes <- subset(melpToMask, melpToMask$V3 == 'mRNA')
# melpToMask_genes2 <- separate(data = melpToMask_genes, col = 'V9', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
# 
# melpToMask_genes2$annot <- gsub("ID=", "", melpToMask_genes2$annot)

# diff_melp_5th_D1D2$rown <- rownames(diff_melp_5th_D1D2)
# results1 = diff_melp_5th_D1D2(diff_melp_5th_D1D2$rown, as.character(melpToMask_genes2$annot))
# diff_melp_5th_D1D2 <- diff_melp_5th_D1D2[results1,]
# 
# diff_melp_D1_5thD2$rown <- rownames(diff_melp_D1_5thD2)
# results1 = diff_melp_D1_5thD2(diff_melp_D1_5thD2$rown, as.character(melpToMask_genes2$annot))
# diff_melp_D1_5thD2 <- diff_melp_D1_5thD2[results1,]
# 
# diff_melp_D2_5thD1$rown <- rownames(diff_melp_D2_5thD1)
# results1 = diff_melp_D2_5thD1(diff_melp_D2_5thD1$rown, as.character(melpToMask_genes2$annot))
# diff_melp_D2_5thD1 <- diff_melp_D2_5thD1[results1,]
##


dev_5thup_erato <- read.table('development/development_5th_p05FC1_eratoUnique_5thup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_5thdown_erato <- read.table('development/development_5th_p05FC1_eratoUnique_5thdown_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_5thup_melp <- read.table('development/development_5th_p05FC1_melpUnique_5thup_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_5thdown_melp <- read.table('development/development_5th_p05FC1_melpUnique_5thdown_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


dev_D1up_erato <- read.table('development/development_D1_p05FC1_eratoUnique_D1up_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D1down_erato <- read.table('development/development_D1_p05FC1_eratoUnique_D1down_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_D1up_melp <- read.table('development/development_D1_p05FC1_melpUnique_D1up_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D1down_melp <- read.table('development/development_D1_p05FC1_melpUnique_D1down_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


dev_D2up_erato <- read.table('development/development_D2_p05FC1_eratoUnique_D2up_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D2down_erato <- read.table('development/development_D2_p05FC1_eratoUnique_D2down_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_D2up_melp <- read.table('development/development_D2_p05FC1_melpUnique_D2up_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D2down_melp <- read.table('development/development_D2_p05FC1_melpUnique_D2down_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


dev_5thup_erato_shared <- read.table('development/development_5th_p05FC1_shared_5thup_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_5thdown_erato_shared <- read.table('development/development_5th_p05FC1_shared_5thdown_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_5thup_melp_shared <- read.table('development/development_5th_p05FC1_shared_5thup_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_5thdown_melp_shared <- read.table('development/development_5th_p05FC1_shared_5thdown_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


dev_D1up_erato_shared <- read.table('development/development_D1_p05FC1_shared_D1up_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D1down_erato_shared <- read.table('development/development_D1_p05FC1_shared_D1down_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_D1up_melp_shared <- read.table('development/development_D1_p05FC1_shared_D1up_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D1down_melp_shared <- read.table('development/development_D1_p05FC1_shared_D1down_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]


dev_D2up_erato_shared <- read.table('development/development_D2_p05FC1_shared_D2up_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D2down_erato_shared <- read.table('development/development_D2_p05FC1_shared_D2down_eratoCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]

dev_D2up_melp_shared <- read.table('development/development_D2_p05FC1_shared_D2up_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]
dev_D2down_melp_shared <- read.table('development/development_D2_p05FC1_shared_D2down_melpCoords_HOMER.txt', sep='\t', skip = 1)[,c(1:11)]



dev_5thup_erato2 <- separate(data = dev_5thup_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thup_erato2$annot <- gsub('ID=', '', dev_5thup_erato2$annot)
dev_5thup_erato2$annot <- gsub('-E\\d+', '', dev_5thup_erato2$annot)
dev_5thup_erato_FC <- diff_erato_5th_D1D2[dev_5thup_erato2$annot,]
dev_5thup_erato_FC$group <- '5thopen'
dev_5thup_erato_FC$species <- 'erato'
dev_5thup_erato_FC$cat <- 'unique'
dev_5thup_erato_FC$stage <- '5th'

dev_5thdown_erato2 <- separate(data = dev_5thdown_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thdown_erato2$annot <- gsub('ID=', '', dev_5thdown_erato2$annot)
dev_5thdown_erato2$annot <- gsub('-E\\d+', '', dev_5thdown_erato2$annot)
dev_5thdown_erato_FC <- diff_erato_5th_D1D2[dev_5thdown_erato2$annot,]
dev_5thdown_erato_FC$group <- '5thclosed'
dev_5thdown_erato_FC$species <- 'erato'
dev_5thdown_erato_FC$cat <- 'unique'
dev_5thdown_erato_FC$stage <- '5th'

dev_5thup_melp2 <- separate(data = dev_5thup_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thup_melp2$annot <- gsub('ID=', '', dev_5thup_melp2$annot)
dev_5thup_melp2$annot <- gsub('-E\\d+', '', dev_5thup_melp2$annot)
dev_5thup_melp_FC <- diff_melp_5th_D1D2[dev_5thup_melp2$annot,]
dev_5thup_melp_FC$group <- '5thopen'
dev_5thup_melp_FC$species <- 'melp'
dev_5thup_melp_FC$cat <- 'unique'
dev_5thup_melp_FC$stage <- '5th'

dev_5thdown_melp2 <- separate(data = dev_5thdown_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thdown_melp2$annot <- gsub('ID=', '', dev_5thdown_melp2$annot)
dev_5thdown_melp2$annot <- gsub('-E\\d+', '', dev_5thdown_melp2$annot)
dev_5thdown_melp_FC <- diff_melp_5th_D1D2[dev_5thdown_melp2$annot,]
dev_5thdown_melp_FC$group <- '5thclosed'
dev_5thdown_melp_FC$species <- 'melp'
dev_5thdown_melp_FC$cat <- 'unique'
dev_5thdown_melp_FC$stage <- '5th'


dev_D1up_erato2 <- separate(data = dev_D1up_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1up_erato2$annot <- gsub('ID=', '', dev_D1up_erato2$annot)
dev_D1up_erato2$annot <- gsub('-E\\d+', '', dev_D1up_erato2$annot)
dev_D1up_erato_FC <- diff_erato_D1_5thD2[dev_D1up_erato2$annot,]
dev_D1up_erato_FC$group <- 'D1open'
dev_D1up_erato_FC$species <- 'erato'
dev_D1up_erato_FC$cat <- 'unique'
dev_D1up_erato_FC$stage <- 'D1'

dev_D1down_erato2 <- separate(data = dev_D1down_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1down_erato2$annot <- gsub('ID=', '', dev_D1down_erato2$annot)
dev_D1down_erato2$annot <- gsub('-E\\d+', '', dev_D1down_erato2$annot)
dev_D1down_erato_FC <- diff_erato_D1_5thD2[dev_D1down_erato2$annot,]
dev_D1down_erato_FC$group <- 'D1closed'
dev_D1down_erato_FC$species <- 'erato'
dev_D1down_erato_FC$cat <- 'unique'
dev_D1down_erato_FC$stage <- 'D1'

dev_D1up_melp2 <- separate(data = dev_D1up_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1up_melp2$annot <- gsub('ID=', '', dev_D1up_melp2$annot)
dev_D1up_melp2$annot <- gsub('-E\\d+', '', dev_D1up_melp2$annot)
dev_D1up_melp_FC <- diff_melp_D1_5thD2[dev_D1up_melp2$annot,]
dev_D1up_melp_FC$group <- 'D1open'
dev_D1up_melp_FC$species <- 'melp'
dev_D1up_melp_FC$cat <- 'unique'
dev_D1up_melp_FC$stage <- 'D1'

dev_D1down_melp2 <- separate(data = dev_D1down_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1down_melp2$annot <- gsub('ID=', '', dev_D1down_melp2$annot)
dev_D1down_melp2$annot <- gsub('-E\\d+', '', dev_D1down_melp2$annot)
dev_D1down_melp_FC <- diff_melp_D1_5thD2[dev_D1down_melp2$annot,]
dev_D1down_melp_FC$group <- 'D1closed'
dev_D1down_melp_FC$species <- 'melp'
dev_D1down_melp_FC$cat <- 'unique'
dev_D1down_melp_FC$stage <- 'D1'


dev_D2up_erato2 <- separate(data = dev_D2up_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2up_erato2$annot <- gsub('ID=', '', dev_D2up_erato2$annot)
dev_D2up_erato2$annot <- gsub('-E\\d+', '', dev_D2up_erato2$annot)
dev_D2up_erato_FC <- diff_erato_D2_5thD1[dev_D2up_erato2$annot,]
dev_D2up_erato_FC$group <- 'D2open'
dev_D2up_erato_FC$species <- 'erato'
dev_D2up_erato_FC$cat <- 'unique'
dev_D2up_erato_FC$stage <- 'D2'

dev_D2down_erato2 <- separate(data = dev_D2down_erato, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2down_erato2$annot <- gsub('ID=', '', dev_D2down_erato2$annot)
dev_D2down_erato2$annot <- gsub('-E\\d+', '', dev_D2down_erato2$annot)
dev_D2down_erato_FC <- diff_erato_D2_5thD1[dev_D2down_erato2$annot,]
dev_D2down_erato_FC$group <- 'D2closed'
dev_D2down_erato_FC$species <- 'erato'
dev_D2down_erato_FC$cat <- 'unique'
dev_D2down_erato_FC$stage <- 'D2'

dev_D2up_melp2 <- separate(data = dev_D2up_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2up_melp2$annot <- gsub('ID=', '', dev_D2up_melp2$annot)
dev_D2up_melp2$annot <- gsub('-E\\d+', '', dev_D2up_melp2$annot)
dev_D2up_melp_FC <- diff_melp_D2_5thD1[dev_D2up_melp2$annot,]
dev_D2up_melp_FC$group <- 'D2open'
dev_D2up_melp_FC$species <- 'melp'
dev_D2up_melp_FC$cat <- 'unique'
dev_D2up_melp_FC$stage <- 'D2'

dev_D2down_melp2 <- separate(data = dev_D2down_melp, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2down_melp2$annot <- gsub('ID=', '', dev_D2down_melp2$annot)
dev_D2down_melp2$annot <- gsub('-E\\d+', '', dev_D2down_melp2$annot)
dev_D2down_melp_FC <- diff_melp_D2_5thD1[dev_D2down_melp2$annot,]
dev_D2down_melp_FC$group <- 'D2closed'
dev_D2down_melp_FC$species <- 'melp'
dev_D2down_melp_FC$cat <- 'unique'
dev_D2down_melp_FC$stage <- 'D2'

###



dev_5thup_erato_shared2 <- separate(data = dev_5thup_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thup_erato_shared2$annot <- gsub('ID=', '', dev_5thup_erato_shared2$annot)
dev_5thup_erato_shared2$annot <- gsub('-E\\d+', '', dev_5thup_erato_shared2$annot)
dev_5thup_erato_shared_FC <- diff_erato_5th_D1D2[dev_5thup_erato_shared2$annot,]
dev_5thup_erato_shared_FC$group <- '5thopen'
dev_5thup_erato_shared_FC$species <- 'erato'
dev_5thup_erato_shared_FC$cat <- 'shared'
dev_5thup_erato_shared_FC$stage <- '5th'

dev_5thdown_erato_shared2 <- separate(data = dev_5thdown_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thdown_erato_shared2$annot <- gsub('ID=', '', dev_5thdown_erato_shared2$annot)
dev_5thdown_erato_shared2$annot <- gsub('-E\\d+', '', dev_5thdown_erato_shared2$annot)
dev_5thdown_erato_shared_FC <- diff_erato_5th_D1D2[dev_5thdown_erato_shared2$annot,]
dev_5thdown_erato_shared_FC$group <- '5thclosed'
dev_5thdown_erato_shared_FC$species <- 'erato'
dev_5thdown_erato_shared_FC$cat <- 'shared'
dev_5thdown_erato_shared_FC$stage <- '5th'

dev_5thup_melp_shared2 <- separate(data = dev_5thup_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thup_melp_shared2$annot <- gsub('ID=', '', dev_5thup_melp_shared2$annot)
dev_5thup_melp_shared2$annot <- gsub('-E\\d+', '', dev_5thup_melp_shared2$annot)
dev_5thup_melp_shared_FC <- diff_melp_5th_D1D2[dev_5thup_melp_shared2$annot,]
dev_5thup_melp_shared_FC$group <- '5thopen'
dev_5thup_melp_shared_FC$species <- 'melp'
dev_5thup_melp_shared_FC$cat <- 'shared'
dev_5thup_melp_shared_FC$stage <- '5th'

dev_5thdown_melp_shared2 <- separate(data = dev_5thdown_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_5thdown_melp_shared2$annot <- gsub('ID=', '', dev_5thdown_melp_shared2$annot)
dev_5thdown_melp_shared2$annot <- gsub('-E\\d+', '', dev_5thdown_melp_shared2$annot)
dev_5thdown_melp_shared_FC <- diff_melp_5th_D1D2[dev_5thdown_melp_shared2$annot,]
dev_5thdown_melp_shared_FC$group <- '5thclosed'
dev_5thdown_melp_shared_FC$species <- 'melp'
dev_5thdown_melp_shared_FC$cat <- 'shared'
dev_5thdown_melp_shared_FC$stage <- '5th'


dev_D1up_erato_shared2 <- separate(data = dev_D1up_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1up_erato_shared2$annot <- gsub('ID=', '', dev_D1up_erato_shared2$annot)
dev_D1up_erato_shared2$annot <- gsub('-E\\d+', '', dev_D1up_erato_shared2$annot)
dev_D1up_erato_shared_FC <- diff_erato_D1_5thD2[dev_D1up_erato_shared2$annot,]
dev_D1up_erato_shared_FC$group <- 'D1open'
dev_D1up_erato_shared_FC$species <- 'erato'
dev_D1up_erato_shared_FC$cat <- 'shared'
dev_D1up_erato_shared_FC$stage <- 'D1'

dev_D1down_erato_shared2 <- separate(data = dev_D1down_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1down_erato_shared2$annot <- gsub('ID=', '', dev_D1down_erato_shared2$annot)
dev_D1down_erato_shared2$annot <- gsub('-E\\d+', '', dev_D1down_erato_shared2$annot)
dev_D1down_erato_shared_FC <- diff_erato_D1_5thD2[dev_D1down_erato_shared2$annot,]
dev_D1down_erato_shared_FC$group <- 'D1closed'
dev_D1down_erato_shared_FC$species <- 'erato'
dev_D1down_erato_shared_FC$cat <- 'shared'
dev_D1down_erato_shared_FC$stage <- 'D1'

dev_D1up_melp_shared2 <- separate(data = dev_D1up_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1up_melp_shared2$annot <- gsub('ID=', '', dev_D1up_melp_shared2$annot)
dev_D1up_melp_shared2$annot <- gsub('-E\\d+', '', dev_D1up_melp_shared2$annot)
dev_D1up_melp_shared_FC <- diff_melp_D1_5thD2[dev_D1up_melp_shared2$annot,]
dev_D1up_melp_shared_FC$group <- 'D1open'
dev_D1up_melp_shared_FC$species <- 'melp'
dev_D1up_melp_shared_FC$cat <- 'shared'
dev_D1up_melp_shared_FC$stage <- 'D1'

dev_D1down_melp_shared2 <- separate(data = dev_D1down_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D1down_melp_shared2$annot <- gsub('ID=', '', dev_D1down_melp_shared2$annot)
dev_D1down_melp_shared2$annot <- gsub('-E\\d+', '', dev_D1down_melp_shared2$annot)
dev_D1down_melp_shared_FC <- diff_melp_D1_5thD2[dev_D1down_melp_shared2$annot,]
dev_D1down_melp_shared_FC$group <- 'D1closed'
dev_D1down_melp_shared_FC$species <- 'melp'
dev_D1down_melp_shared_FC$cat <- 'shared'
dev_D1down_melp_shared_FC$stage <- 'D1'


dev_D2up_erato_shared2 <- separate(data = dev_D2up_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2up_erato_shared2$annot <- gsub('ID=', '', dev_D2up_erato_shared2$annot)
dev_D2up_erato_shared2$annot <- gsub('-E\\d+', '', dev_D2up_erato_shared2$annot)
dev_D2up_erato_shared_FC <- diff_erato_D2_5thD1[dev_D2up_erato_shared2$annot,]
dev_D2up_erato_shared_FC$group <- 'D2open'
dev_D2up_erato_shared_FC$species <- 'erato'
dev_D2up_erato_shared_FC$cat <- 'shared'
dev_D2up_erato_shared_FC$stage <- 'D2'

dev_D2down_erato_shared2 <- separate(data = dev_D2down_erato_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2down_erato_shared2$annot <- gsub('ID=', '', dev_D2down_erato_shared2$annot)
dev_D2down_erato_shared2$annot <- gsub('-E\\d+', '', dev_D2down_erato_shared2$annot)
dev_D2down_erato_shared_FC <- diff_erato_D2_5thD1[dev_D2down_erato_shared2$annot,]
dev_D2down_erato_shared_FC$group <- 'D2closed'
dev_D2down_erato_shared_FC$species <- 'erato'
dev_D2down_erato_shared_FC$cat <- 'shared'
dev_D2down_erato_shared_FC$stage <- 'D2'

dev_D2up_melp_shared2 <- separate(data = dev_D2up_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2up_melp_shared2$annot <- gsub('ID=', '', dev_D2up_melp_shared2$annot)
dev_D2up_melp_shared2$annot <- gsub('-E\\d+', '', dev_D2up_melp_shared2$annot)
dev_D2up_melp_shared_FC <- diff_melp_D2_5thD1[dev_D2up_melp_shared2$annot,]
dev_D2up_melp_shared_FC$group <- 'D2open'
dev_D2up_melp_shared_FC$species <- 'melp'
dev_D2up_melp_shared_FC$cat <- 'shared'
dev_D2up_melp_shared_FC$stage <- 'D2'

dev_D2down_melp_shared2 <- separate(data = dev_D2down_melp_shared, col = 'V11', into = c("annot", "annot1"), sep = ";(?=[^;]+$)")
dev_D2down_melp_shared2$annot <- gsub('ID=', '', dev_D2down_melp_shared2$annot)
dev_D2down_melp_shared2$annot <- gsub('-E\\d+', '', dev_D2down_melp_shared2$annot)
dev_D2down_melp_shared_FC <- diff_melp_D2_5thD1[dev_D2down_melp_shared2$annot,]
dev_D2down_melp_shared_FC$group <- 'D2closed'
dev_D2down_melp_shared_FC$species <- 'melp'
dev_D2down_melp_shared_FC$cat <- 'shared'
dev_D2down_melp_shared_FC$stage <- 'D2'


dev_5thup_erato_FC$log2FoldChange <- dev_5thup_erato_FC$log2FoldChange*-1
dev_5thdown_erato_FC$log2FoldChange <- dev_5thdown_erato_FC$log2FoldChange*-1
dev_5thup_melp_FC$log2FoldChange <- dev_5thup_melp_FC$log2FoldChange*-1
dev_5thdown_melp_FC$log2FoldChange <- dev_5thdown_melp_FC$log2FoldChange*-1
dev_5thup_erato_shared_FC$log2FoldChange <- dev_5thup_erato_shared_FC$log2FoldChange*-1
dev_5thdown_erato_shared_FC$log2FoldChange <- dev_5thdown_erato_shared_FC$log2FoldChange*-1
dev_5thup_melp_shared_FC$log2FoldChange <- dev_5thup_melp_shared_FC$log2FoldChange*-1
dev_5thdown_melp_shared_FC$log2FoldChange <- dev_5thdown_melp_shared_FC$log2FoldChange*-1


FC_all <- rbind(dev_5thup_erato_FC, dev_5thdown_erato_FC, dev_D1up_erato_FC, dev_D1down_erato_FC, dev_D2up_erato_FC, dev_D2down_erato_FC,
                dev_5thup_melp_FC, dev_5thdown_melp_FC, dev_D1up_melp_FC, dev_D1down_melp_FC, dev_D2up_melp_FC, dev_D2down_melp_FC,
                dev_5thup_erato_shared_FC, dev_5thdown_erato_shared_FC, dev_D1up_erato_shared_FC, dev_D1down_erato_shared_FC, dev_D2up_erato_shared_FC, dev_D2down_erato_shared_FC,
                dev_5thup_melp_shared_FC, dev_5thdown_melp_shared_FC, dev_D1up_melp_shared_FC, dev_D1down_melp_shared_FC, dev_D2up_melp_shared_FC, dev_D2down_melp_shared_FC)
FC_all[is.na(FC_all)] <- 0


library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)

ggplot(FC_all, aes(x=group, y=log2FoldChange, group=interaction(cat, group, stage), fill=group)) +
  geom_boxplot(outlier.size=0) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.5, alpha=0.2, width=0.1) +
  theme(
    plot.title = element_text(size=11)
  ) +facet_grid(~cat*species) + ylim(-10,10)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# violin
# FC_all %>%
#   ggplot(aes(x=group, y=log2FoldChange, fill=group, group=interaction(cat, group, stage))) +
#   geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
#   scale_fill_viridis(discrete=T, name="") +
#   geom_jitter(color="black", size=0.5, alpha=0.2, width=0.1) +
#   theme(
#     plot.title = element_text(size=11)
#   ) +facet_grid(~cat*species) + ylim(-10,10)

t.test(dev_5thdown_erato_FC$log2FoldChange, dev_5thup_erato_FC$log2FoldChange)
t.test(dev_D1down_erato_FC$log2FoldChange, dev_D1up_erato_FC$log2FoldChange)
t.test(dev_D2down_erato_FC$log2FoldChange, dev_D2up_erato_FC$log2FoldChange)

t.test(dev_5thdown_melp_FC$log2FoldChange, dev_5thup_melp_FC$log2FoldChange)
t.test(dev_D1down_melp_FC$log2FoldChange, dev_D1up_melp_FC$log2FoldChange)
t.test(dev_D2down_melp_FC$log2FoldChange, dev_D2up_melp_FC$log2FoldChange)

t.test(dev_5thdown_erato_shared_FC$log2FoldChange, dev_5thup_erato_shared_FC$log2FoldChange)
t.test(dev_D1down_erato_shared_FC$log2FoldChange, dev_D1up_erato_shared_FC$log2FoldChange)
t.test(dev_D2down_erato_shared_FC$log2FoldChange, dev_D2up_erato_shared_FC$log2FoldChange)

t.test(dev_5thdown_melp_shared_FC$log2FoldChange, dev_5thup_melp_shared_FC$log2FoldChange)
t.test(dev_D1down_melp_shared_FC$log2FoldChange, dev_D1up_melp_shared_FC$log2FoldChange)
t.test(dev_D2down_melp_shared_FC$log2FoldChange, dev_D2up_melp_shared_FC$log2FoldChange)

t.test(dev_5thdown_erato_FC$log2FoldChange, dev_5thdown_erato_shared_FC$log2FoldChange)
t.test(dev_5thup_erato_FC$log2FoldChange, dev_5thup_erato_shared_FC$log2FoldChange)
