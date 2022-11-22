
# sections

sectionIDY <- read.table("Section_unique_PAN_IDY.txt", h = T)

head(sectionIDY)
nrow(sectionIDY)

mean(sectionIDY$IDY)
sd(sectionIDY$IDY)

nrow(subset(sectionIDY, sectionIDY$IDY == 0))/(nrow(sectionIDY)-103)

# FW/HW

FWHW_shared <- read.table("FWHW_shared_PAN_IDY.txt", h = T)

nrow(FWHW_shared)

mean(FWHW_shared$IDY)
sd(FWHW_shared$IDY)


FWHW_unique <- read.table("FWHW_unique_PAN_IDY.txt", h = T)

FWHW_unique <- FWHW_unique[!FWHW_unique$start %in% FWHW_shared$start,]

nrow(FWHW_unique)

mean(FWHW_unique$IDY)
sd(FWHW_unique$IDY)

nrow(subset(FWHW_unique, FWHW_unique$IDY == 0))/(nrow(FWHW_unique))

