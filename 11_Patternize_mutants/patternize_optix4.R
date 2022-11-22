library(patternize)
library(RColorBrewer)

IDList_D <- c("Optix4-Mutant-Demophoon-XV9-279-D-male",
              "Optix4-Mutant-Demophoon-XVI3-294-D-male",
              "Optix4-Mutant-Demophoon-XIII6-268-D-Female",
              "Optix4-Mutant-Demophoon-XV7-288-D-Female",
              "Optix4-Mutant-Demophoon-XV10-287-D-Female",
              "Optix4-Mutant-Demophoon-XVI2-303-D-Female",
              "Optix4-Mutant-Demophoon-XVII1-301-D-Female",
              "Optix4-Mutant-Demophoon-XVI7-293-D",
              "Optix4-Mutant-Demophoon-XVI9-296-D-Female",
              "Optix4-Mutant-Demophoon-XVII3-292-D",
              "Optix4-Mutant-Demophoon-XVII4-297-D-Female",
              "Optix4-Mutant-Demophoon-XVIII4-299-D-Female",
              "Optix4-Mutant-Demophoon-XVIII14-300-D-Male",
              "Optix4-Mutant-Demophoon-XX2-295-D-Female",
              "Optix4-Mutant-Demophoon-XX5-314-D-Male",
              
              "Optix4-H.e.demophoon-DQ5-Mutant-Dorsal-M",
              "Optix4-H.e.demophoon-DR9-Mutant-Dorsal-M",
              "Optix4-H.e.demophoon-DR11-Mutant-Dorsal-F",
              "Optix4-H.e.demophoon-DR13-Mutant-Dorsal-F",
              "Optix4-H.e.demophoon-DR18-Mutant-Dorsal-M",
              "Optix4-H.e.demophoon-DS7-Mutant-Dorsal-F",
              "Optix4-H.e.demophoon-DS10-Mutant-Dorsal-F")
              # "Optix4-Rosina-XIX2-59-D-male",
              # "Optix4-Rosina-XXI1-D")

IDList_V <- c("Optix4-Mutant-Demophoon-XV9-279-V-male",
             "Optix4-Mutant-Demophoon-XVI3-294-V-male",
             "Optix4-Mutant-Demophoon-XIII6-268-V-Female",
             "Optix4-Mutant-Demophoon-XV7-288-V-Female",
             "Optix4-Mutant-Demophoon-XV10-287-V-Female",
             "Optix4-Mutant-Demophoon-XVI2-303-V-Female",
             "Optix4-Mutant-Demophoon-XVII1-301-V-Female",
             "Optix4-Mutant-Demophoon-XVI7-293-V",
             "Optix4-Mutant-Demophoon-XVI9-296-V-Female",
             "Optix4-Mutant-Demophoon-XVII3-292-V",
             "Optix4-Mutant-Demophoon-XVII4-297-V-Female",
             "Optix4-Mutant-Demophoon-XVIII4-299-V-Female",
             "Optix4-Mutant-Demophoon-XVIII14-300-V-Male",
             "Optix4-Mutant-Demophoon-XX2-295-V-Female",
             "Optix4-Mutant-Demophoon-XX5-314-V-Male")
             # "Optix4-Rosina-XIX2-59-V-male",
             # "Optix4-Rosina-XXI1-V")


prepath <- "Optix4-Mutants" 
extension <- '_landmarks_R.txt' 
landmarkList_D_R <- makeList(IDList_D, 'landmark', prepath, extension)
landmarkList_V_R <- makeList(IDList_V, 'landmark', prepath, extension)

extension <- '_landmarks_L.txt' 
landmarkList_D_L <- makeList(IDList_D, 'landmark', prepath, extension)
landmarkList_V_L <- makeList(IDList_V, 'landmark', prepath, extension)



prepath <-  "Optix4-Mutants" 
extension <- '.jpg' 
imageList_D <- makeList(IDList_D, 'image', prepath, extension)
imageList_V <- makeList(IDList_V, 'image', prepath, extension)


# load in target image
target <- as.matrix(read.table('Cartoons2/BC0004_landmarks_LFW.txt',h = F))

outline_BC0004 <- read.table('Cartoons2/BC0004_outline.txt', h= F)
lines_BC0004 <- list.files(path="Cartoons2", pattern='BC0004_vein', full.names = T)

image_BC0004 <- list()
image_BC0004[['BC0004']] <- raster::stack('Cartoons2/BC0004-D.JPG')


# sampleRGB(imageList_D[[1]]) 

RGB <- c(208,43,41) # red
rasterList_lanRGB_mutant_red_D_R <- patLanRGB(imageList_D, landmarkList_D_R, RGB, transformRef = target, resampleFactor = 2,
                                      colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_red_D_R, file = 'optix4_rda/rasterList_lanRGB_mutant_red_D_R.rda')


rasterList_lanRGB_mutant_red_D_L <- patLanRGB(imageList_D, landmarkList_D_L, RGB, transformRef = target, resampleFactor = 3,
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_red_D_L, file = 'optix4_rda/rasterList_lanRGB_mutant_red_D_L.rda')


# sampleRGB(imageList_V[[1]]) 

RGB <- c(208,93,59) # red 
rasterList_lanRGB_mutant_red_V_R <- patLanRGB(imageList_V, landmarkList_V_R, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_red_V_R, file = 'optix4_rda/rasterList_lanRGB_mutant_red_V_R.rda')


rasterList_lanRGB_mutant_red_V_L <- patLanRGB(imageList_V, landmarkList_V_L, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_red_V_L, file = 'optix4_rda/rasterList_lanRGB_mutant_red_V_L.rda')


# sampleRGB(imageList_D[[1]])

RGB <- c(230,216,141) # yellow
rasterList_lanRGB_mutant_yellow_D_R <- patLanRGB(imageList_D, landmarkList_D_R, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_yellow_D_R, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_D_R.rda')


rasterList_lanRGB_mutant_yellow_D_L <- patLanRGB(imageList_D, landmarkList_D_L, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_yellow_D_L, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_D_L.rda')


# sampleRGB(imageList_V[[3]])

RGB <- c(195,185,123) # yellow
rasterList_lanRGB_mutant_yellow_V_R <- patLanRGB(imageList_V, landmarkList_V_R, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_yellow_V_R, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_V_R.rda')


rasterList_lanRGB_mutant_yellow_V_L <- patLanRGB(imageList_V, landmarkList_V_L, RGB, transformRef = target, resampleFactor = 3, 
                                              colOffset = 0.15, crop = TRUE, res = 200, adjustCoords = TRUE, plot = 'stack', cropOffset = c(1,1,1,1))

save(rasterList_lanRGB_mutant_yellow_V_L, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_V_L.rda')


###



summedRaster_lanRGB_mutant_red_D_R <- sumRaster(rasterList_lanRGB_mutant_red_D_R, IDList_D, type = 'RGB')
summedRaster_lanRGB_mutant_red_D_L <- sumRaster(rasterList_lanRGB_mutant_red_D_L, IDList_D, type = 'RGB')

summedRaster_lanRGB_mutant_red_V_R <- sumRaster(rasterList_lanRGB_mutant_red_V_R, IDList_V, type = 'RGB')
summedRaster_lanRGB_mutant_red_V_L <- sumRaster(rasterList_lanRGB_mutant_red_V_L, IDList_V, type = 'RGB')

summedRaster_lanRGB_mutant_yellow_D_R <- sumRaster(rasterList_lanRGB_mutant_yellow_D_R, IDList_D, type = 'RGB')
summedRaster_lanRGB_mutant_yellow_D_L <- sumRaster(rasterList_lanRGB_mutant_yellow_D_L, IDList_D, type = 'RGB')

summedRaster_lanRGB_mutant_yellow_V_R <- sumRaster(rasterList_lanRGB_mutant_yellow_V_R, IDList_V, type = 'RGB')
summedRaster_lanRGB_mutant_yellow_V_L <- sumRaster(rasterList_lanRGB_mutant_yellow_V_L, IDList_V, type = 'RGB')

setMask(summedRaster_lanRGB_mutant_red_D_R, IDList_D, filename = 'mask_mutant_red_D_R.txt')
setMask(summedRaster_lanRGB_mutant_red_D_L, IDList_D, filename = 'mask_mutant_red_D_L.txt')

setMask(summedRaster_lanRGB_mutant_red_V_R, IDList_V, filename = 'mask_mutant_red_V_R.txt')
setMask(summedRaster_lanRGB_mutant_red_V_L, IDList_V, filename = 'mask_mutant_red_V_L.txt')

setMask(summedRaster_lanRGB_mutant_yellow_D_R, IDList_D, filename = 'mask_mutant_yellow_D_R.txt')
setMask(summedRaster_lanRGB_mutant_yellow_D_L, IDList_D, filename = 'mask_mutant_yellow_D_L.txt')

setMask(summedRaster_lanRGB_mutant_yellow_V_R, IDList_V, filename = 'mask_mutant_yellow_V_R.txt')
setMask(summedRaster_lanRGB_mutant_yellow_V_L, IDList_V, filename = 'mask_mutant_yellow_V_L.txt')

mask_mutant_red_V_R <- read.table("mask_mutant_red_V_R.txt", h=T)
mask_mutant_red_V_L <- read.table("mask_mutant_red_V_L.txt", h=T)

mask_mutant_yellow_D_R <- read.table("mask_mutant_yellow_D_R.txt", h=T)
mask_mutant_yellow_D_L <- read.table("mask_mutant_yellow_D_L.txt", h=T)

mask_mutant_yellow_V_R <- read.table("mask_mutant_yellow_V_R.txt", h=T)
mask_mutant_yellow_V_L <- read.table("mask_mutant_yellow_V_L.txt", h=T)

rasterList_lanRGB_mutant_red_D_R_M <- rasterList_lanRGB_mutant_red_D_R
rasterList_lanRGB_mutant_red_D_L_M <- rasterList_lanRGB_mutant_red_D_L

rasterList_lanRGB_mutant_red_V_R_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_red_V_R)){
ID <- names(rasterList_lanRGB_mutant_red_V_R)[[e]]
rasterList_lanRGB_mutant_red_V_R_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_red_V_R[[ID]], IDlist = IDList_V, mask_mutant_red_V_R,
                                            refShape = 'target', imageList = imageList_V)
}

rasterList_lanRGB_mutant_red_V_L_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_red_V_L)){
  ID <- names(rasterList_lanRGB_mutant_red_V_L)[[e]]
  rasterList_lanRGB_mutant_red_V_L_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_red_V_L[[ID]], IDlist = IDList_V, mask_mutant_red_V_L,
                                                          refShape = 'target', imageList = imageList_V)
}

rasterList_lanRGB_mutant_yellow_D_R_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_yellow_D_R)){
  ID <- names(rasterList_lanRGB_mutant_yellow_D_R)[[e]]
  rasterList_lanRGB_mutant_yellow_D_R_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_yellow_D_R[[ID]], IDlist = IDList_D, mask_mutant_yellow_D_R,
                                                          refShape = 'target', imageList = imageList_D)
}

rasterList_lanRGB_mutant_yellow_D_L_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_yellow_D_L)){
  ID <- names(rasterList_lanRGB_mutant_yellow_D_L)[[e]]
  rasterList_lanRGB_mutant_yellow_D_L_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_yellow_D_L[[ID]], IDlist = IDList_D, mask_mutant_yellow_D_L,
                                                          refShape = 'target', imageList = imageList_D)
}

rasterList_lanRGB_mutant_yellow_V_R_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_yellow_V_R)){
  ID <- names(rasterList_lanRGB_mutant_yellow_V_R)[[e]]
  rasterList_lanRGB_mutant_yellow_V_R_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_yellow_V_R[[ID]], IDlist = IDList_V, mask_mutant_yellow_V_R,
                                                          refShape = 'target', imageList = imageList_V)
}

rasterList_lanRGB_mutant_yellow_V_L_M <-list()
for(e in 1:length(rasterList_lanRGB_mutant_yellow_V_L)){
  ID <- names(rasterList_lanRGB_mutant_yellow_V_L)[[e]]
  rasterList_lanRGB_mutant_yellow_V_L_M[[ID]] <- maskOutline(rasterList_lanRGB_mutant_yellow_V_L[[ID]], IDlist = IDList_V, mask_mutant_yellow_V_L,
                                                          refShape = 'target', imageList = imageList_V)
}

save(rasterList_lanRGB_mutant_red_D_R_M, file = 'optix4_rda/rasterList_lanRGB_mutant_red_D_R_M.rda')
save(rasterList_lanRGB_mutant_red_D_L_M, file = 'optix4_rda/rasterList_lanRGB_mutant_red_D_L_M.rda')

save(rasterList_lanRGB_mutant_red_V_R_M, file = 'optix4_rda/rasterList_lanRGB_mutant_red_V_R_M.rda')
save(rasterList_lanRGB_mutant_red_V_L_M, file = 'optix4_rda/rasterList_lanRGB_mutant_red_V_L_M.rda')

save(rasterList_lanRGB_mutant_yellow_D_R_M, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_D_R_M.rda')
save(rasterList_lanRGB_mutant_yellow_D_L_M, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_D_L_M.rda')

save(rasterList_lanRGB_mutant_yellow_V_R_M, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_V_R_M.rda')
save(rasterList_lanRGB_mutant_yellow_V_L_M, file = 'optix4_rda/rasterList_lanRGB_mutant_yellow_V_L_M.rda')

load('optix4_rda/rasterList_lanRGB_mutant_red_D_R_M.rda')
load('optix4_rda/rasterList_lanRGB_mutant_red_D_L_M.rda')

load('optix4_rda/rasterList_lanRGB_mutant_red_V_R_M.rda')
load('optix4_rda/rasterList_lanRGB_mutant_red_V_L_M.rda')

load('optix4_rda/rasterList_lanRGB_mutant_yellow_D_R_M.rda')
load('optix4_rda/rasterList_lanRGB_mutant_yellow_D_L_M.rda')

load('optix4_rda/rasterList_lanRGB_mutant_yellow_V_R_M.rda')
load('optix4_rda/rasterList_lanRGB_mutant_yellow_V_L_M.rda')

load('Demophoon-wt-cagephotos/rasterList_lanRGB_WT.rda')

summedRaster_lanRGB_mutant_red_D_R_M <- sumRaster(rasterList_lanRGB_mutant_red_D_R_M, IDList_D, type = 'RGB')
summedRaster_lanRGB_mutant_red_D_L_M <- sumRaster(rasterList_lanRGB_mutant_red_D_L_M, IDList_D, type = 'RGB')

summedRaster_lanRGB_mutant_red_D_M <- summedRaster_lanRGB_mutant_red_D_R_M + summedRaster_lanRGB_mutant_red_D_L_M

summedRaster_lanRGB_mutant_red_V_R_M <- sumRaster(rasterList_lanRGB_mutant_red_V_R_M, IDList_V, type = 'RGB')
summedRaster_lanRGB_mutant_red_V_L_M <- sumRaster(rasterList_lanRGB_mutant_red_V_L_M, IDList_V, type = 'RGB')

summedRaster_lanRGB_mutant_red_V_M <- summedRaster_lanRGB_mutant_red_V_R_M + summedRaster_lanRGB_mutant_red_V_L_M

summedRaster_lanRGB_mutant_yellow_D_R_M <- sumRaster(rasterList_lanRGB_mutant_yellow_D_R_M, IDList_D, type = 'RGB')
summedRaster_lanRGB_mutant_yellow_D_L_M <- sumRaster(rasterList_lanRGB_mutant_yellow_D_L_M, IDList_D, type = 'RGB')

summedRaster_lanRGB_mutant_yellow_D_M <- summedRaster_lanRGB_mutant_yellow_D_R_M + summedRaster_lanRGB_mutant_yellow_D_L_M

summedRaster_lanRGB_mutant_yellow_V_R_M <- sumRaster(rasterList_lanRGB_mutant_yellow_V_R_M, IDList_V, type = 'RGB')
summedRaster_lanRGB_mutant_yellow_V_L_M <- sumRaster(rasterList_lanRGB_mutant_yellow_V_L_M, IDList_V, type = 'RGB')

summedRaster_lanRGB_mutant_yellow_V_M <- summedRaster_lanRGB_mutant_yellow_V_R_M + summedRaster_lanRGB_mutant_yellow_V_L_M


IDList_WT_D <- names(rasterList_lanRGB_WT)

prepath <-  "Demophoon-wt-cagephotos" 
extension <- '.jpg' 
imageList_WT_D <- makeList(IDList_WT_D, 'image', prepath, extension)


summedRaster_lanRGB_WT_red_D <- sumRaster(rasterList_lanRGB_WT, IDList_WT_D, type = 'RGB')

setMask(summedRaster_lanRGB_WT_red_D, IDList_WT_D, filename = 'mask_WT_red_D.txt')

mask_mutant_red_D <- read.table("mask_WT_red_D.txt", h=T)

rasterList_lanRGB_WT_M <-list()
for(e in 1:length(rasterList_lanRGB_WT)){
  ID <- names(rasterList_lanRGB_WT)[[e]]
  rasterList_lanRGB_WT_M[[ID]] <- maskOutline(rasterList_lanRGB_WT[[ID]], IDlist = IDList_WT_D, mask_mutant_red_D,
                                                             refShape = 'target', imageList = imageList_WT_D)
}
summedRaster_lanRGB_WT_red_D_M <- sumRaster(rasterList_lanRGB_WT_M, IDList_WT_D, type = 'RGB')

colfunc <- colorRampPalette(brewer.pal(9,'Blues'))(100)[100:1]

plotHeat(summedRaster_lanRGB_WT_red_D_M, IDList_WT_D, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "WT red dorsal")

###

colfunc <- colorRampPalette(brewer.pal(9,'Blues'))(100)[100:1]

par(mfrow=c(1,4))
colfunc <- c('black','red')
plotHeat(summedRaster_lanRGB_mutant_red_D_M, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "mutant red dorsal")

plotHeat(summedRaster_lanRGB_mutant_red_V_M, c(IDList_V, IDList_V), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = TRUE, main = "mutant red ventral")

colfunc <- c('black', 'yellow')
plotHeat(summedRaster_lanRGB_mutant_yellow_D_M, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "mutant yellow dorsal")

plotHeat(summedRaster_lanRGB_mutant_yellow_V_M, c(IDList_V, IDList_V), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = TRUE, main = "mutant yellow ventral")

summedRaster_lanRGB_mutant_yellow_D_M_comp <- summedRaster_lanRGB_mutant_yellow_D_M
summedRaster_lanRGB_mutant_yellow_V_M_comp <- summedRaster_lanRGB_mutant_yellow_V_M

summedRaster_lanRGB_mutant_yellow_D_M_comp[summedRaster_lanRGB_mutant_yellow_D_M_comp > 0] <- 1
summedRaster_lanRGB_mutant_yellow_V_M_comp[summedRaster_lanRGB_mutant_yellow_V_M_comp > 0] <- 1

par(mfrow=c(1,2))
colfunc <- c('black', 'yellow')
plotHeat(summedRaster_lanRGB_mutant_yellow_D_M_comp, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "mutant yellow dorsal", normalized = T)

plotHeat(summedRaster_lanRGB_mutant_yellow_V_M_comp, c(IDList_V, IDList_V), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = TRUE, main = "mutant yellow ventral", normalized = T)


summedRaster_lanRGB_mutant_red_D_M[summedRaster_lanRGB_mutant_red_D_M == 0] <- NA
summedRaster_lanRGB_mutant_red_D_M_invert <- 30-summedRaster_lanRGB_mutant_red_D_M
summedRaster_lanRGB_mutant_red_D_M_invert[summedRaster_lanRGB_mutant_red_D_M_invert > 0] <- 1

# raster::extent(summedRaster_lanRGB_WT_red_D_M) <- raster::extent(summedRaster_lanRGB_mutant_red_D_M)
# subtr_red_D <- summedRaster_lanRGB_WT_red_D_M/21 - summedRaster_lanRGB_mutant_red_D_M/30

par(mfrow=c(2,3), mar = c(1,1,1,2))

colfunc <- c('black','red')

plotHeat(summedRaster_lanRGB_WT_red_D_M, IDList_WT_D, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = F, main = "WT dorsal red frequency", zlim = c(0,1))

colfunc <- c('black','red')
plotHeat(summedRaster_lanRGB_mutant_red_D_M, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = T, main = "mutant dorsal red frequency", zlim = c(0,1))

colfunc <- c('black','red')
plotHeat(summedRaster_lanRGB_mutant_red_D_M_invert, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "mutant dorsal red reduction composit", zlim = c(0,1))

colfunc <- c('black', 'yellow')
plotHeat(summedRaster_lanRGB_mutant_yellow_D_M-30, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = F, main = "WT dorsal yellow frequency", zlim = c(0,0.3))

colfunc <- c('black', 'yellow')
plotHeat(summedRaster_lanRGB_mutant_yellow_D_M, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = T, main = "mutant dorsal yellow frequency", zlim = c(0,0.3))

colfunc <- c('gray35', 'yellow')
plotHeat(summedRaster_lanRGB_mutant_yellow_D_M_comp, c(IDList_D, IDList_D), plotCartoon = TRUE, refShape = 'target', outline = outline_BC0004, 
         lines = lines_BC0004, landList = landmarkList_WT, adjustCoords = TRUE, 
         imageList = image_BC0004, cartoonID = 'BC0004', cartoonFill = 'black', cartoonOrder = 'under', 
         colpalette = colfunc, legend = FALSE, main = "mutant dorsal yellow composit", normalized = T)
