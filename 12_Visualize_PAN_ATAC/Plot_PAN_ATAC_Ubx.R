library(rtracklayer)
library(valr)

# interval
start = 59107173 -20000
end = 60107173 

end-start

startP = 120000
endP = 360000

# startP = 240000
# endP = 300000

startPos <- start

minLemgth = 0

genomes <- c(6,5,3,1,4,2)
cols <- c('cornflowerblue','mediumspringgreen')
names <- c('H. e. dem', 'H. m. ros')

###
# ATAC data
###

range <- paste('pan:',start+startP,'-',start+endP,sep='')

erato_5th_brain <- import.bw("bw/brain_5th_H_erato.bw", which = GRanges(range))
erato_5th_FW <- import.bw("bw/FW_5th_H_erato.bw", which = GRanges(range))
erato_5th_HW <- import.bw("bw/HW_5th_H_erato.bw", which = GRanges(range))

erato_D1_FD <- import.bw("bw/FD_D1_H_erato.bw", which = GRanges(range))
erato_D1_FM <- import.bw("bw/FM_D1_H_erato.bw", which = GRanges(range))
erato_D1_FP <- import.bw("bw/FP_D1_H_erato.bw", which = GRanges(range))
erato_D1_HW <- import.bw("bw/HW_D1_H_erato.bw", which = GRanges(range))

erato_D2_FD <- import.bw("bw/FD_D2_H_erato.bw", which = GRanges(range))
erato_D2_FM <- import.bw("bw/FM_D2_H_erato.bw", which = GRanges(range))
erato_D2_FP <- import.bw("bw/FP_D2_H_erato.bw", which = GRanges(range))
erato_D2_HW <- import.bw("bw/HW_D2_H_erato.bw", which = GRanges(range))

# eratoFiles <- c("xx", "bw/FW_5th_H_erato.bw","bw/HW_5th_H_erato.bw", "bw/FD_D1_H_erato.bw","bw/FM_D1_H_erato.bw","bw/FP_D1_H_erato.bw","bw/HW_D1_H_erato.bw","bw/FD_D2_H_erato.bw","bw/FM_D2_H_erato.bw","bw/FP_D2_H_erato.bw","bw/HW_D2_H_erato.bw")

melp_5th_brain <- import.bw("bw/brain_5th_H_melp.bw", which = GRanges(range))
melp_5th_FW <- import.bw("bw/FW_5th_H_melp.bw", which = GRanges(range))
melp_5th_HW <- import.bw("bw/HW_5th_H_melp.bw", which = GRanges(range))

melp_D1_FD <- import.bw("bw/FD_D1_H_melp.bw", which = GRanges(range))
melp_D1_FM <- import.bw("bw/FM_D1_H_melp.bw", which = GRanges(range))
melp_D1_FP <- import.bw("bw/FP_D1_H_melp.bw", which = GRanges(range))
melp_D1_HW <- import.bw("bw/HW_D1_H_melp.bw", which = GRanges(range))

melp_D2_FD <- import.bw("bw/FD_D2_H_melp.bw", which = GRanges(range))
melp_D2_FM <- import.bw("bw/FM_D2_H_melp.bw", which = GRanges(range))
melp_D2_FP <- import.bw("bw/FP_D2_H_melp.bw", which = GRanges(range))
melp_D2_HW <- import.bw("bw/HW_D2_H_melp.bw", which = GRanges(range))

# melpFiles <- c("xx", "bw/FW_5th_H_melp.bw","bw/HW_5th_H_melp.bw", "bw/FD_D1_H_melp.bw","bw/FM_D1_H_melp.bw","bw/FP_D1_H_melp.bw","bw/HW_D1_H_melp.bw","bw/FD_D2_H_melp.bw","bw/FM_D2_H_melp.bw","bw/FP_D2_H_melp.bw","bw/HW_D2_H_melp.bw")
# 
# 
# erato_list <- list(erato_5th_brain, erato_5th_FW, erato_5th_HW, 
#                    erato_D1_FD, erato_D1_FM, erato_D1_FP, erato_D1_HW,
#                    erato_D2_FD, erato_D2_FM, erato_D2_FP, erato_D2_HW)
# 
# melp_list <- list(melp_5th_brain, melp_5th_FW, melp_5th_HW, 
#                   melp_D1_FD, melp_D1_FM, melp_D1_FP, melp_D1_HW,
#                   melp_D2_FD, melp_D2_FM, melp_D2_FP, melp_D2_HW)
# 
# names <- c("5th brain", "5th FW", "5th HW", 
#            "D1 FD", "D1 FM", "D1 FP", "D1 HW",
#            "D2 FD", "D2 FM", "D2 FP", "D2 HW")

eratoFiles <- c("xx", "bw/FW_5th_H_erato.bw","bw/HW_5th_H_erato.bw")

melpFiles <- c("xx", "bw/FW_5th_H_melp.bw","bw/HW_5th_H_melp.bw")


erato_list <- list(erato_5th_brain, erato_5th_FW, erato_5th_HW)

melp_list <- list(melp_5th_brain, melp_5th_FW, melp_5th_HW)

names <- c("5th brain", "5th FW", "5th HW")


match_ME <- read.table("MACS2/match_melp_to_erato.1bp.txt", h=F)
match_EM <- read.table("MACS2/match_erato_to_melp.1bp.txt", h=F)

match_ME_S <- subset(match_ME, match_ME$V2 > startPos+startP & match_ME$V3 < startPos+endP & match_ME$V4 != '.')[,c(4:6)]
match_EM_S <- subset(match_EM, match_EM$V2 > startPos+startP & match_EM$V3 < startPos+endP & match_EM$V4 != '.')[,c(4:6)]

colnames(match_ME_S) <- c("chrom","start","end")
colnames(match_EM_S) <- c("chrom","start","end")



###
# Diff expr tracks
###

instar5_FW_erato <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_erato_FWup_PAN.txt', h=F)
instar5_HW_erato <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_erato_HWup_PAN.txt', h=F)
instar5_FW_melp <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_melp_FWup_PAN.txt', h=F)
instar5_HW_melp <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_melp_HWup_PAN.txt', h=F)

D1_FW_erato <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_erato_FWup_PAN.txt', h=F)
D1_HW_erato <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_erato_HWup_PAN.txt', h=F)
D1_FW_melp <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_melp_FWup_PAN.txt', h=F)
D1_HW_melp <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_melp_HWup_PAN.txt', h=F)

D2_FW_erato <- read.table('FWHW_pan_coords/FWHW_D2_p05FC1_erato_FWup_PAN.txt', h=F)
D2_HW_erato <- read.table('FWHW_pan_coords/FWHW_D2_p05FC1_erato_HWup_PAN.txt', h=F)
D2_FW_melp <-NULL
D2_HW_melp <- read.table('FWHW_pan_coords/FWHW_D2_p05FC1_melp_HWup_PAN.txt', h=F)

instar5_FW_shared <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_shared_FWup_panCoords.bed', h=F)
instar5_HW_shared <- read.table('FWHW_pan_coords/FWHW_5th_p05FC1_shared_HWup_notStrict.txt', h=F)

D1_FW_shared <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_shared_FWup_panCoords.bed', h=F)
D1_HW_shared <- read.table('FWHW_pan_coords/FWHW_D1_p05FC1_shared_HWup_notStrict.txt', h=F)

D2_FW_shared <- NULL
D2_HW_shared <- read.table('FWHW_pan_coords/FWHW_D2_p05FC1_shared_HWup_notStrict.txt', h=F)

instar5_FW_erato_S <- subset(instar5_FW_erato, instar5_FW_erato$V2 > startPos+startP & instar5_FW_erato$V3 < startPos+endP )
instar5_HW_erato_S <- subset(instar5_HW_erato, instar5_HW_erato$V2 > startPos+startP & instar5_HW_erato$V3 < startPos+endP )

D1_FW_erato_S <- subset(D1_FW_erato, D1_FW_erato$V2 > startPos+startP & D1_FW_erato$V3 < startPos+endP )
D1_HW_erato_S <- subset(D1_HW_erato, D1_HW_erato$V2 > startPos+startP & D1_HW_erato$V3 < startPos+endP )

D2_FW_erato_S <- subset(D2_FW_erato, D2_FW_erato$V2 > startPos+startP & D2_FW_erato$V3 < startPos+endP )
D2_HW_erato_S <- subset(D2_HW_erato, D2_HW_erato$V2 > startPos+startP & D2_HW_erato$V3 < startPos+endP )

instar5_FW_melp_S <- subset(instar5_FW_melp, instar5_FW_melp$V2 > startPos+startP & instar5_FW_melp$V3 < startPos+endP )
instar5_HW_melp_S <- subset(instar5_HW_melp, instar5_HW_melp$V2 > startPos+startP & instar5_HW_melp$V3 < startPos+endP )

D1_FW_melp_S <- subset(D1_FW_melp, D1_FW_melp$V2 > startPos+startP & D1_FW_melp$V3 < startPos+endP )
D1_HW_melp_S <- subset(D1_HW_melp, D1_HW_melp$V2 > startPos+startP & D1_HW_melp$V3 < startPos+endP )

D2_FW_melp_S <- subset(D2_FW_melp, D2_FW_melp$V2 > startPos+startP & D2_FW_melp$V3 < startPos+endP )
D2_HW_melp_S <- subset(D2_HW_melp, D2_HW_melp$V2 > startPos+startP & D2_HW_melp$V3 < startPos+endP )

instar5_FW_shared_S <- subset(instar5_FW_shared, instar5_FW_shared$V2 > startPos+startP & instar5_FW_shared$V3 < startPos+endP )
instar5_HW_shared_S <- subset(instar5_HW_shared, instar5_HW_shared$V2 > startPos+startP & instar5_HW_shared$V3 < startPos+endP )

D1_FW_shared_S <- subset(D1_FW_shared, D1_FW_shared$V2 > startPos+startP & D1_FW_shared$V3 < startPos+endP )
D1_HW_shared_S <- subset(D1_HW_shared, D1_HW_shared$V2 > startPos+startP & D1_HW_shared$V3 < startPos+endP )

D2_FW_shared_S <- subset(D2_FW_shared, D2_FW_shared$V2 > startPos+startP & D2_FW_shared$V3 < startPos+endP )
D2_HW_shared_S <- subset(D2_HW_shared, D2_HW_shared$V2 > startPos+startP & D2_HW_shared$V3 < startPos+endP )

xxx <- c()
# erato_DF_list <- list(xxx,instar5_FW_erato_S,instar5_HW_erato_S,D1_FW_erato_S,D1_FW_erato_S,D1_FW_erato_S,D1_HW_erato_S,D2_FW_erato_S,D2_FW_erato_S,D2_FW_erato_S,D2_HW_erato_S)
# melp_DF_list <- list(xxx,instar5_FW_melp_S,instar5_HW_melp_S,D1_FW_melp_S,D1_FW_melp_S,D1_FW_melp_S,D1_HW_melp_S,D2_FW_melp_S,D2_FW_melp_S,D2_FW_melp_S,D2_HW_melp_S)
# 
# shared_DF_list <- list(xxx,instar5_FW_shared_S,instar5_HW_shared_S,D1_FW_shared_S,D1_FW_shared_S,D1_FW_shared_S,D1_HW_shared_S,D2_FW_shared_S,D2_FW_shared_S,D2_FW_shared_S,D2_HW_shared_S)

erato_DF_list <- list(xxx,instar5_FW_erato_S,instar5_HW_erato_S)
melp_DF_list <- list(xxx,instar5_FW_melp_S,instar5_HW_melp_S)

shared_DF_list <- list(xxx,instar5_FW_shared_S,instar5_HW_shared_S)




alphas <- c(0.3, 0.3, 0.3, 0.3, 0.3)

png('Ubx_PAN_postman_sub_11112022.png', width=5000, height=1000)

layout(matrix(c(1:9), nrow=9, byrow=TRUE), height = c(0.5,0.5,1,1,1,1,1,1,1))
# layout.show(n=25)

par(mar = c(0,4,0,4))

###
# optix
###

ANNOT_start <- read.table("gff/ubx_era_start_pan.txt", sep='\t', h=T)
ANNOT_end <- read.table("gff/ubx_era_end_pan.txt", sep='\t', h=T)

ANNOT <- as.data.frame(cbind(ANNOT_start[,2], ANNOT_end[,2]))
names(ANNOT) <- c("con_start", "con_end")

plot(NULL, xlim=c(startP,endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

for (g in 1:3){
  rect(ANNOT$con_start[g]-start, 0.2, ANNOT$con_end[g]-start, 0.8, col = "black", border = "black", lwd=8)
}
rect(min(ANNOT$con_start[c(1:3)])-start, 0.5, max(ANNOT$con_end[c(1:3)])-start, 0.5, col = "black", border = "black", lwd=8)

text(10000+max(ANNOT$con_end[c(1:3)]-start), 0.5, substitute(paste(italic(nn)), list(nn='Abd-A')), cex=5)

for (g in 4:6){
  rect(ANNOT$con_start[g]-start, 0.2, ANNOT$con_end[g]-start, 0.8, col = "black", border = "black", lwd=8)
}
rect(min(ANNOT$con_start[c(4:6)])-start, 0.5, max(ANNOT$con_end[c(4:6)])-start, 0.5, col = "black", border = "black", lwd=8)


text(10000+max(ANNOT$con_end-start), 0.5, substitute(paste(italic(nn)), list(nn='Ubx')), cex=5)

###
# PAN genome
###

blocks_erato <- read.table("blocks_1_min_2.bed", h=F)
blocks_melp <- read.table("blocks_2_min_1.bed", h=F)

colnames(blocks_erato) <- c('scaf', 'start', 'end')
colnames(blocks_melp) <- c('scaf', 'start', 'end')

plot(NULL, xlim=c(start+startP,start+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

top = 1
bot = 0

rect(start+startP, 0, start+endP, 1, col = "black", border = NA)

# genomes <- c(1,2)
# for(e in genomes){
  

unique1 <- subset(blocks_erato, (blocks_erato$end-blocks_erato$start >= minLemgth) & (blocks_erato$start > startPos+startP) & (blocks_erato$end < startPos+endP))
unique2 <- subset(blocks_melp, (blocks_melp$end-blocks_melp$start >= minLemgth) & (blocks_melp$start > startPos+startP) & (blocks_melp$end < startPos+endP))

for(i in 1:nrow(unique1)){
  rect(unique1$start[i], bot, unique1$end[i], top, col = cols[1], border = NA)
}
for(i in 1:nrow(unique2)){
  rect(unique2$start[i], bot, unique2$end[i], top, col = cols[2], border = NA)
}
# }




###
# ATAC
###

bot2 = 0
top2 = 200

for(e in 1:length(erato_list)){
  
  unique <- blocks_erato
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > startPos+startP) & (unique$end < startPos+endP))
  
  plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
  
  if(e == 3 | e == 4 | e == 5){
    rect(start+startP, bot2, start+endP, top2, border = NA, col = adjustcolor('red', alpha = 0.05))
  }
  for(i in 1:nrow(unique)){
    rect(unique$start[i], bot2, unique$end[i], top2, col = adjustcolor(cols[1], alpha = alphas[e]), border = NA)
  }
  
  par(new=T)
  plot(0.5*(start(erato_list[[e]]) + end(erato_list[[e]])), erato_list[[e]]$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = "black")
  
  if(any(e==c(3,4,5))){
    
    eratoFileName <- eratoFiles[e]
    for(n in 1:nrow(erato_DF_list[[e]])){
      range <- paste('pan:',erato_DF_list[[e]]$V2[n],'-',erato_DF_list[[e]]$V3[n],sep='')
      table <- import.bw(eratoFileName, which = GRanges(range))
      
      par(new=T)
      plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'purple')
    }
    
    # for(n in 1:nrow(shared_DF_list[[e]])){
    #   range <- paste('pan:',shared_DF_list[[e]]$V2[n],'-',shared_DF_list[[e]]$V3[n],sep='')
    #   table <- import.bw(eratoFileName, which = GRanges(range))
    #   
    #   par(new=T)
    #   plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'red')
    # }
    
    # intersect
    colnames(erato_DF_list[[e]]) <- c("chrom","start","end")
    
    res <- bed_intersect(erato_DF_list[[e]], match_ME_S)
    intersected <- dplyr::mutate(res, start = pmax(start.x, start.y),end = pmin(end.x, end.y))[,c(7,8)]
    
    for(n in 1:nrow(intersected)){
      range <- paste('pan:',intersected$start[n],'-',intersected$end[n],sep='')
      table <- import.bw(eratoFileName, which = GRanges(range))

      par(new=T)
      plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'red')
    }
  }
    
  # mtext(names[e], side = 1, cex=5, padj = -2, las = 1, adj=0)
  
}

# plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
# for(i in 1:nrow(match_EM_S)){
#   rect(match_EM_S$start[i], bot2, match_EM_S$end[i], top2, col = adjustcolor('black', alpha = 1), border = NA)
# }
# 
# plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)

bot2 = 0
top2 = 200

for(e in 1:length(melp_list)){
  
  unique <- blocks_melp
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > startPos+startP) & (unique$end < startPos+endP))
  
  plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
  if(e == 3 | e == 4 | e == 5){
    rect(start+startP, bot2, start+endP, top2, border = NA, col = adjustcolor('red', alpha = 0.05))
  }
  for(i in 1:nrow(unique)){
    rect(unique$start[i], bot2, unique$end[i], top2, col = adjustcolor(cols[2], alpha = alphas[e]), border = NA)
  }
  
  par(new=T)
  plot(0.5*(start(melp_list[[e]]) + end(melp_list[[e]])), melp_list[[e]]$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = "black")
  
  if(any(e==c(3,4,5))){
    
    melpFileName <- melpFiles[e]
    for(n in 1:nrow(melp_DF_list[[e]])){
      range <- paste('pan:',melp_DF_list[[e]]$V2[n],'-',melp_DF_list[[e]]$V3[n],sep='')
      table <- import.bw(melpFileName, which = GRanges(range))
      
      par(new=T)
      plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'purple')
    }
    
    # for(n in 1:nrow(shared_DF_list[[e]])){
    #   range <- paste('pan:',shared_DF_list[[e]]$V2[n],'-',shared_DF_list[[e]]$V3[n],sep='')
    #   table <- import.bw(melpFileName, which = GRanges(range))
    #   
    #   par(new=T)
    #   plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'red')
    # }
    
    # intersect
    colnames(melp_DF_list[[e]]) <- c("chrom","start","end")
    
    res <- bed_intersect(melp_DF_list[[e]], match_EM_S)
    intersected <- dplyr::mutate(res, start = pmax(start.x, start.y),end = pmin(end.x, end.y))[,c(7,8)]
    
    for(n in 1:nrow(intersected)){
      range <- paste('pan:',intersected$start[n],'-',intersected$end[n],sep='')
      table <- import.bw(melpFileName, which = GRanges(range))
      
      par(new=T)
      plot(0.5*(start(table) + end(table)), table$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none', col = 'red')
    }
  }
  
  # mtext(names[e], side = 1, cex=5, padj = -2, las = 1, adj=0)
  
} 

# plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
# for(i in 1:nrow(match_ME_S)){
#   rect(match_ME_S$start[i], bot2, match_ME_S$end[i], top2, col = adjustcolor('black', alpha = 1), border = NA)
# }
# 
# plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)

# points(59297502,-5, pch=19, col="red")
# points(59299121,-5, pch=19, col="red")
# points(59205514,-5, pch=19, col="red")
# points(59281076,-5, pch=19, col="red")
# 
# points(59231758,-5, pch=19, col="blue")
# points(59264426,-5, pch=19, col="blue")
# points(59284808,-5, pch=19, col="blue")
# points(59300794,-5, pch=19, col="blue")
# points(59318674,-5, pch=19, col="blue")
# points(59335601,-5, pch=19, col="blue")
# points(59345920,-5, pch=19, col="blue")
# points(59351671,-5, pch=19, col="blue")
# points(59355882,-5, pch=19, col="blue")
# points(59357343,-5, pch=19, col="blue")
# points(59409505,-5, pch=19, col="blue")
# points(59427668,-5, pch=19, col="blue")
# points(59430854,-5, pch=19, col="blue")
# points(59442871,-5, pch=19, col="blue")
# 
# points(59221523,-5, pch=19, col="green")
# points(59233126,-5, pch=19, col="green")
# points(59238885,-5, pch=19, col="green")
# points(59306718,-5, pch=19, col="green")
# points(59308425,-5, pch=19, col="green")
# points(59313272,-5, pch=19, col="green")
# points(59315415,-5, pch=19, col="green")
# points(59316893,-5, pch=19, col="green")
# points(59322289,-5, pch=19, col="green")
# points(59327118,-5, pch=19, col="green")
# points(59330089,-5, pch=19, col="green")
# points(59331060,-5, pch=19, col="green")
# points(59335038,-5, pch=19, col="green")
# points(59343180,-5, pch=19, col="green")
# points(59346975,-5, pch=19, col="green")
# points(59371285,-5, pch=19, col="green")
# points(59372053,-5, pch=19, col="green")
# points(59373317,-5, pch=19, col="green")
# points(59378495,-5, pch=19, col="green")
# points(59392651,-5, pch=19, col="green")
# points(59414998,-5, pch=19, col="green")
# points(59418998,-5, pch=19, col="green")


plot(NULL, xlim=c(start+startP,start+endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1, at = seq(start+startP,start+endP, by=10000), labels = NA, line =-10, lwd = 5, lwd.ticks = 5, tck = -0.2)
axis(1, at = seq(start+startP,start+endP, by=100000), labels = seq(0,(endP-startP)/1000000, by=0.1), line =-3, cex.axis = 6, lwd=0)


dev.off()

