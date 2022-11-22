library(rtracklayer)


# interval
start = 380272889 
end = 382272889 

end-start

startP = 1100000#900000
endP = 1240000#1500000


minLemgth = 200

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

erato_list <- list(erato_5th_brain, erato_5th_FW, erato_5th_HW, 
                   erato_D1_FD, erato_D1_FM, erato_D1_FP, erato_D1_HW,
                   erato_D2_FD, erato_D2_FM, erato_D2_FP, erato_D2_HW)

melp_list <- list(melp_5th_brain, melp_5th_FW, melp_5th_HW, 
                  melp_D1_FD, melp_D1_FM, melp_D1_FP, melp_D1_HW,
                  melp_D2_FD, melp_D2_FM, melp_D2_FP, melp_D2_HW)

names <- c("5th brain", "5th FW", "5th HW", 
           "D1 FD", "D1 FM", "D1 FP", "D1 HW",
           "D2 FD", "D2 FM", "D2 FP", "D2 HW")

alphas <- c(0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 0.3, 0.3, 0.3, 0.3)

png('optix_PAN_postman_FW_melp.png', width=5000, height=5000)

layout(matrix(c(1:25), nrow=25, byrow=TRUE), height = c(0.5,0.5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5))
# layout.show(n=25)

par(mar = c(0,4,0,4))

###
# optix
###

ANNOT_start <- read.table("gff/optix_start_pan.txt", sep='\t', h=T)
ANNOT_end <- read.table("gff/optix_end_pan.txt", sep='\t', h=T)

ANNOT <- as.data.frame(cbind(ANNOT_start[,2], ANNOT_end[,2]))
names(ANNOT) <- c("con_start", "con_end")

plot(NULL, xlim=c(startP,endP), ylim = c(0,1.4), axes=FALSE, ann=FALSE)

for (g in 1:nrow(ANNOT)){
  rect(ANNOT$con_start[g]-start, 0.2, ANNOT$con_end[g]-start, 0.8, col = "black", border = "black", lwd=8)
}
rect(min(ANNOT$con_start)-start, 0.5, max(ANNOT$con_end)-start, 0.5, col = "black", border = "black", lwd=8)

text(10000+(min(ANNOT$con_start-start)+max(ANNOT$con_end-start))/2, 0.5, substitute(paste(italic(nn)), list(nn='optix')), cex=5)

rect(381439684-start, 1.1, 381450070-start, 1.3, col = "red", border = "red", lwd=8)
rect(381475963-start, 1.1, 381490465-start, 1.3, col = "yellow", border = "yellow", lwd=8)
rect(381499756-start, 1.1, 381501318-start, 1.3, col = "yellow", border = "yellow", lwd=8)
rect(381490488-start, 1.1, 381497917-start, 1.3, col = "gray35", border = "gray35", lwd=8)
rect(381501605-start, 1.1, 381509560-start, 1.3, col = "gray35", border = "gray35", lwd=8)

rect(381404732-start, 1.2, 381406466-start, 1.4, col = "yellow", border = "yellow", lwd=8)

Lewis_start <- read.table("Lewis_Hel_optix_peaks_erato_start_pan.txt", h=T, sep='\t')
Lewis_end <- read.table("Lewis_Hel_optix_peaks_erato_end_pan.txt", h=T, sep='\t')

Lewis <- as.data.frame(cbind(Lewis_start[,2], Lewis_end[,2]))
colnames(Lewis) <- c('start','end')
Lewis$y <- 0.8

par(new=T)
plot((Lewis$start-start+Lewis$end-start)/2, Lewis$y, xlim=c(startP,endP), ylim = c(0,1.3), axes=FALSE, ann=FALSE, col = 'blue', cex = 5, pch=19)

###
# PAN genome
###

blocks_erato <- read.table("blocks_unique_1_optix.txt", h=F)
blocks_melp <- read.table("blocks_unique_2_optix.txt", h=F)

plot(NULL, xlim=c(startP,endP), ylim = c(0,1), axes=FALSE, ann=FALSE)

top = 1
bot = 0

rect(startP, 0, endP, 1, col = "black", border = NA)

genomes <- c(1,2)
for(e in genomes){
  
  unique <- read.table(paste("blocks_unique_",e,"_optix.txt",sep=""), h = F)
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > startP) & (unique$end < endP))
  
  for(i in 1:nrow(unique)){
    rect(unique$start[i], bot, unique$end[i], top, col = cols[e], border = NA)
  }
}



###
# diff analysis
###

D1_FP_erato <- read.table('sections_pan_coords/sections_D1_p05FC1_erato_FPup_pan.bed', h=F)
D1_FM_erato <- read.table('sections_pan_coords/sections_D1_p05FC1_erato_FMup_pan.bed', h=F)
D1_FD_erato <- read.table('sections_pan_coords/sections_D1_p05FC1_erato_FDup_pan.bed', h=F)

D2_FP_erato <- read.table('sections_pan_coords/sections_D2_p05FC1_erato_FPup_pan.bed', h=F)
D2_FM_erato <- read.table('sections_pan_coords/sections_D2_p05FC1_erato_FMup_pan.bed', h=F)
D2_FD_erato <- read.table('sections_pan_coords/sections_D2_p05FC1_erato_FDup_pan.bed', h=F)

D1_FP_melp <- read.table('sections_pan_coords/sections_D1_p05FC1_melp_FPup_pan.bed', h=F)
D1_FM_melp <- read.table('sections_pan_coords/sections_D1_p05FC1_melp_FMup_pan.bed', h=F)
D1_FD_melp <- read.table('sections_pan_coords/sections_D1_p05FC1_melp_FDup_pan.bed', h=F)

# D2_FP_melp <- read.table('sections_D2_p05FC1_melp_FPup_pan.bed', h=F)
# D2_FM_melp <- read.table('sections_D2_p05FC1_melp_FMup_pan.bed', h=F)
# D2_FD_melp <- read.table('sections_D2_p05FC1_melp_FDup_pan.bed', h=F)

D1_FP_erato$y <- 90
D1_FM_erato$y <- 90
D1_FD_erato$y <- 90
D2_FP_erato$y <- 90
D2_FM_erato$y <- 90
D2_FD_erato$y <- 90
D1_FP_melp$y <- 90
D1_FM_melp$y <- 90
D1_FD_melp$y <- 90

xxx <- c()
erato_DF_list <- list(xxx,xxx,xxx,D1_FP_erato,D1_FM_erato,D1_FD_erato,xxx,D2_FP_erato,D2_FM_erato,D2_FD_erato,xxx)
melp_DF_list <- list(xxx,xxx,xxx,D1_FP_melp,D1_FM_melp,D1_FD_melp,xxx,xxx,xxx,xxx,xxx)
###
# ATAC
###

bot2 = 0
top2 = 100

for(e in 1:length(erato_list)){
  
  unique <- blocks_erato
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > startP) & (unique$end < endP))
  
  plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
  for(i in 1:nrow(unique)){
    rect(unique$start[i]+start, bot2, unique$end[i]+start, top2, col = adjustcolor(cols[1], alpha = alphas[e]), border = NA)
  }
  
  par(new=T)
  plot(0.5*(start(erato_list[[e]]) + end(erato_list[[e]])), erato_list[[e]]$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none')
  
  if(any(e==c(4,5,6,8,9,10))){
    par(new=T)
    plot((erato_DF_list[[e]]$V2+erato_DF_list[[e]]$V3)/2, erato_DF_list[[e]]$y, xlim = c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE, col = 'red', cex = 5, pch=19)
    
  }
  mtext(names[e], side = 1, cex=5, padj = -2, las = 1, adj=0)
  
}
  

  
bot2 = 0
top2 = 100

for(e in 1:length(melp_list)){
  
  unique <- blocks_melp
  colnames(unique) <- c('scaf', 'start', 'end')
  unique <- subset(unique, (unique$end-unique$start >= minLemgth) & (unique$start > startP) & (unique$end < endP))
  
  plot(NULL, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE)
  for(i in 1:nrow(unique)){
    rect(unique$start[i]+start, bot2, unique$end[i]+start, top2, col = adjustcolor(cols[2], alpha = alphas[e]), border = NA)
  }
  rect(381404732, bot2, 381406466, top2, col = adjustcolor("yellow",0.2), border = NA, lwd=8)
  par(new=T)
  plot(381407644, 95, xlim=c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE, col = adjustcolor('blue',0.2), cex = 5, pch=19)
  
  par(new=T)
  plot(0.5*(start(melp_list[[e]]) + end(melp_list[[e]])), melp_list[[e]]$score, type='l', xlim = c(start+startP,start+endP), ylim = c(bot2,top2), ylab = "", yaxt = "n", lwd = 5, xlab = "", xaxt = "n", main = "", bty='none')
  
  if(any(e==c(4,5,6))){
    par(new=T)
    plot((melp_DF_list[[e]]$V2+melp_DF_list[[e]]$V3)/2, melp_DF_list[[e]]$y, xlim = c(start+startP,start+endP), ylim = c(bot2,top2), axes=FALSE, ann=FALSE, col = 'red', cex = 5, pch=19)
    
  }
  mtext(names[e], side = 1, cex=5, padj = -2, las = 1, adj=0)
  
} 


plot(NULL, xlim=c(startP,endP), ylim = c(0,1), axes=FALSE, ann=FALSE)
axis(1, at = seq(startP,endP, by=10000), labels = NA, line =-10, lwd = 5, lwd.ticks = 5, tck = -0.2)
axis(1, at = seq(startP,endP, by=100000), labels = seq(startP/1000000,endP/1000000, by=0.1), line =-3, cex.axis = 6, lwd=0)


dev.off()
  
