
files <- list.files(path='counts_nonStranded', pattern="*.counts")

files_erato <- files[grepl('erato', files)]

files_melp <- files[grepl('melp', files)]



### erato 
for(e in 1:length(files_erato)){
  
  data <- read.table(paste('counts_nonStranded/', files_erato[e], sep=''), h=F)
  
  if(e==1){
    
    counts_erato <- data
  }
  else{
    counts_erato <- cbind(counts_erato, data[,2])
  }
}


header_erato <- gsub("f_erato.counts", "", files_erato)

colnames(counts_erato) <- c("gene_id", header_erato)

head(counts_erato)

counts_erato <- subset(counts_erato,  counts_erato$gene_id != '__no_feature'	&
                                      counts_erato$gene_id != '__ambiguous'	&
                                      counts_erato$gene_id != '__too_low_aQual'	&
                                      counts_erato$gene_id != '__not_aligned'	&
                                      counts_erato$gene_id != '__alignment_not_unique')	

nrow(counts_erato)


### melp
for(e in 1:length(files_melp)){
  
  data <- read.table(paste('counts_nonStranded/', files_melp[e], sep=''), h=F)
  
  if(e==1){
    
    counts_melp <- data
  }
  else{
    counts_melp <- cbind(counts_melp, data[,2])
  }
}


header_melp <- gsub("f_melp.counts", "", files_melp)

colnames(counts_melp) <- c("gene_id", header_melp)

head(counts_melp)

counts_melp <- subset(counts_melp,  counts_melp$gene_id != '__no_feature'	&
                         counts_melp$gene_id != '__ambiguous'	&
                         counts_melp$gene_id != '__too_low_aQual'	&
                         counts_melp$gene_id != '__not_aligned'	&
                         counts_melp$gene_id != '__alignment_not_unique')	

nrow(counts_melp)

### write

write.table(counts_erato, file = "counts_table_erato_nonStranded.txt", quote = F, row.names = F)
write.table(counts_melp, file = "counts_table_melp_nonStranded.txt", quote = F, row.names = F)
