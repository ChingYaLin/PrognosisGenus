####Package####
library (EDASeq)
library(edgeR)
library(biomaRt)

####Functions####
RenameSample <- function(df){
  name <- names(df)
  for (i in 1:length(name)) {
    temp_len <- nchar(name[i])
    temp <- substr(name[i],temp_len-5,temp_len)
    names(df)[names(df) == name[i]] <- temp
  }
  return(df)
}

####Data preprocessing####
# ----> Set path----
setwd("D:/Research/Data/Raw/LUAD_origin")

# ----> Load data----
## Loading raw count of gene expression
raw_count <-  read.csv("LUAD_gene_expression_raw(138sample).csv", header = TRUE,row.names = 1)

# ----> Clean useless row and column----
## remove normal sample
raw_count <- raw_count[,-grep("N",names(raw_count))]
raw_count <- raw_count[-which(rownames(raw_count)%in%c("__no_feature","__ambiguous",
                                                       "__too_low_aQual","__not_aligned","__alignment_not")),]
## Rename the column name to sample
raw_count <- RenameSample(raw_count)

# ----> Cut low express genes----
cpm_log <- cpm(raw_count, log = TRUE)
mean_log2_cpm <- apply(cpm_log, 1, mean)
expr_cutoff <- 0
## The histogram already plot in preprocessing part
#hist(mean_log2_cpm)
#abline(v = expr_cutoff, col = "red", lwd = 3)
sum(mean_log2_cpm > expr_cutoff)
data_clean <- raw_count[mean_log2_cpm > expr_cutoff, ]

'''
# ----> Search the gene length(for tpm)----
###!!! Need ti take long time!!!!
ensembl_list <- rownames(data_clean)
result <- EDASeq::getGeneLengthAndGCContent(ensembl_list[1:1000], "hsa")
for (i in 2:(length(ensembl_list)%/%1000)) {
  temp <- EDASeq::getGeneLengthAndGCContent(ensembl_list[(i*1000+1):(i+1)*1000])
  result <- rbind(result,temp)
}
temp <- EDASeq::getGeneLengthAndGCContent(ensembl_list[(length(ensembl_list)%/%1000)*1000+1:length(ensembl_list)], "hsa")
result <- rbind(result,temp)

# ----> Remove na rows----
## remove the rows that delete by na
result <- na.omit(result) 

# ----> Save the gene length file----
write.csv(result,file="D:/Research/Data/Redo/Gene_length.csv",row.names = T)
'''

# ----> Load gene length file----
result <- read.csv("D:/Research/Data/Redo/Gene_length.csv",row.names = 1)
result <- na.omit(result)
data_clean <- data_clean[which(rownames(data_clean)%in%rownames(result)),]

# ----> Convert to HGNC gene symbol----
## Transfer the ID name to the HSA ID that TIMER2.0 need
ensembl_list <- rownames(data_clean)
human <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- biomaRt::getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id", values=ensembl_list, mart=human)

# ----> Rename the rownames----
data_clean <- data_clean[which(rownames(data_clean)%in%gene_coords$ensembl_gene_id),]
data_clean$X <- rownames(data_clean)
data_clean$HGNC<-0
for (i in 1:length(data_clean$HGNC)) {
  data_clean$HGNC[i]<-gene_coords[which(gene_coords$ensembl_gene_id==data_clean$X[i]),1]
}
## Remove the rows that not have HGNC symbol
data_clean <- data_clean[-which(data_clean$HGNC%in%c("")),]
result <- result[which(rownames(result)%in%data_clean$X),]

rownames(data_clean) <-data_clean$HGNC
data_clean <- data_clean[,-which(names(data_clean)%in%c("HGNC","X"))]

# ----> Transfer to TPM----
for (i in 1:length(names(data_clean))) {
  countt<-data_clean[names(data_clean)[i]]
  len <- result$length
  data_clean[names(data_clean)[i]] <- countt/len
  data_clean[names(data_clean)[i]] <- data_clean[names(data_clean)[i]]/sum(data_clean[names(data_clean)[i]],na.rm = T)*1e6
}
# ----> Write TIMER2.0 file----
write.csv(data_clean,file="D:/Research/Data/Redo/138sample_tpm_HGNC.csv",row.names = T)



