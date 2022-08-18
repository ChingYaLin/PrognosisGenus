#### Functions ####
PathName <- function(foldername,filename){
  # To read the folder file, set every file's path
  path <- paste("D:/Research/TCGA",foldername,filename,sep = '/')
  return(path)
}
#### Main ####
# ----> Set pathway and package----
setwd("D:/Research/TCGA")
library("dplyr")
# ----> Read data----
## manifest that 
manifest <- read.table("gdc_manifest_20210628_120505.txt",header = T)
allsample <- read.csv("D:/Research/Data/Raw/TCGA.csv",header = T)
## Grasp the information of file path
tumor_dir <- allsample$tumor
normal_dir<- allsample$normal
normal_dir <- normal_dir[!normal_dir %in% c("")]
fil <- allsample %>% dplyr::filter(normal %in% normal_dir)
tumor_c <- fil$tumor

# ----> Merge tumor raw count(138)----
## Merge the tumor sample count together and rename columns
# Read first file to create data.frame to combine with other files
filepath <- PathName(manifest$id[1],manifest$filename[1])
sample_name <- allsample$X[which(allsample$tumor==manifest$id[1])]
com <- read.table(test)
names(com)[names(com) == "V2"] <- paste("TM_1_",sample_name,sep = "")
# Merge the other files into the dataframe
for (i in 2:length(tumor_dir)) {
  match_index <- match(tumor_dir[i],manifest$id)
  filepath <- PathName(manifest$id[match_index],manifest$filename[match_index])
  temp <- read.table(filepath)
  temp$V1<-substr(temp$V1,1,15) #Remove the extra string in gene name
  rownames(temp) <- temp$V1
  #Check whether the sample have the normal data
  if(manifest$id[match_index] %in% tumor_c){
    sample_name <- allsample$X[which(allsample$tumor==manifest$id[match_index])]
    # TM = tumor have matched normal
    name <- paste("TM_",as.character(count),"_",sample_name,sep = "")
    names(temp)[names(temp) == "V2"] <- name
  }
  else{
    sample_name <- allsample$X[which(allsample$tumor==manifest$id[match_index])]
    # TS =  tumor have single sample
    name <- paste("TS_",as.character(count),"_",sample_name,sep = "")
    names(temp)[names(temp) == "V2"] <- name
  }
  com <- cbind(com,temp)
  com<-com[,!(names(com) %in% c("V1"))]
}
# ----> Merge normal raw count(34)----
## Merge the normal sample file into com and rename by N_1,N_2...
for (i in 1:length(normal_dir)) {
  match_index <- match(normal_dir[i],manifest$id)
  filepath <- PathName(manifest$id[match_index],manifest$filename[match_index])
  temp <- read.table(filepath)
  temp$V1<-substr(temp$V1,1,15) #Remove the extra string in gene name
  rownames(temp) <- temp$V1
  name <- paste("NO_",as.character(count),sep = "")
  names(temp)[names(temp) == "V2"] <- name
  com <- cbind(com,temp)
  com<-com[,!(names(com) %in% c("V1"))]
}
# ----> Write the gene raw file----
write.csv(com,file = 'D:/Research/Data/Raw/LUAD_origin/LUAD_gene_expression_raw(138sample).csv',row.names = T)
