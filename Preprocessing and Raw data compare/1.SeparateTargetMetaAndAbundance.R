#### This program want to filter the LUAD data from all kind of cancer data
####Functions part####
Filtermetadata <- function(metadata,cancer_type){
  #Filter metadata that RNA-seq and WGS in target cancer type and Primary tumor
  clean_data <- metadata[metadata$investigation==cancer_type, ]
  clean_data <- clean_data[clean_data$sample_type=='Primary Tumor', ]
  data_wgs <- clean_data[clean_data$experimental_strategy == "WGS",]
  data_rnaseq <- clean_data[clean_data$experimental_strategy == 'RNA-Seq',]
  return(list(wgs = data_wgs,rna = data_rnaseq,all_target = clean_data))
}
ExtractAbundanceData <- function(target_file,sample_list){
  #This function can extract sub-dataframe from whole data frame by sample names
  sub_data <- target_file[target_file$X %in% c(sample_list),]
  return(sub_data)
}
#### Main####
# ---->Set file pathway----
setwd("D:/Research/Data/Raw/Paper origin data")

# ---->Read data----
#Metadata
df_meta <- read.csv("Metadata-TCGA-Kraken-17625-Samples.csv")
#Microbiome abundance
df_kraken <- read.csv("Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")

# ---->Check the data----
head(df_meta)
head(df_kraken)

# ---->Get the target data from whole TCGA sample----
# We filter the LUAD data from raw data
clean_meta <- Filtermetadata(df_meta,cancer_type = "TCGA-LUAD")
#Whole target (We separate the wgs and rna type data)
meta_wgs<- clean_meta[["wgs"]]
meta_rnaseq <- clean_meta[["rna"]]
#Remove replicate sample in RNA-seq and WGS
filter_rna <- meta_rnaseq[!duplicated(meta_rnaseq[ , c("sample_uuid")]), ] 
filter_wgs <- meta_wgs[!duplicated(meta_wgs[ , c("sample_uuid")]), ] 
#Extract the samples that have WGS data and RNA data in the same time
same_UUID <- intersect(filter_rna$sample_uuid,filter_wgs$sample_uuid)
filter_rna <- filter_rna[filter_rna$sample_uuid %in% c(same_UUID),]
filter_wgs <- filter_wgs[filter_wgs$sample_uuid %in% c(same_UUID),]

# ---->Write the WGS and RNA metadata file----
write.csv(meta_wgs,file = 'D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_WGS.csv',row.names = F)
write.csv(meta_rnaseq,file = 'D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_RNAseq.csv',row.names = F)
write.csv(filter_wgs,file = "D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_WGS_138sample.csv",row.names = F)
write.csv(filter_rna,file = "D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_RNA_138sample.csv",row.names = F)

# ---->Filter target cancer type abundance data by sample list----
#Less sample abundance extract
data_filter_rna <- ExtractAbundanceData(df_kraken,filter_rna$X)
data_filter_WGS <- ExtractAbundanceData(df_kraken,filter_wgs$X)
#All target samples
data_filter <- ExtractAbundanceData(df_kraken,clean_meta$all_target$X)
data_filter_rna_all <- ExtractAbundanceData(df_kraken,meta_rnaseq$X)
data_filter_WGS_all <- ExtractAbundanceData(df_kraken,meta_wgs$X)

# ---->Write the WGS and RNA abundance file----
write.csv(data_filter_rna,file = "D:/Research/Data/Raw/LUAD_origin/LUAD_RNA_abundance_138sample.csv",row.names = F)
write.csv(data_filter_WGS,file = "D:/Research/Data/Raw/LUAD_origin/LUAD_WGS_abundance_138sample.csv",row.names = F)
write.csv(data_filter,file="D:/Research/Data/Raw/LUAD_origin/LUAD_kraken_voom_snm.csv",row.names = F)
write.csv(data_filter_rna_all,file = "D:/Research/Data/Raw/LUAD_origin/LUAD_RNA_seq_abundance.csv",row.names = F)
write.csv(data_filter_WGS_all,file = "D:/Research/Data/Raw/LUAD_origin/LUAD_WGS_abundance.csv",row.names = F)
