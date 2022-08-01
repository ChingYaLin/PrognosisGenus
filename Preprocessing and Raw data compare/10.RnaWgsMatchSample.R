####Match ID####
setwd("D:/Research/Data/Raw/LUAD_origin")

# ----> Load data and order----
rna <- read.csv("Metadata_LUAD_RNA_138sample.csv",header = T)
wgs <- read.csv("Metadata_LUAD_WGS_138sample.csv",header = T)
#Order by sample uuid
rna <- rna[order(rna$sample_uuid),]
wgs <- wgs[order(wgs$sample_uuid),]

# ----> Combine RNA and WGS sample ID ----
#Combine RNA and WGS by same sample uuid
match_table <- cbind(rna$X,wgs$X,rna$sample_uuid,wgs$sample_uuid)
match_table <- data.frame(match_table)
#Remove extra column and rename column
match_table <- match_table[,1:3]
names(match_table)<-c("rna_sample","wgs_sample","sample_uuid")

# ----> Write file ----
write.csv(match_table,"Match_rna_wgs_sample.csv",row.names = F)

#### Create the RNA sample name's files####
# ----> Gene expression ----
gene_data <- read.csv("D:/Research/Data/Redo/log2cpm_clean_gene(CutLow).csv",header = T,row.names = 1)
match_table <- match_table[order(match_table$wgs_sample),]
gene_data <- gene_data[,order(names(gene_data))]
names(gene_data) <- match_table$rna_sample
write.csv(gene_data,"D:/Research/Data/Redo/log2cpm_clean_gene_rna.csv",row.names = T)

# ----> Genus abundance ----
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
genus_rna <- read.csv("LUAD_RNA_abundance_138sample.csv", header = TRUE,row.names = 1)
genus_rna <- genus_rna[,which(names(genus_rna)%in%names(genus))]
write.csv(genus_rna,"D:/Research/Data/Redo/voom_svm_clean_genus_rna.csv",row.names = T)
