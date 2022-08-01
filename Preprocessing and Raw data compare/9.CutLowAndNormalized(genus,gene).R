####This Script want to cut low expression of genus normalized data and gene expression data
####Main####
# ----> Import library----
library(edgeR)
library(dplyr)
# ----> Set environment----
setwd("D:/Research/Data/Raw/LUAD_origin")

# ----> Cut low express genus----
## Read genus abundance file
genus_raw <- read.table("LUAD_WGS_abundance_138sample.csv", header = TRUE,row.names = 1,sep = ",")
##Clean the contaminant column
contaminant <- c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
genus_raw <- genus_raw[,-which(names(genus_raw)%in%contaminant)]
##Calculate the mean of normalized data
mean_normalized_data <- apply(genus_raw, 2, mean)
## Draw the histogram and set the cutoff
png(file="D:/Research/PlotNew/24.Genus abundance cut low express.png",
    width=800, height=600)
hist(x=mean_normalized_data,
     main="Histogram of genus abunance",         
     xlab="Mean of relative abundance",                      
     ylab="Count",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = seq(-5,20,1))

expr_cutoff <- 0
abline(v = expr_cutoff, col = "red", lwd = 3)
dev.off()
## Check how many genus would remain and filter the genus data
sum(mean_normalized_data > expr_cutoff)
cut_genus <- genus_raw[, mean_normalized_data > expr_cutoff]

# ----> Cut low express gene----
##deal the gene expression data
gene_raw <- read.table("LUAD_gene_expression_raw(138sample).csv", header = TRUE,row.names = 1,sep = ",")
##Clean the non gene raws
gene_raw <- gene_raw[-c(60485:60488),]
##Remove the normal sample
gene_raw <- dplyr::select(gene_raw,-starts_with("N"))
##Transfer raw count to log cpm
cpm_log <- cpm(gene_raw, log = TRUE)
cpm_log <- data.frame(cpm_log)
mean_log2_cpm <- apply(cpm_log, 1, mean)
## Draw the histogram and set the cutoff
png(file="D:/Research/PlotNew/25.Gene expression cut low express.png",
    width=800, height=600)
hist(x=mean_log2_cpm,
     main="Histogram of gene expression(log2cpm)",         
     xlab="Mean of log2cpm",                      
     ylab="Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2)
expr_cutoff <- 0
abline(v = expr_cutoff, col = "red", lwd = 3)
dev.off()
## Check how many gene would remain and filter the gene data
sum(mean_log2_cpm > expr_cutoff)
cut_gene <- cpm_log[mean_log2_cpm > expr_cutoff,]
cut_gene <- data.frame(cut_gene)
##rename the column by sample
name <- names(cut_gene)
for (i in 1:length(name)) {
  temp_len <- nchar(name[i])
  temp <- substr(name[i],temp_len-5,temp_len)
  names(cut_gene)[names(cut_gene) == name[i]] <- temp
}
# ----> Write genus and gene data----
##write the normalized gene data and genus data into CSV
write.csv(cut_gene,file="D:/Research/Data/Redo/log2cpm_clean_gene(CutLow).csv",row.names = T)
write.csv(cut_genus,file="D:/Research/Data/Redo/voom_svm_clean_genus(CutLow).csv",row.names = T)

