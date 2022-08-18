####Functions####
library(plyr)

RenameGenus <- function(df){
  # Rename the genus name into genus level
  colname_genus <- colnames(df)
  for (j in 1:length(colname_genus)) {
    if (colname_genus[j]!='X'){
      temp <- strsplit(colname_genus[j],split = "\\.")
      names(df)[names(df)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
    }
  }
  temp <- gsub("g_","",names(df))
  temp <- gsub("[_]","",temp)
  names(df) <- temp
  return(df)
}

RenameSample <- function(df){
  # Rename the sample's name to TCGA-Barcode
  temp <- colnames(df)
  temp <- substr(temp,1,12)
  temp <- gsub("[.]","-",temp)
  names(df) <- temp
  return(df)
}

TransferTheCategory <- function(df,category,informtable){
  # Transfer the genus abundance dataframe of each sample to other microbiome category level
  temp <- df
  temp$category <- informtable[[category]][which(informtable$query%in%rownames(temp))]
  temp <- plyr::ddply(temp,"category",numcolwise(sum))
  temp <- temp[is.na(temp$category)==F,]
  rownames(temp) <- temp$category
  temp <- temp[,-1]
  temp <- temp[,order(names(temp))]
  return(temp)
}

TransferProportion <- function(df){
  # Count the genus proportion
  df <- data.frame(t(df))
  df <- df/rowSums(df)
  df <- data.frame(t(df))
  return(df)
}

sccCalculate <- function(df1,df2){
  # To calculate the correlation between two type of data
  result <- names(df1)
  result <- data.frame(result)
  result$pvalue <- 0
  result$SCC <- 0
  #calculate the correlation and fill the matrix
  for (i in 1:length(names(df2))) {
    m<-names(df2)[i]
    g<-names(df1)[i]
    test<-cor.test(df2[[m]], df1[[g]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    result$SCC[i] <- rho
    result$pvalue[i] <- pvalue
  }
  return(result)
}

####Main####
# ----> Set path and load files----
#Set the path
setwd("D:/Research/Data/Redo")
#Load the microbiome taxonomy information
miRNA_bacter <- read.csv("Bacteria_miRNA.csv",row.names = 1)
my_bacter <- read.csv("Bacteria_mydata.csv",row.names = 1)
#Remove the phylum that no have answers
my_bacter <- my_bacter[which(my_bacter$phylum!=""),]
#Load the raw data of RNA, WGS and miRNA
miRNA_raw <- read.csv("D:/Research/Data/Raw/miRNA microbiome/LUAD_relative_abundance_miRNA.txt",sep = "\t")
raw_data <- read.csv("D:/Research/Data/Raw/Paper origin data/Kraken-TCGA-Raw-Data-17625-Samples.csv", header = TRUE,row.names = 1)
#RNA and WGS match sample: Inorder to catch raw microbiome data from whole table.
RNA_WGS_matchsample <- read.csv("D:/Research/Data/Redo/match_rna_wgs_sample.csv", header = TRUE)

# ----> Raw count data preprocessing----
#Rename the column by Genus level
raw_data <- RenameGenus(raw_data)
#To extract the microbiome column from abundance data and remove contaminant
raw_data <- raw_data[,which(names(raw_data)%in%my_bacter$query)]

# ----> miRNA data preprocessing----
#Rename the sample name to barcode
miRNA_raw <- RenameSample(miRNA_raw)
tgca <- read.csv("D:/Research/Data/Raw/TCGA.csv", header = TRUE,row.names = 1)
miRNA_raw <-miRNA_raw[,which(names(miRNA_raw)%in%tgca$barcode)]
for (i in 1:length(names(miRNA_raw))) {
  names(miRNA_raw)[i] <- rownames(tgca)[which(tgca$barcode==names(miRNA_raw)[i])]
}
miRNA_raw<-miRNA_raw[,1:134]
miRNA_genus <- miRNA_raw

# ----> Extract Microbiome and Sample----
#WGS sample
WGS_genus <- raw_data[which(rownames(raw_data)%in%names(miRNA_genus)),]
WGS_genus <- data.frame(t(WGS_genus))
#RNA sample (rename sample as WGS)
RNA_genus <- raw_data[which(rownames(raw_data)%in%RNA_WGS_matchsample$rna),]
for (i in 1:length(rownames(RNA_genus))) {
  rownames(RNA_genus)[i] <- RNA_WGS_matchsample$wgs[which(RNA_WGS_matchsample$rna==rownames(RNA_genus)[i])]
}
RNA_genus <- RNA_genus[which(rownames(RNA_genus)%in%names(WGS_genus)),]
RNA_genus <- data.frame(t(RNA_genus))

# ----> Order the sample----
miRNA_genus <- miRNA_genus[,order(names(miRNA_genus))]
WGS_genus <- WGS_genus[,order(names(WGS_genus))]
RNA_genus <- RNA_genus[,order(names(RNA_genus))]

# ---->Create different compare data----
## Genus (find common genus)
common_genus <- intersect(rownames(miRNA_genus),rownames(WGS_genus))
miRNA_g <- miRNA_genus[which(rownames(miRNA_genus)%in%common_genus),]
WGS_g <- WGS_genus[which(rownames(WGS_genus)%in%common_genus),]
RNA_g <- RNA_genus[which(rownames(RNA_genus)%in%common_genus),]
miRNA_g<-miRNA_g[order(rownames(miRNA_g)),]
WGS_g<-WGS_g[order(rownames(WGS_g)),]
RNA_g<-RNA_g[order(rownames(RNA_g)),]

miRNA_g <- TransferProportion(miRNA_g)
WGS_g <- TransferProportion(WGS_g)
RNA_g <- TransferProportion(RNA_g)
miRNA_g<-miRNA_g[,-which(names(miRNA_g)=="s11217")]
WGS_g<-WGS_g[,-which(names(WGS_g)=="s11217")]
RNA_g<-RNA_g[,-which(names(RNA_g)=="s11217")]
## Family
# Transfer the category into Family
miRNA_f <- TransferTheCategory(miRNA_genus,"family",miRNA_bacter)
WGS_f <- TransferTheCategory(WGS_genus,"family",my_bacter)
RNA_f <- TransferTheCategory(RNA_genus,"family",my_bacter)
# Intersect the common family
common_family <- intersect(rownames(miRNA_f),rownames(WGS_f))
miRNA_f <- miRNA_f[which(rownames(miRNA_f)%in%common_family),]
WGS_f <- WGS_f[which(rownames(WGS_f)%in%common_family),]
RNA_f <- RNA_f[which(rownames(RNA_f)%in%common_family),]
# Transfer data to proportional value
miRNA_f <- TransferProportion(miRNA_f)
WGS_f <- TransferProportion(WGS_f)
RNA_f <- TransferProportion(RNA_f)

## Phylum
# Transfer the category into Phylum
miRNA_p <- TransferTheCategory(miRNA_genus,"phylum",miRNA_bacter)
WGS_p <- TransferTheCategory(WGS_genus,"phylum",my_bacter)
RNA_p <- TransferTheCategory(RNA_genus,"phylum",my_bacter)


# Intersect the common Phylum
common_phylum <- intersect(rownames(miRNA_p),rownames(WGS_p))
miRNA_p <- miRNA_p[which(rownames(miRNA_p)%in%common_phylum),]
WGS_p <- WGS_p[which(rownames(WGS_p)%in%common_phylum),]
RNA_p <- RNA_p[which(rownames(RNA_p)%in%common_phylum),]
# Transfer data to proportional value
miRNA_p <- TransferProportion(miRNA_p)
WGS_p <- TransferProportion(WGS_p)
RNA_p <- TransferProportion(RNA_p)


# ----> Calculate correlations----
break1<-seq(0.5,1,0.01)
break2 <- seq(-0.5,0.5,0.05)
####>> Genus####
# ----~~~ WGS vs miRNA----
#By Sample
genusWvsMsample <- sccCalculate(WGS_g,miRNA_g)
hist(genusWvsMsample$SCC,main = "Compare correlation between miRNA and WGS sample",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
genusWvsMsampleList <- na.omit(genusWvsMsample)
#By Microbiome
genusWvsMmicrobiome <- sccCalculate(data.frame(t(WGS_g)),data.frame(t(miRNA_g)))
hist(genusWvsMmicrobiome$SCC,main = "Compare correlation between miRNA and WGS genus",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
genusWvsMmicrobiomeList <- na.omit(genusWvsMmicrobiome)

# ----~~~ WGS vs RNA----
#By Sample
genusWvsRsample <- sccCalculate(WGS_g,RNA_g)
png(file="D:/Research/PlotNew/6.WGSvsRNAsampleSCC.png",
    width=800, height=600)
hist(genusWvsRsample$SCC,main = "Sample correlation (Genus)" ,
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
dev.off()
genusWvsRsampleList <- na.omit(genusWvsRsample)
#By Microbiome
genusWvsRmicrobiome <- sccCalculate(data.frame(t(WGS_g)),data.frame(t(RNA_g)))
png(file="D:/Research/PlotNew/7.WGSvsRNAgenusSCC.png",
    width=800, height=600)
hist(genusWvsRmicrobiome$SCC,main = "Genus correlation",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break2)
dev.off()
genusWvsRmicrobiomeList <- na.omit(genusWvsRmicrobiome)

# ----~~~ RNA vs miRNA----
#By Sample
genusMvsRsample <- sccCalculate(miRNA_g,RNA_g)
hist(genusMvsRsample$SCC,main = "Compare correlation between RNA-seq and miRNA sample",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
genusMvsRsampleList <- na.omit(genusMvsRsample)
#By Microbiome
genusMvsRmicrobiome <- sccCalculate(data.frame(t(miRNA_g)),data.frame(t(RNA_g)))
hist(genusMvsRmicrobiome$SCC,main = "Compare correlation between RNA-seq and miRNA genus",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
genusMvsRmicrobiomeList <- na.omit(genusMvsRmicrobiome)

####>> Family####
# ----~~~ WGS vs miRNA----
#By Sample
familyWvsMsample <- sccCalculate(WGS_f,miRNA_f)
hist(familyWvsMsample$SCC,main = "Compare correlation between miRNA and WGS sample",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
familyWvsMsampleList <- na.omit(familyWvsMsample)
#By Microbiome
familyWvsMmicrobiome <- sccCalculate(data.frame(t(WGS_f)),data.frame(t(miRNA_f)))
hist(familyWvsMmicrobiome$SCC,main = "Compare correlation between miRNA and WGS family",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
familyWvsMmicrobiomeList <- na.omit(familyWvsMmicrobiome)

# ----~~~ WGS vs RNA----
#By Sample
familyWvsRsample <- sccCalculate(WGS_f,RNA_f)
png(file="D:/Research/PlotNew/8.WGSvsRNAsample_f_SCC.png",
    width=800, height=600)
hist(familyWvsRsample$SCC,main = "Sample correlation (Family)",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
dev.off()
familyWvsRsampleList <- na.omit(familyWvsRsample)
#By Microbiome
familyWvsRmicrobiome <- sccCalculate(data.frame(t(WGS_f)),data.frame(t(RNA_f)))
png(file="D:/Research/PlotNew/9.WGSvsRNAfamilySCC.png",
    width=800, height=600)
hist(familyWvsRmicrobiome$SCC,main = "Family correlation",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break2)
dev.off()
familyWvsRmicrobiomeList <- na.omit(familyWvsRmicrobiome)

# ----~~~ RNA vs miRNA----
#By Sample
familyMvsRsample <- sccCalculate(miRNA_f,RNA_f)
hist(familyMvsRsample$SCC,main = "Compare correlation between RNA-seq and miRNA sample",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
familyMvsRsampleList <- na.omit(familyMvsRsample)
#By Microbiome
familyMvsRmicrobiome <- sccCalculate(data.frame(t(miRNA_f)),data.frame(t(RNA_f)))
hist(familyMvsRmicrobiome$SCC,main = "Compare correlation between RNA-seq and miRNA family",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
familyMvsRmicrobiomeList <- na.omit(familyMvsRmicrobiome)


####>> Phylum####
# ----~~~ WGS vs miRNA----
#By Sample
phylumWvsMsample <- sccCalculate(WGS_p,miRNA_p)
hist(phylumWvsMsample$SCC,main = "Compare correlation between miRNA and WGS sample",
     xlab = "SCC",
     ylab = "Frequency")
draw <- data.frame(cbind(WGS = WGS_p$s11321, miRNA = miRNA_p$s11321))
library(ggpubr)
ggscatter(draw, x = "WGS", y = "miRNA",
          add = "reg.line",  
          conf.int = TRUE,    
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 0.25, label.y = 0.5)+
  labs(title = "Scatter plot of sample s11321",
       x = "WGS relative abundance",
       y = "miRNA relative abundance")

phylumWvsMsampleList <- na.omit(phylumWvsMsample)
#By Microbiome
phylumWvsMmicrobiome <- sccCalculate(data.frame(t(WGS_p)),data.frame(t(miRNA_p)))
hist(phylumWvsMmicrobiome$SCC,main = "Compare correlation between miRNA and WGS phylum",
     xlab = "SCC",
     ylab = "Frequency")
phylumWvsMmicrobiomeList <- na.omit(phylumWvsMmicrobiome)

# ----~~~ WGS vs RNA----
#By Sample
phylumWvsRsample <- sccCalculate(WGS_p,RNA_p)
png(file="D:/Research/PlotNew/10.WGSvsRNAsample_p_SCC.png",
    width=800, height=600)
hist(phylumWvsRsample$SCC,main = "Sample correlation (Phylum)",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break1)
dev.off()
phylumWvsRsampleList <- na.omit(phylumWvsRsample)
draw <- data.frame(cbind(WGS = WGS_p$s10650, RNA = RNA_p$s10650))
ggscatter(draw, x = "WGS", y = "RNA",
          add = "reg.line",  
          conf.int = TRUE,    
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 0.25, label.y = 0.5)+
  labs(title = "Scatter plot of sample s10650",
       x = "WGS relative abundance",
       y = "RNA relative abundance")
#By Microbiome
phylumWvsRmicrobiome <- sccCalculate(data.frame(t(WGS_p)),data.frame(t(RNA_p)))
png(file="D:/Research/PlotNew/11.WGSvsRNAphylumSCC.png",
    width=800, height=600)
hist(phylumWvsRmicrobiome$SCC,main = "Phylum correlation",
     xlab = "SCC",
     ylab = "Frequency",
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,
     breaks = break2)
dev.off()
phylumWvsRmicrobiomeList <- na.omit(phylumWvsRmicrobiome)

# ----~~~ RNA vs miRNA----
#By Sample
phylumMvsRsample <- sccCalculate(miRNA_p,RNA_p)
hist(phylumMvsRsample$SCC,main = "Compare correlation between RNA-seq and miRNA sample",
     xlab = "SCC",
     ylab = "Frequency")
phylumMvsRsampleList <- na.omit(phylumMvsRsample)
draw <- data.frame(cbind(miRNA = miRNA_p$s11447, RNA = RNA_p$s11447))
ggscatter(draw, x = "miRNA", y = "RNA",
          add = "reg.line",  
          conf.int = TRUE,    
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 0.25, label.y = 0.5)+
  labs(title = "Scatter plot of sample s11447",
       x = "miRNA relative abundance",
       y = "RNA relative abundance")

#By Microbiome
phylumMvsRmicrobiome <- sccCalculate(data.frame(t(miRNA_p)),data.frame(t(RNA_p)))
hist(phylumMvsRmicrobiome$SCC,main = "Compare correlation between RNA-seq and miRNA phylum",
     xlab = "SCC",
     ylab = "Frequency")
phylumMvsRmicrobiomeList <- na.omit(phylumMvsRmicrobiome)
