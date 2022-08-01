#### Functions ####
RenameGenus <- function(dataframe){
  #Rename the bacteria's rename into only genus level
  colname <- colnames(dataframe)
  for (j in 1:length(colname)) {
    if (colname[j]!='X'){
      temp <- strsplit(colname[j],split = "\\.") # separate string by [.]
      names(dataframe)[names(dataframe)==colname[j]] <- temp[[1]][length(temp[[1]])]
    }
  }
  return(dataframe)
}

#### Main####
# ---->Set file pathway----
setwd("D:/Research/Data/Raw/LUAD_origin")

# ---->Load Data (WGS and RNA-seq)----
RNA_abundance <- read.csv("LUAD_RNA_abundance_138sample.csv",row.names = 1)
WGS_abundance <- read.csv("LUAD_WGS_abundance_138sample.csv",row.names = 1)
RNA_abundance <- RNA_abundance[,-c(1285:1287)]
WGS_abundance <- WGS_abundance[,-c(1285:1287)]
# ---->Rename the genus name----
RNA_abundance <- RenameGenus(RNA_abundance)
WGS_abundance <- RenameGenus(WGS_abundance)
# ----> Create comparing table----
new_wgs <- colMeans(WGS_abundance)
new_rna <- colMeans(RNA_abundance)
compare_table <- cbind.data.frame(wgs_mean = new_wgs, rna_mean = new_rna)

# ----> sign the color by expression----
# Calculate the mean between WGS_mean and RNA_mean
mean_together <- rowMeans(compare_table)
compare_table$mean_together <- mean_together
# Filter low abundance genus
compare_table <- compare_table[!(compare_table$wgs_mean<=0|compare_table$rna_mean<=0),]
compare_table$wgs_proportion <- compare_table$wgs_mean/sum(compare_table$wgs_mean)
compare_table$rna_proportion <- compare_table$rna_mean/sum(compare_table$rna_mean)
compare_table$color <- ifelse((compare_table$wgs_mean>=15|compare_table$rna_mean>=15),rownames(compare_table),"Abundance<15")
# ----> Combine the WGS and RNA-seq part for draw proportion----
rna <- compare_table[,c(-1,-4)]
wgs <- compare_table[,c(-2,-5)]
rna$type <- "RNA_seq"
rna$genus <- rownames(rna)
wgs$type <- "WGS"
wgs$genus <- rownames(wgs)
names(rna) <- c("mean","mean_togather","proportion","color","type","genus")
names(wgs) <- c("mean","mean_togather","proportion","color","type","genus")
drawdata <- rbind(rna,wgs)

# ----> Drawing proportion----
#Import drawing package
library("dplyr")
library("ggplot2")
#Drawing whole proportion data
ggplot2::ggplot(drawdata,aes(x = type,fill = color)) + 
  ggplot2::geom_bar(position = "fill")+
  ggplot2::labs(title = "Proportion of each genus in different data type",
       x="Data type",
       y="Proportion")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text = element_text(size=10))

ggsave("D:/Research/PlotNew/4.GenusProportion compare(WGS,RNA)_allgenus.png",
       width = 10, height = 7.5)
#Filter abundance < 15 genus
drawdata_remove_low_abundance<- drawdata%>%
  filter(drawdata$color != "Abundance<15")
ggplot2::ggplot(drawdata_remove_low_abundance,aes(x = type,fill = color)) + 
  ggplot2::geom_bar(position = "fill")+
  ggplot2::labs(title = "Proportion of each genus in different data type",
                subtitle = "Remove abundance < 15 genus",
                x="Data type",
                y="Proportion")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text = element_text(size=10))
ggsave("D:/Research/PlotNew/5.GenusProportion compare(WGS,RNA)_bigger15.png",
       width = 10, height = 7.5)
