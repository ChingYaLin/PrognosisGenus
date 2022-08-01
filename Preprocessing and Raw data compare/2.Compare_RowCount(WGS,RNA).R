#### Functions ####
RenameGenus <- function(df) {
  #This function can rename the data frame that the column name is genus
  bacteria_name <- colnames(df)
  for (i in 1:length(bacteria_name)) {
    if(bacteria_name[i]!="X"){
      temp <- strsplit(bacteria_name[i],split = "\\.")
      names(df)[names(df)==bacteria_name[i]] <- temp[[1]][length(temp[[1]])]
    }
  }
  df <- df[,!(names(df) %in% c("X"))]
  return(df)
}

CreateRawCountTable <- function(df){
  sums <- rowSums(df[,1:1284])
  sums <- data.frame(sums)
  sums$X <- df$UUID
  sums<-sums[str_order(sums$X),]
  rownames(sums) <- sums$X
  sums <- sums[,!(names(sums) %in% c("X"))]
  return(sums)
}

####Main####
# ----> Import package and setting ----
library(tidyverse)
setwd("D:/Research/Data/Raw/Paper origin data")

# ----> Load data and clean----
df_raw <- read.csv("Kraken-TCGA-Raw-Data-17625-Samples.csv")
filter_rna <- read.csv("D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_RNA_138sample.csv")
filter_wgs <- read.csv("D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_WGS_138sample.csv")
df_kraken <- read.csv("Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")
df_meta <- read.csv("Metadata-TCGA-Kraken-17625-Samples.csv")

#select 1285 columns form 1994 columns
drop <- c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_kraken = df_kraken[,!(names(df_kraken) %in% drop)]
genus_list <- colnames(df_kraken)
df_raw <- df_raw %>% select(genus_list)
df_raw$UUID <- df_meta$sample_uuid
# Filter raw count data that 138 samples
data_filter_rna <- df_raw[df_raw$X %in% c(filter_rna$X),]
data_filter_wgs <- df_raw[df_raw$X %in% c(filter_wgs$X),]

# ----> Rename genus name of data frame----
#rename the bacteria name into genus
data_filter_wgs <- RenameGenus(data_filter_wgs)
data_filter_rna <- RenameGenus(data_filter_rna)

# ----> Create raw count compare table----
#sum the row by sample
rna_sum <- CreateRawCountTable(data_filter_rna)
wgs_sum <- CreateRawCountTable(data_filter_wgs)

compare_table <- rbind(wgs_sum,rna_sum)
compare_table <- data.frame(t(compare_table))

# ----> Compare raw count of different data----
table_test <- compare_table$wgs_sum/compare_table$rna_sum

#Calculate the number of different data type which has two times counts bigger than other
rna_bigger <- sum(table_test <= 0.5)
wgs_bigger <- sum(table_test >= 2)
equal <- sum(table_test>0.5 & table_test<2)
cat("RNA bigger:",rna_bigger,"WGS bigger:",wgs_bigger,"Equal:",equal)
#Different data type's mean sum raw counts
wgs_mean <- round(mean(compare_table$wgs_sum))
rna_mean <- round(mean(compare_table$rna_sum))
cat("WGS mean:",wgs_mean,"/ RNA-seq mean:",rna_mean)

# ----> Drawing the table----
ggplot2::ggplot(compare_table,aes(x=wgs_sum,y=rna_sum))+
  ggplot2::geom_point(colour = "blue", size = 2,shape = 21)+
  ggplot2::geom_hline(yintercept = 500000,color = "red")+
  ggplot2::geom_vline(xintercept = 500000,color="red")+
  ggplot2::labs(title = "Raw counts compare",x="WGS", y="RNA-seq")+
  ggplot2::theme(plot.title = element_text(size =20),
        axis.title = element_text(size =15),
        axis.text = element_text(size =15))+
  ggplot2::annotate(geom="text", x=1.0e+07, y=1750000, 
                    label=paste("The number of sample:\n",
                                "RNA bigger:",rna_bigger,"/ WGS bigger:",wgs_bigger,"/ Equal:",equal,"\n",
                                "The mean of sample raw counts:\n",
                                "WGS mean:",wgs_mean,"/ RNA-seq mean:",rna_mean),
                    size=5)

ggplot2::ggsave("D:/Research/PlotNew/1.Raw counts compare(WGS,RNA).png",
                width = 10, height = 7.5)

