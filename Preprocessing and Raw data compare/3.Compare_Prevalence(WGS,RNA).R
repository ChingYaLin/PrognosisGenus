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

Prevalence_cal <- function(df,cutoff=0){
  # This function can calculate the genus prevalence of whole data frame
  ## Input: df--data frame that column is genus, raw is sample,
  ##        cutoff--the baseline of genus that seen presented in the sample, default is 0
  ## Output: a data frame that have whole genus prevalence
  prevalence <- colSums(df >= cutoff)
  prevalence <- prevalence/length(df[,1])
  return(data.frame(prevalence))
}

####Main####
# ----> Setting environment and load files----
setwd("D:/Research/Data/Raw/LUAD_origin")
### Whole LUAD sample compare(WGS:162 samples,RNA:519 samples)
#df_wgs <- read.csv("LUAD_WGS_abundance.csv") 
#df_rna <- read.csv("LUAD_RNA_seq_abundance.csv")
### Same LUAD sample compare(WGS:138 samples,RNA:138 samples)
df_wgs <- read.csv("LUAD_WGS_abundance_138sample.csv")
df_rna <- read.csv("LUAD_RNA_abundance_138sample.csv")

# ----> Rename the column's name to genus level----
df_wgs <- RenameGenus(df_wgs)
df_rna <- RenameGenus(df_rna)

# ----> Remove the column of contaminant----
contaminant <- c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_rna = df_rna[,!(names(df_rna) %in% contaminant)]
df_wgs = df_wgs[,!(names(df_wgs) %in% contaminant)]

# ----> Calculate the prevalence----
prevalence_rna <- Prevalence_cal(df_rna,2)
prevalence_wgs <- Prevalence_cal(df_wgs,2)

# ----> Combine result----
compare_table <- data.frame(prevalence_rna,prevalence_wgs)
names(compare_table) <- c("RNA","WGS")

# ----> Visualization----
library(ggplot2)
# Sort the value
compare_table <- compare_table[order(-compare_table$RNA,-compare_table$WGS),]
num <- c(1:length(compare_table$RNA))
compare_table$num <- num
compare_table$compare <- compare_table$RNA-compare_table$WGS

#count the number of > and <
rna_bigger <- sum(compare_table$compare>=0.1)
wgs_bigger <- sum(compare_table$compare<=(-0.1))
equal <- sum(compare_table$compare>(-0.1) & compare_table$compare <0.1)
cat("RNA > 0.1:",rna_bigger,"RNA < -0.1:",wgs_bigger,"Equal:",equal)

# Draw the line
ggplot2::ggplot(compare_table) +
  ggplot2::geom_histogram(aes(compare),fill = "blue",binwidth = 0.01)+
  ggplot2::labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="value", y="count")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15))
ggsave("D:/Research/PlotNew/2.Prevalence compare(WGS,RNA)_histogram.png",
       width = 10, height = 7.5)

new <- compare_table[,c("num","compare")]
new <- new[order(new$compare),]
new$num <- num
ggplot2::ggplot(new,aes(x=num,y=compare))+
  ggplot2::geom_line(color="blue")+
  ggplot2::geom_hline(yintercept = 0,color = "red")+
  ggplot2::labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="microbiome", y="comparison")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15))+
  ggplot2::annotate(geom="text", x=400, y=0.3, 
                    label=paste("The number of sample that RNA minus WGS value :\n",
                                "RNA > 0.1:",rna_bigger,"/ RNA < -0.1:",wgs_bigger,"/ Equal:",equal),
                    size=5)

ggsave("D:/Research/PlotNew/3.Prevalence compare(WGS,RNA)_line.png",
       width = 10, height = 7.5)

