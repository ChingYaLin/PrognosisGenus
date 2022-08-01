#This script is to calculate the co-variance between two prognosis-associated genera
####Package####
library(pheatmap)
####Functions####
BlankDataFrame <- function(){
  # This function can create a blank matrix to store the SCC
  df <- data.frame(matrix(0,length(names(prog_genus)),length(names(prog_genus))))
  rownames(df) <- names(prog_genus)
  names(df) <- names(prog_genus)
  return(df)
}

####Main####
#----> Setting environment and data ----
# Set file path
setwd("D:/Research/Data/Redo")
# Load the abundance matrix of prognosis-associated genus
prog_genus <- read.csv('Five_year_prognosis.csv',row.names = 1)
# Remove the useless column
prog_genus <- prog_genus[,-which(names(prog_genus)%in%c('os.time','OS'))]

# ----> SCC calculation----
#Create a blank dataframe
df_cor <- BlankDataFrame() #To store SCC value
df_pvalue <- BlankDataFrame() #To store p-value of SCC

# ----> Calculation----
for (i in 1:length(names(prog_genus))) {
  for (j in 1:length(names(prog_genus))) {
    m<-names(prog_genus)[i]
    g<-names(prog_genus)[j]
    test<-cor.test(prog_genus[[m]], prog_genus[[g]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    df_cor[which(rownames(df_cor)==g),which(names(df_cor)==m)] <- rho
    df_pvalue[which(rownames(df_pvalue)==g),which(names(df_pvalue)==m)] <- pvalue
  }
}

#----> Heatmap visualization----
png(filename = "D:/Research/PlotNew/Prognosis/PrognosisGenus_CoVariance.png",
    width = 800, height = 750, units = "px")
pheatmap::pheatmap(df_cor,main = "The co-variance of prognosis-associated genera",
                   fontsize = 15)
dev.off()