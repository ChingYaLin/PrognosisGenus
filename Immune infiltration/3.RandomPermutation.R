####Functions####
BlankDataFrame <- function(){
  df <- data.frame(matrix(0,length(names(immune)),length(names(genus))))
  rownames(df) <- names(immune)
  names(df) <- names(genus)
  return(df)
}

####Set path and Load Data####
# ----> Set path----
setwd("D:/Research/Data/Redo")
# ----> Load data and clean----
##Immune infiltration
immune <- read.csv("estimation_matrix.csv", header = TRUE,row.names = 1)
immune <- data.frame(t(immune))
## Genus abundance
genus <- read.csv("voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
##Grep CIBERSORT.ABS column
immune <- immune[,grep("CIBERSORT.ABS",names(immune))]

# ----> Random permutation----
## Random permutation(reorder sample value)
for (i in 1:length(colnames(genus))) {
  temp <- colnames(genus)[i]
  genus[[temp]] <-sample(genus[[temp]])
}
##Sort sample name
genus <- genus[ order(row.names(genus)), ]
immune <- immune[ order(row.names(immune)), ]

####SCC calculation####
# ----> Create a blank dataframe----
df_cor <- BlankDataFrame()
df_pvalue <- BlankDataFrame()

# ----> Calculation----
for (i in 1:length(names(genus))) {
  for (j in 1:length(names(immune))) {
    m<-names(genus)[i]
    g<-names(immune)[j]
    test<-cor.test(genus[[m]], immune[[g]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    df_cor[which(rownames(df_cor)==g),which(names(df_cor)==m)] <- rho
    df_pvalue[which(rownames(df_pvalue)==g),which(names(df_pvalue)==m)] <- pvalue
  }
}

# ----> Write file----
write.csv(df_cor,file="ImmuneCorrelation_SCC_random(prognosis).csv",row.names = T)
write.csv(df_pvalue,file="ImmuneCorrelation_pvalue_random(prognosis).csv",row.names = T)

