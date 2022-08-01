####Functions####
BlankDataFrame <- function(){
  df <- data.frame(matrix(0,length(names(gene)),length(names(genus))))
  rownames(df) <- names(gene)
  names(df) <- names(genus)
  return(df)
}

####Data prepare####
# ----> Set work path----
setwd("D:/Research/Data/Redo")

# ----> Loading data----
## load cut low abundance genus matrix and cut low normalized gene data
gene <- read.csv("log2cpm_clean_gene(CutLow).csv", header = TRUE,row.names = 1)
gene <-data.frame(t(gene))
genus <- read.csv("voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
prognosis_genus <- read.csv("Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv",header = T)
# ----> Filter genus----
## Filter the genus in prognosis
genus <- genus[,which(names(genus)%in%prognosis_genus$Genus)]

# ----> Random permutation----
## Random permutation(reorder sample value)
for (i in 1:length(colnames(gene))) {
  temp <- colnames(gene)[i]
  gene[[temp]] <-sample(gene[[temp]])
}

#sort the row name
genus <- genus[ order(row.names(genus)), ]
gene <- gene[ order(row.names(gene)), ]

# ----> Store data frame----
## Create a blank dataframe
df_cor <- BlankDataFrame()
df_pvalue <- BlankDataFrame()

####Spearman Correlation####
#calculate the correlation and fill the matrix
for (i in 1:length(names(genus))) {
  for (j in 1:length(names(gene))) {
    m<-names(genus)[i]
    g<-names(gene)[j]
    SCC<-cor.test(genus[[m]], gene[[g]], method=c("spearman"),exact=F)
    rho <- SCC$estimate[["rho"]]
    pvalue <- SCC$p.value
    df_cor[which(rownames(df_cor)==g),which(names(df_cor)==m)] <- rho
    df_pvalue[which(rownames(df_pvalue)==g),which(names(df_pvalue)==m)] <- pvalue
  }
}

write.csv(df_cor,file="Correlation_SCC_random(prognosis).csv",row.names = T)
write.csv(df_pvalue,file="Correlation_pvalue_random(prognosis).csv",row.names = T)



