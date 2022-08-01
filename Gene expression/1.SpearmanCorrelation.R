####Package####
library("ggpubr")
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
## Sort the row name
genus <- genus[order(row.names(genus)), ]
gene <-gene[order(row.names(gene)), ]

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

####Write SCC files####
write.csv(df_cor,file="Correlation_SCC_origin(prognosis).csv",row.names = T)
write.csv(df_pvalue,file="Correlation_pvalue_origin(prognosis).csv",row.names = T)

#### Big correlation dotplot####
#----> Find the big correlate ----
temp_m <- c()
temp_g <- c()
for (i in 1:length(names(genus))) {
  for (j in 1:length(names(gene))) {
    if(abs(df_cor[[names(genus)[i]]][j])>0.35){
      temp_m <- c(temp_m,names(genus)[i])
      temp_g <- c(temp_g,names(gene)[j])
    }
  }
}

# ----> Plot SCC dotplot ----
drawcorrelation <- function(num) {
  highscc <- cbind(genus[[temp_m[num]]],gene[[temp_g[num]]])
  highscc <- data.frame(highscc)
  
  
  temp_genus_name <- strsplit(temp_m[num],split = "\\.")
  temp_genus_name <- temp_genus_name[[1]][length(temp_genus_name[[1]])]
  
  ggpubr::ggscatter(highscc, x = "X1", y = "X2", 
            add = "loess", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = temp_genus_name, ylab = temp_g[num])+
    ggplot2::labs(title="Spearman correlation")+
    ggplot2::theme(plot.title = element_text(size = 20),
                   axis.title = element_text(size = 15),
                   axis.text = element_text(size = 15),
                   legend.text = element_text(size = 15),
                   legend.title = element_text(size = 15),
                   plot.subtitle = element_text(size = 15))
}

pdf("D:/Research/PlotNew/SCC_dotplot.pdf") 
for (i in 1:length(temp_g)) {
  plot(drawcorrelation(i))
}
dev.off()

