setwd("D:/Research/Data/GSEA_prerank")

####Data####
# ----> Read data----
## Prognosis genus
select_genus <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv")
select_genus <- select_genus$Genus
## Genus x Gene SCC
genus_gene_SCC <- read.csv("D:/Research/Data/Redo/Correlation_SCC_origin.csv",row.names = 1)
# ----> Filter the potential genus----
selected_SCC <- genus_gene_SCC[,which(names(genus_gene_SCC)%in%select_genus)]
name_r <- rownames(selected_SCC)

####Create rnk file####
for (i in 1:23) {
  table <-name_r
  table <- data.frame(table)
  table$scc <- selected_SCC[,i]
  name <- names(selected_SCC)[i]
  
  table <- table[order(-table$scc),]
  name <- strsplit(name,split = "\\.")
  name <- name[[1]][length(name[[1]])]
  names(table)[1] <- name
  filename <- paste(name,".rnk",sep = "")
  
  write.table(table,file=filename,quote=F,sep="\t",row.names=F,col.names = F)
}