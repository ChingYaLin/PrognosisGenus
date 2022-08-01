####Package####
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
####Main####
# ----> Set path and read data----
## Check species of dataset
msigdbr::msigdbr_species()
all_gene_sets = msigdbr::msigdbr(species = "Homo sapiens")
head(all_gene_sets)
## Check gene set
msigdbr::msigdbr_collections()
msigdbr_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)
## Read genus
select_genus <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv")
select_genus <- select_genus$Genus

# ----> GSEA----
result <- list()
rank_store <- list()
for (i in 1:length(select_genus)) {
  name <- strsplit(select_genus[i],split = "\\.")
  name <- name[[1]][length(name[[1]])]
  filename <- paste("D:/Research/Data/GSEA_prerank/",name,".rnk",sep = "")
  ranks <- read.table(filename,header=FALSE, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$V2, ranks$V1)
  fgseaRes <- fgsea(msigdbr_list, ranks,maxSize=500)
  rank_store[[name]] <- ranks
  result[[name]] <- fgseaRes
}
save(rank_store,result, file = "GSEA.RData")
## Check the result
head(result[[name]])
head(result[[name]][order(pval), ])

# ----> Plot enrichment----
plotEnrichment(msigdbr_list[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]],
               ranks) + labs(title="GOBP_ADAPTIVE_IMMUNE_RESPONSE")

#----> Combine top pathways----
load("GSEA.RData")
for (i in 1:23) {
  target <- i
  topPathwaysUp <- result[[target]][NES > 1][head(order(pval), n=10), pathway]
  topPathwaysDown <- result[[target]][NES < -1][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  #plot.new()
  #plotGseaTable(msigdbr_list[topPathways], rank_store[[target]], result[[target]], 
  #              gseaParam=0.5)
  
  
  # ----> Plot significant pathway----
  plot <- result[[target]]
  plot <- plot[which(plot$pathway%in%topPathways),]
  
  ggplot(plot, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste("Hallmark pathways NES from GSEA (",names(result)[target],")",sep = "")) + 
    theme_minimal()+
    ggplot2::theme(axis.text= element_text(size = 15),
                   axis.title = element_text(size = 15),
                   plot.title = element_text(size = 20),
                   legend.text = element_text(size = 15),
                   legend.title = element_text(size = 15))
  filepath <- paste("D:/Research/PlotNew/Hallmark pathway/Hallmark pathway_",names(result)[target],".png",sep="")
  ggsave(filepath,dpi=800,
         width = 600,
         height = 300,
         units = c("mm"))
}




#### Reference####
# ----> msigdbr data set----
#https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html?fbclid=IwAR3hV7VtVtUjluJyXz_6sJiG_L9A-kroZ2lhbS9psIe7nizu-SVGnVR8j7s
# ----> fgsea tutorial----
#http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html?fbclid=IwAR20mU0BDXHuq-4D-ZbvT3IbEas0port9hw3Le46YlNZKhQJaOjiSyVF3J8
# ----> DESeq analysis by fgsea----
#https://stephenturner.github.io/deseq-to-fgsea/
# ----> RNA-seq analysis----
#https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-gene-set-testing.nb.html

