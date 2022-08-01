setwd("D:/Research/Script_Prognosis_workflow/Gene expression")
load("GSEA.RData")

# ----> Combine result----
result[[1]][["tranferP"]] <- ifelse(result[[1]]$NES<0,log10(result[[1]]$padj),-log10(result[[1]]$padj))
combine <- data.frame(x = result[[1]][["tranferP"]],row.names = result[[1]]$pathway)
names(combine)[which(names(combine)=="x")] <- names(result)[1]
for (i in 2:length(result)) {
  result[[i]][["tranferP"]] <- ifelse(result[[i]]$NES<0,log10(result[[i]]$padj),-log10(result[[i]]$padj))
  combine$x <- result[[i]][["tranferP"]][match(rownames(combine),result[[i]]$pathway)]
  names(combine)[which(names(combine)=="x")] <- names(result)[i]
}

# ----> Remove non significant----
## abs(p value) combine
combinep <- data.frame(x = abs(result[[1]]$padj),row.names = result[[1]]$pathway)
names(combinep)[which(names(combinep)=="x")] <- names(result)[1]
for (i in 2:length(result)) {
  combinep$x <- result[[i]]$pval[match(rownames(combinep),result[[i]]$pathway)]
  names(combinep)[which(names(combinep)=="x")] <- names(result)[i]
}
check_rm <- data.frame(p_maen = rowMeans(combinep,na.rm =T))
check_rm$check <- check_rm$p_maen<0.05
# ----> Visualization----
library("pheatmap")
draw <- combine
pheatmap(draw,cutree_rows = 20,show_colnames = T,show_rownames = F,clustering_distance_cols="euclidean",
         main = "GSEA functional annotation about prognosis genus",
         legend_breaks = c((-4), (-2), 0, 2,4, max(draw)),
         legend_labels = c("-4", "-2", "0", "2","4", "-log10 p-value\n"))


out <- pheatmap(draw,cutree_rows = 20,show_colnames = T,show_rownames = F,fontsize_row = 3,clustering_distance_cols="euclidean")
rownames(draw[out$tree_row[["order"]],])
colnames(draw[,out$tree_col[["order"]]])
sort(cutree(out$tree_row, k=20))

cluster = which(sort(cutree(out$tree_row, k=20))==1)
cluster_rm <- names(cluster)

cut_heat_map <- draw[-which(rownames(draw)%in%cluster_rm),]
png(file = "D:/Research/PlotNew/FunctionAnnotationAndGenus.png",   # The directory you want to save the file in
    width = 600, # The width of the plot in inches
    height = 600, units = "px")
par(mar=c(10,10,3,3))
pheatmap(cut_heat_map,cutree_rows = 5,show_rownames = F,fontsize_col = 15 ,show_colnames = T,
         main = "Genus associated functions")
dev.off()

out1 <- pheatmap(cut_heat_map,cutree_rows = 5,show_rownames = F,fontsize_row =3 ,show_colnames = T,
         main = "Genus associated functions")

cluster_new <- names(which(sort(cutree(out1$tree_row, k=5))==5))

cluster1 <- draw[which(rownames(draw)%in%cluster_new),]
cluster2 <- draw[which(rownames(draw)%in%cluster_new),]
cluster3 <- draw[which(rownames(draw)%in%cluster_new),]
cluster4 <- draw[which(rownames(draw)%in%cluster_new),]
cluster5 <- draw[which(rownames(draw)%in%cluster_new),]

png("D:/Research/PlotNew/GSEA_Cluster5.png",width = 800,height = 800)
pheatmap(cluster5,fontsize_col = 15,show_rownames = F,cluster_rows = T)
dev.off()

rownames(cluster4)
write(rownames(cluster1),"cluster1.txt")
