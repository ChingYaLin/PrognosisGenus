####Package####
library("pheatmap")
library(tibble)
library(tidyr)
library("ggpubr")
library(igraph)
####Functions####
RenameGenus <- function(df){
  colname_genus <- colnames(df)
  for (j in 1:length(colname_genus)) {
    temp <- strsplit(colname_genus[j],split = "\\.")
    names(df)[names(df)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
  }
  return(df)
}

RenameGenussingle <- function(vector){
  for (j in 1:length(vector)) {
    temp <- strsplit(vector[j],split = "\\__")
    vector[j] <- temp[[1]][length(temp[[1]])]
  }
  return(vector)
}

RenameImmune <- function(df){
  rowname_I <- rownames(df)
  for (i in 1:length(rowname_I)) {
    temp <- strsplit(rowname_I[i],split = "_")
    rownames(df)[rownames(df)==rowname_I[i]] <- temp[[1]][1]
  }
  return(df)
}

RenameImmunesingle <- function(vector){
  for (i in 1:length(vector)) {
    temp <- strsplit(vector[i],split = "_")
    vector[i] <- temp[[1]][1]
  }
  return(vector)
}
####Main####
# ----> Set path and read files----
# Set path
setwd("D:/Research/Data/Redo")
# SCC table (immune vs genus)
immune_genus_scc <- read.csv("ImmuneCorrelation_SCC_origin(prognosis).csv", header = TRUE,row.names = 1)
## Prognosis genus
prognosis_genus <- read.csv("Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv")

# ----> Rename the table ----
##Filter the prognosis genus 
plot_table <- immune_genus_scc[,which(names(immune_genus_scc)%in%prognosis_genus$Genus)]
## Rename the genus name
plot_table <- RenameGenus(plot_table)
## Rename the immune name
plot_table <- RenameImmune(plot_table)
plot_table <- na.omit(plot_table)


#### Visualization####
# ----> Heatmap----
png(file = "D:/Research/PlotNew/Immune/GenusImmuneHeatmap.png",   # The directory you want to save the file in
    width = 800, height = 800, units = "px")
pheatmap::pheatmap(plot_table,cutree_cols = 2,
                   main = "The correlation of single genus and immune infiltration",
                   fontsize = 15)
dev.off()

# ----> Density----
data <- read.csv("ImmuneCorrelation_SCC_origin(prognosis).csv", header = TRUE,row.names = 1)
ran_data <- read.csv("ImmuneCorrelation_SCC_random(prognosis).csv", header = TRUE,row.names = 1)
DF <- rownames_to_column(data, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")
DF_ran <- rownames_to_column(ran_data, var = "Xval")
DF_ran <- DF_ran %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")
## wo distribution compare
DF$class <- "Origin"
DF_ran$class <- "Random permutation"
density_plot <- rbind(DF,DF_ran)

ggplot(density_plot)+
  geom_density(aes(x=Value,color = class,fill = class),alpha = 0.2)+
  labs(title = "SCC density")+
  ggplot2::theme_minimal()+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank())
ggplot2::ggsave("D:/Research/PlotNew/Immune/SCC Density_Immune(Compare).png",dpi = 500)

# ----> Network----
DF <- DF[which(DF$Source%in%prognosis_genus$Genus),]

select_combination <- DF[which(abs(DF$Value)>=0.15),]
select_combination <- select_combination[order(select_combination$Source),]
select_combination$Source <- RenameGenussingle(select_combination$Source)
select_combination$Xval <- RenameImmunesingle(select_combination$Xval)

nodes <- unique(c(select_combination$Xval,select_combination$Source))
nodes <- data.frame(id = nodes)
nodes$degree <- 0 
test <- rbind(data.frame(table(select_combination$Xval)),data.frame(table(select_combination$Source)))
matchh <-match(nodes$id,test$Var1)
nodes$degree <- test$Freq[matchh]
nodes$class <- 1
nodes$class[15:31] <-2
nodes$module_num <- 0
nodes$name <- nodes$id


net <- graph_from_data_frame(d=select_combination, vertices=nodes, directed=F) 
plot(net, edge.arrow.size=.4,vertex.label=nodes$id,layout=layout_components)


net <- simplify(net, remove.multiple = F, remove.loops = T) 
colrs <- c("#FF9933", "gold")
V(net)$color <- colrs[V(net)$class]
# Set node size based on audience size:
V(net)$size <- V(net)$degree*5

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"

# Set edge width based on weight:
E(net)$width <- abs(E(net)$Value)*10

#change arrow size and edge color:
edge.color <- ifelse(E(net)$Value>0,"#FF6666","#66B2FF")


pdf(file = "D:/Research/PlotNew/Immune/Network.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)

plot(net,vertex.label.cex=.8,layout=layout.auto,edge.color=edge.color,
     vertex.label=nodes$name) 
legend(x=-1.1, y=-1.1, c("Immune cell","Genus","Positive correlation","Negative correlation"), pch=c(21, 21, NA,NA),
       lty = c(NA,NA,1,1),
       col=c("#777777","#777777","#FF6666","#66B2FF"), pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)

dev.off()


# ----> Scatter plot----
genus <- select_combination$Source

test <- table(genus)
test <- data.frame(test)
## Scatter plot
plot_data1 <- rbind(immune_data,genus_abundance)
plot_data1 <- data.frame(t(plot_data1))

ggscatter(plot_data, x = "T.cell.NK_XCELL", y = "k__Viruses.f__Baculoviridae.g__Betabaculovirus", 
          add = "loess", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab ="T.cell.NK_XCELL", ylab = "g__Betabaculovirus")+
  labs(title="Spearman correlation")
