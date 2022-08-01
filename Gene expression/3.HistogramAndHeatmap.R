####Package####
library(tibble)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(reshape2)
####Functions####
RenameGenus <- function(df){
  # The df column name is genus
  genus_name <- colnames(df)
  for (j in 1:length(genus_name)) {
    temp <- strsplit(genus_name[j],split = "\\.")
    names(df)[names(df)==genus_name[j]] <- temp[[1]][length(temp[[1]])]
  }
  return(df)
}

####Set path and Load data####
# ----> Set work path----
setwd("D:/Research/Data/Redo")
# ----> Load data----
## Correlation data of gene and prognosis genus
data <- read.csv("Correlation_SCC_origin(prognosis).csv", header = TRUE,row.names = 1)
ran_data <- read.csv("Correlation_SCC_random(prognosis).csv", header = TRUE,row.names = 1)

# ----> Rename genus----
data <- RenameGenus(data)
ran_data <- RenameGenus(ran_data)
####Prepare plot data####
## Construct the value to one column that can plot by histogram
DF <- tibble::rownames_to_column(data, var = "Xval")
DF <- DF %>% tidyr::pivot_longer(cols = 2:length(names(DF)), names_to = "Source", values_to = "Value")
DF_ran <- tibble::rownames_to_column(ran_data, var = "Xval")
DF_ran <- DF_ran %>% tidyr::pivot_longer(cols = 2:length(names(DF_ran)), names_to = "Source", values_to = "Value")

####Histogram####
# ----> Origin data----
ggplot(DF)+
  ggplot2::geom_histogram(aes(x=Value),binwidth = 0.05,color = "royalblue",fill="slategray1")+
  ggplot2::labs(title = "Correlation distribution",subtitle = "15323 genes with 23 genus",
       x="Spearman Correlation coefficient")+
  ggplot2::theme_classic()+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/SCC distribution(origin).png",dpi = 500)
# ----> Random permutation----
ggplot(DF_ran)+
  ggplot2::geom_histogram(aes(x=Value),binwidth = 0.05,color = "seagreen",fill="lightgreen")+
  ggplot2::labs(title = "Correlation distribution",subtitle = "15323 genes with 23 genus",
                x="Spearman Correlation coefficient")+
  ggplot2::theme_classic()+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/SCC distribution(random).png",dpi = 500)

####Heatmap####
# ----> Origin data----
origin <- scale(data)
ord <- hclust(dist(origin, method = "euclidean"), method = "ward.D" )$order
pd <- as.data.frame(origin)
pd$Gene <- rownames(pd)
pd.m <- reshape2::melt( pd, id.vars = "Gene", variable.name = "Genus" )
pd.m$Genus <- factor( pd.m$Genus, levels = colnames(origin), labels = colnames(origin))
pd.m$Gene <- factor( pd.m$Gene, levels = rownames(origin)[ord],  labels =rownames(origin)[ord])

ggplot( pd.m, aes(Gene, Genus) ) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue", high="red")+
  ggplot2::labs(title ="SCC Heatmap(origin)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank())

ggplot2::ggsave("D:/Research/PlotNew/SCC Heatmap(origin).png",dpi = 500)  

#heatmap(t(data),main = 'Heatmap(origin)')
#heatmap(t(ran_data),main = 'Heatmap(random)')
# ----> Origin data----
random <- scale(ran_data)
ord <- hclust(dist(random, method = "euclidean"), method = "ward.D" )$order
pd <- as.data.frame(random)
pd$Gene <- rownames(pd)
pd.m <- reshape2::melt( pd, id.vars = "Gene", variable.name = "Genus" )
pd.m$Genus <- factor( pd.m$Genus, levels = colnames(random), labels = colnames(random))
pd.m$Gene <- factor( pd.m$Gene, levels = rownames(random)[ord],  labels =rownames(random)[ord])

ggplot( pd.m, aes(Gene, Genus) ) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue", high="red")+
  ggplot2::labs(title ="SCC Heatmap(random)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank())

ggplot2::ggsave("D:/Research/PlotNew/SCC Heatmap(random).png",dpi = 500)  

####Density####
# ----> Two distribution compare----
DF$class <- "Origin"
DF_ran$class <- "Random permutation"
density_plot <- rbind(DF,DF_ran)

ggplot(density_plot)+
  geom_density(aes(x=Value,color = class,fill = class),alpha = 0.2)+
  labs(title = "SCC Density")+
  ggplot2::theme_minimal()+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank())
ggplot2::ggsave("D:/Research/PlotNew/SCC Density(Compare).png",dpi = 500)
  

## The another histogram
#count_origin <- table(cut(DF$Value,breaks = c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)))
#barplot(count_origin)
#count_random <- table(cut(DF_ran$Value,breaks = c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)))
#barplot(count_random)

# ----> Minus result----
#plot the difference count between the origin and random group
count_dif <- count_origin-count_random
png("D:/Research/PlotNew/SCC Histogram(Minus).png",
    width = 800, height = 600, units = "px")

barplot(count_dif,col = "tan1",border="darkorange",
        main="The barplot of origin minus random",
        xlab="SCC interval",ylab="count",
        cex.axis = 1.2,cex.names = 1.2,cex.lab = 1.2,cex.main=1.5)
dev.off()
summary_count <- data.frame(cbind(count_origin,count_random,count_dif))
write.csv(summary_count,"SCC distribution compare.csv")
