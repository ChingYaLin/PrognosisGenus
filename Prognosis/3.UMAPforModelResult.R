####Package####
library(plyr)
library(umap)
library(survival)
library(survminer)
####Function####
RenameGenus <- function(df){
  colname_genus <- colnames(df)
  for (j in 1:length(colname_genus)) {
    if (colname_genus[j]!='X'){
      temp <- strsplit(colname_genus[j],split = "\\.")
      names(df)[names(df)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
    }
  }
  temp <- gsub("g_","",names(df))
  temp <- gsub("[_]","",temp)
  names(df) <- temp
  return(df)
}

TransferTheCategory <- function(df,category,informtable){
  temp <- df
  temp$category <- informtable[[category]][which(informtable$query%in%rownames(temp))]
  temp <- plyr::ddply(temp,"category",numcolwise(sum))
  temp <- temp[-which(is.na(temp$category)),]
  rownames(temp) <- temp$category
  temp <- temp[,-1]
  temp <- temp[,order(names(temp))]
}

TransferProportion <- function(df){
  df <- data.frame(t(df))
  df <- df/rowSums(df)
  df <- data.frame(t(df))
  return(df)
}
####UMAP####
# ----> Set path and read data----
setwd("D:/Research/Data/Redo")
## Reading data
umap_test_all <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
prognosis_genus <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv",header = T)
five_year <- read.csv("Five_year_prognosis.csv",header = T,row.names = 1)
## Grab the 23 prognosis genus
prognosis_genus <- prognosis_genus$Genus
umap_test <- umap_test_all[,which(names(umap_test_all)%in%prognosis_genus)]
## Rename the genus column
umap_test_all <- RenameGenus(umap_test_all)
umap_test<-RenameGenus(umap_test)

# ----> umap data frame----
genus.umap = umap::umap(umap_test, )
genus.umap_all = umap::umap(umap_test_all, )

# ----> Plot preprocess----
## 23 genus
umap_plot_df <- data.frame(genus.umap$layout)
umap_plot_df<-umap_plot_df[order(rownames(umap_plot_df)),]
umap_plot_df$vital <- five_year$OS
umap_plot_df$cluster <- ifelse(umap_plot_df$X2>2.5,"cluster 1","cluster 2")
umap_plot_df$OS.time <- five_year$os.time
## 1146 genus
umap_plot_df_all <- data.frame(genus.umap_all$layout)
umap_plot_df_all<-umap_plot_df_all[order(rownames(umap_plot_df_all)),]
umap_plot_df_all$vital <- five_year$OS
umap_plot_df_all$cluster <- ifelse(umap_plot_df$X2>2.5,"cluster 1","cluster 2")
umap_plot_df_all$x1_order23 <- umap_plot_df$X1[match(rownames(umap_plot_df_all),rownames(umap_plot_df))]

####Visualization (23 genus)####
# ----> KM-plot for cluster 1 and 2----
par(mar=c(10,10,3,3))
new <- umap_plot_df
new$OS.time <- ifelse(new$OS.time >= 2000, 2000, new$OS.time)
new$OS <- ifelse(new$vital=='Dead',1,0)
new$OS <- ifelse(new$OS.time >= 2000, 0,new$OS)
new$cluster <- factor(new$cluster)

fit_t <- survfit(Surv(OS.time, OS) ~ cluster, data = new)
ggsurvplot(fit_t, data = new, pval = TRUE,title = 'KM-plot for different clusters',
           conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           surv.median.line = "hv",
           risk.table.col = "cluster",
           risk.table.y.text.col = T,
           palette = c("#E7B800", "#2E9FDF"),
           xlab = "Time in days",
           ggtheme = theme_light())

ggsave("D:/Research/PlotNew/Prognosis/KMplotForCluster.pdf", width=10, height=10)


# ----> Vital----
ggplot2::ggplot(umap_plot_df,aes(x = X1,y = X2,color = vital)) +
  ggplot2::geom_point(size=3,alpha = 0.5)+
  ggplot2::labs(title = "Umap (23 genus)",
       subtitle = "After five year",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_23_fiveYear(vital).png",dpi = 450)

# ----> Cluster----
ggplot2::ggplot(umap_plot_df,aes(x = X1,y = X2,color = cluster)) +
  ggplot2::geom_point(size=3,alpha = 0.5)+
  ggplot2::labs(title = "Umap (23 genus)",
                subtitle = "After five year",
                x="X1(UMAP)",
                y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_23_fiveYear(cluster).png",dpi = 450)

# ----> x1 of 23 genus----
ggplot2::ggplot(umap_plot_df,aes(x = X1,y = X2,color = X1)) +
  ggplot2::geom_point(size=4,alpha = 0.5)+
  ggplot2::scale_colour_gradient2(low = "blue", mid = "snow3",high = "red", na.value = NA)+
  ggplot2::labs(title = "Umap (23 genus)",
                subtitle = "After five year",
                x="X1(UMAP)",
                y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_23_fiveYear(X1).png",dpi = 450)

# ----> Odds ratio----
class_matrix <- data.frame("Dead" = c(sum(umap_plot_df$cluster=="cluster 1" & umap_plot_df$vital =="Dead"),
                                      sum(umap_plot_df$cluster=="cluster 2" & umap_plot_df$vital =="Dead")),
                           "Alive"= c(sum(umap_plot_df$cluster=="cluster 1" & umap_plot_df$vital =="Alive"),
                                      sum(umap_plot_df$cluster=="cluster 2" & umap_plot_df$vital =="Alive")),
                           row.names = c("cluster1","cluster2"))
png("D:/Research/PlotNew/Cluster distribution.png",
    width = 600, height = 600, units = "px")
par( cex.lab=1.5, cex.main=2, cex.sub=1.5)
mosaicplot(class_matrix, color = TRUE,
           main = "Five year after",cex.axis=1.5)  
dev.off()
fisher.test(class_matrix)


####Visualization (1146 genus)####
# ----> Vital----
ggplot2::ggplot(umap_plot_df_all,aes(x = X1,y = X2,color = vital)) +
  ggplot2::geom_point(size=3,alpha = 0.5)+
  ggplot2::labs(title = "Umap (1146 genus)",
       subtitle = "After fives year",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_1146_fiveYear(vital).png",dpi = 450)

# ----> Cluster----
ggplot2::ggplot(umap_plot_df_all,aes(x = X1,y = X2,color = cluster)) +
  ggplot2::geom_point(size=3,alpha = 0.5)+
  ggplot2::labs(title = "Umap (1146 genus)",
                subtitle = "After fives year",
                x="X1(UMAP)",
                y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_1146_fiveYear(cluster).png",dpi = 450)

# ----> X1 of 23 genus----
ggplot2::ggplot(umap_plot_df_all,aes(x = X1,y = X2,color = x1_order23)) +
  ggplot2::geom_point(size=3,alpha = 0.6)+
  ggplot2::scale_colour_gradient2(low = "blue", mid = "snow3",high = "red", na.value = NA)+
  ggplot2::labs(title = "Umap (1146 genus)",
                subtitle = "After five year",
                x="X1(UMAP)",
                y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/UMAP_1146_fiveYear(X1).png",dpi = 450)


####Label by phylum####
# ----> Read data and preprocess----
raw_data <- read.csv("D:/Research/Data/Raw/Paper origin data/Kraken-TCGA-Raw-Data-17625-Samples.csv", header = TRUE,row.names = 1)
raw_data <- RenameGenus(raw_data)
raw_data <- raw_data[,which(names(raw_data)%in%my_bacter$query)]
my_bacter <- read.csv("Bacteria_mydata.csv",row.names = 1)
WGS_genus <- raw_data[which(rownames(raw_data)%in%rownames(genus_matrix)),]
WGS_genus <- data.frame(t(WGS_genus))
WGS_genus <- WGS_genus[,order(names(WGS_genus))]
WGS_p <- TransferTheCategory(WGS_genus,"phylum",my_bacter)
WGS_p <- TransferProportion(WGS_p)
WGS_p <- data.frame(t(WGS_p))
genus_mean <- rowMeans(WGS_p)
genus_mean <- sort(genus_mean)

# ----> Add phylum to umap information----
## 23 genus
umap_plot_df$p_Proteobacteria <- WGS_p$Proteobacteria[match(rownames(umap_plot_df),rownames(WGS_p))]
umap_plot_df$p_Actinobacteria <- WGS_p$Actinobacteria[match(rownames(umap_plot_df),rownames(WGS_p))]
umap_plot_df$p_Chlamydiae <- WGS_p$Chlamydiae[match(rownames(umap_plot_df),rownames(WGS_p))]
umap_plot_df$p_Firmicutes <- WGS_p$Firmicutes[match(rownames(umap_plot_df),rownames(WGS_p))]
umap_plot_df$p_Bacteroidetes <- WGS_p$Bacteroidetes[match(rownames(umap_plot_df),rownames(WGS_p))]
## all genus
umap_plot_df_all$p_Proteobacteria <- WGS_p$Proteobacteria[match(rownames(umap_plot_df_all),rownames(WGS_p))]
umap_plot_df_all$p_Actinobacteria <- WGS_p$Actinobacteria[match(rownames(umap_plot_df_all),rownames(WGS_p))]
umap_plot_df_all$p_Chlamydiae <- WGS_p$Chlamydiae[match(rownames(umap_plot_df_all),rownames(WGS_p))]
umap_plot_df_all$p_Firmicutes <- WGS_p$Firmicutes[match(rownames(umap_plot_df_all),rownames(WGS_p))]
umap_plot_df_all$p_Bacteroidetes <- WGS_p$Bacteroidetes[match(rownames(umap_plot_df_all),rownames(WGS_p))]

#----> Visualization----
## 23 genus
ggplot(umap_plot_df,aes(x = X1,y = X2,color = p_Bacteroidetes))+
  geom_point(size=4,alpha = 0.8)+
  scale_colour_gradient(low = "steelblue1",high = "tomato", na.value = NA)+
  labs(title = "Umap (23 genus)",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))

## all genus
ggplot(umap_plot_df_all,aes(x = X1,y = X2,color = p_Bacteroidetes))+
  geom_point(size=4,alpha = 0.8)+
  scale_colour_gradient(low = "steelblue1",high = "tomato", na.value = NA)+
  labs(title = "Umap (1146 genus)",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))

####Label by genus####
# ----> Add genus adundance----
umap_plot_df <- cbind(umap_plot_df,umap_test[match(rownames(umap_plot_df), rownames(umap_test)),])
umap_plot_df_all <- cbind(umap_plot_df_all,umap_test[match(rownames(umap_plot_df_all), rownames(umap_test)),])

#----> Visualization----
## 23 genus
ggplot(umap_plot_df,aes(x = X1,y = X2,color = g__Pyramidobacter))+
  geom_point(size=4,alpha = 0.8)+
  scale_colour_gradient(low = "steelblue1",high = "tomato", na.value = NA)+
  labs(title = "Umap (23 genus)",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))

## all genus
ggplot(umap_plot_df_all,aes(x = X1,y = X2,color = g__Hypovirus))+
  geom_point(size=4,alpha = 0.8)+
  scale_colour_gradient(low = "steelblue1",high = "tomato", na.value = NA)+
  labs(title = "Umap (1146 genus)",
       x="X1(UMAP)",
       y="X2(UMAP)")+
  ggplot2::theme(plot.title = element_text(size = 20),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 plot.subtitle = element_text(size = 15))
