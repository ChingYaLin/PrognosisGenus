####UMAP####
# ----> Set file path and import library----
library(umap)
library(ggplot2)
setwd("D:/Research/Data/Redo")
# ----> Read data----
## Clinical data
clinical_trait <- read.csv("Clinical Traits_AfterClean.csv",
                           header = T,row.names = 1)
## Genus abundance matrix
genus_abundance <- read.csv("voom_svm_clean_genus(CutLow).csv",
                            header = T,row.names = 1)
# ----> UMAP----
## Dimensional reduction by UMAP
genus.umap = umap::umap(genus_abundance, )
## Create the plot matrix
umap_plot_df <- data.frame(genus.umap$layout)
match_index <- match(rownames(umap_plot_df),rownames(clinical_trait))
## Add clinical label to plot matrix
umap_plot_df$gender <- clinical_trait$gender[match_index]
umap_plot_df$race <- clinical_trait$race[match_index]
umap_plot_df$status <- clinical_trait$vital_status[match_index]
umap_plot_df$new_tumor_event <- clinical_trait$new_tumor_event_type[match_index]
umap_plot_df$treatment_first_course <- clinical_trait$treatment_outcome_first_course[match_index]
umap_plot_df$T_label <- clinical_trait$pathologic_t_label[match_index]
umap_plot_df$N_label <- clinical_trait$pathologic_t_label[match_index]
umap_plot_df$stage <- clinical_trait$pathologic_stage_label[match_index]
umap_plot_df$servere <- clinical_trait$stage_severe[match_index]
  
# ----> Visualization----
## gender
ggplot(umap_plot_df,aes(x = X1,y = X2,color = as.factor(gender))) +
  ggplot2::geom_point(size=4,alpha = 0.5)+
  labs(title = "Umap (gender)",subtitle = "whole genus (1146)",x="X1(UMAP)",y="X2(UMAP)")+
  scale_color_manual(label=c("Female","Male"),
                     values = c("#FF6666","#3399FF"),
                     name = "Gender")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text =element_text(size = 15),
                 legend.title = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/28.UMAP(gender)_1146.png",dpi=500)

## status
ggplot(umap_plot_df,aes(x = X1,y = X2,color = as.factor(status))) +
  ggplot2::geom_point(size=4,alpha = 0.5)+
  labs(title = "Umap (status)",subtitle = "whole genus (1146)",x="X1(UMAP)",y="X2(UMAP)")+
  scale_color_manual(label=c("Alive","Dead"),
                     values = c("#FF6666","#3399FF"),
                     name = "Status")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text =element_text(size = 15),
                 legend.title = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/29.UMAP(status)_1146.png",dpi=500)

## new tumor event
ggplot(umap_plot_df,aes(x = X1,y = X2,color = as.factor(new_tumor_event))) +
  ggplot2::geom_point(size=4,alpha = 0.5)+
  labs(title = "Umap (new tumor event)",subtitle = "whole genus (1146)",x="X1(UMAP)",y="X2(UMAP)")+
  scale_color_manual(label=c("NA","Distant Metastasis","Distant Metastasis|New Primary Tumor",
                             "Locoregioal Recurrence","Locoregioal Recurrence|Distant Metastasis"),
                     values = c("#FF6666","#3399FF","#FF9933","#00CC00","#B266FF"),
                     name = "Event")+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text =element_text(size = 15),
                 legend.title = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/30.UMAP(new tumor event)_1146.png",dpi=500)


