####Function####
AssignChrToNum <- function(df,column){
  # Transfer the different chracters of clinical traits to differet number
  temp <- unique(df[[column]])
  temp <- temp[!is.na(temp)]
  temp <- data.frame(type=sort(temp),num=c(0:(length(temp)-1)))
  df[[column]] <- ifelse(is.na(df[[column]]),(-1),temp$num[match(df[[column]],temp$type)])
  result = list("DF"=df,"AssignTable"=temp)
  return(result)
}

####Main####
# ----> Import library----
library(readxl)
library(ggplot2)
# ----> Set work path and Load file----
setwd("D:/Research/Data/Raw")
# The file content the UUID and matched Barcode
tcga <-  read.csv("TCGA.csv", header = TRUE)
tcga_cdr <-  readxl::read_excel("D:/Research/Data/Raw/Clinical meta from CDR/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
metadata <- read.csv("D:/Research/Data/Raw/LUAD_origin/Metadata_LUAD_WGS_138sample.csv",header = T,row.names = 1)

# ----> Filter the sample----
#Grab the 138 sample from TCGA barcode and give sample name
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%tcga$barcode),]
match_barcode <- match(tcga_select$bcr_patient_barcode,tcga$barcode)
tcga_select$X <- tcga$X[match_barcode]

# ----> Filter the column about available clinical traits----
# Combine Metadata and clinical data by sample name
tcga_select <- tcga_select[order(tcga_select$X),]
metadata <- metadata[order(rownames(metadata)),]
clinical_trait <- cbind(tcga_select,metadata)
# Remain column set:
remain <- c("age_at_diagnosis","days_to_death","gender",
            "vital_status","pathologic_t_label","pathologic_n_label",
            "pathologic_stage_label","protion_weight","stage_severe",
            "age_at_initial_pathologic_diagnosis","race","ajcc_pathologic_tumor_stage",
            "last_contact_days_to","death_days_to","new_tumor_event_type",
            "new_tumor_event_dx_days_to","treatment_outcome_first_course")
clinical_trait <- clinical_trait[,which(names(clinical_trait)%in%remain)]
## Compare the replicate column(Vital status,stage label)
#Vital status
status <- data.frame(table(clinical_trait$vital_status))
status2 <- data.frame(table(clinical_trait$vital_status.1))
status$type <- "TCGA-CDR"
status2$type <- "Meta"
status_plot <- rbind(status,status2)
names(status_plot) <- c("status","Freq","type")
# Plot colum
ggplot2::ggplot(data = status_plot,aes(x=type,
                                       y=Freq))+
  geom_col(aes(fill=status),width=0.3)+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15))+
  ggplot2::labs(title = "Two clinical data compare",
       subtitle = "The difference of Vital status",
       x="Clinical data source",
       y="Count")
ggplot2::ggsave("D:/Research/PlotNew/26.Status_compare.png",dpi=300)
# Stage label
stage <- data.frame(table(clinical_trait$ajcc_pathologic_tumor_stage))
stage2 <- data.frame(table(clinical_trait$pathologic_stage_label))
stage$type <- "TCGA-CDR"
stage2$type <- "Meta"
stage_plot <- rbind(stage,stage2)
names(stage_plot) <- c("stage","Freq","type")
# Plot colum
ggplot2::ggplot(data = stage_plot,aes(x=type,
                                       y=Freq))+
  geom_col(aes(fill=stage),width=0.3)+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15))+
  ggplot2::labs(title = "Two clinical data compare",
       subtitle = "The difference of stage level",
       x="Clinical data source",
       y="Count")
ggplot2::ggsave("D:/Research/PlotNew/27.stage_compare.png",dpi=300)
# Remove replicate column
clinical_trait <- clinical_trait[,-which(names(clinical_trait)%in%c("gender.1","vital_status.1","race.1","age_at_diagnosis"))]
# Remain this object

# ----> Transfer the chr column to num----
transfer_trait <- clinical_trait
summary(transfer_trait)
# ---->> Gender----
result <- AssignChrToNum(transfer_trait,"gender")
transfer_trait <- result$DF
gender_a <- result$AssignTable
# ---->> Race----
result <- AssignChrToNum(transfer_trait,"race")
transfer_trait <- result$DF
rece_a <- result$AssignTable
# ---->> AJCC stage----
result <- AssignChrToNum(transfer_trait,"ajcc_pathologic_tumor_stage")
transfer_trait <- result$DF
ajccStage_a <- result$AssignTable
# ---->> Status----
result <- AssignChrToNum(transfer_trait,"vital_status")
transfer_trait <- result$DF
status_a <- result$AssignTable
#---->> New tumor event----
result <- AssignChrToNum(transfer_trait,"new_tumor_event_type")
transfer_trait <- result$DF
newTumorevent_a <- result$AssignTable
#---->> Treatment outcome first course----
result <- AssignChrToNum(transfer_trait,"treatment_outcome_first_course")
transfer_trait <- result$DF
treatment_a <- result$AssignTable
# ---->> Pathologic T label----
result <- AssignChrToNum(transfer_trait,"pathologic_t_label")
transfer_trait <- result$DF
Tlabel_a <- result$AssignTable
# ---->> Pathologic N label----
result <- AssignChrToNum(transfer_trait,"pathologic_n_label")
transfer_trait <- result$DF
Nlabel_a <- result$AssignTable
# ---->> Pathologic stage label----
result <- AssignChrToNum(transfer_trait,"pathologic_stage_label")
transfer_trait <- result$DF
Stage_a <- result$AssignTable

transfer_trait <- transfer_trait[,-which(names(clinical_trait)%in%c("ajcc_pathologic_tumor_stage","death_days_to","last_contact_days_to","new_tumor_event_dx_days_to","days_to_death"))]

write.csv(transfer_trait,"D:/Research/Data/Redo/Clinical Traits_AfterClean.csv")
