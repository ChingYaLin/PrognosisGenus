#### Functions####
#Count the prevalence of each genus
Prevalence_genus <- function(cutoff){
  #The cutoff is the setting how much value is non-appear in human body
  result <- data.frame(names(raw_abundance))
  result$prevalence <- 0
  for (i in 1:length(names(raw_abundance))) {
    count_have_genus <- sum(raw_abundance[i]>cutoff) #count the number of sample that value > cutoff
    prevalence_temp <- count_have_genus/length(rownames(raw_abundance)) #calculate the prevalence
    result$prevalence[which(result$names.raw_abundance.==names(raw_abundance)[i])] <- prevalence_temp #fill the prevalence of the genus into reult
  }
  return(result)
}
PrevalenceExceed20 <- function(df){
  # Checking whether the provalence of genus is exceed 20% of total samples, if yes, it would return 'Ture', otherwise, it would return 'False'
  df <- df[order(df$prevalence),]
  df$num <- c(1:1284)
  df$prevalence_exceed_20 <- ifelse(df$prevalence>0.2,T,F)
  return(df)
}
#### Main ####
# ----> Import library and set work path----
setwd("D:/Research/Data/Raw/LUAD_origin")
library(tibble)
library(tidyr)
library(ggplot2)

# ----> Load data----
## Genus abundance
raw_abundance <-  read.csv("LUAD_WGS_abundance_138sample.csv", header = TRUE,row.names = 1)
# ----> Clean data and \combine all values----
## Remove the contaminant column
raw_abundance <- raw_abundance[,-which(names(raw_abundance)%in%c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters"))]
## Whole matrix value
all_genus_value <- tibble::rownames_to_column(raw_abundance, var = "Xval")
all_genus_value <- all_genus_value %>% tidyr::pivot_longer(cols = 2:1285, names_to = "Source", values_to = "Value")

# ----> Histogram of abundance----
## The histogram of all relative abundance value, and set the cutoff represent the genus exist in the sample
png(file="D:/Research/PlotNew/18.genus_exist_cutoff2.png",
    width=800, height=600)
bin <- seq(-10,25)
hist(all_genus_value$Value,
     main = "Relative abundace values",
     xlab = "Values (Relative abundance)"
     ,ylab = "counts",breaks = bin,
     cex.lab = 1.5,cex.axis=1.5,cex.main=2)
abline(v = 2, col = "red", lwd = 3)
dev.off()
# ----> Calculate the prevalence----
# Setting the different genus abundance cutoff to calculate the prevalence of each genus
zero <- Prevalence_genus(0)
one <- Prevalence_genus(1)
two <- Prevalence_genus(2)
three <- Prevalence_genus(3)

# ----> Histogram of prevalence----
## Can change any baseline of files
png(file="D:/Research/PlotNew/19.Prevalence distribution Histogram_cutoff2.png",
    width=800, height=600)
bin <- seq(0,1,0.05)
hist(two$prevalence,
     main = "Histogram of prevalence distribution",
     xlab = "Prevalence (%)",ylab = "Counts",
     sub = "cutoff = 2",
     breaks = bin,
     cex.lab = 1.5,cex.axis=1.5,cex.main=2,cex.sub = 1.3)
abline(v = 0.2, col = "red", lwd = 3)
dev.off()
# ----> Check Prevalence > 0.2 and sort----
zero <- PrevalenceExceed20(zero)
one <- PrevalenceExceed20(one)
two <- PrevalenceExceed20(two)
three <- PrevalenceExceed20(three)

# ----> Plot the prevalence of genus----
ggplot2::ggplot(data = two,aes(x=num,y=prevalence)) +
  ggplot2::geom_line(size = 1)+
  ggplot2::geom_hline(yintercept =0.2,color = "red")+
  ggplot2::labs(title = "Prevalence of each genus",
       subtitle = paste("Cutoff=2",",the number of genus that prevalence>20%:",sum(two$prevalence_exceed_20 == T)),
       x="Genus",y="Prevalence")+
  theme_light()+
  ggplot2::theme(plot.title = element_text(size =20),
                 axis.title = element_text(size =15),
                 axis.text = element_text(size =15),
                 plot.subtitle = element_text(size = 15),
                 legend.text = element_text(size=10))
  
ggsave("D:/Research/PlotNew/20.Each genus Prevalence(Cutoff=2).png",
       width = 8, height = 6)

# ----> Write data----
write.csv(zero,file="D:/Research/Data/Redo/Cutoff_0_genus.csv",row.names = F)
write.csv(one,file="D:/Research/Data/Redo/Cutoff_1_genus.csv",row.names = F)
write.csv(two,file="D:/Research/Data/Redo/Cutoff_2_genus.csv",row.names = F)
write.csv(three,file="D:/Research/Data/Redo/Cutoff_3_genus.csv",row.names = F)
