#### Main####
# ----> Import library and work path----
setwd("D:/Research/Data")
library(taxize)
# ----> Read files ----
bacteria <- read.csv("D:/Research/Data/Raw/COVID microbiome/relab_taxonomy_bacteria_metagenome_S.txt",sep = "\t",row.names = 1)

# ----> COVID metagenome----
## Transfer to propotional
scaled.bacteria <- bacteria
for(i in 1:length(names(bacteria))){
  for (j in 1:length(rownames(bacteria))) {
    scaled.bacteria[j,i] <- bacteria[j,i]/sum(bacteria[,i])
  }
}
## Calculate the prevalence
bacteria_select <- data.frame(name = rownames(bacteria),
                              prevelance = rep(0,4424))
for (i in 1:length(bacteria_select$prevelance)) {
  count <- 0
  for(j in 1:length(scaled.bacteria[i,])){
    if(scaled.bacteria[i,j] > 0){
      count <- count + 1
    }
  }
  prevelance <- count/189
  bacteria_select$prevelance[i] <- prevelance
}
bacteria_select$prevelance_20 <- ifelse(bacteria_select$prevelance > 0.2 ,T,F)
bacteriafilter <- bacteria_select[which(bacteria_select$prevelance_20 == T),]
bacteriafilter <- bacteriafilter$name
bacteria_list_metagenome <- taxize::tax_name(sci = bacteriafilter, get = c("genus","family","phylum"), db = "ncbi")
write.csv(bacteria_list_metagenome,"D:/Research/Data/Redo/Covid_bacteria_metagenome.csv")

# ----> COVID transcriptome----
bacteria2 <- read.csv("D:/Research/Data/Raw/COVID microbiome/BAL_bacteria_relative_abundance.csv")
bacterianame2 <- rownames(bacteria2)
bacteria_list_metatranscriptome <- tax_name(sci = bacterianame2, get = c("genus","family","phylum"), db = "ncbi")

write.csv(bacteria_list_metatranscriptome,"D:/Research/Data/Redo/Covid_bacteria_metatranscriptome.csv")

# ----> My Prognosis genus ----
prognosis <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv")
## Rename genus
for (j in 1:length(prognosis$Genus)) {
    temp <- strsplit(prognosis$Genus[j],split = "_") # separate string by [.]
    prognosis$Genus[j] <- temp[[1]][length(temp[[1]])]
}
bacteria_list_prognosis <- tax_name(sci = prognosis$Genus, get = c("genus","family","phylum"), db = "ncbi")

# ----> miRNA genus----
miRNA <- read.csv("D:/Research/Data/Raw/miRNA microbiome/LUAD_relative_abundance_miRNA.txt",sep = "\t")
miRNA_bacteria <- tax_name(sci = rownames(miRNA), get = c("genus","family","phylum"), db = "ncbi")
write.csv(miRNA_bacteria,"D:/Research/Data/Redo/Bacteria_miRNA.csv")

# ----> My cutoff genus----
my_data <- read.csv("D:/Research/Data/Raw/LUAD_origin/LUAD_WGS_abundance_138sample.csv")
my_data <- colnames(my_data)
my_data <- my_data[2:1285]
for (j in 1:length(my_data)) {
  temp <- strsplit(my_data[j],split = "_") # separate string by [.]
  my_data[j] <- temp[[1]][length(temp[[1]])]
}
mydata_bacteria <- tax_name(sci = my_data, get = c("genus","family","phylum"), db = "ncbi")
write.csv(mydata_bacteria,"D:/Research/Data/Redo/Bacteria_mydata.csv")

# ----> Colon genus----
colon_microbiome <- read.csv("D:/Research/Data/Raw/Colon microbiome/GMREPO_relative_abundance_of_all_species_genus_in_all_phenotypes_summary.tsv",sep = "\t")
colon_microbiome <- colon_microbiome[which(colon_microbiome$taxonomic.rank=="genus"),]
colon_microbiome <- colon_microbiome$scientific.name
colon_microbiome <- unique(colon_microbiome)

# ----> Create check exist matrix----
prognosis$COVID_metagenome <- ifelse(prognosis$Genus %in% bacteria_list_metagenome$genus,T,F)
prognosis$COVID_metatranscriptome <- ifelse(prognosis$Genus%in%bacteria_list_metatranscriptome$genus,T,F)
prognosis$miRNA <- ifelse(prognosis$Genus%in%miRNA_bacteria$genus,T,F)
prognosis$colon <- ifelse(prognosis$Genus%in%colon_microbiome,T,F)
prognosis$check_T <- rowSums(prognosis == T)
prognosis$check_exist <- ifelse(prognosis$check_T >= 3,T,F)
write.csv(prognosis,"D:/Research/Data/Redo/Prognosis_check_genus_exist.csv")
