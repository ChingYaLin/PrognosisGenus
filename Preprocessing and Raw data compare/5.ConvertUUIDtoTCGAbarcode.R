#### Main####
# ----> Set file path and import package----
setwd("D:/Research/Data/Raw/LUAD_origin")
library(TCGAutils)
# ----> Read data ----
data <- read.csv("Metadata_LUAD_RNA_138sample.csv")
uuids <- data$case_uuid
barcode <- UUIDtoBarcode(uuids, from_type = "case_id")
barcode$X <- data$X[match(barcode$case_id,data$case_uuid)]
write.csv(barcode,file = 'D:/Research/Data/Raw/TCGA_Barcode_rna.csv',row.names = F)

#Consider about how to create the TCGA file with this!!!!!!!!!!!!!!
