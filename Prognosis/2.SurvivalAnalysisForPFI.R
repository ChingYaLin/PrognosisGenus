####Package####
library("survival")
library("survminer")
library("readxl")
library('dplyr')
library('pheatmap')
library('ggvenn')
####Function####
RenameGenus <- function(df){
  #Rename the microbiome name into genus
  colname_genus <- colnames(df)
  for (j in 1:length(colname_genus)) {
    if (colname_genus[j]!='X'){
      temp <- strsplit(colname_genus[j],split = "\\.")
      names(df)[names(df)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
    }
  }
  return(df)
}

Seperate3express <- function(df,genus,express){
  df[[express]] <- "Median"
  for (i in 1:length(df[[genus]])) {
    if(df[[genus]][i]<=quantile(df[[genus]])[2]){
      df[[express]][i] <- "Low"
    }
    else if(df[[genus]][i]>=quantile(df[[genus]])[4]){
      df[[express]][i] <-"High"
    }
  }
  return(df)
}

####Preprocess####
setwd("D:/Research/Data/Redo")
# ----> Read data----
tcga <-  read.csv("D:/Research/Data/Raw/TCGA.csv", header = TRUE)
tcga_cdr <-  readxl::read_excel("D:/Research/Data/Raw/Clinical meta from CDR/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
genus <- read.csv("voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
clinical_clean <- read.csv("D:/Research/Data/Redo/Clinical Traits_AfterClean.csv", header = TRUE,row.names = 1)
# ----> Filter clinical data----
barcode <- tcga$barcode
## Filter the 138sample of clinical data
clinical <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%barcode),c("bcr_patient_barcode","OS","PFI","OS.time","PFI.time")]
## Assign sample name to clinical data
clinical$X <- tcga$X[match(clinical$bcr_patient_barcode,tcga$barcode)]
clinical <- clinical[order(clinical$X),]
clinical_clean <- clinical_clean[order(rownames(clinical_clean)),]
clinical_clean <- cbind(clinical_clean,clinical)

# ----> Rename genus----
## Rename the genus name by genus level
genus <- RenameGenus(genus)

# ----> Combine analysis table----
## Combine clinical data and genus data
clinical <- clinical[order(clinical$X),]
genus <- genus[order(rownames(genus)),]
## Scale the genus abindance
genus_scale <- data.frame(scale(genus))
ana_table <- cbind(clinical,genus_scale)

# ----> Reset the PFI status----
## Reset status by 2000 days to be end point
ana_table$PFI <- ifelse(ana_table$PFI.time >= 2000, 0,ana_table$PFI)
ana_table$PFI.time <- ifelse(ana_table$PFI.time >= 2000, 2000, ana_table$PFI.time)

# ----> Clinical data setting ----
clinical_clean$PFI <- ifelse(clinical_clean$PFI.time >= 2000, 0,clinical_clean$PFI)
clinical_clean$PFI.time <- ifelse(clinical_clean$PFI.time >= 2000, 2000, clinical_clean$PFI.time)


####Survival analysis####
# ----> Univariate----
covariates <- names(genus)
## Set formulas
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(PFI.time, PFI)~', x)))
## Fit model
univ_models <- lapply( univ_formulas, function(x){survival::coxph(x, data = ana_table)})
## Result
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- data.frame(res)
## Filter the p.value < 0.05 genus
summ <- res[which(res$p.value<0.05),]
summary(summ)

# ----> Multivariate----
## For 54 genus find by PFI index
select_genus <- rownames(summ)
select_genus <- paste(select_genus,collapse = "+")
formular <- as.formula(paste('Surv(PFI.time, PFI)~', select_genus))

res.cox <- survival::coxph(formular, data =  ana_table)
multiple_result <- summary(res.cox)
multiple_coeffecient <- data.frame(multiple_result[["coefficients"]])


# ----> Combine uni and multi result----
genus_result <- cbind(summ,multiple_coeffecient)
genus_result$mullow95 <- multiple_result$conf.int[,3]
genus_result$mulhigh95 <- multiple_result$conf.int[,4]
remove_c <- c("wald.test","se.coef.","z")
genus_result <- genus_result[,-which(names(genus_result)%in%remove_c)]
genus_result$mullow95 <- round(genus_result$mullow95,digits = 2)
genus_result$mulhigh95 <- round(genus_result$mulhigh95,digits = 2)
genus_result$exp.coef. <- round(genus_result$exp.coef.,digits=2)
genus_result$exp.coef. <- paste(genus_result$exp.coef.," (",genus_result$mullow95,
                                   "-",genus_result$mulhigh95,")",sep = "")
genus_result$coef <- round(genus_result$coef,digits = 2)

genus_result <- genus_result[,1:6]
names(genus_result) <- c("coef_uni","HR(95CI)_uni","p.value_uni",
                            "coef_mul","HR(95CI)_mul","p.value_mul")

## Filter the genus by p-value of multivariate <= 0.1
genus_result <- genus_result[which(genus_result$p.value_mul <= 0.1),]
write.csv(genus_result,"UniAndMulti_survivalanalysis_PFI.csv")

# ----> Compare OS and PFI genus result----
os_genus <- read.csv("UniAndMulti_trainResult.csv")
os_genus <- os_genus$X[2:12]
list_venn <-list('OS_genus'=os_genus,'PFI_genus'=rownames(genus_result))
ggvenn::ggvenn(list_venn,text_size = 8,set_name_size = 8) 
ggsave("D:/Research/PlotNew/Prognosis/CompareTwoSurvivalIndexGenus.png",dpi = 500)

dupicate_genus <- intersect(rownames(summ),rownames(summ_pfi))
