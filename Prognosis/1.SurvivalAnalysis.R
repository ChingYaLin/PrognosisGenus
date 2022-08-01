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
##Risk score
riskscore <- read.csv('Survival_riskscore.csv')

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

# ----> Reset the OS status----
## Reset status by 2000 days to be end point
ana_table$OS <- ifelse(ana_table$OS.time >= 2000, 0,ana_table$OS)
ana_table$PFI <- ifelse(ana_table$PFI.time >= 2000, 0,ana_table$PFI)
ana_table$OS.time <- ifelse(ana_table$OS.time >= 2000, 2000, ana_table$OS.time)
ana_table$PFI.time <- ifelse(ana_table$PFI.time >= 2000, 2000, ana_table$PFI.time)

# ----> Clinical data setting ----
clinical_clean$OS <- ifelse(clinical_clean$OS.time >= 2000, 0,clinical_clean$OS)
clinical_clean$PFI <- ifelse(clinical_clean$PFI.time >= 2000, 0,clinical_clean$PFI)
clinical_clean$OS.time <- ifelse(clinical_clean$OS.time >= 2000, 2000, clinical_clean$OS.time)
clinical_clean$PFI.time <- ifelse(clinical_clean$PFI.time >= 2000, 2000, clinical_clean$PFI.time)
## risk score setting
clinical_clean$riskscore_uni <- riskscore$risk_score_uni
clinical_clean$riskscore_mul <- riskscore$risk_score_mul

####Survival analysis####
# ----> Clinical trait analysis----
clinic_covariates <- names(clinical_clean)[1:9]
clinic_covariates <- clinic_covariates[-4]
clinic_covariates <- c(clinic_covariates,'riskscore_uni','riskscore_mul')
## Set formulas
univ_formulas <- sapply(clinic_covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
## Fit model
univ_models <- lapply( univ_formulas, function(x){survival::coxph(x, data = clinical_clean)})
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
## Multivariate
select_feature <- paste(clinic_covariates,collapse = "+")
formular <- as.formula(paste('Surv(OS.time, OS)~', select_feature))

res.cox <- survival::coxph(formular, data =  clinical_clean)
multiple_result <- summary(res.cox)
multiple_coeffecient <- data.frame(multiple_result[["coefficients"]])

clinical_result <- cbind(res,multiple_coeffecient)
clinical_result$mullow95 <- multiple_result$conf.int[,3]
clinical_result$mulhigh95 <- multiple_result$conf.int[,4]
remove_c <- c("wald.test","se.coef.","z")
clinical_result <- clinical_result[,-which(names(clinical_result)%in%remove_c)]
clinical_result$mullow95 <- round(clinical_result$mullow95,digits = 2)
clinical_result$mulhigh95 <- round(clinical_result$mulhigh95,digits = 2)
clinical_result$exp.coef. <- round(clinical_result$exp.coef.,digits=2)
clinical_result$exp.coef. <- paste(clinical_result$exp.coef.,"(",clinical_result$mullow95,
                                   "-",clinical_result$mulhigh95,")",sep = "")
clinical_result$coef <- round(clinical_result$coef,digits = 2)
clinical_result$Pr...z.. <- round(clinical_result$Pr...z..,digits = 2)
clinical_result <- clinical_result[,1:6]
names(clinical_result) <- c("coef_uni","HR(95CI)_uni","p.value_uni",
                            "coef_mul","HR(95CI)_mul","p.value_mul")
write.csv(clinical_result,"clinical_survivalAnalysis_addRiskScore.csv")

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
summ_pfi <- res[which(res$p.value<0.05),]
summ <- res[which(res$p.value<0.05),]
summary(summ)
## Compare OS and PFI genus result
list_venn <-list('OS_genus'=rownames(summ),'PFI_genus'=rownames(summ_pfi))
ggvenn::ggvenn(list_venn,text_size = 20,set_name_size = 20) 
ggsave("D:/Research/PlotNew/31.CompareTwoSurvivalIndexGenus.png",dpi = 500)

dupicate_genus <- intersect(rownames(summ),rownames(summ_pfi))
dupicate <- cbind(summ[which(rownames(summ)%in%dupicate_genus),],summ_pfi[which(rownames(summ_pfi)%in%dupicate_genus),])
names(dupicate) <- c("OS.beta","OS.HR","OS.Wald.test","OS.Pvalue",
                  "PFI.beta","PFI.HR","PFI.Wald.test","PFI.Pvalue")

write.csv(dupicate,"OSandPFIoverlapGenus.csv")
# ----> Multivariate----
## For 23 genus find by OS index
select_genus <- rownames(summ)
select_genus <- paste(select_genus,collapse = "+")
formular <- as.formula(paste('Surv(OS.time, OS)~', select_genus))

res.cox <- survival::coxph(formular, data =  ana_table)
multiple_result <- summary(res.cox)
multiple_coeffecient <- data.frame(multiple_result[["coefficients"]])

## For duplicate of OS and PFI index genus
dup_genus <- rownames(dupicate)
dup_genus <- paste(dup_genus,collapse = "+")
formular_dup <- as.formula(paste('Surv(OS.time, OS)~', dup_genus))

res.cox.dup <- survival::coxph(formular_dup, data =  ana_table)
multiple_result_dup <- summary(res.cox.dup)
multiple_coeffecient_dup <- data.frame(multiple_result_dup[["coefficients"]])

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
write.csv(genus_result,"UniAndMulti_survivalanalysis.csv")

####Risk score####
prognosis_genus_table <- ana_table[,which(names(ana_table)%in%c(rownames(summ),"OS","OS.time"))]

# ----> calculate risk score----
## Univariate coefficient
prognosis_genus_table$risk_score_uni <- 0
for (i in 1:length(prognosis_genus_table$risk_score_uni)) {
  riskscore <- 0
  for (j in 1:length(rownames(summ))) {
    riskscore <- riskscore+as.numeric(summ$beta[j])*prognosis_genus_table[[rownames(summ)[j]]][i]
  }
  prognosis_genus_table$risk_score_uni[i] <- riskscore
}
## Multivariate coefficient
prognosis_genus_table$risk_score_mul <- 0
for (i in 1:length(prognosis_genus_table$risk_score_mul)) {
  riskscore <- 0
  for (j in 1:length(rownames(multiple_coeffecient))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient$coef[j])*prognosis_genus_table[[rownames(multiple_coeffecient)[j]]][i]
  }
  prognosis_genus_table$risk_score_mul[i] <- riskscore
}
##Duplicate 10 genus univariate
prognosis_genus_table$risk_score_uni10 <- 0
for (i in 1:length(prognosis_genus_table$risk_score_uni10)) {
  riskscore <- 0
  for (j in 1:length(rownames(dupicate))) {
    riskscore <- riskscore+as.numeric(dupicate$OS.beta[j])*prognosis_genus_table[[rownames(dupicate)[j]]][i]
  }
  prognosis_genus_table$risk_score_uni10[i] <- riskscore
}
##Duplicate 10 genus mulivariate
prognosis_genus_table$risk_score_mul10 <- 0
for (i in 1:length(prognosis_genus_table$risk_score_mul10)) {
  riskscore <- 0
  for (j in 1:length(rownames(dupicate))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient_dup$coef[j])*prognosis_genus_table[[rownames(dupicate)[j]]][i]
  }
  prognosis_genus_table$risk_score_mul10[i] <- riskscore
}


# ----> Save the risk score file----
## Save the file to doing classification in model training
write.csv(prognosis_genus_table,"Survival_riskscore.csv")
write.csv(summ,"Univarite_survivalResult.csv")
write.csv(multiple_coeffecient,"Multivarite_survivalResult.csv")
write.csv(dupicate,"Univariate_duplicate10survival.csv")
####Visualization####
# ----> Survival probability----
## Survival probability in all prognosis genus
test <- survival::survfit(res.cox)
survminer::ggsurvplot(test,data=ana_table, palette= '#2E9FDF',
           ggtheme = theme_minimal())

## Individual genus survival probability plot
survminer::ggsurvplot(survfit(univ_models$g__Simplexvirus),data=ana_table,
           ggtheme = theme_minimal(),title = "g__Simplexvirus")


# ----> High and low abundance----
## Add High and low group of significant multivariate result
new <- ana_table %>% dplyr::mutate(Mamastro_exp = ifelse(g__Mamastrovirus >=median(ana_table$g__Mamastrovirus), "High", "Low"),
                            Nitro_exp = ifelse(g__Candidatus_Nitrosopelagicus >=median(ana_table$g__Candidatus_Nitrosopelagicus), "High", "Low"),
                            Brachy_exp = ifelse(g__Brachyspira >=median(ana_table$g__Brachyspira), "High", "Low"),
                            Terrim_exp = ifelse(g__Terrimonas >=median(ana_table$g__Terrimonas), "High", "Low"),
                            Thermo_exp = ifelse(g__Thermocrispum >=median(ana_table$g__Thermocrispum), "High", "Low"),
                            Ornith_exp = ifelse(g__Ornithinibacillus >=median(ana_table$g__Ornithinibacillus), "High", "Low"))

significant_genus <- c('Nitro_exp','Brachy_exp','Terrim_exp',
                       'Mamastro_exp','Thermo_exp','Ornith_exp')
## Transfer the charactor into factor
for (i in significant_genus) {
  new[[i]] <- factor(new[[i]])
}


fit1 <- survfit(Surv(OS.time, DSS) ~ Mamastro_exp, data = new)
survminer::ggsurvplot(fit1, data = new, pval = TRUE)

# ----> 25 percent plot----
percent_25 <- new
percent_25 <- Seperate3express(percent_25,"g__Mamastrovirus",'Mamastro_exp')
percent_25 <- Seperate3express(percent_25,'g__Candidatus_Nitrosopelagicus','Nitro_exp')
percent_25 <- Seperate3express(percent_25,'g__Brachyspira','Brachy_exp')
percent_25 <- Seperate3express(percent_25,'g__Terrimonas','Terrim_exp')
percent_25 <- Seperate3express(percent_25,'g__Thermocrispum','Thermo_exp')
percent_25 <- Seperate3express(percent_25,'g__Ornithinibacillus','Ornith_exp')


fit1 <- survfit(Surv(OS.time, DSS) ~ Mamastro_exp, data = percent_25)
survminer::ggsurvplot(fit1, data = percent_25, pval = TRUE)


# ----> Heatmap for specific genus----
## setting filter
genus_prognosis <- rownames(summ)
filter <- multiple_result$coefficients[,5] <= 0.05
## Plot the genus that Pr(>|z|) <= 0.05
genus_abundance_whole <- genus[,filter]

test <- 200:-200
test <- test*0.01
draw <- data.frame(t(genus_abundance_whole))
pheatmap::pheatmap(draw,scale = "row",
                   color = hcl.colors(400, "RdBu"),
                   breaks = test,
                   show_colnames = F,
                   annotation_col = subset(ana_table, select = c("OS")))

## Plot the genus that select by survival analysis
genus_abundance_prognosis <- ana_table[,which(names(ana_table)%in%genus_prognosis)]

test <- 200:-200
test <- test*0.01
draw <- data.frame(t(genus_abundance_prognosis))
#>>> The os label
pdf("D:/Research/PlotNew/33.Heatmap of OS (23genus).pdf")
par(mar=c(8,8,8,8))
pheatmap::pheatmap(draw,scale = "row",
                   color = hcl.colors(400, "RdBu"),
                   breaks = test,
                   show_colnames = F,
                   annotation_col = subset(ana_table, select = c("OS")))
dev.off()

pdf("D:/Research/PlotNew/34.Heatmap of PFI (23genus).pdf")
par(mar=c(8,8,8,8))
pheatmap::pheatmap(draw,scale = "row",
                   color = hcl.colors(400, "RdBu"),
                   breaks = test,
                   show_colnames = F,
                   annotation_col = subset(ana_table, select = c("PFI")))
dev.off()
# ----> KM-plot----
splots <- list()
for (i in 1:length(rownames(summ))) {
  par(mar=c(10,10,3,3))
  new <- ana_table %>% mutate(exp = ifelse(ana_table[[rownames(summ)[i]]] >=median(ana_table[[rownames(summ)[i]]]), "High", "Low"))
  new$exp <- factor(new$exp)
  fit_t <- survfit(Surv(OS.time, OS) ~ exp, data = new)
  splots[[i]] <- ggsurvplot(fit_t, data = new, pval = TRUE,title = rownames(summ)[i],
             conf.int = TRUE,
             risk.table = TRUE, # Add risk table
             surv.median.line = "hv",
             risk.table.col = "exp",
             risk.table.y.text.col = T,
             palette = c("#E7B800", "#2E9FDF"),
             xlab = "Time in days",
             ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res <- arrange_ggsurvplots(splots, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/PlotNew/32.SurvivalKMplot.pdf", res, width=6, height=6)
}

####Write result####
summary(fit_t)$table
####store the results of multiple analysis
results <- multiple_result[["coefficients"]]
write.csv(results,file = "D:/Research/Data/Redo/Cox_surv_multipleResult.csv")

####Referance####
## Cox regression model
#http://www.sthda.com/english/wiki/cox-proportional-hazards-model

