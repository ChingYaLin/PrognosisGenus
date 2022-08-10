####Function####
RenameGenus <- function(df){
  colname <- colnames(df)
  for (j in 1:length(colname)) {
    temp <- strsplit(colname[j],split = "\\.")
    names(df)[names(df)==colname[j]] <- temp[[1]][length(temp[[1]])]
  }
  return(df)
}
# 1=dead,0=alive

DiffYearData <- function(df,year){
  temp <- df
  days <- year*365
  for (i in 1:length(df$OS.time)) {
    if(df$OS.time[i]>days & df$OS[i]=="Dead"){
      temp$OS[i] <- "Alive"
    }
    if(df$PFI.time[i]>days & df$PFI[i]==1){
      temp$PFI[i] <- 0
    }
  }
  return(temp)
}

PredictToResultOS <- function(fit_model,type,testdata){
  pred_class <- stats::predict(fit_model,
                               new_data = testdata,
                               type = "class")
  
  ##Prediction Probabilities
  pred_proba <- stats::predict(fit_model,
                               new_data = testdata,
                               type = "prob")
  ##Model Evaluation
  data_results <- testdata %>%
    dplyr::select(OS) %>%
    dplyr::bind_cols(pred_class, pred_proba,model = type)
  return(data_results)
}

PredictToResultPFI <- function(fit_model,type,testdata){
  pred_class <- stats::predict(fit_model,
                               new_data = testdata,
                               type = "class")
  
  ##Prediction Probabilities
  pred_proba <- stats::predict(fit_model,
                               new_data = testdata,
                               type = "prob")
  ##Model Evaluation
  data_results <- testdata %>%
    dplyr::select(PFI) %>%
    dplyr::bind_cols(pred_class, pred_proba,model = type)
  return(data_results)
}


####Preprocessing####
# ----> Set path and library----
setwd("D:/Research/Data/Raw")
library("readxl")
# ----> Read data----
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus(CutLow).csv", header = TRUE,row.names = 1)
tcga <-  read.csv("TCGA.csv", header = TRUE)
tcga_cdr <-  readxl::read_excel("Clinical meta from CDR/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
prognosis_genus <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv",header = T)
## Read risk score data
prognosis_genus_table <-  read.csv("D:/Research/Data/Redo/Survival_riskscore.csv",header = T,row.names = 1)
surv_result <- read.csv("D:/Research/Data/Redo/UniAndMulti_survivalanalysis.csv",header = T,row.names = 1)
surv_result <- rownames(surv_result)[which(as.numeric(surv_result$X.4) <= 0.1)]
# ----> Filter the sample----
##Grab the 138 sample from TCGA barcode and give sample name
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%tcga$barcode),]
match_barcode <- match(tcga_select$bcr_patient_barcode,tcga$barcode)
tcga_select$X <- tcga$X[match_barcode]
## Grab the 23 prognosis genus
prognosis_genus <- prognosis_genus$Genus
genus_abundance <- genus[,which(names(genus)%in%prognosis_genus)]
## Rename the genus column and sort the sample
genus_abundance<-RenameGenus(genus_abundance)
genus_matrix <- genus_abundance[order(row.names(genus_abundance)),]

## for 11 significant genus
genus_abundance <- genus_abundance[,which(names(genus_abundance)%in%surv_result)]
genus_matrix <- genus_abundance[order(row.names(genus_abundance)),]
# ---> Combine analysis data----
#combine data to analysis
match_sample <- match(rownames(genus_matrix),tcga_select$X)
genus_matrix$OS.time <- tcga_select$OS.time[match_sample]
genus_matrix$OS <- tcga_select$vital_status[match_sample]
genus_matrix$PFI <- tcga_select$PFI[match_sample]
genus_matrix$PFI.time <- tcga_select$PFI.time[match_sample]
genus_matrix$risk_score_uni <- prognosis_genus_table$risk_score_uni[match(rownames(genus_matrix),rownames(prognosis_genus_table))]
genus_matrix$risk_score_mul <- prognosis_genus_table$risk_score_mul[match(rownames(genus_matrix),rownames(prognosis_genus_table))]
genus_matrix$risk_score_uni10 <- prognosis_genus_table$risk_score_uni10[match(rownames(genus_matrix),rownames(prognosis_genus_table))]
genus_matrix$risk_score_mul10 <- prognosis_genus_table$risk_score_mul10[match(rownames(genus_matrix),rownames(prognosis_genus_table))]


# ----> Seperate diff prognosis data----
## After one years
one_year <- DiffYearData(genus_matrix,1)
## After three years
three_year <- DiffYearData(genus_matrix,3)
## After five Years
five_year <- DiffYearData(genus_matrix,5)

write.csv(five_year,"Five_year_prognosis.csv")
####Model building (OS)####
library(mlbench)
library(tidymodels)
library(grid)
library(pheatmap)
library(xgboost)
library(randomForest)
library(kernlab)
library("survival")
library("survminer")
library("readxl")
library('dplyr')
# ----> Separate training data and testing data----
five_year$OS <- as.factor(five_year$OS)
data_split <- rsample::initial_split(five_year,
                                     prop = 0.8,
                                     strata = OS)
train.data <- data_split %>%
  training()
test.data <- data_split %>%
  testing()

# ----> Fit cox model and riskscore by train data----

## Univariate
train <- train.data[,1:13]
train$OS <- ifelse(train$OS=="Alive",0,1)
covariates <- names(genus_abundance)
## Set formulas
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
## Fit model
univ_models <- lapply( univ_formulas, function(x){survival::coxph(x, data = train)})
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
select_genus <- paste(covariates,collapse = "+")
formular <- as.formula(paste('Surv(OS.time, OS)~', select_genus))

res.cox <- survival::coxph(formular, data =  train)
multiple_result <- summary(res.cox)
multiple_coeffecient <- data.frame(multiple_result[["coefficients"]])
## Combine uni and Multi result
genus_result <- cbind(res,multiple_coeffecient)
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
#write.csv(genus_result,"D:/Research/Data/Redo/UniAndMulti_trainResult.csv")

###Calculate risk score
##Univariate coefficient (Train.data)
train.data$riskscore_uni <- 0
for (i in 1:length(train.data$riskscore_uni)) {
  riskscore <- 0
  for (j in 1:length(rownames(res))) {
    riskscore <- riskscore+as.numeric(res$beta[j])*train.data[[rownames(res)[j]]][i]
  }
  train.data$riskscore_uni[i] <- riskscore
}
## Multivariate coefficient (Train.data)
train.data$riskscore_mul <- 0
for (i in 1:length(train.data$riskscore_mul)) {
  riskscore <- 0
  for (j in 1:length(rownames(multiple_coeffecient))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient$coef[j])*train.data[[rownames(multiple_coeffecient)[j]]][i]
  }
  train.data$riskscore_mul[i] <- riskscore
}

##Univariate coefficient (Test.data)
test.data$riskscore_uni <- 0
for (i in 1:length(test.data$riskscore_uni)) {
  riskscore <- 0
  for (j in 1:length(rownames(res))) {
    riskscore <- riskscore+as.numeric(res$beta[j])*test.data[[rownames(res)[j]]][i]
  }
  test.data$riskscore_uni[i] <- riskscore
}
## Multivariate coefficient (Test.data)
test.data$riskscore_mul <- 0
for (i in 1:length(test.data$riskscore_mul)) {
  riskscore <- 0
  for (j in 1:length(rownames(multiple_coeffecient))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient$coef[j])*test.data[[rownames(multiple_coeffecient)[j]]][i]
  }
  test.data$riskscore_mul[i] <- riskscore
}



# ----> Setting the formular for training model----
select_genus <- names(genus_abundance)
select_genus <- paste(select_genus,collapse = "+")
formular <- as.formula(paste('OS~', select_genus))

# ----> Logistic regression----
log_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(formular, data = train.data)

# ----> Random forest----
rf_spec <- rand_forest()%>%
  set_engine(engine = "randomForest") %>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> XGBoost----
xgb_spec <- boost_tree()%>%
  set_engine(engine = "xgboost") %>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> SVM----
svm_spec <- svm_poly()%>%
  set_engine(engine = "kernlab")%>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> Risk score----
test <- test.data
test$OS <- ifelse(test$OS=="Alive",0,1)
## Reset status by 2000 days to be end point
test$OS <- ifelse(test$OS.time >= 2000, 0,test$OS)
test$OS.time <- ifelse(test$OS.time >= 2000, 2000, test$OS.time)


test <- DiffYearData(test,5)
##UNI
'
quartilelist <- quantile(test$riskscore_uni)[2:4]
threshold <- c("third_q","median_q","first_q")
splots <- list()
for (i in 1:length(threshold)) {
  par(mar=c(10,10,3,3))
  new <- test %>% mutate(risk = ifelse(test$riskscore_uni >=quartilelist[i], "High", "Low"))
  new$risk <- factor(new$risk)
  fit_t <- survfit(Surv(OS.time, OS) ~ risk, data = new)
  splots[[i]] <- ggsurvplot(fit_t, data = new, pval = TRUE,title = threshold[i],
                            conf.int = TRUE,
                            risk.table = TRUE, # Add risk table
                            surv.median.line = "hv",
                            risk.table.col = "risk",
                            risk.table.y.text.col = T,
                            palette = c("#E7B800", "#2E9FDF"),
                            xlab = "Time in days",
                            ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res <- arrange_ggsurvplots(splots, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_OSuni_test.pdf", res, width=6, height=6)
}


## MUL
quartilelist2 <- quantile(test$riskscore_mul)[2:4]
threshold <- c("third_q","median_q","first_q")
splots1 <- list()
for (i in 1:length(threshold)) {
  par(mar=c(10,10,3,3))
  new1 <- test %>% mutate(risk = ifelse(test$riskscore_mul >=quartilelist2[i], "High", "Low"))
  new1$risk <- factor(new1$risk)
  fit_t1 <- survfit(Surv(OS.time, OS) ~ risk, data = new1)
  splots1[[i]] <- ggsurvplot(fit_t1, data = new1, pval = TRUE,title = threshold[i],
                            conf.int = TRUE,
                            risk.table = TRUE, # Add risk table
                            surv.median.line = "hv",
                            risk.table.col = "risk",
                            risk.table.y.text.col = T,
                            palette = c("#E7B800", "#2E9FDF"),
                            xlab = "Time in days",
                            ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res1 <- arrange_ggsurvplots(splots1, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_OSmul_test.pdf", res1, width=6, height=6)
}
'
## Last I decide median to be threshol
riskuni_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(OS ~ riskscore_uni, data = train.data)

riskmul_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(OS ~ riskscore_mul, data = train.data)

#### Predict data ####
log_result <- PredictToResultOS(log_spec,"logistic",test.data)
rf_result <- PredictToResultOS(rf_spec,"random forest",test.data)
xgb_result <- PredictToResultOS(xgb_spec,"XGBoost",test.data)
svm_result <- PredictToResultOS(svm_spec,"svm",test.data)
riskuni_result <- PredictToResultOS(riskuni_spec ,"risk score uni",test.data)
riskmul_result <- PredictToResultOS(riskmul_spec ,"risk score mul",test.data)

out <- bind_rows(log_result, rf_result,xgb_result,svm_result,
                 riskuni_result,riskmul_result)
# ----> Confusion matrix----
confusion_matrix <- out %>%
  group_by(model)%>%
  yardstick::conf_mat(truth = OS,estimate = .pred_class)
  
'
out %>%
  filter(model=="risk score mul")%>%
  yardstick::conf_mat(truth = OS,estimate = .pred_class)%>%
  autoplot(type="heatmap")+
  ggplot2::labs(title = "Confusion matrix (Risk score_multivariate)")+
  theme(plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        text= element_text(size = 20))
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/Confusion_matrix(RS_M).png",dpi = 450)
'
# ----> Custom matrix ----
custom_metrics <- yardstick::metric_set(accuracy, sens, spec, precision, recall, f_meas, kap, mcc)
cus <- out %>%
  group_by(model)%>%
  custom_metrics(truth = OS,estimate = .pred_class)

# ----> Write custom matrix to excel----
# Write the first data set in a new workbook
#library(openxlsx)

dataset_names <- list('XGBoost' = cus[which(cus$model=="XGBoost"),],
                      'logistic' = cus[which(cus$model=="logistic"),],
                      'RF' = cus[which(cus$model=="random forest"),],
                      'svm' = cus[which(cus$model=="svm"),],
                      'RS_U' = cus[which(cus$model=="risk score uni"),],
                      'RS_M' = cus[which(cus$model=="risk score mul"),])
#openxlsx::write.xlsx(dataset_names, file = "D:/Research/Data/Redo/Prognosis_custommatrix_OS.xlsx")
# ----> AUC-ROC ----
auc9 <- out %>%
  group_by(model) %>% yardstick::roc_auc(truth = OS,
                     .pred_Alive)

# ----> 10 times auc----
auc10time <- data.frame(model = auc$model,auc1 = auc$.estimate)
auc10time$auc10 <- auc9$.estimate

rownames(auc10time) <- auc10time$model
auc10time <- auc10time[,-1]
auc10time$mean <- rowMeans(auc10time[,1:10])
auc10time$max <- 0
auc10time$min <- 0
for (i in 1:length(rownames(auc10time))) {
  auc10time$max[i] <- max(auc10time[i,1:10])
  auc10time$min[i] <- min(auc10time[i,1:10])
}
auc10time_result <- auc10time[,which(names(auc10time)%in%c('mean','max','min'))]
names(auc10time)[11] <- "auc_mean"
names(auc10time)[12] <- 'auc_max'
names(auc10time)[13] <- 'auc_min'
write.csv(auc10time,'C:/Users/user/Desktop/10timesAUC_multisignificant.csv')
## plot roc
out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = OS, .pred_Alive) %>% # get values to plot an ROC curve
  ggplot2::autoplot()+
  geom_abline(slope = 1, intercept = 0, size = 0.4)+
  ggplot2::labs(title = "ROC Curve (After three year)",
                subtitle = "Testing data")+
  ggplot2::annotate(geom="text", x=0.75, y=0.2, 
                    label=paste("AUC.log =",round(auc$.estimate[which(auc$model=="logistic")],digits = 4),"\n",
                                "AUC.rf =",round(auc$.estimate[which(auc$model=="random forest")],digits = 4),"\n",
                                "AUC.svm =",round(auc$.estimate[which(auc$model=="svm")],digits = 4),"\n",
                                "AUC.xgb =",round(auc$.estimate[which(auc$model=="XGBoost")],digits = 4),"\n",
                                "AUC.risk_uni =",round(auc$.estimate[which(auc$model=="risk score uni")],digits = 4),"\n",
                                "AUC.risk_mul =",round(auc$.estimate[which(auc$model=="risk score mul")],digits = 4)),
                    size=4)+
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/ROC_threeYear_OS_test.png",dpi = 450)

####Model building (PFI)####
overlap_genus <- read.csv('D:/Research/Data/Redo/OSandPFIoverlapGenus.csv')
overlap_genus <- overlap_genus$X
# ----> Separate training data and testing data----
five_year$PFI <- as.factor(five_year$PFI)
data_split <- rsample::initial_split(five_year,
                                     prop = 0.8,
                                     strata = PFI)
train.data <- data_split %>%
  training()
test.data <- data_split %>%
  testing()

# ----> Fit cox model and riskscore by train data----
stay_column <- c(overlap_genus, 'PFI','PFI.time')
## Univariate
train <- train.data[,which(names(train.data)%in%stay_column)]
train$PFI <- ifelse(train$PFI==0,0,1)

covariates <- overlap_genus
## Set formulas
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(PFI.time, PFI)~', x)))
## Fit model
univ_models <- lapply( univ_formulas, function(x){survival::coxph(x, data = train)})
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
select_genus <- paste(covariates,collapse = "+")
formular <- as.formula(paste('Surv(PFI.time, PFI)~', select_genus))

res.cox <- survival::coxph(formular, data =  train)
multiple_result <- summary(res.cox)
multiple_coeffecient <- data.frame(multiple_result[["coefficients"]])
## Combine uni and Multi result
genus_result <- cbind(res,multiple_coeffecient)
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
write.csv(genus_result,"D:/Research/Data/Redo/UniAndMulti_trainResultForPFI.csv")

###Calculate risk score
##Univariate coefficient (Train.data)
train.data$riskscore_uni <- 0
for (i in 1:length(train.data$riskscore_uni)) {
  riskscore <- 0
  for (j in 1:length(rownames(res))) {
    riskscore <- riskscore+as.numeric(res$beta[j])*train.data[[rownames(res)[j]]][i]
  }
  train.data$riskscore_uni[i] <- riskscore
}
## Multivariate coefficient (Train.data)
train.data$riskscore_mul <- 0
for (i in 1:length(train.data$riskscore_mul)) {
  riskscore <- 0
  for (j in 1:length(rownames(multiple_coeffecient))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient$coef[j])*train.data[[rownames(multiple_coeffecient)[j]]][i]
  }
  train.data$riskscore_mul[i] <- riskscore
}

##Univariate coefficient (Test.data)
test.data$riskscore_uni <- 0
for (i in 1:length(test.data$riskscore_uni)) {
  riskscore <- 0
  for (j in 1:length(rownames(res))) {
    riskscore <- riskscore+as.numeric(res$beta[j])*test.data[[rownames(res)[j]]][i]
  }
  test.data$riskscore_uni[i] <- riskscore
}
## Multivariate coefficient (Test.data)
test.data$riskscore_mul <- 0
for (i in 1:length(test.data$riskscore_mul)) {
  riskscore <- 0
  for (j in 1:length(rownames(multiple_coeffecient))) {
    riskscore <- riskscore+as.numeric(multiple_coeffecient$coef[j])*test.data[[rownames(multiple_coeffecient)[j]]][i]
  }
  test.data$riskscore_mul[i] <- riskscore
}



# ----> Setting the formular for training model----
select_genus <- overlap_genus
select_genus <- paste(select_genus,collapse = "+")
formular <- as.formula(paste('PFI~', select_genus))

# ----> Logistic regression----
log_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(formular, data = train.data)

# ----> Random forest----
rf_spec <- rand_forest()%>%
  set_engine(engine = "randomForest") %>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> XGBoost----
xgb_spec <- boost_tree()%>%
  set_engine(engine = "xgboost") %>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> SVM----
svm_spec <- svm_poly()%>%
  set_engine(engine = "kernlab")%>%
  set_mode("classification")%>%
  fit(formular, data = train.data)

# ----> Risk score----
test <- test.data
## Reset status by 2000 days to be end point
test$PFI <- ifelse(test$PFI.time >= 2000, 0,1)
test$PFI.time <- ifelse(test$PFI.time >= 2000, 2000, test$PFI.time)

##UNI
quartilelist <- quantile(test$riskscore_uni)[2:4]
threshold <- c("third_q","median_q","first_q")
splots <- list()
for (i in 1:length(threshold)) {
  par(mar=c(10,10,3,3))
  new <- test %>% mutate(risk = ifelse(test$riskscore_uni >=quartilelist[i], "High", "Low"))
  new$risk <- factor(new$risk)
  fit_t <- survfit(Surv(PFI.time, PFI) ~ risk, data = new)
  splots[[i]] <- ggsurvplot(fit_t, data = new, pval = TRUE,title = threshold[i],
                            conf.int = TRUE,
                            risk.table = TRUE, # Add risk table
                            surv.median.line = "hv",
                            risk.table.col = "risk",
                            risk.table.y.text.col = T,
                            palette = c("#E7B800", "#2E9FDF"),
                            xlab = "Time in days",
                            ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res <- arrange_ggsurvplots(splots, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_PFIuni_test.pdf", res, width=6, height=6)
}

## MUL
quartilelist2 <- quantile(test$riskscore_mul)[2:4]
threshold <- c("third_q","median_q","first_q")
splots1 <- list()
for (i in 1:length(threshold)) {
  par(mar=c(10,10,3,3))
  new1 <- test %>% mutate(risk = ifelse(test$riskscore_mul >=quartilelist2[i], "High", "Low"))
  new1$risk <- factor(new1$risk)
  fit_t1 <- survfit(Surv(PFI.time, PFI) ~ risk, data = new1)
  splots1[[i]] <- ggsurvplot(fit_t1, data = new1, pval = TRUE,title = threshold[i],
                             conf.int = TRUE,
                             risk.table = TRUE, # Add risk table
                             surv.median.line = "hv",
                             risk.table.col = "risk",
                             risk.table.y.text.col = T,
                             palette = c("#E7B800", "#2E9FDF"),
                             xlab = "Time in days",
                             ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res1 <- arrange_ggsurvplots(splots1, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_PFImul_test.pdf", res1, width=6, height=6)
}

## Last I decide median to be threshol
riskuni_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(PFI ~ riskscore_uni, data = train.data)

riskmul_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(PFI ~ riskscore_mul, data = train.data)

#### Predict data ####
log_result <- PredictToResultPFI(log_spec,"logistic",test.data)
rf_result <- PredictToResultPFI(rf_spec,"random forest",test.data)
xgb_result <- PredictToResultPFI(xgb_spec,"XGBoost",test.data)
svm_result <- PredictToResultPFI(svm_spec,"svm",test.data)
riskuni_result <- PredictToResultPFI(riskuni_spec ,"risk score uni",test.data)
riskmul_result <- PredictToResultPFI(riskmul_spec ,"risk score mul",test.data)

out <- bind_rows(log_result, rf_result,xgb_result,svm_result,
                 riskuni_result,riskmul_result)
# ----> Confusion matrix----
confusion_matrix <- out %>%
  group_by(model)%>%
  yardstick::conf_mat(truth = PFI,estimate = .pred_class)


out %>%
  filter(model=="svm")%>%
  yardstick::conf_mat(truth = PFI,estimate = .pred_class)%>%
  autoplot(type="heatmap")+
  ggplot2::labs(title = "Confusion matrix (svm)")+
  theme(plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        text= element_text(size = 20))
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/Confusion_matrix(svm).png",dpi = 450)

# ----> Custom matrix ----
custom_metrics <- yardstick::metric_set(accuracy, sens, spec, precision, recall, f_meas, kap, mcc)
cus <- out %>%
  group_by(model)%>%
  custom_metrics(truth = PFI,estimate = .pred_class)

# ----> Write custom matrix to excel----
# Write the first data set in a new workbook
library(openxlsx)

dataset_names <- list('XGBoost' = cus[which(cus$model=="XGBoost"),],
                      'logistic' = cus[which(cus$model=="logistic"),],
                      'RF' = cus[which(cus$model=="random forest"),],
                      'svm' = cus[which(cus$model=="svm"),],
                      'RS_U' = cus[which(cus$model=="risk score uni"),],
                      'RS_M' = cus[which(cus$model=="risk score mul"),])
openxlsx::write.xlsx(dataset_names, file = "D:/Research/Data/Redo/Prognosis_custommatrix_PFI.xlsx")
# ----> AUC-ROC ----
auc <- out %>%
  group_by(model) %>% yardstick::roc_auc(truth = PFI,
                                         .pred_0)

## plot roc
out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = PFI, .pred_0) %>% # get values to plot an ROC curve
  ggplot2::autoplot()+
  geom_abline(slope = 1, intercept = 0, size = 0.4)+
  ggplot2::labs(title = "ROC Curve (After five year)",
                subtitle = "Tresting data")+
  ggplot2::annotate(geom="text", x=0.75, y=0.2, 
                    label=paste("AUC.log =",round(auc$.estimate[which(auc$model=="logistic")],digits = 4),"\n",
                                "AUC.rf =",round(auc$.estimate[which(auc$model=="random forest")],digits = 4),"\n",
                                "AUC.svm =",round(auc$.estimate[which(auc$model=="svm")],digits = 4),"\n",
                                "AUC.xgb =",round(auc$.estimate[which(auc$model=="XGBoost")],digits = 4),"\n",
                                "AUC.risk_uni =",round(auc$.estimate[which(auc$model=="risk score uni")],digits = 4),"\n",
                                "AUC.risk_mul =",round(auc$.estimate[which(auc$model=="risk score mul")],digits = 4)),
                    size=4)+
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/ROC_threeYear_PFI_test.png",dpi = 450)


