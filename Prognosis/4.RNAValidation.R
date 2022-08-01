####Package####
# ----> Survival analysis----
library(survival)
library(survminer)
library(readxl)
library(dplyr)
# ----> Model building----
library(mlbench)
library(tidymodels)
library(grid)
library(pheatmap)
library(xgboost)
library(randomForest)
library(kernlab)


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

DiffYearData <- function(df,year){
  temp <- df
  days <- year*365
  for (i in 1:length(df$OS.time)) {
    if(df$OS.time[i]>days & df$OS[i]==1){
      temp$OS[i] <- 0
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

####Preprocess####
setwd("D:/Research/Data/Redo")

# ----> Read data----
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus_rna.csv", header = TRUE,row.names = 1)
tcga_rna <-  read.csv("D:/Research/Data/Raw/TCGA_Barcode_rna.csv", header = TRUE)
tcga_cdr <-  readxl::read_excel("D:/Research/Data/Raw/Clinical meta from CDR/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
prognosis_genus <- read.csv('D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv',
                            header = T)
genus10 <- read.csv('D:/Research/Data/Redo/Univariate_duplicate10survival.csv',header = T)

# ----> Filter the sample----
##Grab the 138 sample from TCGA barcode and give sample name
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%tcga_rna$submitter_id),]
tcga_select$X <- tcga_rna$X[match(tcga_select$bcr_patient_barcode,tcga_rna$submitter_id)]
## Grab the 23 prognosis genus
genus_abundance <- genus[,which(names(genus)%in%prognosis_genus$Genus)]

# ----> Rename genus----
## Rename the genus column and sort the sample
genus_abundance<-RenameGenus(genus_abundance)
genus_matrix <- genus_abundance[order(row.names(genus_abundance)),]


# ----> Filter clinical data----
## Filter the 138sample of clinical data
clinical <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%tcga_rna$submitter_id),c("bcr_patient_barcode","OS","PFI","OS.time","PFI.time")]
## Assign sample name to clinical data
clinical$X <- tcga_rna$X[match(clinical$bcr_patient_barcode,tcga_rna$submitter_id)]


# ----> Combine analysis table----
## Combine clinical data and genus data
clinical <- clinical[order(clinical$X),]
genus_matrix <- genus_matrix[order(rownames(genus_matrix)),]
## Scale the genus abindance
genus_scale <- data.frame(scale(genus_matrix))
ana_table <- cbind(clinical,genus_scale)

# ----> Reset the OS status----
## Reset status by 2000 days to be end point
ana_table$OS <- ifelse(ana_table$OS.time >= 2000, 0,ana_table$OS)
ana_table$PFI <- ifelse(ana_table$PFI.time >= 2000, 0,ana_table$PFI)
ana_table$OS.time <- ifelse(ana_table$OS.time >= 2000, 2000, ana_table$OS.time)
ana_table$PFI.time <- ifelse(ana_table$PFI.time >= 2000, 2000, ana_table$PFI.time)
five_year <- DiffYearData(ana_table,5)
####Model building (OS)####
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
train <- train.data
train$OS <- ifelse(train$OS==0,0,1)
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
write.csv(genus_result,"D:/Research/Data/Redo/UniAndMulti_trainResult_RNA.csv")

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
riskuni_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(OS ~ riskscore_uni, data = train.data)

riskmul_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(OS ~ riskscore_mul, data = train.data)

#### Predict data ####
log_result <- PredictToResultOS(log_spec,"logistic",train.data)
rf_result <- PredictToResultOS(rf_spec,"random forest",train.data)
xgb_result <- PredictToResultOS(xgb_spec,"XGBoost",train.data)
svm_result <- PredictToResultOS(svm_spec,"svm",train.data)
riskuni_result <- PredictToResultOS(riskuni_spec ,"risk score uni",train.data)
riskmul_result <- PredictToResultOS(riskmul_spec ,"risk score mul",train.data)

out <- bind_rows(log_result, rf_result,xgb_result,svm_result,
                 riskuni_result,riskmul_result)
# ----> AUC-ROC ----
auc <- out %>%
  group_by(model) %>% yardstick::roc_auc(truth = OS,
                                         .pred_0)

## plot roc
out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = OS, .pred_0) %>% # get values to plot an ROC curve
  ggplot2::autoplot()+
  geom_abline(slope = 1, intercept = 0, size = 0.4)+
  ggplot2::labs(title = "ROC Curve (After five year)",
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
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/ROC_fiveYear_OS_train.png",dpi = 450)

# ----> Risk score KM plot----
test <- test.data
test$OS <- ifelse(test$OS==0,0,1)
## Reset status by 2000 days to be end point
test$OS <- ifelse(test$OS.time >= 2000, 0,test$OS)
test$OS.time <- ifelse(test$OS.time >= 2000, 2000, test$OS.time)


test <- DiffYearData(test,5)
##UNI
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
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_OSuni_test_RNA.pdf", res, width=6, height=6)
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
  ggsave("D:/Research/PlotNew/Prognosis/SurvivalKMplot_OSmul_test_RNA.pdf", res1, width=6, height=6)
}



####Model building (PFI)####
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
stay_column <- c(genus10$X, 'PFI','PFI.time')
## Univariate
train <- train.data[,which(names(train.data)%in%stay_column)]
train$PFI <- ifelse(train$PFI==0,0,1)

covariates <- genus10$X
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
write.csv(genus_result,"D:/Research/Data/Redo/UniAndMulti_trainResultForPFI_RNA.csv")

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
select_genus <- genus10$X
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
riskuni_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(PFI ~ riskscore_uni, data = train.data)

riskmul_spec <- logistic_reg() %>%
  set_engine(engine = "glm") %>%
  set_mode("classification") %>%
  fit(PFI ~ riskscore_mul, data = train.data)

#### Predict data ####
log_result <- PredictToResultPFI(log_spec,"logistic",train.data)
rf_result <- PredictToResultPFI(rf_spec,"random forest",train.data)
xgb_result <- PredictToResultPFI(xgb_spec,"XGBoost",train.data)
svm_result <- PredictToResultPFI(svm_spec,"svm",train.data)
riskuni_result <- PredictToResultPFI(riskuni_spec ,"risk score uni",train.data)
riskmul_result <- PredictToResultPFI(riskmul_spec ,"risk score mul",train.data)

out <- bind_rows(log_result, rf_result,xgb_result,svm_result,
                 riskuni_result,riskmul_result)

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
                subtitle = "Training data")+
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
ggplot2::ggsave("D:/Research/PlotNew/Prognosis/ROC_fiveYear_PFI_train.png",dpi = 450)

# ---->Risk score KM plot----
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
