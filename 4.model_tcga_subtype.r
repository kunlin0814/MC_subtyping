library(randomForest)
library(xgboost)
library(DMwR)
library(grid)
library(Boruta)
library(vars)
library(dplyr)
library(varSelRF)
library(caret)
library(e1071)
library(ranger)
library(smotefamily)
library(mlr3measures)
library(kernlab)
source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')

comparison <- c("LumA","LumB")

base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
#"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
minor_group <- each_pheno_num[each_pheno_num$n ==min(each_pheno_num$n),]$SUBTYPE


# We will use both train and test data

comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
load(paste(prev_results_base,'feature_selection_tcga_subtype.rdata',sep="/"))


train_control <- trainControl(method = "cv" , number = 5,
                              classProbs = TRUE)
# data <- data_model
# nround <- 5
Model_50_RF <- function(data,nround) {
  PC0_rf <- 0 ; PC1_rf <- 0;RC0_rf <- 0;RC1_rf <- 0;F10_rf <- 0;F11_rf <- 0;MF1_rf<-0 ; PR_AUC_rf <- 0;AUC_rf <- 0;MPC_rf <-0 ; MRC_rf <-0 ; MCC_rf <-0
  nround = nround
  for (i in 1:nround) {
    set.seed(i)
    print(i)
    train_indices <- createDataPartition(data_model$SUBTYPE, p = 0.8, list = FALSE, times = 1)
    train_data <- data_model[train_indices, ]
    test_data <- data_model[-train_indices, ]
    
    # We use ADASYN (ADASYN is an method to oversampling the minority class)
    genData_ADAS <- ADAS(X = train_data[, -1], target = train_data$SUBTYPE, K = 5)
    train_df_adas <- genData_ADAS[["data"]]
    class_df <- train_df_adas[ncol(train_df_adas)]
    train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
    train_data_balanced <- train_df_adas
    names(train_data_balanced)[1] <- "SUBTYPE"
    
    train_data_balanced <- train_data_balanced %>%
      mutate(across(.cols = 1, .fns = factor)) %>%
      mutate(across(.cols = 2:length(.), .fns = as.numeric))
    
    train_control <- trainControl(
      method = "cv", number = 5, classProbs = TRUE)
    # Random Forest
    rf_model <- caret::train(SUBTYPE ~ ., data = train_data_balanced, method = "ranger",
                             verbose = FALSE, trControl = train_control, metric = 'Accuracy',
                             tuneLength = 10)
    
    # Prediction
    res_rf <- Res_CMT(prob = predict(rf_model, test_data, type = "prob"),
                      pred = predict(rf_model, test_data),
                      test_data = test_data, model =rf_model )
    
    # Store results in vectors
    PC0_rf[i] <- res_rf[[1]]
    PC1_rf[i] <- res_rf[[2]]
    RC0_rf[i] <- res_rf[[3]]
    RC1_rf[i] <- res_rf[[4]]
    F10_rf[i] <- res_rf[[5]]
    F11_rf[i] <- res_rf[[6]]
    MPC_rf[i] <- res_rf[[7]]
    MRC_rf[i] <- res_rf[[8]]
    MF1_rf[i] <- res_rf[[9]]
    PR_AUC_rf[i] <- res_rf[[10]]
    AUC_rf[i] <- res_rf[[11]]
    MCC_rf[i] <- res_rf[[12]]
  }
  
  df_res_rf <- data.frame(PC0 = PC0_rf, PC1 = PC1_rf, RC0 = RC0_rf, RC1 = RC1_rf,
                          F10 = F10_rf, F11 = F11_rf, MPC = MPC_rf, MRC = MRC_rf,
                          MF1 = MF1_rf, ROC_AUC = AUC_rf, PR_AUC = PR_AUC_rf, MCC = MCC_rf)
  
  df_res_rf <- as.data.frame(t(df_res_rf))
  df_res_rf[is.na(df_res_rf)] <- 0
  df_res_rf$mean <- apply(df_res_rf, 1, mean)
  df_res_rf$sd <- apply(df_res_rf, 1, sd)
  df_res_rf$CI95 <- 2 * sqrt((df_res_rf$mean * (1 - df_res_rf$mean)) / nround)
  df_res_rf <- round(df_res_rf, 3)
  
  return(df_res_rf)
}


Model_50_SVM <- function(data,nround) {
  PC0_svm <- 0 ; PC1_svm <- 0;RC0_svm <- 0;RC1_svm <- 0;F10_svm <- 0;F11_svm <- 0;MF1_svm<-0 ; PR_AUC_svm <- 0;AUC_svm <- 0;MPC_svm <-0 ; MRC_svm <-0 ; MCC_svm <-0
  
  nround = nround
  for (i in 1:nround) {
    print(i)
    set.seed(i)
    train_indices <- createDataPartition(data$SUBTYPE, p = 0.8, list = FALSE, times = 1)
    train_data <- data[train_indices, ]
    test_data <- data[-train_indices, ]
    
    # We use ADASYN
    genData_ADAS <- ADAS(X = train_data[, -1], target = train_data$SUBTYPE, K = 5)
    train_df_adas <- genData_ADAS[["data"]]
    class_df <- train_df_adas[ncol(train_df_adas)]
    train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
    train_data_balanced <- train_df_adas
    names(train_data_balanced)[1] <- "SUBTYPE"
    
    train_data_balanced <- train_data_balanced %>%
      mutate(across(.cols = 1, .fns = factor)) %>%
      mutate(across(.cols = 2:length(.), .fns = as.numeric))
    
    # SVM
    svm_fit <- caret::train(SUBTYPE ~ ., data = train_data_balanced, method = "svmRadial",
                            trControl = train_control, preProcess = c("center", "scale"),
                            tuneLength = 10)
    
    # Prediction
    predictions <- predict(svm_fit, newdata = test_data)
    probabilities <- predict(svm_fit, newdata = test_data, type = "prob")
    res_svm <- Res_CMT(prob = probabilities, pred = predictions, 
                       test_data = test_data,model= svm_fit)
    
    # Store results in vectors
    PC0_svm[i] <- res_svm[[1]]
    PC1_svm[i] <- res_svm[[2]]
    RC0_svm[i] <- res_svm[[3]]
    RC1_svm[i] <- res_svm[[4]]
    F10_svm[i] <- res_svm[[5]]
    F11_svm[i] <- res_svm[[6]]
    MPC_svm[i] <- res_svm[[7]]
    MRC_svm[i] <- res_svm[[8]]
    MF1_svm[i] <- res_svm[[9]]
    PR_AUC_svm[i] <- res_svm[[10]]
    AUC_svm[i] <- res_svm[[11]]
    MCC_svm[i] <- res_svm[[12]]
  }
  
  df_res_svm <- data.frame(PC0 = PC0_svm, PC1 = PC1_svm, RC0 = RC0_svm, RC1 = RC1_svm,
                           F10 = F10_svm, F11 = F11_svm, MPC = MPC_svm, MRC = MRC_svm,
                           MF1 = MF1_svm, ROC_AUC = AUC_svm, PR_AUC = PR_AUC_svm, MCC = MCC_svm)
  
  df_res_svm <- as.data.frame(t(df_res_svm))
  df_res_svm[is.na(df_res_svm)] <- 0
  df_res_svm$mean <- apply(df_res_svm, 1, mean)
  df_res_svm$sd <- apply(df_res_svm, 1, sd)
  df_res_svm$CI95 <- 2 * sqrt((df_res_svm$mean * (1 - df_res_svm$mean)) / nround)
  df_res_svm <- round(df_res_svm, 3)
  
  return(df_res_svm)
}


Res_CMT <- function(prob,pred,test_data,model) {
  #AUC 
  prob = predict(model, test_data, type = "prob")
  pred = predict(model, test_data)
  #AUC 
  roc_c <- pROC::auc(test_data$SUBTYPE, prob[,2])
  ROC_AUC <- roc_c[1]
  
  PR_AUC <- prauc(truth = as.factor(test_data$SUBTYPE) , prob=prob[, 2]
                  ,positive= minor_group) ## positive, we can choose minority group as the positive
  
  #confusion matrix
  actual = as.factor(test_data$SUBTYPE)
  predicted = as.factor(pred)
  cm = as.matrix(table(Actual = actual, Predicted = predicted))
  
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of SUBTYPEes
  diag = diag(cm) # number of correctly SUBTYPEified instances per SUBTYPE 
  rowsums = apply(cm, 1, sum) # number of instances per SUBTYPE
  colsums = apply(cm, 2, sum) # number of predictions per SUBTYPE
  p = rowsums / n # distribution of instances over the actual SUBTYPEes
  q = colsums / n # distribution of instances over the predicted SUBTYPEes
  
  precision = diag / colsums 
  recall = diag / rowsums 
  f1 = 2 * precision * recall / (precision + recall) 
  #data.frame(precision, recall, f1) 
  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1)
  mcc = mltools::mcc(preds = predicted,actuals = actual)
  
  
  #All results 
  Res <- matrix(c(precision,recall,f1,macroPrecision,macroRecall, macroF1, PR_AUC, ROC_AUC,mcc))
  return(Res) #we tell the function what to output
}

############

# comparison <- c("Basal","LumA")
# # We will use both train and test data 
# voom_tt0.05 <- read.csv('BasalvsLumA_DEG_train_tgca_subtype.csv',header = T, row.names = 1)
# comparison_header  <- paste(comparison, collapse = 'vs')
results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)

# tcga_data <- fread("all_tcga_combat_corrected.csv",header = T)
# setDF(tcga_data)
# row.names(tcga_data) <- tcga_data$V1
# tcga_data <- tcga_data[,-1]

# pheno_tcga <- read.csv("phenotype_all_tcga.csv",header = T)
# row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
# pheno_tcga$X <- NULL
# pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison,]
# tcga_data <- tcga_data[,colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

#data_model <- read.csv( "data_model_tcga_subtype.csv",header = T,row.names = 1)

#Freq = 25 times 
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq25)
dim(data_model) #filter = 25 times, 240 genes
#664 241

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq25 <- Model_50_RF(data = data_model,1)
#View(res_tcga_subtype_rf_freq25)

res_tcga_subtype_svm_freq25 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq25)


#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq25,paste(results_base,"res_tcga_subtype_rf_freq25.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq25,paste(results_base,"res_tcga_subtype_svm_freq25.csv",sep="/"))

##########

#Freq freq40 times
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq40)

dim(data_model) #freq = 40 times, 506 genes
#664 207

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq40 <- Model_50_RF(data = data_model,5)
#View(res_tcga_subtype_rf_freq40)

res_tcga_subtype_svm_freq40 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq40)

#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq40,paste(results_base,"res_tcga_subtype_rf_freq40.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq40,paste(results_base,"res_tcga_subtype_svm_freq40.csv",sep="/"))

###### 
#Freq freq50 times 
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq50)

dim(data_model) #freq = 15 times, 7 genes
#664 138

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq50 <- Model_50_RF(data = data_model,5)
#View(res_tcga_subtype_rf_freq50)

res_tcga_subtype_svm_freq50 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq50)

#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq50,paste(results_base,"res_tcga_subtype_rf_freq50.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq50,paste(results_base,"res_tcga_subtype_svm_freq50.csv",sep="/"))

sum_model_ALL <- data.frame(metrics = rownames(res_tcga_subtype_rf_freq50),
                        rf_freq25= res_tcga_subtype_rf_freq25$mean,
                        svm_freq25= res_tcga_subtype_svm_freq25$mean,
                        rf_freq40= res_tcga_subtype_rf_freq40$mean,
                        svm_freq40= res_tcga_subtype_svm_freq40$mean,
                        rf_freq50= res_tcga_subtype_rf_freq50$mean,
                        svm_freq50= res_tcga_subtype_svm_freq50$mean)
#View(sum_model_ALL)

write.csv(sum_model_ALL,paste(results_base,"sum_model_tcga_subtype.csv",sep="/"))

save.image(paste(results_base,"model_tcga_subtype.rdata",sep="/"))

df_genes <- data.frame(freq50 = boruta_features50)
write.csv(df_genes,paste(results_base,"df_genes_freq50_subtypes.csv",sep="/"))
