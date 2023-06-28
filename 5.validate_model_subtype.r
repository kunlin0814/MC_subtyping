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
library(sva)
library(smotefamily)
prev_results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
load(paste(prev_results_base,'model_tcga_subtype.rdata',sep="/"))
Res_CMT <- function(prob,pred,test_data,model) {
  #AUC 
  prob = predict(model, test_data, type = "prob")
  pred = predict(model, test_data)
  #AUC 
  roc_c <- pROC::auc(test_data$SUBTYPE, prob[,2])
  ROC_AUC <- roc_c[1]
  
  PR_AUC <- prauc(truth = as.factor(test_data$SUBTYPE) , prob=prob[, 2]
                  ,positive= comparison[2])
  
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

data_test<- read.csv("cmt_combat_corrected.csv",header = T, row.names = 1)
pheno_test <- read.csv("phenotype_cmt.csv",header = T)
pheno_train <- pheno_tcga
dim(data_test)
#  11856    78

data_test <- data.frame(data_test)

train_freq25 <- data.frame(t(pca_data_freq25))
train_freq40 <-data.frame(t(pca_data_freq40))
train_freq50 <- data.frame(t(pca_data_freq50))

dim(train_freq25)

test_freq25 <- data_test[ row.names(data_test) %in% row.names(train_freq25),]
test_freq40 <- data_test[ row.names(data_test) %in% row.names(train_freq40),]
test_freq50 <- data_test[ row.names(data_test) %in% row.names(train_freq50),]
dim(test_freq25)

abc <- colnames(test_freq25)== pheno_test$PATIENT_ID # ALL TRUE


RF_validate <- function(train_freq, test_freq) {
  # freq25 times
  train_data <- data.frame(SUBTYPE= factor(pheno_train$SUBTYPE), t(train_freq))
  test_data <- data.frame(SUBTYPE= factor(pheno_test$SUBTYPE), t(test_freq))
  
  
  if (sum(table(train_data$SUBTYPE)) >= 0.8 * sum(table(train_data$SUBTYPE)) &&
      sum(table(train_data$SUBTYPE)) <= 1.2 * sum(table(train_data$SUBTYPE))) {
    train_data_balanced <- train_data
  } else {
    # We use ADASYN
    genData_ADAS <- ADAS(X = train_data[, -1], target = train_data$SUBTYPE, K = 5)
    train_df_adas <- genData_ADAS[["data"]]
    class_df <- train_df_adas[ncol(train_df_adas)]
    train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
    train_data_balanced <- train_df_adas
    names(train_data_balanced)[1] <- "SUBTYPE"
  }
  
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
  res_rf_freq <- Res_CMT(prob = predict(rf_model, test_data, type = "prob"),
                         pred = predict(rf_model, test_data),
                         test_data = test_data, model = rf_model)
  print(res_rf_freq)
}

res_rf_freq25 = RF_validate(train_freq25, test_freq25)
res_rf_freq40 = RF_validate(train_freq40, test_freq40)
res_rf_freq50 = RF_validate(train_freq50, test_freq50)

test_model_sum <- data.frame(freq25 = res_rf_freq25,
                             freq40 = res_rf_freq40,
                             freq50 = res_rf_freq50)

View(test_model_sum)

write.csv(test_model_sum, "model_sum_validate_subtype.csv")
save.image("validate_model_subtype.rdata")

#PCA


PCA_plot2(t(test_freq25),pheno_test$SUBTYPE, title ="Selected genes-cmt_subtype(Freq=25,n=240)")
PCA_plot2(t(test_freq25),pheno_test$SUBTYPE, title ="Selected genes-cmt_subtype(Freq=40,n=206)")
PCA_plot2(t(test_freq25),pheno_test$SUBTYPE, title ="Selected genes-cmt_subtype(Freq=50,n=137)")


