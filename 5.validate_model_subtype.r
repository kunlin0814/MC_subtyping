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
library(mlr3measures)
library(kernlab)

base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
#"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

comparison <- c("Basal","LumA")
other_group <- comparison[!comparison %in% "Basal"]
comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
load(paste(prev_results_base,'model_tcga_subtype.rdata',sep="/"))

each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
minor_group <- each_pheno_num[each_pheno_num$n ==min(each_pheno_num$n),]$SUBTYPE

data_test<- read.csv(paste(base,"all_cmt_combat_corrected.csv",sep="/"),header = T, row.names = 1)
pheno_test <- read.csv(paste(base,"phenotype_cmt.csv",sep="/"),header = T)
pheno_test[pheno_test$SUBTYPE!="basal",]$SUBTYPE <- other_group
pheno_test[pheno_test$SUBTYPE=="basal",]$SUBTYPE <- "Basal"
pheno_train <- pheno_tcga
dim(data_test)
#  11856    78


Res_CMT <- function(prob,pred,test_data,model) {
  #AUC 
  prob = predict(model, test_data, type = "prob")
  pred = predict(model, test_data)
  #AUC 
  roc_c <- pROC::auc(test_data$SUBTYPE, prob[,2])
  ROC_AUC <- roc_c[1]
  
  PR_AUC <- prauc(truth = as.factor(test_data$SUBTYPE) , prob=prob[, 2]
                  ,positive= minor_group)
  
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
results_base <- paste(base,'Step5_model_validation',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)




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

#View(test_model_sum)

write.csv(test_model_sum, paste(results_base,"model_sum_validate_subtype.csv",sep="/"))
save.image(paste(results_base,"validate_model_subtype.rdata",sep="/"))

#PCA

pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_25.pdf",sep=""),
               sep="/") 
              ,height = 4.5,width = 6)
p <- PCA_plot2(t(test_freq25),pheno_test$SUBTYPE, 
               title =paste("Selected genes-cmt_subtype Freq=25, n= ", nrow(train_freq25),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_40.pdf",sep=""),sep="/")
                ,height = 4.5,width = 6)

p <- PCA_plot2(t(test_freq40),pheno_test$SUBTYPE, 
               title =paste("Selected genes-cmt_subtype Freq=40, n= ", nrow(train_freq40),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_50.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot2(t(test_freq50),pheno_test$SUBTYPE, 
               title =paste("Selected genes-cmt_subtype Freq=25, n= ", nrow(train_freq50),sep=""))

print(p)
dev.off()
