library(data.table)
library(MultiBaC)
library(MultiAssayExperiment)
library(sva)
library(pvca)
#library(DMwR)
library(grid)
library(Boruta)
library(vars)
library(varSelRF)
library(caret)
library(e1071)
library(ranger)
library(sva)
library(smotefamily)
library(mlr3measures)
library(kernlab)
library(edgeR)
library(limma)
library(tmaptools)
library(tidyverse)
library(reshape2)
#library(biomaRt)
library(randomForest)

PCA_plot <- function(data,group,scale_status=T,label= "", title = "") {
  Data_pca <- prcomp(scale(data,center = T,scale = scale_status))
  contributions <- summary(Data_pca)$importance[2,]*100
  dataScores12 <-data.frame(Data_pca$x[, 1:2])
  
  if (label == "") {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ 
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16) + ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ geom_text(aes(label=label))+
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16)+ ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  return(p)
  
}


## phenotype table has the same sample order as the sample order in the colname of the data
CBFpvcaFunction <- function(data,phenotypedata){
  
  #packages <- c("pvca", "ggplot2", "Biobase") 
  #lapply(packages, library, character.only = TRUE, quietly = TRUE)
  
  if(length(seq_len(nrow(phenotypedata))) != length(seq_len(nrow(data)))) {
    warning(paste0("Data needs to be a wide", "\n", "Transposing data", "\n"))
    data <- t(data)
  }
  
  
  if(!identical(rownames(phenotypedata), rownames(data))) {
    warning(paste0("Rownames in phenotypedata do not match the rownames of data", "\n"))
    warning(paste0("Matching rownames in phenotypedata to data", "\n"))
    
    if(nrow(phenotypedata) != nrow(data)) {
      stop("Phenotypedata and data have different amounts of rows")
    } else {
      rownames(phenotypedata) <- rownames(data)
    }
  }
  
  data1_pheno_formatted<-new("AnnotatedDataFrame",data=phenotypedata)
  
  data_forPVCA<-ExpressionSet(assayData=t(data),phenoData=data1_pheno_formatted)
  
  pvcaObj_data <- pvcaBatchAssess (abatch = data_forPVCA, batch.factors = colnames(data1_pheno_formatted@data), threshold = 0.7)
  
  data_pvca_res<-data.frame(as.data.frame(pvcaObj_data$label),t(as.data.frame(pvcaObj_data$dat)))
  
  colnames(data_pvca_res)<-c("effect","variance")
  
  
  p<-ggplot((data_pvca_res[-nrow(data_pvca_res),]), aes(x= effect, y = variance)) + 
    geom_bar(stat = 'identity',  position = 'dodge', col ='transparent')+ 
    scale_fill_discrete(guide = 'none') + 
    theme_bw(base_size = 16)+
    theme(plot.title =  element_text(hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))+
    labs(x = 'Effects', y = 'Weighted average proportion variance')+
    ggtitle("PVCA estimation bar chart corrected data")
  
  return(p)
  
}

### We use ADAS algorithm to overampling the minority class to deal with imbalance dataset
Balance_data <- function(data){
  train_indices <- createDataPartition(data$SUBTYPE, p = 0.8, list = FALSE, times = 1)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  genData_ADAS <- ADAS(X = train_data[, -1], target = train_data$SUBTYPE, K = 5)
  train_df_adas <- genData_ADAS[["data"]]
  class_df <- train_df_adas[ncol(train_df_adas)]
  train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
  train_data_balanced <- train_df_adas
  names(train_data_balanced)[1] <- "SUBTYPE"
  
  train_data_balanced <- train_data_balanced %>%
    mutate(across(.cols = 1, .fns = factor)) %>%
    mutate(across(.cols = 2:length(.), .fns = as.numeric))
  
  return(train_data_balanced)
}

Model_50_RF <- function(data,nround) {
  PC0_rf <- 0 ; PC1_rf <- 0;RC0_rf <- 0;RC1_rf <- 0;F10_rf <- 0;F11_rf <- 0;MF1_rf<-0 ; PR_AUC_rf <- 0;AUC_rf <- 0;MPC_rf <-0 ; MRC_rf <-0 ; MCC_rf <-0
  nround = nround
  for (i in 1:nround) {
    set.seed(i)
    print(i)
    train_data_balanced <-Balance_data(data)
    
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
    train_data_balanced <-Balance_data(data)
    
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

