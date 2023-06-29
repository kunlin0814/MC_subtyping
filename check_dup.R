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