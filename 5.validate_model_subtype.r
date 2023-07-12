# This script utilizes the genes selected by the Boruta algorithm (selected 25 or 50 times) as features to create and evaluate the performance of a random forest model using the dog dataset.
# Note: 
# 1. The script automatically identifies the group with the minimum number of samples between two classes as the minority group.
# 2. Prior to running this script, make sure to execute the "3.feature_selection_tcga_subtype.r" script to obtain the Boruta-selected genes data.

# The main objective of this script is to train a random forest model using the selected genes from the human dataset as input features and test its performance on the dog dataset.

source('C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')

base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
#"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

comparison <- c("Basal", "LumA")
comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <-
  paste(base, 'Step3_feature_selection', comparison_header, sep = "/")
results_base <-
  paste(base, 'Step5_model_validation', comparison_header, sep = "/")

Res_CMT <- function(prob, pred, test_data, model) {
  prob <-  predict(model, test_data, type = "prob")
  pred <-  predict(model, test_data)
  #AUC
  roc_c <- pROC::auc(test_data$SUBTYPE, prob[, 2])
  ROC_AUC <- roc_c[1]
  
  PR_AUC <-
    prauc(
      truth = as.factor(test_data$SUBTYPE) ,
      prob = prob[, 2]
      ,
      positive = minor_group
    )
  
  #confusion matrix
  actual <-  as.factor(test_data$SUBTYPE)
  predicted <-  as.factor(pred)
  cm <-  as.matrix(table(Actual = actual, Predicted = predicted))
  
  n <-  sum(cm) # number of instances
  nc <-  nrow(cm) # number of SUBTYPEes
  diag <-  diag(cm) # number of correctly SUBTYPEified instances per SUBTYPE
  rowsums <-  apply(cm, 1, sum) # number of instances per SUBTYPE
  colsums <-  apply(cm, 2, sum) # number of predictions per SUBTYPE
  p = rowsums / n # distribution of instances over the actual SUBTYPEes
  q = colsums / n # distribution of instances over the predicted SUBTYPEes
  
  precision <-  diag / colsums
  recall <-  diag / rowsums
  f1 <-  2 * precision * recall / (precision + recall)
  #data.frame(precision, recall, f1)
  macroPrecision <-  mean(precision)
  macroRecall <-  mean(recall)
  macroF1 <-  mean(f1)
  mcc = mltools::mcc(preds = predicted, actuals = actual)
  
  
  #All results
  Res <-
    matrix(c(
      precision,
      recall,
      f1,
      macroPrecision,
      macroRecall,
      macroF1,
      PR_AUC,
      ROC_AUC,
      mcc
    ))
  return(Res) 
}
RF_validate <- function(train_freq, test_freq) {
  train_data <-
    data.frame(SUBTYPE = factor(pheno_train$SUBTYPE), t(train_freq))
  test_data <-
    data.frame(SUBTYPE = factor(pheno_test$SUBTYPE), t(test_freq))
  
  
  if (sum(table(train_data$SUBTYPE)) >= 0.8 * sum(table(train_data$SUBTYPE)) &&
      sum(table(train_data$SUBTYPE)) <= 1.2 * sum(table(train_data$SUBTYPE))) {
    train_data_balanced <- train_data
  } else {
    # We use ADASYN (ADASYN is an method to oversampling the minority class)
    genData_ADAS <-
      ADAS(X = train_data[, -1],
           target = train_data$SUBTYPE,
           K = 5)
    train_df_adas <- genData_ADAS[["data"]]
    class_df <- train_df_adas[ncol(train_df_adas)]
    train_df_adas <-
      cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
    train_data_balanced <- train_df_adas
    names(train_data_balanced)[1] <- "SUBTYPE"
  }
  
  train_data_balanced <- train_data_balanced %>%
    mutate(across(.cols = 1, .fns = factor)) %>%
    mutate(across(.cols = 2:length(.), .fns = as.numeric))
  train_control <- trainControl(method = "cv",
                                number = 5,
                                classProbs = TRUE)
  
  # Random Forest
  rf_model <-
    caret::train(
      SUBTYPE ~ .,
      data = train_data_balanced,
      method = "ranger",
      verbose = FALSE,
      trControl = train_control,
      metric = 'Accuracy',
      tuneLength = 10
    )
  # Prediction
  res_rf_freq <-
    Res_CMT(
      prob = predict(rf_model, test_data, type = "prob"),
      pred = predict(rf_model, test_data),
      test_data = test_data,
      model = rf_model
    )
  #print(res_rf_freq)
  return (res_rf_freq)
}

dir.create(results_base, recursive = TRUE)
data_test <-
  read.csv(
    paste(base, "all_cmt_combat_corrected.csv", sep = "/"),
    header = T,
    row.names = 1
  )
pheno_test <-
  read.csv(paste(base, "phenotype_cmt.csv", sep = "/"), header = T)

pheno_tcga <-
  read.csv(paste(base, "phenotype_all_tcga.csv", sep = '/'), header = T)
row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga$X <- NULL
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison, ]
pheno_train <- pheno_tcga
each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
minor_group <-
  each_pheno_num[each_pheno_num$n == min(each_pheno_num$n),]$SUBTYPE
pca_data_freq25 <- read.table(
  paste(
    prev_results_base,
    paste(comparison_header, 'PCA_data_gene_list_freq25.txt', sep =
            "_")
    ,
    sep = '/'
  ),
  stringsAsFactors = F,
  row.names = 1
)

pca_data_freq50 <- read.table(
  paste(
    prev_results_base,
    paste(comparison_header, 'PCA_data_gene_list_freq50.txt', sep =
            "_")
    ,
    sep = '/'
  ),
  stringsAsFactors = F,
  row.names = 1
)


if (toupper('basal') %in% toupper(comparison)) {
  other_group <- comparison[!comparison %in% "Basal"]
  pheno_test[pheno_test$SUBTYPE != "basal",]$SUBTYPE <- other_group
  pheno_test[pheno_test$SUBTYPE == "basal",]$SUBTYPE <- "Basal"
} else if (!toupper('basal') %in% toupper(comparison) &
           toupper('LumA') %in% toupper(comparison)) {
  other_group <- comparison[!comparison %in% "LumA"]
  pheno_test[pheno_test$SUBTYPE != "basal",]$SUBTYPE <- "LumA"
  pheno_test[pheno_test$SUBTYPE == "basal",]$SUBTYPE <- other_group
}

data_test <- data.frame(data_test)

train_freq25 <- data.frame(t(pca_data_freq25))
train_freq50 <- data.frame(t(pca_data_freq50))


test_freq25 <-
  data_test[row.names(data_test) %in% row.names(train_freq25),]
test_freq50 <-
  data_test[row.names(data_test) %in% row.names(train_freq50),]


res_rf_freq25 <-  RF_validate(train_freq25, test_freq25)
res_rf_freq50 <-  RF_validate(train_freq50, test_freq50)

test_model_sum <- data.frame(freq25 = res_rf_freq25,
                             freq50 = res_rf_freq50)


rownames(test_model_sum) <-
  c('PC0',
    'PC1',
    'RC0',
    'RC1',
    'F10',
    'F11',
    'MPC',
    'MRC',
    'MF1',
    'ROC_AUC',
    'PR_AUC',
    'MCC')

write.csv(test_model_sum,
          paste(results_base, "model_sum_validate_subtype.csv", sep = "/"))

#Use PCA to show the separation between two subtypes in the test dataset

pdf(
  file = paste(
    results_base,
    paste(
      comparison_header,
      "Selected genes-tcga_subtype_Freq_25.pdf",
      sep = ""
    ),
    sep = "/"
  )
  ,
  height = 4.5,
  width = 6
)
p <- PCA_plot(
  t(test_freq25),
  pheno_test$SUBTYPE,
  scale_status = F,
  title = paste(
    "Selected genes-cmt_subtype Freq=25, n= ",
    nrow(train_freq25),
    sep = ""
  )
)
print(p)
dev.off()

pdf(
  file = paste(
    results_base,
    paste(
      comparison_header,
      "Selected genes-tcga_subtype_Freq_50.pdf",
      sep = ""
    ),
    sep = "/"
  )
  ,
  height = 4.5,
  width = 6
)
p <- PCA_plot(
  t(test_freq50),
  pheno_test$SUBTYPE,
  scale_status = F,
  title = paste(
    "Selected genes-cmt_subtype Freq=50, n= ",
    nrow(train_freq50),
    sep = ""
  )
)

print(p)
dev.off()
