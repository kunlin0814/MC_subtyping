

# feature_selection_tcga_tcga
library(Boruta)
library(vars)
library(dplyr)
library(varSelRF)
library(xgboost)
library(caret)
library(e1071)
library(ranger)

# We will use both train and test data 
data <- data.frame(t(tcga_data[rownames(tcga_data)%in% rownames(voom_tt0.05),]))
dim(data)

data$SUBTYPE <- factor(pheno_tcga$SUBTYPE)
table(data$SUBTYPE )

#basal nonbasal 
# 171      493 

# Perform Boruta feature selection with 50 runs

boruta_results <- vector("list", 50)
#varsel_results <- vector("list", 50)

for (i in 1:50) {
  set.seed(i)  # Set seed for reproducibility
  
  train_indices <- createDataPartition(data$SUBTYPE,
                                       p = 0.8,
                                       list = FALSE,
                                       times = 1)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Boruta
  boruta_genes <- Boruta(SUBTYPE ~ ., train_data)
  boruta_fixed <- TentativeRoughFix(boruta_genes)
  boruta_final <- getSelectedAttributes(boruta_fixed, withTentative = FALSE)
  print(paste("boruta_round:", i))
  print(boruta_final)
  boruta_results[[i]] <- boruta_final
  
  # VarselRF
  #varsel_genes <- varSelRF(train_data[, -ncol(train_data)], as.factor(train_data$SUBTYPE))
  #varsel_genes <- varsel_genes[["selected.vars"]]
  #print(paste("varsel_round:", i))
  #print(varsel_genes)
  #varsel_results[[i]] <- varsel_genes
}

# Create frequency tables for Boruta and VarselRF results
boruta_freq <- table(unlist(boruta_results))
#varsel_freq <- table(unlist(varsel_results))

# Convert frequency tables to data frames
boruta_df <- data.frame(feature = names(boruta_freq), frequency = as.numeric(boruta_freq))
#varsel_df <- data.frame(feature = names(varsel_freq), frequency = as.numeric(varsel_freq))

# SOrt 
boruta_df <- boruta_df[order(boruta_df$frequency, decreasing = TRUE), ]
View(boruta_df)
#varsel_df <- varsel_df[order(varsel_df$frequency, decreasing = TRUE), ]
#View(varsel_df)

# Apply cutoffs to filter features
boruta_filtered25 <- boruta_df[boruta_df$frequency >= 25, , drop = FALSE]
#varsel_filtered5 <- varsel_df[varsel_df$frequency >= 5, , drop = FALSE]

boruta_filtered40<- boruta_df[boruta_df$frequency >= 40, , drop = FALSE]
#varsel_filtered10 <- varsel_df[varsel_df$frequency >= 10, , drop = FALSE]

boruta_filtered50 <- boruta_df[boruta_df$frequency >= 50, , drop = FALSE]
#varsel_filtered15 <- varsel_df[varsel_df$frequency >= 15, , drop = FALSE]



# Extract filtered features
boruta_features25 <- boruta_filtered25$feature
#varsel_features5 <- varsel_filtered5$feature
boruta_features40 <- boruta_filtered40$feature
boruta_features50 <- boruta_filtered50$feature


tcga_data_t <- data.frame(t(tcga_data))

#PCA
pca_data_freq25 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features25]
pca_data_freq40 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features40]
pca_data_freq50 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features50]


PCA_plot(pca_data_freq25,pheno_tcga$SUBTYPE, title ="Selected genes-tcga_subtype(Freq=25,n=240)")
PCA_plot(pca_data_freq40,pheno_tcga$SUBTYPE, title ="Selected genes-tcga_subtype(Freq=40,n=206)")
PCA_plot(pca_data_freq50,pheno_tcga$SUBTYPE, title ="Selected genes-tcga_subtype(Freq=50,n=137)")


save.image("feature_selection_tcga_subtype.rdata")

Saveplots(getwd())
