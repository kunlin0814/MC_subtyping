# feature_selection_tcga_tcga
library(Boruta)
library(vars)
library(dplyr)
library(varSelRF)
library(xgboost)
library(caret)
library(e1071)
library(ranger)
#load("E:/My Drive/Josh_MC_Paper_data/ML_gene_set/orig_results/3.feature_selection_tcga_subtype.rdata")

source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
base <- #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
  "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

comparison <- c("Basal","LumB")
comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <- paste(base,'Step2DEG',comparison_header,sep="/")
load(paste(prev_results_base,paste(comparison_header,"_tcga_DEG_subtype.rdata",sep=""),sep='/'))

# We will use both train and test data 

# comparison_header  <- paste(comparison, collapse = 'vs')
results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)


# tcga_data <- fread("all_tcga_combat_corrected.csv",header = T)
# setDF(tcga_data)
# row.names(tcga_data) <- tcga_data$V1
# tcga_data <- tcga_data[,-1]
# 
# pheno_tcga <- read.csv("phenotype_all_tcga.csv",header = T)
# row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
# pheno_tcga$X <- NULL
# pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison,]
# tcga_data <- tcga_data[,colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

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
  #print(boruta_final)
  boruta_results[[i]] <- boruta_final
  
  # VarselRF
  #varsel_genes <- varSelRF(train_data[, -ncol(train_data)], as.factor(train_data$SUBTYPE))
  #varsel_genes <- varsel_genes[["selected.vars"]]
  #print(paste("varsel_round:", i))
  #print(varsel_genes)
  #varsel_results[[i]] <- varsel_genes
}

# Create frequency tables for Boruta and VarselRF results
gene_list <- unlist(boruta_results)
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

pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_25.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq25,pheno_tcga$SUBTYPE, 
              title =paste("Selected genes-tcga_subtype Freq=25"," n=",length(boruta_features25),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_40.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq40,pheno_tcga$SUBTYPE, 
              title =paste("Selected genes-tcga_subtype Freq=40"," n=",length(boruta_features40),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_50.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq50,pheno_tcga$SUBTYPE, 
              title =paste("Selected genes-tcga_subtype Freq=50,","n=",length(boruta_features50),sep=""))
print(p)
dev.off()

save.image(paste(results_base,"feature_selection_tcga_subtype.rdata",sep ="/"))

gene_list_freq25 <-data.table(gene=colnames(pca_data_freq25))
gene_list_freq40 <-data.table(gene=colnames(pca_data_freq40))
gene_list_freq50 <-data.table(gene=colnames(pca_data_freq50))
fwrite(gene_list_freq25, file =paste(results_base,
                        paste(comparison_header,"gene_list_freq25.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)

fwrite(gene_list_freq40, file =paste(results_base,
                                     paste(comparison_header,"gene_list_freq40.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)
fwrite(gene_list_freq50, file =paste(results_base,
                                     paste(comparison_header,"gene_list_freq50.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)
#Saveplots(getwd())
