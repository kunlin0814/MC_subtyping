source('C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')

base <- #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
  "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

comparison <- c("Basal", "LumA")
comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <-
  paste(base, 'Step2DEG', comparison_header, sep = "/")

## gene file from limma with p<0.05
voom_tt0.05 <- read.csv(paste(
  prev_results_base,
  paste(comparison_header, "DEG_train_tgca_subtype.csv", sep = "_")
  ,
  sep = "/"
), row.names = 1)


results_base <-
  paste(base, 'Step3_feature_selection', comparison_header, sep = "/")
dir.create(results_base, recursive = TRUE)

# We will use both train and test data to do the feature selection
tcga_data <-
  fread(paste(base, "all_tcga_combat_corrected.csv", sep = '/'),
        header = T)
setDF(tcga_data)
row.names(tcga_data) <- tcga_data$V1
tcga_data <- tcga_data[, -1]

pheno_tcga <-
  read.csv(paste(base, "phenotype_all_tcga.csv", sep = '/'), header = T)
row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID

pheno_tcga$X <- NULL
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison, ]
tcga_data <-
  tcga_data[, colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

data <-
  data.frame(t(tcga_data[rownames(tcga_data) %in% rownames(voom_tt0.05), ]))
data$SUBTYPE <- factor(pheno_tcga$SUBTYPE)

# Perform Boruta feature selection with 50 runs
boruta_results <- list()

for (i in 1:50) {
  set.seed(i)  # Set seed for reproducibility
  
  train_indices <- createDataPartition(data$SUBTYPE,
                                       p = 0.8,
                                       list = FALSE,
                                       times = 1)
  train_data <- data[train_indices,]
  test_data <- data[-train_indices,]
  
  boruta_genes <- Boruta(SUBTYPE ~ ., train_data)
  boruta_fixed <- TentativeRoughFix(boruta_genes)
  boruta_final <-
    getSelectedAttributes(boruta_fixed, withTentative = FALSE)
  print(paste("boruta_round:", i))
  boruta_results[[i]] <- boruta_final
  
}

# Create frequency tables for Boruta results
gene_list <- unlist(boruta_results)
boruta_freq <- table(unlist(boruta_results))

boruta_df <-
  data.frame(feature = names(boruta_freq),
             frequency = as.numeric(boruta_freq))

# Sort
boruta_df <-
  boruta_df[order(boruta_df$frequency, decreasing = TRUE),]

# Apply cutoffs to filter features
boruta_filtered25 <-
  boruta_df[boruta_df$frequency >= 25, , drop = FALSE]
boruta_filtered50 <-
  boruta_df[boruta_df$frequency >= 50, , drop = FALSE]


# Extract filtered features
boruta_features25 <- boruta_filtered25$feature
boruta_features50 <- boruta_filtered50$feature


tcga_data_t <- data.frame(t(tcga_data))

#PCA
pca_data_freq25 <-
  tcga_data_t[, colnames(tcga_data_t) %in% boruta_features25]
pca_data_freq50 <-
  tcga_data_t[, colnames(tcga_data_t) %in% boruta_features50]

gene_list_freq25 <- data.table(gene = colnames(pca_data_freq25))
gene_list_freq50 <- data.table(gene = colnames(pca_data_freq50))

fwrite(
  gene_list_freq25,
  file = paste(
    results_base,
    paste(comparison_header, "gene_list_freq25.txt", sep = "_"),
    sep = "/"
  )
  ,
  eol = "\n",
  row.names = F,
  col.names = T
)

fwrite(
  gene_list_freq50,
  file = paste(
    results_base,
    paste(comparison_header, "gene_list_freq50.txt", sep = "_"),
    sep = "/"
  )
  ,
  eol = "\n",
  row.names = F,
  col.names = T
)

## Visualization of the PCA from selected genes
# pdf(file=paste(results_base,
#                paste(comparison_header,"Selected genes-tcga_subtype_Freq_25.pdf",sep=""),
#                sep="/") ,height = 4.5,width = 6)
# p <- PCA_plot(pca_data_freq25,pheno_tcga$SUBTYPE,
#               title =paste("Selected genes-tcga_subtype Freq=25"," n=",length(boruta_features25),sep=""))
# print(p)
# dev.off()
# pdf(file=paste(results_base,
#                paste(comparison_header,"Selected genes-tcga_subtype_Freq_50.pdf",sep=""),
#                sep="/")
#     ,height = 4.5,width = 6)
# p <- PCA_plot(pca_data_freq50,pheno_tcga$SUBTYPE,
#               title =paste("Selected genes-tcga_subtype Freq=50,","n=",length(boruta_features50),sep=""))
# print(p)
# dev.off()

# save.image(paste(results_base,"feature_selection_tcga_subtype.rdata",sep ="/"))
