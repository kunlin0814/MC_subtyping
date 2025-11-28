# Script: 04_model_training.R
# Description: Trains RF and SVM models using selected features.

source("R/utils.R")
config <- load_config()
load_dependencies(c("data.table", "caret", "randomForest", "e1071", "ranger", "kernlab", "smotefamily", "mlr3measures", "pROC", "MLmetrics", "mltools"))
source("R/analysis.R")

# --- Parameters ---
comparison <- c("Basal", "LumA")
comparison_header <- paste(comparison, collapse = "vs")
feature_selection_path <- file.path(config$paths$results_dir, paste0("Step3_feature_selection_", comparison_header))
results_path <- file.path(config$paths$results_dir, paste0("Step4_model_training_", comparison_header))
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
nround <- config$parameters$model_rounds
cv_folds <- config$parameters$model_cv_folds
tune_length <- config$parameters$model_tune_length

# --- Load Data ---
message("Loading data...")
tcga_data_path <- file.path(config$paths$processed_data_dir, "all_tcga_combat_corrected.csv")
pheno_path <- file.path(config$paths$processed_data_dir, "phenotype_all_tcga.csv")

tcga_data <- fread(tcga_data_path, header = TRUE)
setDF(tcga_data)
rownames(tcga_data) <- tcga_data$V1
tcga_data <- tcga_data[, -1]

pheno_tcga <- read.csv(pheno_path, header = TRUE)
rownames(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison, ]
tcga_data <- tcga_data[, colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

# Load Gene Lists
gene_list_25 <- read.table(file.path(feature_selection_path, "gene_list_freq25.txt"), header = TRUE, stringsAsFactors = FALSE)$gene
gene_list_50 <- read.table(file.path(feature_selection_path, "gene_list_freq50.txt"), header = TRUE, stringsAsFactors = FALSE)$gene

# --- Model Training Function ---
run_models <- function(genes, freq_label) {
    message("Training models for ", freq_label, " (n=", length(genes), ")...")

    data_subset <- t(tcga_data[rownames(tcga_data) %in% genes, ])
    data_subset <- as.data.frame(data_subset)
    data_subset$SUBTYPE <- factor(pheno_tcga$SUBTYPE)

    # Random Forest
    message("  Running Random Forest...")
    rf_res <- Model_50_RF(data_subset, nround = nround, cv_folds = cv_folds, tune_length = tune_length)
    write.csv(rf_res, file.path(results_path, paste0("RF_results_", freq_label, ".csv")))

    # SVM
    message("  Running SVM...")
    svm_res <- Model_50_SVM(data_subset, nround = nround, cv_folds = cv_folds, tune_length = tune_length)
    write.csv(svm_res, file.path(results_path, paste0("SVM_results_", freq_label, ".csv")))
}

# Run for both gene lists
if (length(gene_list_25) > 0) run_models(gene_list_25, "freq25")
if (length(gene_list_50) > 0) run_models(gene_list_50, "freq50")

message("Step 04 complete.")
