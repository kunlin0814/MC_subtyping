# Script: 05_validation.R
# Description: Validates the model on Dog (CMT) data.

source("R/utils.R")
config <- load_config()
load_dependencies(c("data.table", "caret", "randomForest", "ranger", "pROC", "MLmetrics", "mltools", "dplyr"))
source("R/analysis.R")
source("R/plotting.R")

# --- Parameters ---
comparison <- c("Basal", "LumA")
comparison_header <- paste(comparison, collapse = "vs")
feature_selection_path <- file.path(config$paths$results_dir, paste0("Step3_feature_selection_", comparison_header))
results_path <- file.path(config$paths$results_dir, paste0("Step5_model_validation_", comparison_header))
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
cv_folds <- config$parameters$model_cv_folds
tune_length <- config$parameters$model_tune_length

# --- Load Data ---
message("Loading data...")
# Human Data (Training)
tcga_data_path <- file.path(config$paths$processed_data_dir, "all_tcga_combat_corrected.csv")
pheno_tcga_path <- file.path(config$paths$processed_data_dir, "phenotype_all_tcga.csv")

tcga_data <- fread(tcga_data_path, header = TRUE)
setDF(tcga_data)
rownames(tcga_data) <- tcga_data$V1
tcga_data <- tcga_data[, -1]

pheno_tcga <- read.csv(pheno_tcga_path, header = TRUE)
rownames(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison, ]
tcga_data <- tcga_data[, colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

# Dog Data (Testing)
cmt_data_path <- file.path(config$paths$processed_data_dir, "all_cmt_combat_corrected.csv")
pheno_cmt_path <- file.path(config$paths$processed_data_dir, "phenotype_all_cmt.csv")

cmt_data <- fread(cmt_data_path, header = TRUE)
setDF(cmt_data)
rownames(cmt_data) <- cmt_data$V1
cmt_data <- cmt_data[, -1]

pheno_cmt <- read.csv(pheno_cmt_path, header = TRUE)
rownames(pheno_cmt) <- pheno_cmt$PATIENT_ID # Assuming PATIENT_ID exists

# Map Dog Subtypes
# Original logic:
# if 'Basal' in comparison, map non-Basal to 'Other' (or second group)?
# The original script logic is a bit specific:
# if Basal in comparison: Basal -> Basal, others -> comparison[!Basal]
# if LumA in comparison (and not Basal): LumA -> LumA, others -> comparison[!LumA]
# This assumes binary classification where one class is of primary interest.

# Let's replicate this logic safely.
target_subtypes <- comparison
pheno_cmt_filtered <- pheno_cmt

if ("Basal" %in% target_subtypes) {
    other_group <- setdiff(target_subtypes, "Basal")
    pheno_cmt_filtered$SUBTYPE <- ifelse(pheno_cmt_filtered$SUBTYPE == "Basal", "Basal", other_group)
} else if ("LumA" %in% target_subtypes) {
    other_group <- setdiff(target_subtypes, "LumA")
    pheno_cmt_filtered$SUBTYPE <- ifelse(pheno_cmt_filtered$SUBTYPE == "LumA", "LumA", other_group)
} else {
    warning("Neither Basal nor LumA found in comparison. Subtype mapping might be incorrect.")
}

# Filter cmt data to match phenotype
cmt_data <- cmt_data[, colnames(cmt_data) %in% pheno_cmt_filtered$PATIENT_ID]

# Load Gene Lists
gene_list_25 <- read.table(file.path(feature_selection_path, "gene_list_freq25.txt"), header = TRUE, stringsAsFactors = FALSE)$gene
gene_list_50 <- read.table(file.path(feature_selection_path, "gene_list_freq50.txt"), header = TRUE, stringsAsFactors = FALSE)$gene

# --- Validation Function ---
validate_model <- function(genes, freq_label) {
    message("Validating model for ", freq_label, " (n=", length(genes), ")...")

    # Prepare Train Data (Human)
    train_data <- t(tcga_data[rownames(tcga_data) %in% genes, ])
    train_data <- as.data.frame(train_data)
    train_data$SUBTYPE <- factor(pheno_tcga$SUBTYPE)

    # Prepare Test Data (Dog)
    test_data <- t(cmt_data[rownames(cmt_data) %in% genes, ])
    test_data <- as.data.frame(test_data)
    test_data$SUBTYPE <- factor(pheno_cmt_filtered$SUBTYPE, levels = levels(train_data$SUBTYPE))

    # Balance Training Data
    train_data_balanced <- Balance_data(train_data)

    # Train RF
    train_control <- trainControl(method = "cv", number = cv_folds, classProbs = TRUE)
    rf_model <- caret::train(SUBTYPE ~ .,
        data = train_data_balanced, method = "ranger",
        verbose = FALSE, trControl = train_control, metric = "Accuracy",
        tuneLength = tune_length
    )

    # Predict on Dog Data
    res <- Res_CMT(
        prob = predict(rf_model, test_data, type = "prob"),
        pred = predict(rf_model, test_data),
        test_data = test_data, model = rf_model
    )

    return(res)
}

# Run Validation
results_list <- list()

if (length(gene_list_25) > 0) {
    res25 <- validate_model(gene_list_25, "freq25")
    results_list$freq25 <- res25
}

if (length(gene_list_50) > 0) {
    res50 <- validate_model(gene_list_50, "freq50")
    results_list$freq50 <- res50
}

# Combine and Save
if (length(results_list) > 0) {
    results_df <- do.call(cbind, results_list)
    rownames(results_df) <- c(
        "PC0", "PC1", "RC0", "RC1", "F10", "F11",
        "MPC", "MRC", "MF1", "PR_AUC", "ROC_AUC", "MCC"
    )
    write.csv(results_df, file.path(results_path, "model_sum_validate_subtype.csv"))
}

# --- Visualization (PCA on Dog Data) ---
if (length(gene_list_25) > 0) {
    pca_data_dog <- t(cmt_data[rownames(cmt_data) %in% gene_list_25, ])
    p <- PCA_plot(pca_data_dog, pheno_cmt_filtered$SUBTYPE, title = paste("Dog Data PCA (Freq25, n=", length(gene_list_25), ")"))
    ggsave(file.path(results_path, "pca_dog_freq25.pdf"), p, width = 6, height = 4)
}

message("Step 05 complete.")
