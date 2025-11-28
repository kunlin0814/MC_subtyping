# Script: 03_feature_selection.R
# Description: Selects features using Boruta algorithm.

source("R/utils.R")
config <- load_config()
load_dependencies(c("data.table", "Boruta", "caret", "dplyr"))
source("R/plotting.R")

# --- Parameters ---
comparison <- c("Basal", "LumA")
comparison_header <- paste(comparison, collapse = "vs")
deg_results_path <- file.path(config$paths$results_dir, paste0("Step2DEG_", comparison_header))
results_path <- file.path(config$paths$results_dir, paste0("Step3_feature_selection_", comparison_header))
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
message("Loading data...")
deg_file <- file.path(deg_results_path, paste0(comparison_header, "_DEG_train_tcga_subtype.csv"))
if (!file.exists(deg_file)) stop("DEG file not found: ", deg_file)

deg_list <- read.csv(deg_file, row.names = 1)
deg_genes <- rownames(deg_list)

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

# Subset to DEG genes
data_subset <- t(tcga_data[rownames(tcga_data) %in% deg_genes, ])
data_subset <- as.data.frame(data_subset)
data_subset$SUBTYPE <- factor(pheno_tcga$SUBTYPE)

# --- Boruta Feature Selection ---
message("Running Boruta Feature Selection (", config$parameters$boruta_runs, " runs)...")

boruta_results <- list()
n_runs <- config$parameters$boruta_runs

for (i in 1:n_runs) {
    set.seed(i)
    # Split logic inside loop as per original script?
    # Original script: "We will use both train and test data to do the feature selection" -> Wait.
    # Original script loop:
    # train_indices <- createDataPartition...
    # train_data <- data[train_indices, ]
    # boruta_genes <- Boruta(SUBTYPE ~ ., train_data)

    train_indices <- createDataPartition(data_subset$SUBTYPE, p = config$parameters$boruta_train_split, list = FALSE)
    train_data <- data_subset[train_indices, ]

    boruta_out <- Boruta(SUBTYPE ~ ., data = train_data, doTrace = 0)
    boruta_fixed <- TentativeRoughFix(boruta_out)
    boruta_final <- getSelectedAttributes(boruta_fixed, withTentative = FALSE)

    boruta_results[[i]] <- boruta_final
    if (i %% 10 == 0) message("Round ", i, " complete.")
}

# Frequency Analysis
boruta_freq <- table(unlist(boruta_results))
boruta_df <- data.frame(feature = names(boruta_freq), frequency = as.numeric(boruta_freq)) %>%
    arrange(desc(frequency))

# Filter
boruta_filtered25 <- boruta_df[boruta_df$frequency >= 25, ]
boruta_filtered50 <- boruta_df[boruta_df$frequency >= 50, ]

message("Selected ", nrow(boruta_filtered25), " features (freq >= 25)")
message("Selected ", nrow(boruta_filtered50), " features (freq >= 50)")

# --- Save Results ---
write.csv(boruta_df, file.path(results_path, "boruta_frequency.csv"), row.names = FALSE)
write.table(boruta_filtered25$feature, file.path(results_path, "gene_list_freq25.txt"),
    row.names = FALSE, col.names = "gene", quote = FALSE
)
write.table(boruta_filtered50$feature, file.path(results_path, "gene_list_freq50.txt"),
    row.names = FALSE, col.names = "gene", quote = FALSE
)

# Save PCA data for selected genes
pca_data_freq25 <- t(tcga_data[rownames(tcga_data) %in% boruta_filtered25$feature, ])
pca_data_freq50 <- t(tcga_data[rownames(tcga_data) %in% boruta_filtered50$feature, ])

write.csv(pca_data_freq25, file.path(results_path, "PCA_data_gene_list_freq25.csv"))
write.csv(pca_data_freq50, file.path(results_path, "PCA_data_gene_list_freq50.csv"))

message("Step 03 complete.")
