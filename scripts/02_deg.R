# Script: 02_deg.R
# Description: Identifies differentially expressed genes between subtypes using limma/voom.

source("R/utils.R")
config <- load_config()
load_dependencies(c("data.table", "limma", "edgeR", "ggplot2", "EnhancedVolcano", "caret"))
source("R/plotting.R")

# --- Parameters ---
comparison <- c("Basal", "LumA") # Could be moved to config if needed
comparison_header <- paste(comparison, collapse = "vs")
results_path <- file.path(config$paths$results_dir, paste0("Step2DEG_", comparison_header))
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
message("Loading data...")
tcga_data_path <- file.path(config$paths$processed_data_dir, "all_tcga_combat_corrected.csv")
pheno_path <- file.path(config$paths$processed_data_dir, "phenotype_all_tcga.csv")

if (!file.exists(tcga_data_path)) stop("TCGA data not found: ", tcga_data_path)
if (!file.exists(pheno_path)) stop("Phenotype data not found: ", pheno_path)

tcga_data <- fread(tcga_data_path, header = TRUE)
setDF(tcga_data)
rownames(tcga_data) <- tcga_data$V1
tcga_data <- tcga_data[, -1]

pheno_tcga <- read.csv(pheno_path, header = TRUE)
rownames(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga$X <- NULL

# Filter for comparison groups
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison, ]
tcga_data <- tcga_data[, colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

# Transpose for processing (samples as rows for some steps, but DGEList needs genes as rows)
# The original script transposes back and forth. Let's be careful.
# DGEList expects counts (or TPMs) with genes in rows and samples in columns.
# tcga_data is currently genes x samples.

# Remove negative counts (log2TPM shouldn't be negative unless < 0, but log2(TPM+1) is >= 0)
# Check for negative values
genes_with_negative_counts <- rownames(tcga_data)[apply(tcga_data, 1, function(x) any(x < 0))]
if (length(genes_with_negative_counts) > 0) {
    message("Removing ", length(genes_with_negative_counts), " genes with negative counts.")
    tcga_data <- tcga_data[!rownames(tcga_data) %in% genes_with_negative_counts, ]
}

# Split Train/Test (80/20)
set.seed(config$parameters$seed)
train_indices <- createDataPartition(pheno_tcga$SUBTYPE, p = 0.8, list = FALSE, times = 1)

# We need to subset the expression data columns based on the phenotype rows
train_samples <- rownames(pheno_tcga)[train_indices]
test_samples <- rownames(pheno_tcga)[-train_indices]

tcga_data_train <- tcga_data[, train_samples]
pheno_train <- pheno_tcga[train_samples, ]

# --- Differential Expression Analysis ---
message("Running Differential Expression Analysis...")
dge <- DGEList(counts = tcga_data_train) # Note: These are log2TPM, not raw counts. Voom handles this?
# Voom is designed for counts, but can work with log-counts if configured?
# Original script used voom on log2TPM. This is technically not ideal (voom expects counts),
# but we are replicating the workflow.
# Actually, if input is log2TPM, we might just use limma directly without voom, or assume it's fine.
# The original script used `voom(tcga_data_dge, design, plot = TRUE)`.

group <- factor(pheno_train$SUBTYPE)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

voom_data <- voom(dge, design, plot = FALSE)
fit <- lmFit(voom_data, design)

# Contrast
contrast_formula <- paste(comparison, collapse = "-")
cont_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# Results
results <- topTable(fit2, number = Inf)
results_sig <- results[results$adj.P.Val <= config$parameters$p_value_cutoff, ]

message("Found ", nrow(results_sig), " significant genes.")

# --- Save Results ---
write.csv(results_sig, file.path(results_path, paste0(comparison_header, "_DEG_train_tcga_subtype.csv")))

# --- Visualization ---
# Volcano Plot
p_volcano <- EnhancedVolcano(results,
    lab = rownames(results),
    x = "logFC",
    y = "adj.P.Val",
    title = paste(comparison_header, "Volcano Plot"),
    pCutoff = config$parameters$p_value_cutoff,
    FCcutoff = 0.5
)
ggsave(file.path(results_path, "volcano_plot.pdf"), p_volcano, width = 8, height = 6)

# PCA of DEGs
pca_data <- t(tcga_data_train[rownames(tcga_data_train) %in% rownames(results_sig), ])
p_pca <- PCA_plot(pca_data, pheno_train$SUBTYPE, title = paste("DEGs PCA (n=", nrow(results_sig), ")"))
ggsave(file.path(results_path, "pca_plot_deg.pdf"), p_pca, width = 6, height = 4)

message("Step 02 complete.")
