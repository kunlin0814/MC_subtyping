# Script: 01_combat.R
# Description: Performs combat correction to address batch effects between dog and human data.

# Load utilities and config
source("R/utils.R")
config <- load_config()

# Load dependencies
load_dependencies(c("data.table", "sva", "ggplot2", "Biobase", "pvca"))
source("R/plotting.R")

# --- Load Data ---
message("Loading data...")
human_tpm_path <- file.path(config$paths$raw_data_dir, config$files$human_tpm)
dog_tpm_path <- file.path(config$paths$raw_data_dir, config$files$dog_tpm)
metadata_path <- file.path(config$paths$raw_data_dir, config$files$metadata)

if (!file.exists(human_tpm_path)) stop("Human TPM file not found: ", human_tpm_path)
if (!file.exists(dog_tpm_path)) stop("Dog TPM file not found: ", dog_tpm_path)
if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)

A.rna <- fread(human_tpm_path, header = TRUE)
colnames(A.rna)[1] <- "GeneName"

B.rna <- fread(dog_tpm_path, header = TRUE)

# --- Preprocessing ---
message("Preprocessing data...")
overlap_gene_list <- intersect(B.rna$GeneName, A.rna$GeneName)
B.rna <- B.rna[GeneName %in% overlap_gene_list, ]
A.rna <- A.rna[GeneName %in% overlap_gene_list, ]

B.rna <- as.data.frame(B.rna)
A.rna <- as.data.frame(A.rna)
B.rna <- na.omit(B.rna)

rownames(B.rna) <- B.rna$GeneName
rownames(A.rna) <- A.rna$GeneName
A.rna <- A.rna[, -1]
B.rna <- B.rna[, -1]

AB.rna <- cbind(A.rna, B.rna)

phenotype_AB <- read.csv(file = metadata_path, header = TRUE, sep = "\t")
AB.rna <- AB.rna[, colnames(AB.rna) %in% phenotype_AB$PATIENT_ID]

# Filter phenotype to match data
phenotype_AB <- phenotype_AB[match(colnames(AB.rna), phenotype_AB$PATIENT_ID), ]
rownames(phenotype_AB) <- phenotype_AB$PATIENT_ID

A_pheno <- phenotype_AB[phenotype_AB$SOURCE == "tcga", ]
B_pheno <- phenotype_AB[phenotype_AB$SOURCE == "cmt", ]

data_combined <- data.frame(gene_id = rownames(AB.rna), AB.rna)

# --- PVCA (Before Combat) ---
# message("Running PVCA before combat...")
# p <- CBFpvcaFunction(t(data_combined[,-1]), phenotypedata = phenotype_AB[,c("SOURCE","SUBTYPE")])
# ggsave(file.path(config$paths$figures_dir, "pvca_before_combat.pdf"), p, width = 6, height = 4)

# --- ComBat Correction ---
message("Running ComBat correction...")
# modcombat = model.matrix( ~ 1 + factor(SUBTYPE), data = phenotype_AB) # Optional covariate

combat_batch_corrected <- ComBat(
    dat = data_combined[, -1],
    batch = phenotype_AB$SOURCE,
    par.prior = TRUE,
    prior.plots = FALSE # Set to TRUE if you want plots
)

rownames(combat_batch_corrected) <- data_combined$gene_id

tcga_corrected <- data.frame(combat_batch_corrected[, colnames(combat_batch_corrected) %in% A_pheno$PATIENT_ID])
cmt_corrected <- data.frame(combat_batch_corrected[, colnames(combat_batch_corrected) %in% B_pheno$PATIENT_ID])

# --- Save Results ---
message("Saving results...")
write.csv(A_pheno, file.path(config$paths$processed_data_dir, "phenotype_all_tcga.csv"))
write.csv(B_pheno, file.path(config$paths$processed_data_dir, "phenotype_all_cmt.csv"))

write.csv(combat_batch_corrected, file.path(config$paths$processed_data_dir, "all_tcga_cmt_combat_corrected.csv"))
write.csv(tcga_corrected, file.path(config$paths$processed_data_dir, "all_tcga_combat_corrected.csv"))
write.csv(cmt_corrected, file.path(config$paths$processed_data_dir, "all_cmt_combat_corrected.csv"))

# --- PCA Plots ---
message("Generating PCA plots...")
# Before
p1 <- PCA_plot(t(AB.rna), phenotype_AB$SOURCE, title = "TCGA+CMT before combat")
ggsave(file.path(config$paths$figures_dir, "pca_before_combat_source.pdf"), p1, width = 6, height = 4)

# After
data2 <- data.frame(t(combat_batch_corrected))
p2 <- PCA_plot(data = data2, group = phenotype_AB$SOURCE, title = "TCGA+CMT after combat")
ggsave(file.path(config$paths$figures_dir, "pca_after_combat_source.pdf"), p2, width = 6, height = 4)

message("Step 01 complete.")
