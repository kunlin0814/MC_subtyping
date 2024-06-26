# This script performs combat correction to address batch effects between dog data and human data.
# It requires three input files:
# 1. The human TPM file.
# 2. The dog TPM values file.
# 3. The subtype metadata for both human and dog, including sample names and subtype information.

# The goal of this script is to select genes that are shared between human and dog for further analysis.
# These genes will be tested in the dog dataset using the genes selected from the model.

# Note: Combat correction is applied to mitigate the impact of batch effects, enabling more accurate comparison and analysis between the datasets.
# Batch effects can arise due to technical variations in data collection or processing, leading to undesired variations in the results.

# Please make sure to have the necessary input files and the required R packages installed before running this script.

source('C:/github/MC_subtyping/MC_subtyping_module.R')

base <- #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
  "C:/Users/kho/Downloads/Josh_MC_Paper_data/ML_gene_set"
#"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)
A.rna  <- fread("TCGA_BRCA_log2TPM+1.csv", header = T)
colnames(A.rna)[1] <- 'GeneName'

results_base <- paste(base, 'Step1_combat', sep = "/")
dir.create(results_base)

B.rna  <- fread("cmt_all.csv", header = T)
overlap_gene_list <- intersect(B.rna$GeneName, A.rna$GeneName)
B.rna <- B.rna[GeneName %in% overlap_gene_list,]
A.rna <- A.rna[GeneName %in% overlap_gene_list,]
B.rna <- data.frame(B.rna)
A.rna <- data.frame(A.rna)
B.rna <- na.omit(B.rna)
row.names(B.rna) <- B.rna$GeneName
row.names(A.rna) <- A.rna$GeneName
A.rna <- A.rna[, -1]
B.rna <- B.rna[, -1]

AB.rna <- cbind(A.rna, B.rna)

phenotype_AB <-
  read.csv(file = "All_dog_TCGA_subtype_meta.txt", header = T, sep = "\t")
AB.rna <- AB.rna[, colnames(AB.rna) %in% phenotype_AB$PATIENT_ID]
# CMT dataset has basal and nonbasal subtype
B_pheno <- phenotype_AB[phenotype_AB$SOURCE == 'cmt', ]
B.rna <- AB.rna[, colnames(AB.rna) %in% B_pheno$PATIENT_ID]
phenotype_AB <-
  phenotype_AB[match(colnames(AB.rna), phenotype_AB$PATIENT_ID), ]
# human dataset
A_pheno <- phenotype_AB[phenotype_AB$SOURCE == "tcga", ]
# dog dataset
B_pheno <- phenotype_AB[phenotype_AB$SOURCE == "cmt", ]
data <- data.frame(gene_id = rownames(AB.rna), AB.rna)
rownames(phenotype_AB) <- phenotype_AB$PATIENT_ID

## pass data[,-1] (without the gene_id column) to the CBFpcvaFunction.
## PCVA function looks for grouping variables to assess contribution of variance,
## so it is important not to pass any columns that contain unique identifiers.
# pdf(file="CBFpvcaFunction.pdf", width = 6,height = 4)
# p <- CBFpvcaFunction(t(data[,-1]), phenotypedata = phenotype_AB[,c("SOURCE","SUBTYPE")])
# print(p)
# dev.off()


# Here we build a model matrix with  covariates = phenotype_AB$CLASS
modcombat = model.matrix( ~ 1 + factor(SUBTYPE), data = phenotype_AB)

#### Correct SAMPLE
# In case there are covariate confunding with batch, we can ignore mod = modcombat for the ComBat correction
combat_batch_corrected = ComBat(
  dat = data[, -1],
  batch = phenotype_AB$SOURCE,
  #mod = modcombat,
  par.prior = TRUE,
  prior.plots = TRUE
)

rownames(combat_batch_corrected) <- data$gene_id

tcga_corrected <-
  data.frame(combat_batch_corrected[, colnames(combat_batch_corrected) %in%
                                      A_pheno$PATIENT_ID])

cmt_corrected <-
  data.frame(combat_batch_corrected[, colnames(combat_batch_corrected) %in%
                                      B_pheno$PATIENT_ID])

write.csv(A_pheno, paste(base, "phenotype_all_tcga.csv", sep = "/"))
write.csv(B_pheno, paste(base, "phenotype_all_cmt.csv", sep = '/'))

write.csv(combat_batch_corrected,
          paste(base, "all_tcga_cmt_combat_corrected.csv", sep = '/'))
write.csv(tcga_corrected,
          paste(base, "all_tcga_combat_corrected.csv", sep = "/"))
write.csv(cmt_corrected,
          paste(base, "all_cmt_combat_corrected.csv", sep = '/'))

#save.image("all_tcga_combat_2data.rdata")

## plot the PCA cluster before combat
# pdf(file="pca_plot_source.pdf", width = 6,height = 4)
# p <- PCA_plot(t(AB.rna),phenotype_AB$SOURCE ,scale_status = T, title = "TCGA+CMT before combat")
# print(p)
# dev.off()
# pdf(file="pca_plot_sub.pdf", width = 6,height = 4)
# p <- PCA_plot(t(AB.rna),phenotype_AB$SUBTYPE ,scale_status = T, title = "TCGA+CMT before combat")
# print(p)
# dev.off()
# pdf(file="TCGA before combat.pdf", width = 6,height = 4)
# p <- PCA_plot(t(A.rna),A_pheno$SUBTYPE ,scale_status = F, title = "TCGA before combat")
# print(p)
# dev.off()
# pdf(file="CMT before combat.pdf", width = 6,height = 4)
# p <- PCA_plot(t(B.rna),B_pheno$SUBTYPE , scale_status = F, title = "CMT before combat")
# print(p)
# dev.off()
#Plot PCA again after combat
# data2  <- data.frame(t(combat_batch_corrected))
# pdf('TCGA_CMT after combat source.pdf', width = 6,height = 4)
# p <- PCA_plot(data = data2 , group = phenotype_AB$SOURCE , title = "TCGA+CMT after combat"  )
# print(p)
# dev.off()
# pdf('TCGA_CMT after combat subtype.pdf', width = 6,height = 4)
# p <- PCA_plot(data = data2 , group = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )
# print(p)
# dev.off()
# pdf('TCGA_CMT after combat source_subtype.pdf', width = 6,height = 4)
# p <- PCA_plot(data = data2 , group = phenotype_AB$SOURCE ,
#           label = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )
# print(p)
# dev.off()

# pdf('TCGA after combat.pdf', width = 6,height = 4)
# p <- PCA_plot2(data = t(tcga_corrected) , group = A_pheno$SUBTYPE , title = "TCGA after combat"  )
# print(p)
# dev.off()
#
# pdf('CMT after combat.pdf', width = 6,height = 4)
# p <- PCA_plot2(data = t(cmt_corrected) , group = B_pheno$SUBTYPE , title = "CMT after combat"  )
# print(p)
# dev.off()

# pdf(file="CBFpvcaFunction_combat_batch_corrected.pdf", width = 6,height = 4)
# CBFpvcaFunction(t(combat_batch_corrected), phenotypedata = phenotype_AB[,c("SOURCE","SUBTYPE")])
# dev.off()
