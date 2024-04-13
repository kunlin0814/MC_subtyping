# This script identifies differentially expressed genes between subtypes using the limma package.
# It requires two input files:
# 1. The combat-corrected human TPM file.
# 2. The subtype metadata for human, including sample names and subtype information.

# The goal of this script is to create a gene list that shows differential expression with a p-value less than 0.05.
# By comparing subtypes, we can identify genes that are significantly differentially expressed between them.

# Note: The combat correction step is assumed to have been performed prior to running this script.
# Combat correction helps to address batch effects, improving the accuracy of differential expression analysis.

# Please ensure that the required input files are available and the necessary R packages are installed before running this script.


source('C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')

base <- "D:/§Úªº¶³ºÝµwºÐ/Josh_MC_Paper_data/ML_gene_set"
setwd(base)

## TCGA breast cancer subtype
# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
comparison <- c("Basal", "LumA")
comparison_header  <- paste(comparison, collapse = 'vs')
results_base <- paste(base, 'Step2DEG', comparison_header, sep = "/")
dir.create(results_base)

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

tcga_data <- data.frame(t(tcga_data))

tcga_data <- tcga_data %>%
  mutate(across(.cols = 1:length(.), .fns = as.numeric))

tcga_data <- data.frame(t(tcga_data))

### If there are negative count in the log2TPM, DGEList will create error, and we need to remove gene with negative counts
# Identify genes with negative counts
genes_with_negative_counts <-
  rownames(tcga_data)[apply(tcga_data, 1, function(x)
    any(x < 0))]

# Exclude genes with negative counts from the analysis
tcga_data_filtered <-
  tcga_data[!rownames(tcga_data) %in% genes_with_negative_counts,]

tcga_data_filtered <-
  data.frame(SUBTYPE = pheno_tcga$SUBTYPE, t(tcga_data_filtered))

set.seed(1)  # Set seed for reproducibility
train_indices <- createDataPartition(
  tcga_data_filtered$SUBTYPE,
  p = 0.8,
  list = FALSE,
  times = 1
)
train_data <- tcga_data_filtered[train_indices,]
test_data <- tcga_data_filtered[-train_indices,]



tcga_data_filtered_train <- data.frame(t(train_data[, -1]))

tcga_data_dge <- DGEList(tcga_data_filtered_train)

#We can make a vector of factors from our phenotype table that contains sample group information
samp_groups <- factor(train_data$SUBTYPE)

#Lets now reassign the group values in the samples data.frame and check the result.
tcga_data_dge[["samples"]]$group <- samp_groups


#### Testing for Differential Expression ####

##### Limma #######

#Set up our basic design matrix with 1 as intercept

design <- model.matrix( ~ group, data = tcga_data_dge[["samples"]])
colnames(design) <- comparison


pdf(
  file = paste(
    results_base,
    paste(comparison_header, "voom_mean_variance_plot.pdf", sep = "")
    ,
    sep = "/"
  ),
  height = 4.5,
  width = 6
)
## voom : Transform count data to log2-counts per million (logCPM) for Linear modelling
voom_data <- voom(tcga_data_dge, design, plot = TRUE)
# dev.off()
voom_fit <- lmFit(object = voom_data, design = design)

## use variable to execute the following command in R
comparison_command <- paste(comparison, collapse = '-')
prestr <- "makeContrasts("
poststr <- ",levels=design)"
commandstr <-
  paste(prestr,
        comparison_header,
        "=",
        comparison_command,
        poststr,
        sep = "")
cont.matrix <- eval(parse(text = commandstr))
voom_fit <- contrasts.fit(fit = voom_fit, contrasts = cont.matrix)
voom_fit <- eBayes(voom_fit)

## We can run a quick diagnostic plot toplot our mean variances after fitting our linear model and estimating.
## The blue line here represents our residual standard deviation estimated by eBayes.
pdf(
  file = paste(
    results_base,
    paste(comparison_header, "voom_fit_plot.pdf", sep = ""),
    sep = "/"
  ) ,
  height = 4.5,
  width = 6
)
plotSA(voom_fit)
dev.off()

voom_tt <- topTable(voom_fit, p.value = 1, number = Inf)
## select only differential expressed genes with p<=0.05
voom_tt0.05 <- voom_tt[voom_tt$adj.P.Val <= 0.05, ]
pca_data <-
  tcga_data_filtered_train[rownames(tcga_data_filtered_train) %in% rownames(voom_tt0.05), ]

pdf(
  file = paste(
    results_base,
    paste(comparison_header, "PCA_plot.pdf", sep = ""),
    sep = "/"
  )
  ,
  height = 4.5,
  width = 6
)
p <- PCA_plot(t(pca_data),
              pheno_train_data,
              title = paste("DEGs of train tcga data ", paste("n=", nrow(voom_tt0.05))))
print(p)
dev.off()

write.csv(voom_tt0.05, paste(
  results_base,
  paste(comparison_header, "DEG_train_tgca_subtype.csv", sep = "_")
  ,
  sep = "/"
))

# save.image(paste(results_base,paste(comparison_header,"_tcga_DEG_subtype.rdata",sep="")
#                  ,sep="/"))
