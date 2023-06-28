library(edgeR)
library(limma)
library(tmaptools)
library(tidyverse)
library(reshape2)
library(biomaRt)
library(caret)

#load("E:/My Drive/Josh_MC_Paper_data/ML_gene_set/orig_results/2.tcga_DEG_subtype.rdata")

source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
#'/Volumes/Research/GitHub/R_util/my_util.R')
#'
# package_location <- as.character(args[1])
# new_data_file <- as.character(args[2])
# out_file_name <- as.character(args[3])

base <- #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
  "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"
#"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)
# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
comparison <- c("Basal","LumB")
comparison_header  <- paste(comparison, collapse = 'vs')

results_base <- paste(base,'Step2DEG',comparison_header,sep="/")
dir.create(results_base)
tcga_data <- fread("all_tcga_combat_corrected.csv",header = T)
setDF(tcga_data)
row.names(tcga_data) <- tcga_data$V1
tcga_data <- tcga_data[,-1]

pheno_tcga <- read.csv("phenotype_all_tcga.csv",header = T)
row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga$X <- NULL
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison,]
tcga_data <- tcga_data[,colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

length(colnames(tcga_data))==length(pheno_tcga$PATIENT_ID)

tcga_data <- data.frame(t(tcga_data))

tcga_data <- tcga_data %>%
  mutate(across(.cols = 1:length(.), .fns = as.numeric))
class(tcga_data$A1BG) # "numeric"


aaa = pheno_tcga$PATIENT_ID == row.names(tcga_data)#TRUE 


#tcga_data_dge <- DGEList(t(tcga_data))
#Error: Negative counts not allowed


tcga_data <- data.frame(t(tcga_data))

# Identify genes with negative counts
genes_with_negative_counts <- rownames(tcga_data)[apply(tcga_data, 1, function(x) any(x < 0))]

# Exclude genes with negative counts from the analysis
tcga_data_filtered <- tcga_data[!rownames(tcga_data) %in% genes_with_negative_counts, ]

tcga_data_filtered <- data.frame( SUBTYPE= pheno_tcga$SUBTYPE, t(tcga_data_filtered))
dim(tcga_data_filtered)
#  664 9384
#tcga_data_filtered <- data.frame( SUBTYPE= pheno_tcga$SUBTYPE, tcga_data_filtered[,-1])


set.seed(1)  # Set seed for reproducibility
train_indices <- createDataPartition(tcga_data_filtered$SUBTYPE,
                                     p = 0.8,
                                     list = FALSE,
                                     times = 1)
train_data <- tcga_data_filtered[train_indices, ]
test_data <- tcga_data_filtered[-train_indices, ]

pheno_train_data <- factor(train_data$SUBTYPE)

tcga_data_filtered_train <- data.frame(t(train_data[,-1]))

tcga_data_dge <- DGEList(tcga_data_filtered_train)

class(tcga_data_dge)
head(tcga_data_dge[["counts"]])

#rownames(tcga_data_dge[["counts"]]) <- colnames(tcga_data)

#We can check our rownames were sucsessfully assigned by viewing our data.
head(rownames(tcga_data_dge[["counts"]]))
head(tcga_data_dge[["counts"]])

dim(tcga_data_dge[["counts"]])
# 9559  138
head(tcga_data_dge[["samples"]])

#We can make a vector of factors from our phenotype table that contains sample group information
samp_groups <- pheno_train_data
#samp_groups

#Lets now reassign the group values in the samples data.frame and check the result.
tcga_data_dge[["samples"]]$group <- samp_groups
tcga_data_dge[["samples"]]


#### Testing for Differential Expression ####

##### Limma ####### 

#Set up our basic design matrix with no intercept
#design <- model.matrix(~group, data = tcga_data_dge[["samples"]])

design <- model.matrix(~ 0 + group, data = tcga_data_dge[["samples"]])
colnames(design) <- comparison
#View(design)
pdf(file=paste(results_base,
               paste(comparison_header,"voom_mean_variance_plot.pdf",sep="")
               ,sep="/")
               ,height = 4.5,width = 6)
voom_data <- voom(tcga_data_dge, design, plot = TRUE)
dev.off()
voom_fit <- lmFit(object = voom_data, design = design)

## use the following 

comparison_command <- paste(comparison, collapse = '-')
prestr <- "makeContrasts("
poststr <- ",levels=design)"
commandstr=paste(prestr,comparison_header,"=",comparison_command,poststr,sep="")
# commandstr
cont.matrix <- eval(parse(text=commandstr))
#cont.matrix <- makeContrasts(LumAvsLumB = LumA - LumB, levels = design)

voom_fit <- contrasts.fit(fit = voom_fit, contrasts = cont.matrix)

voom_fit <- eBayes(voom_fit)

#We can run a quick diagnostic plot to again plot our mean variances after fitting our linear model and estimating.
#The blue line here represents our residual standard deviation estimated by eBayes.
pdf(file=paste(results_base,
               paste(comparison_header,"voom_fit_plot.pdf",sep=""),
               sep="/") ,height = 4.5,width = 6)
plotSA(voom_fit)
dev.off()

voom_tt <- topTable(voom_fit, p.value = 1, number = Inf)
voom_tt0.05 <- voom_tt[voom_tt$adj.P.Val <=0.05,]

#7624 genes 

pca_data <- tcga_data_filtered_train[rownames(tcga_data_filtered_train)%in% rownames(voom_tt0.05),]
dim(pca_data)

pdf(file=paste(results_base,
               paste(comparison_header,"PCA_plot.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(t(pca_data), pheno_train_data, 
              title = paste("DEGs of train tcga data ", paste("n=",nrow(voom_tt0.05))))
print(p)
dev.off()
#tcga_data_df <-  data.frame(tcga_data_dge[["counts"]])

write.csv(voom_tt0.05,paste(results_base,paste(comparison_header,"DEG_train_tgca_subtype.csv",sep = "_")
                            ,sep="/"))

save.image(paste(results_base,paste(comparison_header,"_tcga_DEG_subtype.rdata",sep="")
                 ,sep="/"))

