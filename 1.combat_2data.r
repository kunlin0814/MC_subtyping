library(MultiBaC)
#load packages 
library(edgeR)
library(limma)
library(MultiAssayExperiment)
library(sva)
library(pvca)
## 06/26/23
## we can just use human TCGA data and ignore dogs
#load("E:/My Drive/Josh_MC_Paper_data/ML_gene_set/orig_results/1.combat_2data.rdata")
source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
#'/Volumes/Research/GitHub/R_util/my_util.R')
#'
#'
base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)

A.rna  <- fread("TCGA_BRCA_log2TPM+1.csv", header = T)
colnames(A.rna)[1] <- 'GeneName'
dim(A.rna) # 11856   664 samples

B.rna  <- fread("cmt_all.csv", header = T)
### Mar-01 duplicated
#B.rna <- B.rna[,-1]
#setDF(B.rna)
dim(B.rna) # 13220   144
overlap_gene_list <- intersect(B.rna$GeneName,A.rna$GeneName)
B.rna <- B.rna[GeneName %in% overlap_gene_list, ]
A.rna <- A.rna[GeneName %in% overlap_gene_list, ]
setDF(B.rna)
setDF(A.rna)
#B.rna <- na.omit(B.rna)
row.names(B.rna) <- B.rna$GeneName
row.names(A.rna) <- A.rna$GeneName
A.rna <- A.rna[,-1]
B.rna <- B.rna[,-1]
abc = row.names(A.rna)== row.names(B.rna) #TRUE
AB.rna <- cbind(A.rna,B.rna)

phenotype_AB <- read.csv(file="All_dog_TCGA_subtype_meta.txt",header=T,sep="\t")
AB.rna <- AB.rna[, colnames(AB.rna)%in% phenotype_AB$PATIENT_ID]
#CMT has only simple cancer = 78 sam
B_pheno <-phenotype_AB[phenotype_AB$SOURCE=='cmt',] 
B.rna <- AB.rna[, colnames(AB.rna)%in% B_pheno$PATIENT_ID]

phenotype_AB <- phenotype_AB[ match(colnames(AB.rna), phenotype_AB$PATIENT_ID ),]

#phenotype_AB<- phenotype_AB[,-1]
#phenotype_AB$SOURCE <- as.factor(phenotype_AB$SOURCE)
#phenotype_AB$LN_status <- as.factor(phenotype_AB$LN_status)
#phenotype_AB$SUBTYPE <- as.factor(phenotype_AB$SUBTYPE)

A_pheno <- phenotype_AB[phenotype_AB$SOURCE== "tcga",]
B_pheno <- phenotype_AB[phenotype_AB$SOURCE== "cmt",]

pdf(file="pca_plot_source.pdf", width = 6,height = 4)
p <- PCA_plot(t(AB.rna),phenotype_AB$SOURCE , title = "TCGA+CMT before combat")
print(p)
dev.off()
pdf(file="pca_plot_sub.pdf", width = 6,height = 4)
p <- PCA_plot(t(AB.rna),phenotype_AB$SUBTYPE , title = "TCGA+CMT before combat")
print(p)
dev.off()
pdf(file="TCGA before combat.pdf", width = 6,height = 4)
p <- PCA_plot2(t(A.rna),A_pheno$SUBTYPE , title = "TCGA before combat")
print(p)
dev.off()
pdf(file="CMT before combat.pdf", width = 6,height = 4)
p <- PCA_plot2(t(B.rna),B_pheno$SUBTYPE , title = "CMT before combat")
print(p)
dev.off()

data <- data.frame(gene_id= rownames(AB.rna),AB.rna )

#First we need to put our Samples into our rownames - because the function
#uses rownames of our phenotype_AB table to match which samples belong to which column
#from the expression data.

rownames(phenotype_AB) <- phenotype_AB$PATIENT_ID
#Lets pass data[,-1] (without the gene_id column) to the CBFpcvaFunction. We will need
#transpose the data using t() to a wide format. We can also pass phenotype_AB[,-1], this
#is necessary because the PCVA function looks for grouping variables to assess contribution
#of variance so it is important not to pass any columns that contain unique identifiers.
pdf(file="CBFpvcaFunction.pdf", width = 6,height = 4)
p <- CBFpvcaFunction(t(data[,-1]), phenotypedata = phenotype_AB[,c("SOURCE","SUBTYPE")])
print(p)
dev.off()


# Here we build a model matrix with  covariates = phenotype_AB$CLASS

modcombat = model.matrix(~1+ phenotype_AB$SUBTYPE, data = phenotype_AB)


#### Correct SAMPLE
combat_batch_corrected = ComBat(dat = as.matrix(data[, -1]), batch = phenotype_AB$SOURCE,
                                mod = modcombat, par.prior = TRUE, prior.plots = TRUE)

rownames(combat_batch_corrected) <- data$gene_id
head(combat_batch_corrected)
pdf(file="CBFpvcaFunction_combat_batch_corrected.pdf", width = 6,height = 4)
p <- CBFpvcaFunction(t(combat_batch_corrected), phenotypedata = phenotype_AB[,c("SOURCE","SUBTYPE")])
print(p)
dev.off()


#Plot PCA again after combat 
data2  <- data.frame(t(combat_batch_corrected))

pdf('TCGA_CMT after combat source.pdf', width = 6,height = 4)
p <- PCA_plot(data = data2 , group = phenotype_AB$SOURCE , title = "TCGA+CMT after combat"  )
print(p)
dev.off()

pdf('TCGA_CMT after combat subtype.pdf', width = 6,height = 4)
p <- PCA_plot(data = data2 , group = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )
print(p)
dev.off()
pdf('TCGA_CMT after combat source_subtype.pdf', width = 6,height = 4)
p <- PCA_plot(data = data2 , group = phenotype_AB$SOURCE , 
          label = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )
print(p)
dev.off()

tcga_corrected <- data.frame(combat_batch_corrected[,colnames(combat_batch_corrected)%in%A_pheno$PATIENT_ID ])
dim(tcga_corrected) # 11856   664

cmt_corrected <- data.frame(combat_batch_corrected[,colnames(combat_batch_corrected)%in%B_pheno$PATIENT_ID ])
dim(cmt_corrected) #11856    78

pdf('TCGA after combat.pdf', width = 6,height = 4)
p <- PCA_plot2(data = t(tcga_corrected) , group = A_pheno$SUBTYPE , title = "TCGA after combat"  )
print(p)
dev.off()

pdf('CMT after combat.pdf', width = 6,height = 4)
p <- PCA_plot2(data = t(cmt_corrected) , group = B_pheno$SUBTYPE , title = "CMT after combat"  )
print(p)
dev.off()


write.csv(A_pheno,"phenotype_all_tcga.csv" )
write.csv(B_pheno,"phenotype_all_cmt.csv" )

write.csv(combat_batch_corrected,"all_tcga_cmt_combat_corrected.csv")
write.csv(tcga_corrected,"all_tcga_combat_corrected.csv")
write.csv(cmt_corrected,"all_cmt_combat_corrected.csv")

save.image("all_tcga_combat_2data.rdata")

