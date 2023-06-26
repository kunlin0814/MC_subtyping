library(MultiBaC)
#load packages 
library(edgeR)
library(limma)
library(MultiAssayExperiment)
library(sva)

## 06/26/23
## we can just use human TCGA data and ignore dogs
load("1.combat_2data.rdata")
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
#'/Volumes/Research/GitHub/R_util/my_util.R')
#'
#'
base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)

A.rna  <- read.csv("TCGA_BRCA_rawTPM.csv", row.names = 1, header = T)
dim(A.rna) # 11856   664 samples

phenotype_human <- read.csv(file="All_TCGA_subtype.txt",header=T)
A.rna <- A.rna[, colnames(A.rna)%in% phenotype_human$PATIENT_ID]

phenotype_human <- phenotype_human[ match(colnames(A.rna), phenotype_human$PATIENT_ID ),]

phenotype_human<- phenotype_human[,-1]
#phenotype_human$SOURCE <- as.factor(phenotype_human$SOURCE)
#phenotype_human$LN_status <- as.factor(phenotype_human$LN_status)
#phenotype_human$SUBTYPE <- as.factor(phenotype_human$SUBTYPE)



A_pheno <- phenotype_human

PCA_plot2(t(A.rna),A_pheno$SUBTYPE , title = "TCGA before combat")
# PCA_plot2(t(B.rna),B_pheno$SUBTYPE , title = "CMT before combat")


data <- data.frame(gene_id= rownames(A.rna),A.rna )

#First we need to put our Samples into our rownames - because the function
#uses rownames of our phenotype_human table to match which samples belong to which column
#from the expression data.

rownames(phenotype_human) <- phenotype_human$PATIENT_ID
#Lets pass data[,-1] (without the gene_id column) to the CBFpcvaFunction. We will need
#transpose the data using t() to a wide format. We can also pass phenotype_human[,-1], this
#is necessary because the PCVA function looks for grouping variables to assess contribution
#of variance so it is important not to pass any columns that contain unique identifiers.

CBFpvcaFunction(t(data[,-1]), phenotypedata = phenotype_human[,c(2,3,4)])



# Here we build a model matrix with  covariates = phenotype_human$CLASS

modcombat = model.matrix(~1+ phenotype_human$SUBTYPE, data = phenotype_human)


#### Correct SAMPLE
combat_batch_corrected = ComBat(dat = as.matrix(data[, -1]), batch = phenotype_human$SOURCE,
                                mod = modcombat, par.prior = TRUE, prior.plots = TRUE)

rownames(combat_batch_corrected) <- data$gene_id
head(combat_batch_corrected)

CBFpvcaFunction(t(combat_batch_corrected), phenotypedata = phenotype_human[,c(2,3,4)])

#Plot PCA again after combat 
data2  <- data.frame(t(combat_batch_corrected))

PCA_plot(data = data2 , group = phenotype_human$SOURCE , title = "TCGA after combat"  )

PCA_plot(data = data2 , group = phenotype_human$SUBTYPE , title = "TCGA after combat"  )

PCA_plot(data = data2 , group = phenotype_human$SOURCE , 
         label = phenotype_human$SUBTYPE , title = "TCGA after combat"  )


tcga_corrected <- data.frame(combat_batch_corrected[,colnames(combat_batch_corrected)%in%A_pheno$PATIENT_ID ])
dim(tcga_corrected) # 11856   664

cmt_corrected <- data.frame(combat_batch_corrected[,colnames(combat_batch_corrected)%in%B_pheno$PATIENT_ID ])
dim(cmt_corrected) #11856    78

PCA_plot2(data = t(tcga_corrected) , group = A_pheno$SUBTYPE , title = "TCGA after combat"  )
PCA_plot2(data = t(cmt_corrected) , group = B_pheno$SUBTYPE , title = "CMT after combat"  )

write.csv(A_pheno,"phenotype_tcga.csv" )
write.csv(B_pheno,"phenotype_cmt.csv" )

write.csv(combat_batch_corrected,"AB_combat_corrected.csv")
write.csv(tcga_corrected,"tcga_combat_corrected.csv")
write.csv(cmt_corrected,"cmt_combat_corrected.csv")

save.image("combat_2data.rdata")

