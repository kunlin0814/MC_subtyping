
library(MultiBaC)
#load packages 
library(edgeR)
library(limma)
library(MultiAssayExperiment)
library(sva)


A.rna  <- read.csv("tcga_all.csv", row.names = 1, header = T)
dim(A.rna) # 11856   664

B.rna  <- read.csv("cmt_all.csv", header = T)
dim(B.rna) # 13220   144

B.rna <- B.rna[ B.rna$GeneName %in% row.names(A.rna), ]
B.rna <- B.rna[ match( row.names(A.rna),B.rna$GeneName ),]
row.names(B.rna) <- B.rna$GeneName
B.rna <- B.rna[,-1]
abc = row.names(A.rna)== row.names(B.rna) #TRUE

AB.rna <- cbind(A.rna,B.rna)
phenotype_AB <- read.csv(file="phenotype_tcga_cmt.csv",header=T)
AB.rna <- AB.rna[, colnames(AB.rna)%in% phenotype_AB$PATIENT_ID]
#CMT has only simple cancer = 78 sam
B.rna <- AB.rna[, colnames(AB.rna)%in% B_pheno$PATIENT_ID]

phenotype_AB <- phenotype_AB[ match(colnames(AB.rna), phenotype_AB$PATIENT_ID ),]

phenotype_AB<- phenotype_AB[,-1]
#phenotype_AB$SOURCE <- as.factor(phenotype_AB$SOURCE)
#phenotype_AB$LN_status <- as.factor(phenotype_AB$LN_status)
#phenotype_AB$SUBTYPE <- as.factor(phenotype_AB$SUBTYPE)



A_pheno <- phenotype_AB[phenotype_AB$SOURCE== "tcga",]
B_pheno <- phenotype_AB[phenotype_AB$SOURCE== "cmt",]

PCA_plot(t(AB.rna),phenotype_AB$SOURCE , title = "TCGA+CMT before combat")

PCA_plot(t(AB.rna),phenotype_AB$SUBTYPE , title = "TCGA+CMT before combat")

PCA_plot2(t(A.rna),A_pheno$SUBTYPE , title = "TCGA before combat")
PCA_plot2(t(B.rna),B_pheno$SUBTYPE , title = "CMT before combat")


data <- data.frame(gene_id= rownames(AB.rna),AB.rna )

#First we need to put our Samples into our rownames - because the function
#uses rownames of our phenotype_AB table to match which samples belong to which column
#from the expression data.

rownames(phenotype_AB) <- phenotype_AB$PATIENT_ID
#Lets pass data[,-1] (without the gene_id column) to the CBFpcvaFunction. We will need
#transpose the data using t() to a wide format. We can also pass phenotype_AB[,-1], this
#is necessary because the PCVA function looks for grouping variables to assess contribution
#of variance so it is important not to pass any columns that contain unique identifiers.

CBFpvcaFunction(t(data[,-1]), phenotypedata = phenotype_AB[,c(2,3,4)])



# Here we build a model matrix with  covariates = phenotype_AB$CLASS

modcombat = model.matrix(~1+ phenotype_AB$SUBTYPE, data = phenotype_AB)


#### Correct SAMPLE
combat_batch_corrected = ComBat(dat = as.matrix(data[, -1]), batch = phenotype_AB$SOURCE,
                                mod = modcombat, par.prior = TRUE, prior.plots = TRUE)

rownames(combat_batch_corrected) <- data$gene_id
head(combat_batch_corrected)

CBFpvcaFunction(t(combat_batch_corrected), phenotypedata = phenotype_AB[,c(2,3,4)])

#Plot PCA again after combat 
data2  <- data.frame(t(combat_batch_corrected))

PCA_plot(data = data2 , group = phenotype_AB$SOURCE , title = "TCGA+CMT after combat"  )

PCA_plot(data = data2 , group = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )

PCA_plot(data = data2 , group = phenotype_AB$SOURCE , 
          label = phenotype_AB$SUBTYPE , title = "TCGA+CMT after combat"  )


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

