library(randomForest)
library(xgboost)
library(DMwR)
library(grid)
library(Boruta)
library(vars)
library(dplyr)
library(varSelRF)
library(caret)
library(e1071)
library(ranger)
library(sva)
library(smotefamily)
library(mlr3measures)
library(kernlab)

source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
base <- #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
  "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"

cat <- "Basal"

# prev_results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
# load(paste(prev_results_base,'model_tcga_subtype.rdata',sep="/"))


#each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
#minor_group <- each_pheno_num[each_pheno_num$n ==min(each_pheno_num$n),]$SUBTYPE

data_test<- read.csv(paste(base,"all_cmt_combat_corrected.csv",sep="/"),header = T, row.names = 1)
pheno_test <- read.csv(paste(base,"phenotype_cmt.csv",sep="/"),header = T)
pheno_test[pheno_test$SUBTYPE!="basal",]$SUBTYPE <- "LumA"
pheno_test[pheno_test$SUBTYPE=="basal",]$SUBTYPE <- "Basal"

pheno_tcga <- read.csv(paste(base,"phenotype_all_tcga.csv",sep="/"),header = T)
row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
pheno_tcga$X <- NULL
#pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison,]
pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE!="Normal",]
tcga_data <- tcga_data[,colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

pheno_train <- pheno_tcga
dim(data_test)
#  11856    78



############
results_base <- paste(base,'Gene_overlap_list',cat,sep="/")
dir.create(results_base,recursive = TRUE)

data_test <- data.frame(data_test)

overlap_freq25 <- fread(paste(base,'Gene_overlap_list',cat,paste(cat,'_overlap_25.txt',sep=""),sep="/"))
overlap_freq50 <- fread(paste(base,'Gene_overlap_list',cat, paste(cat,'_overlap_50.txt',sep=""),sep="/"))
union_freq25 <- fread(paste(base,'Gene_overlap_list',cat,paste(cat,'_union_25.txt',sep=""),sep="/"))
union_freq50 <- fread(paste(base,'Gene_overlap_list',cat, paste(cat,'_union_50.txt',sep=""),sep="/"))

union_freq25$freq25
union_freq50$freq50


test_freq25_overlap <- data_test[ row.names(data_test) %in% overlap_freq25$freq25,]
test_freq50_overlap <- data_test[ row.names(data_test) %in% overlap_freq50$freq50,]
test_freq25_union <- data_test[ row.names(data_test) %in% union_freq25$freq25,]
test_freq50_union <- data_test[ row.names(data_test) %in% union_freq50$freq50,]



#PCA

pdf(file=paste(results_base,
               paste(cat,"_overlap_Freq_25.pdf",sep=""),
               sep="/") 
    ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq25_overlap),pheno_test$SUBTYPE, scale_status = F,
               title =paste(paste(cat, " Freq=25, n= ", nrow(test_freq25_overlap),sep=""),sep=""))
print(p)
dev.off()

pdf(file=paste(results_base,
               paste(cat,"_overlap_Freq_50.pdf",sep=""),
               sep="/") ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq50_overlap),pheno_test$SUBTYPE,scale_status = F,
              title =paste(paste(cat, " Freq=25, n= ", nrow(test_freq50_overlap),sep=""),sep=""))

print(p)
dev.off()

pdf(file=paste(results_base,
               paste(cat,"_union_Freq_25.pdf",sep=""),
               sep="/") 
    ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq25_union),pheno_test$SUBTYPE, scale_status = F,
              title =paste(paste(cat, " Freq=25, n= ", nrow(test_freq25_union),sep=""),sep=""))
print(p)
dev.off()

pdf(file=paste(results_base,
               paste(cat,"_union_Freq_50.pdf",sep=""),
               sep="/") ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq50_union),pheno_test$SUBTYPE,scale_status = F,
              title =paste(paste(cat, " Freq=25, n= ", nrow(test_freq50_union),sep=""),sep=""))

print(p)
dev.off()

