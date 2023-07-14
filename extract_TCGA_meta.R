library(edgeR)
library(limma)
library(tmaptools)
library(tidyverse)
library(reshape2)
library(biomaRt)

source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
#'/Volumes/Research/GitHub/R_util/my_util.R')
base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)
all_TCGA_data <- fread('TCGA_BRCA_rawTPM.csv')
all_dog_meta <- fread("phenotype_cmt.csv")
colnames(all_dog_meta) <- c("Sample_name","PATIENT_ID","SOURCE","LN_status","SUBTYPE")
all_dog_meta <- all_dog_meta[-1,]
all_dog_meta <- unique(all_dog_meta[,.(PATIENT_ID,SUBTYPE)])
all_dog_meta$SOURCE <- 'cmt'
all_human_meta <- all_TCGA_data[1,]
all_human_meta <- all_human_meta[,-1]
sample_list <- c(colnames(all_human_meta))
all_human_meta <- data.frame(t(all_human_meta))
colnames(all_human_meta) <- "SUBTYPE"
all_human_meta$SOURCE <- 'tcga'
all_human_meta$PATIENT_ID <- sample_list
setDT(all_human_meta)
human_dog_meta <- rbind(all_human_meta, all_dog_meta)
fwrite(human_dog_meta, file ="All_dog_TCGA_subtype_meta.txt", sep ="\t", col.names = T, row.names = F, eol = "\n")




