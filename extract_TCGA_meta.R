library(edgeR)
library(limma)
library(tmaptools)
library(tidyverse)
library(reshape2)
library(biomaRt)

source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')
#'/Volumes/Research/GitHub/R_util/my_util.R')

all_data <- fread('TCGA_BRCA_rawTPM.csv')
all_meta <- all_data[1,]
all_meta <- all_meta[,-1]
sample_list <- c(colnames(all_meta))
all_meta <- data.frame(t(all_meta))
colnames(all_meta) <- "SUBTYPE"

all_meta$PATIENT_ID <- sample_list
setDT(all_meta)

fwrite(all_meta, file ="All_TCGA_subtype.txt", sep ="\t", col.names = T, row.names = F, eol = "\n")




