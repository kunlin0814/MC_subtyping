library(ggplot2)
library(data.table)

source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')

base <- #"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"  
"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"

# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
main <- c("Basal","LumA","LumB","Her2")
combine <- combn(main,2)
results_base <-  #"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set/Gene_overlap_list/LumA"
  "E:/My Drive/Josh_MC_Paper_data/ML_gene_set/Step2DEG"
combine[,6] <- c("LumA","Basal")

total_overlap_list <- list()
total_union_list <- list()
# Basal 1:3, LumA 4:6
for ( i in 4:6){
  comparison <- combine[,i]
  comparison_header  <- paste(comparison, collapse = 'vs')
  prev_results_base <- paste(base,'Step2DEG',comparison_header,sep="/")
  setwd(prev_results_base)
  deg_list <- fread(paste(comparison_header,"_DEG_train_tgca_subtype.csv", sep=""))

  genelist = deg_list$V1
  total_union_list[[comparison_header]] <- genelist
} 
union_gene <- c()
for (each in names(total_union_list)){
  
  union_gene <- c(union_gene,total_union_list[[each]])
 
}

union_gene <- unique(union_gene)

overlap_gene <- intersect(total_union_list[[1]],total_union_list[[2]])
overlap_gene <- intersect(overlap_gene,total_union_list[[3]])

