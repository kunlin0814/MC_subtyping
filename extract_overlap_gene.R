library(ggplot2)
library(data.table)

source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')

base <- "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"  
  #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"

# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
main <- c("Basal","LumA","LumB","Her2")
combine <- combn(main,2)
results_base <-  "G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set/Gene_overlap_list/LumA"
  #"E:/My Drive/Josh_MC_Paper_data/ML_gene_set/Gene_overlap_list/Basal"
combine[,6] <- c("LumA","Basal")

total_25_list <- list()
total_50_list <- list()
# Basal 1:3, LumA 4:6
for ( i in 4:6){
  comparison <- combine[,i]
  comparison_header  <- paste(comparison, collapse = 'vs')
  prev_results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
  setwd(prev_results_base)
  freq25 <- fread(paste(comparison_header,"_gene_list_freq25.txt", sep=""))
  freq50 <- fread(paste(comparison_header,"_gene_list_freq50.txt", sep=""))
  freq25_gene = freq25$gene
  freq50_gene = freq50$gene
  total_25_list[[comparison_header]] <- freq25_gene
  total_50_list[[comparison_header]] <- freq50_gene
}
union_25_gene <- c()
union_50_gene <- c()
for (each in names(total_25_list)){

  union_25_gene <- c(union_25_gene,total_25_list[[each]])
  union_50_gene <- c(union_50_gene,total_50_list[[each]])
}

union_25_gene <- unique(union_25_gene)
union_50_gene <- unique(union_50_gene)


overlap_25_gene <- intersect(total_25_list[[1]],total_25_list[[2]])
overlap_25_gene <- intersect(overlap_25_gene,total_25_list[[3]])
  
overlap_50_gene <- intersect(total_50_list[[1]],total_50_list[[2]])
overlap_50_gene <- intersect(overlap_50_gene,total_50_list[[3]])


overlap_25_table <- data.table(freq25=overlap_25_gene)
overlap_50_table <- data.table(freq50=overlap_50_gene)
union_25_table <- data.table(freq25=union_25_gene)
union_50_table <- data.table(freq50=union_50_gene)

fwrite(overlap_25_table, file=paste(results_base,"LumA_overlap_25.txt",sep="/"),
       col.names = T)
fwrite(overlap_50_table, file=paste(results_base,"LumA_overlap_50.txt",sep="/"),
       col.names = T)
fwrite(union_25_table, file=paste(results_base,"LumA_union_25.txt",sep="/"),
       col.names = T)
fwrite(union_50_table, file=paste(results_base,"LumA_union_50.txt",sep="/"),
       col.names = T)


# basalvslumA_50 <- total_50_list$BasalvsLumA
# basalvslumA_25 <- total_25_list$BasalvsLumA
  
  