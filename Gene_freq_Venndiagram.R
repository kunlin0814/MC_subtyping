# ---- New allele identification (Known allele distribution) ----
library(ggplot2)
library(data.table)
library(VennDiagram)
library(futile.logger)


source(
  'C:/Users/abc73/Documents/GitHub/MC_subtyping/MC_subtyping_module.R')
source(
  'C:/Users/abc73/Documents/GitHub/R_util/my_util.R')

base <-  "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"

# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
main <- c("Basal","LumA","LumB","Her2")
combine <- combn(main,2)


comparison <- c("Basal","Her2")

comparison_header  <- paste(comparison, collapse = 'vs')
prev_results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
setwd(prev_results_base)

freq25 <- fread(paste(comparison_header,"_gene_list_freq25.txt", sep=""))
freq40 <- fread(paste(comparison_header,"_gene_list_freq40.txt", sep=""))
freq50 <- fread(paste(comparison_header,"_gene_list_freq50.txt", sep=""))


set1 = freq25$gene
set2 = freq40$gene
set3 = freq50$gene

result = list(set1,set2,set3)

pdf(file=paste(comparison_header,'_gene_freq_overlap.pdf',sep=""),height = 4.5,width =6)
temp <- venn.diagram(
  x = result,
  category.names = c("Freq 25" , "Freq 40" , "Freq 50"),
  filename = NULL,
  #paste(result_base,paste('AlleleValidation_venn_diagramm_',date,'.png',sep="")
  #,sep='/'),
  cex = 2, # font size for the number
  output=TRUE,
  cat.cex = 1.5) # font size for category
  #cat.just=list(c(0,1) , c(0,1) , c(0,1)))
# height = 2500, width = 3000, resolution =
#   500, imagetype = "png", units = "px", compression =
#   "lzw")
# temp <- venn.diagram(list(B = 1:1800, A = 1571:2020),
#                      fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,
#                      lty =2, fontfamily =3, filename = NULL)

grid.draw(temp)
dev.off()  