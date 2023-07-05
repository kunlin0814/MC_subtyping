## The whole pipeline involved 5 steps, but this script included the last four steps in the whole pipelie
# Step1 : Use Combat to perform batch effect correction on the expression data  
# Step2 : Identify Differential Gene expression 
# Step3 : Feature selection in the random forest model with Boruta method and repeat 50 times 
# Step4 : Create a random forest model and evaluate the performance in a training set
# Step5 : Use the feature genes identified in Step3 to evaluate the model performance

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

base <- "E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
#"G:/MAC_Research_Data/Josh_MC_Paper_data/ML_gene_set"
#"E:/My Drive/Josh_MC_Paper_data/ML_gene_set"
setwd(base)
# "LumA"   "LumB"   "Basal"  "Her2"   "Normal"
comparison <- c("LumB","Her2")
comparison_header  <- paste(comparison, collapse = 'vs')

#### DEG
## Step2
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

#We can make a vector of factors from our phenotype table that contains sample group information
samp_groups <- pheno_train_data

#Lets now reassign the group values in the samples data.frame and check the result.
tcga_data_dge[["samples"]]$group <- samp_groups


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

## use eval to execuate the following with variables  
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
p <- PCA_plot(t(pca_data), pheno_train_data, scale_status = T,
              title = paste("DEGs of train tcga data ", paste("n=",nrow(voom_tt0.05))))
print(p)
dev.off()
#tcga_data_df <-  data.frame(tcga_data_dge[["counts"]])

write.csv(voom_tt0.05,paste(results_base,paste(comparison_header,"DEG_train_tgca_subtype.csv",sep = "_")
                            ,sep="/"))

save.image(paste(results_base,paste(comparison_header,"_tcga_DEG_subtype.rdata",sep="")
                 ,sep="/"))

##### Feature selection
## Step3
prev_results_base <- paste(base,'Step2DEG',comparison_header,sep="/")
load(paste(prev_results_base,paste(comparison_header,"_tcga_DEG_subtype.rdata",sep=""),sep='/'))

# We will use both train and test data 

# comparison_header  <- paste(comparison, collapse = 'vs')
results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)


# tcga_data <- fread("all_tcga_combat_corrected.csv",header = T)
# setDF(tcga_data)
# row.names(tcga_data) <- tcga_data$V1
# tcga_data <- tcga_data[,-1]
# 
# pheno_tcga <- read.csv("phenotype_all_tcga.csv",header = T)
# row.names(pheno_tcga) <- pheno_tcga$PATIENT_ID
# pheno_tcga$X <- NULL
# pheno_tcga <- pheno_tcga[pheno_tcga$SUBTYPE %in% comparison,]
# tcga_data <- tcga_data[,colnames(tcga_data) %in% pheno_tcga$PATIENT_ID]

data <- data.frame(t(tcga_data[rownames(tcga_data)%in% rownames(voom_tt0.05),]))
dim(data)

data$SUBTYPE <- factor(pheno_tcga$SUBTYPE)
table(data$SUBTYPE )

#basal nonbasal 
# 171      493 

# Perform Boruta feature selection with 50 runs

boruta_results <- vector("list", 50)
#varsel_results <- vector("list", 50)

for (i in 1:50) {
  set.seed(i)  # Set seed for reproducibility
  
  train_indices <- createDataPartition(data$SUBTYPE,
                                       p = 0.8,
                                       list = FALSE,
                                       times = 1)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Boruta
  boruta_genes <- Boruta(SUBTYPE ~ ., train_data)
  boruta_fixed <- TentativeRoughFix(boruta_genes)
  boruta_final <- getSelectedAttributes(boruta_fixed, withTentative = FALSE)
  print(paste("boruta_round:", i))
  #print(boruta_final)
  boruta_results[[i]] <- boruta_final
  
  # VarselRF
  #varsel_genes <- varSelRF(train_data[, -ncol(train_data)], as.factor(train_data$SUBTYPE))
  #varsel_genes <- varsel_genes[["selected.vars"]]
  #print(paste("varsel_round:", i))
  #print(varsel_genes)
  #varsel_results[[i]] <- varsel_genes
}

# Create frequency tables for Boruta and VarselRF results
gene_list <- unlist(boruta_results)
boruta_freq <- table(unlist(boruta_results))
#varsel_freq <- table(unlist(varsel_results))

# Convert frequency tables to data frames
boruta_df <- data.frame(feature = names(boruta_freq), frequency = as.numeric(boruta_freq))
#varsel_df <- data.frame(feature = names(varsel_freq), frequency = as.numeric(varsel_freq))

# SOrt 
boruta_df <- boruta_df[order(boruta_df$frequency, decreasing = TRUE), ]
# View(boruta_df)
#varsel_df <- varsel_df[order(varsel_df$frequency, decreasing = TRUE), ]
#View(varsel_df)

# Apply cutoffs to filter features
boruta_filtered25 <- boruta_df[boruta_df$frequency >= 25, , drop = FALSE]
#varsel_filtered5 <- varsel_df[varsel_df$frequency >= 5, , drop = FALSE]

boruta_filtered40<- boruta_df[boruta_df$frequency >= 40, , drop = FALSE]
#varsel_filtered10 <- varsel_df[varsel_df$frequency >= 10, , drop = FALSE]

boruta_filtered50 <- boruta_df[boruta_df$frequency >= 50, , drop = FALSE]
#varsel_filtered15 <- varsel_df[varsel_df$frequency >= 15, , drop = FALSE]



# Extract filtered features
boruta_features25 <- boruta_filtered25$feature
#varsel_features5 <- varsel_filtered5$feature
boruta_features40 <- boruta_filtered40$feature
boruta_features50 <- boruta_filtered50$feature


tcga_data_t <- data.frame(t(tcga_data))

#PCA
pca_data_freq25 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features25]
pca_data_freq40 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features40]
pca_data_freq50 <- tcga_data_t[,colnames(tcga_data_t)%in% boruta_features50]

pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_25.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq25,pheno_tcga$SUBTYPE,scale_status = T, 
              title =paste("Selected genes-tcga_subtype Freq=25"," n=",length(boruta_features25),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_40.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq40,pheno_tcga$SUBTYPE, scale_status = T,
              title =paste("Selected genes-tcga_subtype Freq=40"," n=",length(boruta_features40),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_50.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(pca_data_freq50,pheno_tcga$SUBTYPE, scale_status = T,
              title =paste("Selected genes-tcga_subtype Freq=50,","n=",length(boruta_features50),sep=""))
print(p)
dev.off()

save.image(paste(results_base,"feature_selection_tcga_subtype.rdata",sep ="/"))

gene_list_freq25 <-data.table(gene=colnames(pca_data_freq25))
gene_list_freq40 <-data.table(gene=colnames(pca_data_freq40))
gene_list_freq50 <-data.table(gene=colnames(pca_data_freq50))
fwrite(gene_list_freq25, file =paste(results_base,
                                     paste(comparison_header,"gene_list_freq25.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)

fwrite(gene_list_freq40, file =paste(results_base,
                                     paste(comparison_header,"gene_list_freq40.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)
fwrite(gene_list_freq50, file =paste(results_base,
                                     paste(comparison_header,"gene_list_freq50.txt",sep="_"), sep="/")
       ,eol = "\n", row.names = F, col.names = T)
#Saveplots(getwd())

############
### Step4 create a model  

prev_results_base <- paste(base,'Step3_feature_selection',comparison_header,sep="/")
load(paste(prev_results_base,'feature_selection_tcga_subtype.rdata',sep="/"))
each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
minor_group <- each_pheno_num[each_pheno_num$n ==min(each_pheno_num$n),]$SUBTYPE

train_control <- trainControl(method = "cv" , number = 5,
                              classProbs = TRUE)

results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)
#data_model <- read.csv( "data_model_tcga_subtype.csv",header = T,row.names = 1)

#Freq = 25 times 
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq25)
dim(data_model) #filter = 25 times, 240 genes
#664 241

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq25 <- Model_50_RF(data = data_model,5)
#View(res_tcga_subtype_rf_freq25)

res_tcga_subtype_svm_freq25 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq25)


#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq25,paste(results_base,"res_tcga_subtype_rf_freq25.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq25,paste(results_base,"res_tcga_subtype_svm_freq25.csv",sep="/"))

##########

#Freq freq40 times
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq40)

dim(data_model) #freq = 40 times, 506 genes
#664 207

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq40 <- Model_50_RF(data = data_model,5)
#View(res_tcga_subtype_rf_freq40)

res_tcga_subtype_svm_freq40 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq40)

#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq40,paste(results_base,"res_tcga_subtype_rf_freq40.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq40,paste(results_base,"res_tcga_subtype_svm_freq40.csv",sep="/"))

###### 
#Freq freq50 times 
data_model<- data.frame(SUBTYPE= factor(pheno_tcga$SUBTYPE), pca_data_freq50)

dim(data_model) #freq = 15 times, 7 genes
#664 138

#check if class is fortor
class(data_model$SUBTYPE)

res_tcga_subtype_rf_freq50 <- Model_50_RF(data = data_model,5)
#View(res_tcga_subtype_rf_freq50)

res_tcga_subtype_svm_freq50 <- Model_50_SVM(data = data_model,5)
#View(res_tcga_subtype_svm_freq50)

#PCA_plot(data_model[,-1],data_model$SUBTYPE,title = "PCA of tcga data,subtype(n=12)")
write.csv(res_tcga_subtype_rf_freq50,paste(results_base,"res_tcga_subtype_rf_freq50.csv",sep="/"))
write.csv(res_tcga_subtype_svm_freq50,paste(results_base,"res_tcga_subtype_svm_freq50.csv",sep="/"))

sum_model_ALL <- data.frame(metrics = rownames(res_tcga_subtype_rf_freq50),
                            rf_freq25= res_tcga_subtype_rf_freq25$mean,
                            svm_freq25= res_tcga_subtype_svm_freq25$mean,
                            rf_freq40= res_tcga_subtype_rf_freq40$mean,
                            svm_freq40= res_tcga_subtype_svm_freq40$mean,
                            rf_freq50= res_tcga_subtype_rf_freq50$mean,
                            svm_freq50= res_tcga_subtype_svm_freq50$mean)
#View(sum_model_ALL)

write.csv(sum_model_ALL,paste(results_base,"sum_model_tcga_subtype.csv",sep="/"))

save.image(paste(results_base,"model_tcga_subtype.rdata",sep="/"))

df_genes <- data.frame(freq50 = boruta_features50)
write.csv(df_genes,paste(results_base,"df_genes_freq50_subtypes.csv",sep="/"))

############
## step5 Validation in Dogs
prev_results_base <- paste(base,'Step4_model_create',comparison_header,sep="/")
load(paste(prev_results_base,'model_tcga_subtype.rdata',sep="/"))

each_pheno_num <- pheno_tcga %>% count(SUBTYPE)
minor_group <- each_pheno_num[each_pheno_num$n ==min(each_pheno_num$n),]$SUBTYPE

data_test<- read.csv(paste(base,"all_cmt_combat_corrected.csv",sep="/"),header = T, row.names = 1)
pheno_test <- read.csv(paste(base,"phenotype_cmt.csv",sep="/"),header = T)

if (toupper('basal') %in% toupper(comparison)){
  other_group <- comparison[!comparison %in% "Basal"]
  pheno_test[pheno_test$SUBTYPE!="basal",]$SUBTYPE <- other_group
  pheno_test[pheno_test$SUBTYPE=="basal",]$SUBTYPE <- "Basal"
}else if (!toupper('basal') %in% toupper(comparison) & toupper('LumA') %in% toupper(comparison)){
  other_group <- comparison[!comparison %in% "LumA"]
  pheno_test[pheno_test$SUBTYPE!="basal",]$SUBTYPE <- "LumA"
  pheno_test[pheno_test$SUBTYPE=="basal",]$SUBTYPE <- other_group
}


pheno_train <- pheno_tcga
dim(data_test)
#  11856    78



results_base <- paste(base,'Step5_model_validation',comparison_header,sep="/")
dir.create(results_base,recursive = TRUE)

data_test <- data.frame(data_test)

train_freq25 <- data.frame(t(pca_data_freq25))
train_freq40 <-data.frame(t(pca_data_freq40))
train_freq50 <- data.frame(t(pca_data_freq50))

dim(train_freq25)

test_freq25 <- data_test[ row.names(data_test) %in% row.names(train_freq25),]
test_freq40 <- data_test[ row.names(data_test) %in% row.names(train_freq40),]
test_freq50 <- data_test[ row.names(data_test) %in% row.names(train_freq50),]
dim(test_freq25)

abc <- colnames(test_freq25)== pheno_test$PATIENT_ID # ALL TRUE


RF_validate <- function(train_freq, test_freq) {
  # freq25 times
  train_data <- data.frame(SUBTYPE= factor(pheno_train$SUBTYPE), t(train_freq))
  test_data <- data.frame(SUBTYPE= factor(pheno_test$SUBTYPE), t(test_freq))
  
  
  if (sum(table(train_data$SUBTYPE)) >= 0.8 * sum(table(train_data$SUBTYPE)) &&
      sum(table(train_data$SUBTYPE)) <= 1.2 * sum(table(train_data$SUBTYPE))) {
    train_data_balanced <- train_data
  } else {
    # We use ADASYN
    genData_ADAS <- ADAS(X = train_data[, -1], target = train_data$SUBTYPE, K = 5)
    train_df_adas <- genData_ADAS[["data"]]
    class_df <- train_df_adas[ncol(train_df_adas)]
    train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[, -ncol(train_df_adas)])
    train_data_balanced <- train_df_adas
    names(train_data_balanced)[1] <- "SUBTYPE"
  }
  
  train_data_balanced <- train_data_balanced %>%
    mutate(across(.cols = 1, .fns = factor)) %>%
    mutate(across(.cols = 2:length(.), .fns = as.numeric))
  train_control <- trainControl(
    method = "cv", number = 5, classProbs = TRUE)
  
  # Random Forest
  rf_model <- caret::train(SUBTYPE ~ ., data = train_data_balanced, method = "ranger",
                           verbose = FALSE, trControl = train_control, metric = 'Accuracy',
                           tuneLength = 10)
  # Prediction
  res_rf_freq <- Res_CMT(prob = predict(rf_model, test_data, type = "prob"),
                         pred = predict(rf_model, test_data),
                         test_data = test_data, model = rf_model)
  print(res_rf_freq)
}

res_rf_freq25 = RF_validate(train_freq25, test_freq25)
res_rf_freq40 = RF_validate(train_freq40, test_freq40)
res_rf_freq50 = RF_validate(train_freq50, test_freq50)

test_model_sum <- data.frame(freq25 = res_rf_freq25,
                             freq40 = res_rf_freq40,
                             freq50 = res_rf_freq50)

#View(test_model_sum)

write.csv(test_model_sum, paste(results_base,"model_sum_validate_subtype.csv",sep="/"))
save.image(paste(results_base,"validate_model_subtype.rdata",sep="/"))

#PCA

pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_25.pdf",sep=""),
               sep="/") 
    ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq25),pheno_test$SUBTYPE, scale_status = F,
               title =paste("Selected genes-cmt_subtype Freq=25, n= ", nrow(train_freq25),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_40.pdf",sep=""),sep="/")
    ,height = 4.5,width = 6)

p <- PCA_plot(t(test_freq40),pheno_test$SUBTYPE, scale_status = F,
               title =paste("Selected genes-cmt_subtype Freq=40, n= ", nrow(train_freq40),sep=""))
print(p)
dev.off()
pdf(file=paste(results_base,
               paste(comparison_header,"Selected genes-tcga_subtype_Freq_50.pdf",sep=""),
               sep="/")
    ,height = 4.5,width = 6)
p <- PCA_plot(t(test_freq50),pheno_test$SUBTYPE, scale_status = F,
               title =paste("Selected genes-cmt_subtype Freq=50, n= ", nrow(train_freq50),sep=""))

print(p)
dev.off()