PCA_plot <- function(data,group,label= "", title = "") {
  Data_pca <- prcomp(scale(data,center = T,scale = T))
  contributions <- summary(Data_pca)$importance[2,]*100
  dataScores12 <-data.frame(Data_pca$x[, 1:2])
  
  if (label == "") {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ 
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16) + ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ geom_text(aes(label=label))+
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16)+ ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  return(p)
  
}

PCA_plot2 <- function(data,group,label= "", title = "") {
  Data_pca <- prcomp(scale(data,center = T,scale = F))
  contributions <- summary(Data_pca)$importance[2,]*100
  dataScores12 <-data.frame(Data_pca$x[, 1:2])
  
  if (label == "") {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ 
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16) + ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    p <- ggplot(dataScores12,aes(x=PC1,y=PC2))+
      geom_point(aes(col=group),size=2)+ geom_text(aes(label=label))+
      xlab(paste("PC1","(",round(contributions[1],digits=2),"%)",sep=""))+
      ylab(paste("PC2","(",round(contributions[2],digits=2),"%)",sep=""))+
      theme_bw(base_size = 16)+ ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  
}

CBFpvcaFunction <- function(data,phenotypedata){
  
  #The only thing the function now assumes is that your phenotype table 
  #is in the same order as your columns of your data.
  
  #packages <- c("pvca", "ggplot2", "Biobase") 
  #lapply(packages, library, character.only = TRUE, quietly = TRUE)
  
  if(length(seq_len(nrow(phenotypedata))) != length(seq_len(nrow(data)))) {
    warning(paste0("Data needs to be a wide", "\n", "Transposing data", "\n"))
    data <- t(data)
  }
  
  
  if(!identical(rownames(phenotypedata), rownames(data))) {
    warning(paste0("Rownames in phenotypedata do not match the rownames of data", "\n"))
    warning(paste0("Matching rownames in phenotypedata to data", "\n"))
    
    if(nrow(phenotypedata) != nrow(data)) {
      stop("Phenotypedata and data have different amounts of rows")
    } else {
      rownames(phenotypedata) <- rownames(data)
    }
  }
  
  data1_pheno_formatted<-new("AnnotatedDataFrame",data=phenotypedata)
  
  data_forPVCA<-ExpressionSet(assayData=t(data),phenoData=data1_pheno_formatted)
  
  pvcaObj_data <- pvcaBatchAssess (abatch = data_forPVCA, batch.factors = colnames(data1_pheno_formatted@data), threshold = 0.7)
  
  data_pvca_res<-data.frame(as.data.frame(pvcaObj_data$label),t(as.data.frame(pvcaObj_data$dat)))
  
  colnames(data_pvca_res)<-c("effect","variance")
  
  
  p<-ggplot((data_pvca_res[-nrow(data_pvca_res),]), aes(x= effect, y = variance)) + 
    geom_bar(stat = 'identity',  position = 'dodge', col ='transparent')+ 
    scale_fill_discrete(guide = 'none') + 
    theme_bw(base_size = 16)+
    theme(plot.title =  element_text(hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))+
    labs(x = 'Effects', y = 'Weighted average proportion variance')+
    ggtitle("PVCA estimation bar chart corrected data")
  
  return(p)
  
}
