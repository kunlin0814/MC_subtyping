# Dependencies

The following R packages are required to run the analysis pipeline.  

## Core Machine Learning & Modeling
- Boruta – Feature selection with random forests  
- caret – Model training and evaluation framework  
- ranger – Fast implementation of random forests  
- randomForest – Classic random forest implementation  
- mlr3measures – Performance metrics for machine learning  
- varSelRF – Variable selection using random forests  
- e1071 – SVMs and other machine learning algorithms  
- kernlab – Kernel-based machine learning methods  

## Data Preprocessing & Batch Correction
- sva – Surrogate variable analysis / batch correction  
- DMwR – Data mining with R (includes SMOTE)  
- smotefamily – Advanced oversampling methods (imbalanced data)  

## Statistical Analysis
- vars – Vector autoregressive models  
- edgeR – Differential expression analysis of RNA-seq data  
- limma – Linear models for microarray/RNA-seq  

## Data Handling & Visualization
- tidyverse – Collection of packages for data science  
- reshape2 – Reshaping data frames  
- tmaptools – Utility functions for spatial data (if needed)  
- biomaRt – BioMart database queries for gene annotation  
- grid – Base graphics system for advanced plotting  

## Utilities
- data.table – Fast data manipulation  


### Installation

You can install all required packages in R using:

```r
install.packages(c(
  "Boruta", "caret", "ranger", "randomForest", "mlr3measures",
  "varSelRF", "e1071", "kernlab", "sva", "DMwR", "smotefamily",
  "vars", "edgeR", "limma", "tidyverse", "reshape2",
  "tmaptools", "biomaRt", "grid"
))