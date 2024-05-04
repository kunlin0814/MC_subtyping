# Subtyping Canine Mammary Tumors with Selected Genes using Machine Learning Models

## Description

The goal of this project is to utilize machine learning methods to identify the most relevant genes for subtyping human breast cancer. These relevant genes can also be employed to subtype canine mammary tumors, enabling a comparative oncology analysis.

The project implements a pipeline consisting of five steps to analyze RNA-seq expression data from the TCGA dataset. The pipeline aims to select the most relevant genes and evaluate their effectiveness in subtyping mammary tumors in dogs.

The subtyping results were used in the research conducted by Watson, J. et al.,
Titled 'Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs,' published in Breast Cancer Res 25, 114 (2023)

## Pipeline Steps

The pipeline consists of the following five steps:

1. **Step 1**: Batch Effect Correction: Apply Combat to perform batch effect correction on the expression data.

2. **Step 2**: Differential Gene Expression: Identify differentially expressed genes.

3. **Step 3**: Feature Selection: Perform gene (feature) selection using the Boruta method with differentially expressed genes in a random forest model. The process is repeated 50 times.

4. **Step 4**: Model Creation: Build a random forest model and evaluate its performance on a training set.

5. **Step 5**: Model Validation: Evaluate the performance of the model using the feature genes identified in Step 3.

## Usage

To run the script and execute the pipeline, follow these steps:

1. Clone the repository: `git clone https://github.com/your-username/your-repo.git`

2. Run the following scripts in order:

   - `1.combat_2data.R`
   - `2.tcga_DEG_subtype.R`
   - `3.feature_selection_tcga_subtype.R`
   - `4.model_tcga_subtype.R`
   - `5.validate_model_subtype.R`

   Note: Ensure that the necessary variables are loaded and the required R packages are installed before running each script.

## Dependencies

The following R packages are required to run the scripts:

- datable
- DMwR
- grid
- Boruta
- vars
- varSelRF
- caret
- e1071
- ranger
- sva
- smotefamily
- mlr3measures
- kernlab
- edgeR
- limma
- tmaptools
- tidyverse
- reshape2
- biomaRt
- randomForest

## Results

After executing the pipeline, the script will produce the following results:

- Batch-corrected expression data
- Differential gene expression results
- Selected feature genes
- Random forest model
- Model performance evaluation

The overall directory structure should resemble the following:

```
base_folder/
base_folder/Step1_combat/
base_folder/Step2DEG/
base_folder/Step3_feature_selection/
base_folder/Step4_model_create/
base_folder/Step5_model_validation/
```

## Publication

Watson, J., Wang, T., Ho, KL. et al. Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs. Breast Cancer Res 25, 114 (2023). https://doi.org/10.1186/s13058-023-01705-5

## Acknowledgements

We extend our heartfelt gratitude to Tanakamol Mahawan from The University of Liverpool for his invaluable contributions in designing the prototype. His expertise and insights played a pivotal role in shaping the project and ensuring its success.
