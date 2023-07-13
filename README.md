# Subtyping Canine Mammary Tumors with Selected Genes using Machine Learning Models

## Description

The goal of this project is to utilize machine learning methods to identify the most relevant genes for subtyping human breast cancer. These relevant genes can also be employed to subtype canine mammary tumors, enabling a comparative oncology analysis.

The project implements a pipeline consisting of five steps to analyze RNA-seq expression data from the TCGA dataset. The pipeline aims to select the most relevant genes and evaluate their effectiveness in subtyping mammary tumors in dogs.

## Pipeline Steps

The pipeline consists of the following five steps:

1. **Step 1**: Batch Effect Correction: Apply Combat to perform batch effect correction on the expression data.

2. **Step 2**: Differential Gene Expression: Identify differentially expressed genes.

3. **Step 3**: Feature Selection: Perform gene (feature) selection using the Boruta method with differentially expressed genes in a random forest model. The process is repeated 50 times.

4. **Step 4**: Model Creation: Build a random forest model and evaluate its performance on a training set.

5. **Step 5**: Model Validation: Evaluate the performance of the model using the feature genes identified in Step 3.

## Usage

To run the script and execute the pipeline, follow these steps:

1. Clone the repository: `git clone https://github.com/kunlin0814/MC_subtyping`

2. Install the required dependencies:

3. Run the following scripts in order:

   - `1.combat_2data.R`
   - `2.tcga_DEG_subtype.R`
   - `3.feature_selection_tcga_subtype.R`
   - `4.model_tcga_subtype.R`
   - `5.validate_model_subtype.R`

Note:

- Make sure to set the `base` path consistently across all scripts and ensure the comparison used is the same. For example, if you set `comparison <- c("Basal", "LumA")`, it means you are comparing the Basal subtype with the LumA subtype.
- Before running each script, ensure that the necessary variables are loaded and the required R packages are installed in your R environment. You can use the `install.packages()` function to install any missing packages.

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
- Model performance evaluation
- Subypting results in canine mammary tumor dataset

The overall directory structure should resemble the following:

```
base_folder/
base_folder/Step1_combat/
base_folder/Step2DEG/
base_folder/Step3_feature_selection/
base_folder/Step4_model_create/
base_folder/Step5_model_validation/
```

## Acknowledgements

We would like to express our gratitude to Tanakamol Mahawan for his invaluable contribution to designing the pipeline and the prototype design. His expertise and insights were instrumental in shaping the project and its success.
