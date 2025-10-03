# Subtyping Canine Mammary Tumors with Selected Genes using Machine Learning Models

## Description

The goal of this project is to utilize machine learning methods to identify the most relevant genes for subtyping human breast cancer. These relevant genes can also be employed to subtype canine mammary tumors, enabling a comparative oncology analysis.

The project implements a pipeline consisting of five steps to analyze RNA-seq expression data from the TCGA dataset. The pipeline aims to select the most relevant genes and evaluate their effectiveness in subtyping mammary tumors in dogs.

The subtyping results were used in the research conducted by Watson, J. et al.,
Titled 'Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs,' published in Breast Cancer Res 25, 114 (2023)

## Pipeline Overview
1. **Batch Effect Correction** – Adjust expression data with ComBat.  
2. **Differential Expression** – Identify subtype-specific differentially expressed genes (DEGs).  
3. **Feature Selection** – Apply Boruta with random forests (50 iterations) to select robust features.  
4. **Model Training** – Train and evaluate a random forest classifier on TCGA data.  
5. **Validation** – Test selected features and models on canine mammary tumor data.  


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
 ├─ Step1_combat/
 ├─ Step2_DEG/
 ├─ Step3_feature_selection/
 ├─ Step4_model_create/
 └─ Step5_model_validation/
```

## Publication

Watson J, Wang T, Ho KL, et al. Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs.
Breast Cancer Research. 2023;25:114. https://doi.org/10.1186/s13058-023-01705-5

## Acknowledgements

We extend our heartfelt gratitude to Tanakamol Mahawan from The University of Liverpool for his invaluable contributions in designing the prototype. His expertise and insights played a pivotal role in shaping the project and ensuring its success.
