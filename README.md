# Mammary tumor subtyping

## Description

This project implements a pipeline consisting of five steps to analyze RNA-seq expression data and build a random forest model.

## Pipeline Steps

The whole pipeline involved 5 steps, but this script included the last four steps in the whole pipeline.

1. **Step1**: Use Combat to perform batch effect correction on the expression data.

2. **Step2**: Identify Differential Gene expression.

3. **Step3**: Perform feature selection using the Boruta method in a random forest model. The process is repeated 50 times.

4. **Step4**: Create a random forest model and evaluate its performance in a training set.

5. **Step5**: Evaluate the model performance using the feature genes identified in Step3.

## Usage

To run the script and execute the pipeline, follow these steps:

1. Clone the repository: `git clone https://github.com/your-username/your-repo.git`

2. Install the required dependencies:

## Dependencies

The following R packages are required to run the script:

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

3. Run the script: `python script.py`

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

## License

[Specify the license under which the project is released, if applicable.]

## Acknowledgements

[Mention any acknowledgements or credits for external contributions, libraries, or resources used in the project.]
