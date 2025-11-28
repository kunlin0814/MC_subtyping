# Subtyping Canine Mammary Tumors with Selected Genes using Machine Learning Models

This pipeline selects subtype-informative genes from TCGA breast cancer RNA-seq data and evaluates how well those features classify canine mammary tumors.

The work supports the findings in Watson, J. et al., *Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs*, **Breast Cancer Research** (2023).

## Pipeline Overview
1. **Batch Effect Correction** – Adjust human and canine expression data with ComBat.  
2. **Differential Expression** – Identify subtype-specific DEGs in TCGA.  
3. **Feature Selection** – Apply Boruta (random forest) to select robust genes.  
4. **Model Training** – Train and evaluate a random forest classifier on TCGA.  
5. **Validation** – Test features and the model on canine mammary tumor data.  

## Prerequisites
- R (tested with 4.3)  
- Suggested packages: `data.table`, `sva`, `pvca`, `ggplot2`, `limma`, `edgeR`, `EnhancedVolcano`, `Boruta`, `randomForest`, `caret`, `ComplexHeatmap`, plus Bioconductor manager `BiocManager`.  
- Scripts auto-install missing packages via `load_dependencies()` (see `R/utils.R`), but pre-installing can save time.
- Pipeline parameters live in `config/config.yml` (e.g., `boruta_runs`, `model_rounds`, CV folds, tune length, and file names).
- TCGA expression file `TCGA_BRCA_log2TPM+1.csv` is large; provide it on request if not already available locally.

## Data
Expected file names and locations are configurable in `config/config.yml`:
- `data/raw/TCGA_BRCA_log2TPM+1.csv` – human log2(TPM+1) expression  
- `data/raw/cmt_all.csv` – canine expression  
- `data/raw/All_dog_TCGA_subtype_meta.txt` – metadata with `PATIENT_ID`, `SOURCE`, and `SUBTYPE`  

Outputs are written under `data/processed/` and `results/` (figures in `results/figures/`).

## Quickstart
1. Place input files in `data/raw/` (or adjust paths in `config/config.yml`).  
2. From the repo root, run the steps in order:
   ```bash
   Rscript scripts/01_combat.R
   Rscript scripts/02_deg.R
   Rscript scripts/03_feature_selection.R
   Rscript scripts/04_model_training.R
   Rscript scripts/05_validation.R
   ```
   Each script will create its required output directories automatically.
   Or run everything at once: `Rscript main.R`.

## Repository Layout
- `scripts/` – step-wise pipeline scripts (`01_combat.R` … `05_validation.R`).  
- `R/` – helpers (`utils.R`, plotting functions).  
- `config/config.yml` – file paths and model parameters (e.g., Boruta runs, CV folds).  
- `data/` – raw inputs in `data/raw/`, processed data in `data/processed/`.  
- `results/` – intermediate results, figures, and trained models.

## Publication
Watson J, Wang T, Ho KL, et al. Human basal-like breast cancer is represented by one of the two mammary tumor subtypes in dogs. *Breast Cancer Research*. 2023;25:114. https://doi.org/10.1186/s13058-023-01705-5

## Acknowledgements
We extend our heartfelt gratitude to Tanakamol Mahawan from The University of Liverpool for his invaluable contributions in designing the prototype.
