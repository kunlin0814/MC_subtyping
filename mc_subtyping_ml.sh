#!/bin/bash
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=50G                   # Job memory limit (10 GB)
#SBATCH --time=10:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=Step1_extract_somatic.%j.out    # Standard output log
#SBATCH --error=Step1_extract_somatic.%j.err     # Standard error log

scripts_location='/scratch/kh31516/Josh_MC_ML_gene/scripts'

ml R/4.0.0-foss-2019b

# package_location <- as.character(args[1])
# new_data_file <- as.character(args[2])
# out_file_name <- as.character(args[3])

Rscript --vanilla "${scripts_location}/2.tcga_DEG_subtype.r" \


Rscript --vanilla "${scripts_location}/3.feature_selection_tcga_subtype.r" \

Rscript --vanilla "${scripts_location}/4.model_tcga_subtype.r" \

Rscript --vanilla "${scripts_location}/5.validate_model_subtype.r" \
