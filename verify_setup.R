# Verification Script
# Description: Checks if the project structure is correct and functions can be loaded.

message("Verifying project setup...")

# 1. Check Directories
dirs <- c("R", "scripts", "config", "data", "results")
missing_dirs <- dirs[!dir.exists(dirs)]
if (length(missing_dirs) > 0) {
    stop("Missing directories: ", paste(missing_dirs, collapse = ", "))
}
message("Directories OK.")

# 2. Check Files
files <- c(
    "config/config.yml",
    "R/utils.R",
    "R/plotting.R",
    "R/analysis.R",
    "scripts/01_combat.R",
    "scripts/02_deg.R",
    "scripts/03_feature_selection.R",
    "scripts/04_model_training.R",
    "scripts/05_validation.R",
    "main.R"
)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
    stop("Missing files: ", paste(missing_files, collapse = ", "))
}
message("Files OK.")

# 3. Test Loading Functions
tryCatch(
    {
        source("R/utils.R")
        source("R/plotting.R")
        source("R/analysis.R")
        config <- load_config()
        message("Functions and Config loaded OK.")
    },
    error = function(e) {
        stop("Failed to load functions/config: ", e$message)
    }
)

message("Verification successful! Please place your data files in 'data/raw' before running main.R.")
