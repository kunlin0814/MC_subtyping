#' Load Configuration
#'
#' This function loads the configuration from config/config.yml
#' @return A list containing the configuration
#' @export
load_config <- function(config_file = "config/config.yml") {
    if (!file.exists(config_file)) {
        stop("Configuration file not found: ", config_file)
    }
    config <- config::get(file = config_file)

    # Ensure directories exist
    dir.create(config$paths$results_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(config$paths$figures_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(config$paths$processed_data_dir, showWarnings = FALSE, recursive = TRUE)

    return(config)
}

#' Install and Load Dependencies
#'
#' @param packages Character vector of package names
#' @export
load_dependencies <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            message(paste("Installing", pkg, "..."))
            if (pkg %in% c("MultiBaC", "MultiAssayExperiment", "sva", "pvca", "limma", "edgeR", "ComplexHeatmap")) {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager")
                }
                BiocManager::install(pkg)
            } else {
                install.packages(pkg)
            }
        }
        library(pkg, character.only = TRUE)
    }
}
