# Master Workflow Script
# Description: Runs the entire analysis pipeline.

# Load config and utils
source("R/utils.R")
config <- load_config()

# Define steps
steps <- c(
    "scripts/01_combat.R",
    "scripts/02_deg.R",
    "scripts/03_feature_selection.R",
    "scripts/04_model_training.R",
    "scripts/05_validation.R"
)

# Run steps
for (step in steps) {
    message("\n========================================================")
    message("Running step: ", step)
    message("========================================================\n")

    tryCatch(
        {
            source(step)
        },
        error = function(e) {
            message("Error in step ", step, ": ", e$message)
            stop("Pipeline failed at step: ", step)
        }
    )
}

message("\nPipeline completed successfully!")
