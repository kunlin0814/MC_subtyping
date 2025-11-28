#' Balance Data using ADAS
#'
#' @param data Data frame with SUBTYPE column
#' @return Balanced data frame
#' @export
Balance_data <- function(data) {
    # Ensure SUBTYPE is a factor
    data$SUBTYPE <- as.factor(data$SUBTYPE)

    # ADAS requires numeric features
    # Assuming data has SUBTYPE as first column
    genData_ADAS <- ADAS(X = data[, -1], target = data$SUBTYPE, K = 5)
    train_df_adas <- genData_ADAS[["data"]]

    names(train_df_adas)[ncol(train_df_adas)] <- "SUBTYPE"

    # Move SUBTYPE to first column
    train_data_balanced <- train_df_adas %>%
        select(SUBTYPE, everything()) %>%
        mutate(across(.cols = 1, .fns = factor)) %>%
        mutate(across(.cols = 2:length(.), .fns = as.numeric))

    return(train_data_balanced)
}

#' Calculate Model Performance Metrics
#'
#' @param prob Prediction probabilities
#' @param pred Predicted classes
#' @param test_data Test data frame
#' @param model Trained model object
#' @return A vector of performance metrics
#' @export
Res_CMT <- function(prob, pred, test_data, model) {
    levels_subtype <- levels(test_data$SUBTYPE)

    # AUC
    roc_c <- tryCatch(pROC::auc(test_data$SUBTYPE, prob[, 2]), error = function(e) NA)
    ROC_AUC <- if (is.na(roc_c)) 0 else roc_c[1]

    PR_AUC <- tryCatch(
        {
            MLmetrics::PRAUC(y_pred = prob[, 2], y_true = ifelse(test_data$SUBTYPE == levels_subtype[2], 1, 0))
        },
        error = function(e) NA
    )
    if (is.na(PR_AUC)) PR_AUC <- 0

    # Confusion Matrix
    actual <- as.factor(test_data$SUBTYPE)
    predicted <- as.factor(pred)

    # Ensure factor levels match
    predicted <- factor(predicted, levels = levels(actual))

    cm <- as.matrix(table(Actual = actual, Predicted = predicted))

    rowsums <- apply(cm, 1, sum)
    colsums <- apply(cm, 2, sum)

    precision <- diag / colsums
    recall <- diag / rowsums
    f1 <- 2 * precision * recall / (precision + recall)

    # Handle NaN
    precision[is.na(precision)] <- 0
    recall[is.na(recall)] <- 0
    f1[is.na(f1)] <- 0

    macroPrecision <- mean(precision)
    macroRecall <- mean(recall)
    macroF1 <- mean(f1)
    mcc <- mltools::mcc(preds = predicted, actuals = actual)

    # Return vector matching original order
    # PC0, PC1, RC0, RC1, F10, F11, MPC, MRC, MF1, PR_AUC, ROC_AUC, MCC
    Res <- c(
        precision[1], precision[2], recall[1], recall[2], f1[1], f1[2],
        macroPrecision, macroRecall, macroF1, PR_AUC, ROC_AUC, mcc
    )
    return(Res)
}

#' Train Random Forest Model (50 runs)
#'
#' @param data Data frame
#' @param nround Number of rounds
#' @return Data frame of results
#' @export
Model_50_RF <- function(data, nround = 50, cv_folds = 5, tune_length = 10) {
    results_matrix <- matrix(0, nrow = nround, ncol = 12)
    colnames(results_matrix) <- c(
        "PC0", "PC1", "RC0", "RC1", "F10", "F11",
        "MPC", "MRC", "MF1", "PR_AUC", "ROC_AUC", "MCC"
    )

    for (i in 1:nround) {
        set.seed(i)

        train_indices <- createDataPartition(data$SUBTYPE, p = 0.8, list = FALSE, times = 1)
        train_data_raw <- data[train_indices, ]
        test_data <- data[-train_indices, ]

        train_data_balanced <- Balance_data(train_data_raw)

        train_control <- trainControl(method = "cv", number = cv_folds, classProbs = TRUE)

        rf_model <- caret::train(SUBTYPE ~ .,
            data = train_data_balanced, method = "ranger",
            verbose = FALSE, trControl = train_control, metric = "Accuracy",
            tuneLength = tune_length
        )

        res_rf <- Res_CMT(
            prob = predict(rf_model, test_data, type = "prob"),
            pred = predict(rf_model, test_data),
            test_data = test_data, model = rf_model
        )

        results_matrix[i, ] <- res_rf
    }

    df_res <- as.data.frame(results_matrix)

    # Summary statistics
    summary_df <- data.frame(
        mean = apply(df_res, 2, mean),
        sd = apply(df_res, 2, sd)
    )
    summary_df$CI95 <- 2 * sqrt((summary_df$mean * (1 - summary_df$mean)) / nround) # Approx CI
    summary_df <- round(summary_df, 3)

    return(summary_df)
}

#' Train SVM Model (50 runs)
#'
#' @param data Data frame
#' @param nround Number of rounds
#' @return Data frame of results
#' @export
Model_50_SVM <- function(data, nround = 50, cv_folds = 5, tune_length = 10) {
    results_matrix <- matrix(0, nrow = nround, ncol = 12)
    colnames(results_matrix) <- c(
        "PC0", "PC1", "RC0", "RC1", "F10", "F11",
        "MPC", "MRC", "MF1", "PR_AUC", "ROC_AUC", "MCC"
    )

    for (i in 1:nround) {
        set.seed(i)

        train_indices <- createDataPartition(data$SUBTYPE, p = 0.8, list = FALSE, times = 1)
        train_data_raw <- data[train_indices, ]
        test_data <- data[-train_indices, ]

        train_data_balanced <- Balance_data(train_data_raw)

        train_control <- trainControl(method = "cv", number = cv_folds, classProbs = TRUE)

        svm_fit <- caret::train(SUBTYPE ~ .,
            data = train_data_balanced, method = "svmRadial",
            trControl = train_control, preProcess = c("center", "scale"),
            tuneLength = tune_length
        )

        res_svm <- Res_CMT(
            prob = predict(svm_fit, test_data, type = "prob"),
            pred = predict(svm_fit, test_data),
            test_data = test_data, model = svm_fit
        )

        results_matrix[i, ] <- res_svm
    }

    df_res <- as.data.frame(results_matrix)

    summary_df <- data.frame(
        mean = apply(df_res, 2, mean),
        sd = apply(df_res, 2, sd)
    )
    summary_df$CI95 <- 2 * sqrt((summary_df$mean * (1 - summary_df$mean)) / nround)
    summary_df <- round(summary_df, 3)

    return(summary_df)
}
