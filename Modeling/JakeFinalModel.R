## SETUP #######################################################################
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(ranger)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)

## PARAMETERS ##################################################################
params <- list(
  positive_class = "POSITIVE_POLYP",
  negative_class = "NEGATIVE_POLYP",
  tune_num_trees = 500,
  final_num_trees = 1000,
  num_threads = 8,
  jitter_amount = 1e-9,
  min_hc_controls = 10,
  prevalence_grid = c(0.10, 0.20, 0.30),
  mtry_fraction_grid = c(0.001, 0.003, 0.01),
  min_node_size_grid = c(20, 50, 100),
  max_depth_grid = c(10, 15, NA),
  top_n_features = 2000,
  min_feature_variance = 1e-8
)

## HELPERS #####################################################################
build_combined_features <- function(study) {
  tax_mat <- assay(study$relative_abundance, 1)
  pab_mat <- assay(study$pathway_abundance, 1)
  pcov_mat <- assay(study$pathway_coverage, 1)

  rownames(tax_mat) <- paste0("TAX_", rownames(tax_mat))
  rownames(pab_mat) <- paste0("PAB_", rownames(pab_mat))
  rownames(pcov_mat) <- paste0("PCOV_", rownames(pcov_mat))

  common_samples <- Reduce(intersect, list(
    colnames(tax_mat),
    colnames(pab_mat),
    colnames(pcov_mat)
  ))

  tax_mat <- tax_mat[, common_samples, drop = FALSE]
  pab_mat <- pab_mat[, common_samples, drop = FALSE]
  pcov_mat <- pcov_mat[, common_samples, drop = FALSE]

  rbind(tax_mat, pab_mat, pcov_mat)
}

prevalence_filter <- function(mat, threshold) {
  prev <- colSums(mat > 0) / nrow(mat)
  names(prev[prev >= threshold])
}

add_jitter <- function(mat, amount) {
  mat + matrix(
    runif(length(mat), min = 0, max = amount),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
}

percentile_norm <- function(values, controls) {
  if (length(controls) == 0) {
    return(rep(NA_real_, length(values)))
  }

  vapply(values, function(x) mean(controls <= x) * 100, numeric(1))
}

normalize_train_test <- function(
  X_train,
  X_test,
  metadata_train,
  threshold,
  jitter_amount,
  min_hc_controls,
  top_n_features,
  min_feature_variance
) {
  keep_features <- prevalence_filter(X_train, threshold)

  if (length(keep_features) == 0) {
    return(NULL)
  }

  filtered_train <- X_train[, keep_features, drop = FALSE]
  filtered_test <- X_test[, keep_features, drop = FALSE]

  filtered_train <- add_jitter(filtered_train, jitter_amount)
  filtered_test <- add_jitter(filtered_test, jitter_amount)

  control_idx <- which(metadata_train$disease_class == "HC")

  if (length(control_idx) < min_hc_controls) {
    return(NULL)
  }

  pnorm_train <- apply(filtered_train, 2, function(feature_values) {
    percentile_norm(feature_values, feature_values[control_idx])
  })

  pnorm_test <- sapply(seq_len(ncol(filtered_test)), function(j) {
    percentile_norm(filtered_test[, j], filtered_train[control_idx, j])
  })

  pnorm_train <- as.matrix(pnorm_train)
  pnorm_test <- as.matrix(pnorm_test)

  colnames(pnorm_train) <- colnames(filtered_train)
  rownames(pnorm_train) <- rownames(filtered_train)
  colnames(pnorm_test) <- colnames(filtered_test)
  rownames(pnorm_test) <- rownames(filtered_test)

  if (anyNA(pnorm_train) || anyNA(pnorm_test)) {
    return(NULL)
  }

  feature_variances <- apply(pnorm_train, 2, var)
  variance_features <- names(feature_variances[is.finite(feature_variances) & feature_variances > min_feature_variance])

  if (length(variance_features) == 0) {
    return(NULL)
  }

  ranked_features <- names(sort(feature_variances[variance_features], decreasing = TRUE))
  selected_features <- head(ranked_features, top_n_features)

  pnorm_train <- pnorm_train[, selected_features, drop = FALSE]
  pnorm_test <- pnorm_test[, selected_features, drop = FALSE]

  list(
    train_matrix = pnorm_train,
    test_matrix = pnorm_test,
    keep_features = selected_features,
    control_count = length(control_idx)
  )
}

fit_rf_model <- function(train_df, num_trees, max_depth, min_node_size, mtry_fraction, num_threads) {
  class_counts <- table(train_df$RF_Class)
  class_weights <- sum(class_counts) / (length(class_counts) * class_counts)
  names(class_weights) <- names(class_counts)

  p <- ncol(train_df) - 1
  mtry_val <- max(1, ceiling(mtry_fraction * p))

  ranger(
    dependent.variable.name = "RF_Class",
    data = train_df,
    num.trees = num_trees,
    max.depth = if (is.na(max_depth)) NULL else max_depth,
    min.node.size = min_node_size,
    mtry = mtry_val,
    class.weights = class_weights,
    probability = TRUE,
    importance = "impurity",
    num.threads = num_threads
  )
}

f1_at_threshold <- function(y_true, positive_probs, threshold, positive_class, negative_class) {
  pred_class <- ifelse(positive_probs >= threshold, positive_class, negative_class)
  pred_factor <- factor(pred_class, levels = c(positive_class, negative_class))
  truth_factor <- factor(y_true, levels = c(positive_class, negative_class))

  cm <- confusionMatrix(pred_factor, truth_factor, positive = positive_class)
  as.numeric(cm$byClass["F1"])
}

find_best_threshold <- function(y_true, positive_probs, positive_class, negative_class) {
  thresholds <- unique(c(seq(0.10, 0.90, by = 0.01), positive_probs))
  thresholds <- thresholds[thresholds > 0 & thresholds < 1]

  scores <- vapply(
    thresholds,
    function(threshold) {
      f1_at_threshold(y_true, positive_probs, threshold, positive_class, negative_class)
    },
    numeric(1)
  )

  best_idx <- which.max(scores)

  list(
    threshold = thresholds[best_idx],
    f1 = scores[best_idx],
    tuning_curve = data.frame(threshold = thresholds, f1 = scores)
  )
}

score_predictions <- function(y_true, positive_probs, threshold, positive_class, negative_class) {
  pred_class <- ifelse(positive_probs >= threshold, positive_class, negative_class)
  pred_factor <- factor(pred_class, levels = c(positive_class, negative_class))
  truth_factor <- factor(y_true, levels = c(positive_class, negative_class))

  cm <- confusionMatrix(pred_factor, truth_factor, positive = positive_class)
  roc_obj <- roc(
    response = ifelse(truth_factor == positive_class, 1, 0),
    predictor = positive_probs,
    quiet = TRUE
  )

  results_summary <- data.frame(
    accuracy = as.numeric(cm$overall["Accuracy"]),
    kappa = as.numeric(cm$overall["Kappa"]),
    sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    specificity = as.numeric(cm$byClass["Specificity"]),
    pos_pred_value = as.numeric(cm$byClass["Pos Pred Value"]),
    neg_pred_value = as.numeric(cm$byClass["Neg Pred Value"]),
    balanced_accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    f1 = as.numeric(cm$byClass["F1"]),
    auc = as.numeric(auc(roc_obj)),
    threshold = threshold
  )

  list(
    confusion_matrix = cm,
    pred_class = pred_factor,
    pred_probs = positive_probs,
    roc_obj = roc_obj,
    summary = results_summary
  )
}

run_inner_lodo <- function(X_train_all, y_train_all, metadata_train_all, tune_grid, params) {
  inner_studies <- unique(metadata_train_all$study_name)
  tuning_results <- vector("list", nrow(tune_grid))

  for (grid_idx in seq_len(nrow(tune_grid))) {
    grid_row <- tune_grid[grid_idx, , drop = FALSE]
    oof_truth <- character(0)
    oof_probs <- numeric(0)
    fold_details <- list()
    skipped_folds <- character(0)

    for (validation_study in inner_studies) {
      
      message("Inner loop validation study: ", validation_study)
      inner_train_idx <- metadata_train_all$study_name != validation_study
      inner_valid_idx <- metadata_train_all$study_name == validation_study

      X_inner_train <- X_train_all[inner_train_idx, , drop = FALSE]
      X_inner_valid <- X_train_all[inner_valid_idx, , drop = FALSE]
      y_inner_train <- droplevels(y_train_all[inner_train_idx])
      y_inner_valid <- droplevels(y_train_all[inner_valid_idx])
      metadata_inner_train <- metadata_train_all[inner_train_idx, , drop = FALSE]

      if (length(unique(y_inner_train)) < 2 || length(unique(y_inner_valid)) < 2) {
        skipped_folds <- c(skipped_folds, validation_study)
        next
      }

      norm_data <- normalize_train_test(
        X_train = X_inner_train,
        X_test = X_inner_valid,
        metadata_train = metadata_inner_train,
        threshold = grid_row$prevalence_threshold,
        jitter_amount = params$jitter_amount,
        min_hc_controls = params$min_hc_controls,
        top_n_features = params$top_n_features,
        min_feature_variance = params$min_feature_variance
      )
      
      message("Normalization complete.")
      
      if (is.null(norm_data)) {
        skipped_folds <- c(skipped_folds, validation_study)
        next
      }

      train_df <- data.frame(RF_Class = y_inner_train, norm_data$train_matrix, check.names = FALSE)
      valid_df <- data.frame(RF_Class = y_inner_valid, norm_data$test_matrix, check.names = FALSE)

      rf_model <- fit_rf_model(
        train_df = train_df,
        num_trees = params$tune_num_trees,
        max_depth = grid_row$max_depth,
        min_node_size = grid_row$min_node_size,
        mtry_fraction = grid_row$mtry_fraction,
        num_threads = params$num_threads
      )
      
      message("Model generation complete.")
      
      valid_probs <- predict(rf_model, data = valid_df[, -1, drop = FALSE])$predictions[, params$positive_class]

      oof_truth <- c(oof_truth, as.character(y_inner_valid))
      oof_probs <- c(oof_probs, valid_probs)
      fold_details[[validation_study]] <- data.frame(
        sample_id = rownames(valid_df),
        study_name = validation_study,
        truth = as.character(y_inner_valid),
        positive_prob = valid_probs
      )
    }

    if (length(oof_truth) == 0) {
      tuning_results[[grid_idx]] <- data.frame(
        grid_id = grid_idx,
        prevalence_threshold = grid_row$prevalence_threshold,
        mtry_fraction = grid_row$mtry_fraction,
        min_node_size = grid_row$min_node_size,
        max_depth = grid_row$max_depth,
        inner_f1 = NA_real_,
        inner_threshold = NA_real_,
        completed_inner_folds = 0
      )
      next
    }

    threshold_fit <- find_best_threshold(
      y_true = oof_truth,
      positive_probs = oof_probs,
      positive_class = params$positive_class,
      negative_class = params$negative_class
    )

    tuning_results[[grid_idx]] <- data.frame(
      grid_id = grid_idx,
      prevalence_threshold = grid_row$prevalence_threshold,
      mtry_fraction = grid_row$mtry_fraction,
      min_node_size = grid_row$min_node_size,
      max_depth = grid_row$max_depth,
      inner_f1 = threshold_fit$f1,
      inner_threshold = threshold_fit$threshold,
      completed_inner_folds = length(fold_details)
    )

    message("Added to tuning results list.")  }

  
  tuning_results_df <- bind_rows(tuning_results) %>%
    arrange(desc(inner_f1), desc(completed_inner_folds))

  best_row <- tuning_results_df %>%
    filter(!is.na(inner_f1)) %>%
    slice(1)

  if (nrow(best_row) == 0) {
    return(NULL)
  }

  best_row
}

## DATA PREP ###################################################################
combined_study_list <- lapply(CRC_progression_studies, build_combined_features)
common_features <- Reduce(intersect, lapply(combined_study_list, rownames))

combined_study_list_intersection <- lapply(combined_study_list, function(mat) {
  mat[common_features, , drop = FALSE]
})

combined_matrix <- do.call(cbind, combined_study_list_intersection)

metadata_list <- lapply(CRC_progression_studies, function(study) {
  as.data.frame(colData(study$relative_abundance))
})

metadata <- bind_rows(metadata_list)

metadata_filtered <- metadata[
  metadata$study_name != "HMP_2012" &
    !metadata$disease_class %in% c("CRC-M", "PA-M", "CRC-H"),
  ,
  drop = FALSE
]

combined_matrix_filtered <- combined_matrix[, rownames(metadata_filtered), drop = FALSE]

metadata_filtered$RF_Class <- dplyr::case_when(
  metadata_filtered$disease_class %in% c("HC", "Other") ~ params$negative_class,
  metadata_filtered$disease_class %in% c("CRC", "CRC+", "PA", "PA+") ~ params$positive_class,
  TRUE ~ NA_character_
)

metadata_filtered <- metadata_filtered[!is.na(metadata_filtered$RF_Class), , drop = FALSE]
metadata_filtered$RF_Class <- factor(
  metadata_filtered$RF_Class,
  levels = c(params$positive_class, params$negative_class)
)

combined_matrix_filtered <- combined_matrix_filtered[, rownames(metadata_filtered), drop = FALSE]

X <- t(combined_matrix_filtered)
y <- droplevels(metadata_filtered[rownames(X), "RF_Class"])

stopifnot(identical(rownames(X), rownames(metadata_filtered)))

message("Data prep complete.")

## TUNING GRID #################################################################
tune_grid <- expand.grid(
  prevalence_threshold = params$prevalence_grid,
  mtry_fraction = params$mtry_fraction_grid,
  min_node_size = params$min_node_size_grid,
  max_depth = params$max_depth_grid,
  stringsAsFactors = FALSE
)

## OUTER LODO ##################################################################
study_ids <- unique(metadata_filtered$study_name)

lodo_results <- list()
lodo_metrics <- list()
tuning_tables <- list()
cm_plots <- list()
roc_plots <- list()

for (heldout_study in study_ids) {
  cat("\nHolding out study:", heldout_study, "\n")

  outer_train_idx <- metadata_filtered$study_name != heldout_study
  outer_test_idx <- metadata_filtered$study_name == heldout_study

  X_outer_train <- X[outer_train_idx, , drop = FALSE]
  X_outer_test <- X[outer_test_idx, , drop = FALSE]
  y_outer_train <- droplevels(y[outer_train_idx])
  y_outer_test <- droplevels(y[outer_test_idx])
  metadata_outer_train <- metadata_filtered[outer_train_idx, , drop = FALSE]
  metadata_outer_test <- metadata_filtered[outer_test_idx, , drop = FALSE]

  if (length(unique(y_outer_train)) < 2 || length(unique(y_outer_test)) < 2) {
    cat("Skipping", heldout_study, "- outer fold lacks both classes\n")
    next
  }

  best_params <- run_inner_lodo(
    X_train_all = X_outer_train,
    y_train_all = y_outer_train,
    metadata_train_all = metadata_outer_train,
    tune_grid = tune_grid,
    params = params
  )

  if (is.null(best_params)) {
    cat("Skipping", heldout_study, "- no valid inner folds completed\n")
    next
  }

  tuning_tables[[heldout_study]] <- best_params

  norm_data <- normalize_train_test(
    X_train = X_outer_train,
    X_test = X_outer_test,
    metadata_train = metadata_outer_train,
    threshold = best_params$prevalence_threshold,
    jitter_amount = params$jitter_amount,
    min_hc_controls = params$min_hc_controls,
    top_n_features = params$top_n_features,
    min_feature_variance = params$min_feature_variance
  )

  if (is.null(norm_data)) {
    cat("Skipping", heldout_study, "- outer preprocessing failed\n")
    next
  }

  train_df <- data.frame(RF_Class = y_outer_train, norm_data$train_matrix, check.names = FALSE)
  test_df <- data.frame(RF_Class = y_outer_test, norm_data$test_matrix, check.names = FALSE)

  rf_model <- fit_rf_model(
    train_df = train_df,
    num_trees = params$final_num_trees,
    max_depth = best_params$max_depth,
    min_node_size = best_params$min_node_size,
    mtry_fraction = best_params$mtry_fraction,
    num_threads = params$num_threads
  )
  
  message("Model Generation Complete.")

  positive_probs <- predict(rf_model, data = test_df[, -1, drop = FALSE])$predictions[, params$positive_class]

  scored <- score_predictions(
    y_true = y_outer_test,
    positive_probs = positive_probs,
    threshold = best_params$inner_threshold,
    positive_class = params$positive_class,
    negative_class = params$negative_class
  )

  cm_df <- as.data.frame(scored$confusion_matrix$table) %>%
    group_by(Reference) %>%
    mutate(Percent = Freq / sum(Freq) * 100) %>%
    ungroup()

  roc_df <- data.frame(
    TPR = scored$roc_obj$sensitivities,
    FPR = 1 - scored$roc_obj$specificities
  )

  cm_plots[[heldout_study]] <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = paste0(Freq, "\n(", round(Percent, 1), "%)")), size = 4) +
    scale_fill_gradient(low = "white", high = "darkseagreen") +
    theme_bw() +
    labs(
      subtitle = paste("Held-out study:", heldout_study),
      x = "Actual",
      y = "Predicted"
    )

  roc_plots[[heldout_study]] <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(linewidth = 1, color = "steelblue4") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      subtitle = paste("Held-out study:", heldout_study),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )
  
  message("Plotting complete.")

  lodo_results[[heldout_study]] <- list(
    heldout_study = heldout_study,
    best_params = best_params,
    features_used = norm_data$keep_features,
    hc_controls = norm_data$control_count,
    model = rf_model,
    confusion_matrix = scored$confusion_matrix,
    pred_class = scored$pred_class,
    pred_probs = scored$pred_probs,
    y_test = y_outer_test,
    roc_obj = scored$roc_obj,
    summary = scored$summary
  )

  lodo_metrics[[heldout_study]] <- cbind(
    data.frame(
      heldout_study = heldout_study,
      n_test = length(y_outer_test),
      n_features = length(norm_data$keep_features),
      n_hc_controls = norm_data$control_count,
      stringsAsFactors = FALSE
    ),
    best_params,
    scored$summary
  )

  print(lodo_metrics[[heldout_study]])
}

lodo_metrics_df <- bind_rows(lodo_metrics)
tuning_summary_df <- bind_rows(tuning_tables, .id = "heldout_study")

overall_summary <- lodo_metrics_df %>%
  summarise(
    mean_f1 = mean(f1, na.rm = TRUE),
    median_f1 = median(f1, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE)
  )

print(lodo_metrics_df)
print(overall_summary)
