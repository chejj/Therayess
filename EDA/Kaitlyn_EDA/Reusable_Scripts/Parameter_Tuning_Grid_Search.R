## Parameters ##################################################################
params <- list(
  num_trees = c(1000, 3000, 5000),
  max_depth = c(8, 10, 15, 20),
  mtry_fraction = c(0.01, 0.03, 0.05), # % of features considered at each split
  min_node_size = c(10, 15, 20, 25),   # minimum samples per terminal node
  num_threads = 8,      # Set to match HPC core
  splitrule = "gini",   # options: "gini", "extratrees" or "hellinger"
  prev_filter = 0.01,    # proportion filtered out in prevalence filtering step (i.e. 0.1 is 10%)
  train_split = 0.8     # proportion of the split to delegate for training (i.e. 0.8 is 80/20 split)
)

## TUNING GRID #################################################################
tune_grid <- expand.grid(
  num_trees = params$num_trees,
  mtry_fraction = params$mtry_fraction,
  min_node_size = params$min_node_size,
  max_depth = params$max_depth,
  stringsAsFactors = FALSE
)

## Helper ######################################################################

library(ranger)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

## Small bracket around the best region
tune_grid <- expand.grid(
  num_trees      = c(1000, 3000, 5000),
  max_depth      = c(8, 10, 12),
  mtry_fraction  = c(0.02, 0.03, 0.05),
  min_node_size  = c(10, 15, 20),
  stringsAsFactors = FALSE
)

run_ranger_grid <- function(train_df, test_df, tune_grid,
                            positive_class = "POSITIVE_POLYP",
                            splitrule = "gini",
                            num_threads = 8,
                            seed = 123) {
  
  # Class weights from training fold only
  class_counts <- table(train_df$RF_Class)
  class_weights <- sum(class_counts) / (length(class_counts) * class_counts)
  names(class_weights) <- names(class_counts)
  
  p <- ncol(train_df) - 1  # predictors only
  
  results <- vector("list", nrow(tune_grid))
  models  <- vector("list", nrow(tune_grid))
  plots   <- vector("list", nrow(tune_grid))
  
  for (i in seq_len(nrow(tune_grid))) {
    g <- tune_grid[i, ]
    
    set.seed(seed + i)
    
    mtry_val <- max(1, ceiling(g$mtry_fraction * p))
    
    rf_model <- ranger(
      dependent.variable.name = "RF_Class",
      data = train_df,
      num.trees = g$num_trees,
      max.depth = g$max_depth,
      class.weights = class_weights,
      mtry = mtry_val,
      min.node.size = g$min_node_size,
      probability = TRUE,
      importance = "impurity",
      splitrule = splitrule,
      num.threads = num_threads
    )
    
    pred_obj <- predict(rf_model, data = test_df)
    prob_mat  <- pred_obj$predictions
    pred_class <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
    pred_class <- factor(pred_class, levels = levels(test_df$RF_Class))
    truth      <- factor(test_df$RF_Class, levels = levels(test_df$RF_Class))
    
    cm <- caret::confusionMatrix(
      data = pred_class,
      reference = truth,
      positive = positive_class
    )
    
    # Positive-class probability
    pos_prob <- prob_mat[, positive_class]
    
    roc_obj <- pROC::roc(
      response = truth,
      predictor = pos_prob,
      levels = rev(levels(truth)),
      quiet = TRUE
    )
    
    # Manual F1 for the positive class
    tp <- sum(pred_class == positive_class & truth == positive_class)
    fp <- sum(pred_class == positive_class & truth != positive_class)
    fn <- sum(pred_class != positive_class & truth == positive_class)
    
    precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
    recall    <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
    f1 <- ifelse(is.na(precision) || is.na(recall) || (precision + recall) == 0,
                 NA, 2 * precision * recall / (precision + recall))
    
    results[[i]] <- data.frame(
      num_trees = g$num_trees,
      max_depth = g$max_depth,
      mtry_fraction = g$mtry_fraction,
      mtry = mtry_val,
      min_node_size = g$min_node_size,
      oob_error = rf_model$prediction.error,
      accuracy = unname(cm$overall["Accuracy"]),
      kappa = unname(cm$overall["Kappa"]),
      sensitivity = unname(cm$byClass["Sensitivity"]),
      specificity = unname(cm$byClass["Specificity"]),
      balanced_accuracy = unname(cm$byClass["Balanced Accuracy"]),
      f1 = f1,
      auc = as.numeric(pROC::auc(roc_obj)),
      stringsAsFactors = FALSE
    )
    
    models[[i]] <- list(
      model = rf_model,
      cm = cm,
      roc = roc_obj,
      pred_class = pred_class,
      prob_mat = prob_mat,
      truth = truth,
      class_counts = class_counts,
      class_weights = class_weights
    )
    
    cm_df_plot <- as.data.frame(cm$table) %>%
      dplyr::group_by(Reference) %>%
      dplyr::mutate(Percent = Freq / sum(Freq) * 100) %>%
      dplyr::ungroup()
    
    cm_plot <- ggplot(cm_df_plot, aes(x = Reference, y = Prediction, fill = Freq)) +
      geom_tile() +
      geom_text(aes(label = paste0(Freq, "\n(", round(Percent, 1), "%)")), size = 4) +
      scale_fill_gradient(low = "white", high = "darkseagreen") +
      theme_bw() +
      labs(
        subtitle = paste(
          "trees =", g$num_trees,
          "| depth =", g$max_depth,
          "| mtry =", mtry_val,
          "| min.node.size =", g$min_node_size
        ),
        x = "Actual",
        y = "Predicted"
      )
    
    roc_df <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    )
    
    roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line(linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      theme_bw() +
      labs(
        subtitle = paste(
          "trees =", g$num_trees,
          "| depth =", g$max_depth,
          "| mtry =", mtry_val,
          "| min.node.size =", g$min_node_size
        ),
        x = "False Positive Rate",
        y = "True Positive Rate"
      )
    
    plots[[i]] <- list(cm_plot = cm_plot, roc_plot = roc_plot)
  }
  
  metrics_df <- bind_rows(results) %>%
    arrange(desc(balanced_accuracy), desc(f1), desc(auc))
  
  list(
    metrics = metrics_df,
    models = models,
    plots = plots
  )
}

grid_out <- run_ranger_grid(
  train_df = train_df,
  test_df = test_df,
  tune_grid = tune_grid,
  positive_class = "POSITIVE_POLYP",
  splitrule = "gini",
  num_threads = params$num_threads,
  seed = 123
)

# Top parameter sets
grid_out$metrics %>% 
  arrange(-sensitivity) %>% 
  head(10)

best_idx <- which.max(grid_out$metrics$balanced_accuracy)
best_row <- grid_out$metrics[best_idx, ]
best_model <- grid_out$models[[best_idx]]$model

best_row
best_model