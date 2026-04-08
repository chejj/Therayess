## SETUP #######################################################################
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(ranger)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)
######PARAMETERS################################################################
params <- list(
  num_trees = 500,
  max_depth = 15,
  mtry_fraction = 0.01, # % of features considered at each split
  min_node_size = 20,   # minimum samples per terminal node
  num_threads = 8,      # Set to match HPC core
  prev_filter = 0.01,    # proportion filtered out in prevalence filtering step (i.e. 0.1 is 10%)
  train_split = 0.8     # proportion of the split to delegate for training (i.e. 0.8 is 80/20 split)
)

# Step 1: Create a matrix for each study ######################################
build_combined_features <- function(study) {
  # Extract matrices
  tax_mat  <- assay(study$relative_abundance, 1)
  pab_mat  <- assay(study$pathway_abundance, 1)
  pcov_mat <- assay(study$pathway_coverage, 1)
  
  # Prefix feature names so they stay unique
  rownames(tax_mat)  <- paste0("TAX_", rownames(tax_mat))
  rownames(pab_mat)  <- paste0("PAB_", rownames(pab_mat))
  rownames(pcov_mat) <- paste0("PCOV_", rownames(pcov_mat))
  
  # Keep only shared samples across the three data types
  common_samples <- Reduce(intersect, list(
    colnames(tax_mat),
    colnames(pab_mat),
    colnames(pcov_mat)
  ))
  
  tax_mat  <- tax_mat[, common_samples, drop = FALSE]
  pab_mat  <- pab_mat[, common_samples, drop = FALSE]
  pcov_mat <- pcov_mat[, common_samples, drop = FALSE]
  
  # Combine features vertically
  combined_mat <- rbind(tax_mat, pab_mat, pcov_mat)
  
  return(combined_mat)
}

# Apply to all studies 
combined_study_list <- lapply(CRC_progression_studies, build_combined_features)

# Step 2: Find all the shared features ########################################
common_features <- Reduce(intersect, lapply(combined_study_list, rownames))

length(common_features)


combined_study_list_intersection <- lapply(combined_study_list, function(mat) {
  mat[common_features, , drop = FALSE]
})

combined_matrix <- do.call(cbind, combined_study_list_intersection)
dim(combined_matrix)

# Build metadata table 
metadata_list <- lapply(CRC_progression_studies, function(study) {
  as.data.frame(colData(study$relative_abundance))
})

metadata <- bind_rows(metadata_list)
dim(metadata)

# optional alignment check
identical(colnames(combined_matrix), rownames(metadata))

# Step 3: REMOVE HMP-2012 and unwanted classes EARLY ##########################
metadata_filtered <- metadata[
  metadata$study_name != "HMP_2012" &
    !metadata$disease_class %in% c("CRC-M", "PA-M", "CRC-H"),
  ,
  drop = FALSE
]

table(metadata_filtered$disease_class)

# subset combined matrix BEFORE filtering/normalization
combined_matrix_filtered <- combined_matrix[, rownames(metadata_filtered), drop = FALSE]

# alignment check
identical(colnames(combined_matrix_filtered), rownames(metadata_filtered))
dim(combined_matrix_filtered)

# Step 4: Create final grouped labels #########################################
metadata_filtered$RF_Class <- dplyr::case_when(
  metadata_filtered$disease_class %in% c("HC", "Other") ~ "NEGATIVE_POLYP",
  metadata_filtered$disease_class %in% c("CRC", "CRC+", "CRC-M", "PA", "PA+", "PA-M") ~ "POSITIVE_POLYP",
  TRUE ~ NA_character_
) %>% as.factor() %>% factor(levels = c("POSITIVE_POLYP", "NEGATIVE_POLYP")  # 1st = positive class, 2nd = reference
)


y2 <- droplevels(metadata_filtered$RF_Class)
table(y2)

# Step 5: Fit 80/20 Random Forest ###################################################
#### A: Prep data for model ----
X <- t(combined_matrix_filtered)   # rows = samples, cols = features
y <- metadata_filtered[rownames(X), "RF_Class"] %>% droplevels() #hard locks labels to x row order

dim(X)
length(y)

###### i: Run train/test split ----
set.seed(42)

trainIndex <- createDataPartition(y, p = params$train_split, list = FALSE)

X_train <- X[trainIndex, , drop = FALSE]
X_test  <- X[-trainIndex, , drop = FALSE]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]

# keep metadata separate
metadata_train <- metadata_filtered[rownames(X_train), , drop = FALSE]
metadata_test  <- metadata_filtered[rownames(X_test), , drop = FALSE]

identical(rownames(X_train), rownames(metadata_train))
identical(rownames(X_test), rownames(metadata_test))

prop.table(table(droplevels(y_train)))
prop.table(table(droplevels(y_test)))

###### ii: Feature Transformations ----
######## A2: Prevalence filtering ----
prevalence_filter <- function(mat, threshold = params$prev_filter) {
  prev <- colSums(mat > 0) / nrow(mat)   # because features are columns
  names(prev[prev >= threshold])
}

keep_prev_features <- prevalence_filter(X_train, params$prev_filter)

filtered_TRAIN_matrix <- X_train[, keep_prev_features, drop = FALSE]
filtered_TEST_matrix  <- X_test[, keep_prev_features, drop = FALSE]

dim(filtered_TRAIN_matrix)
dim(filtered_TEST_matrix)

# check metadata still aligned
identical(rownames(filtered_TRAIN_matrix), rownames(metadata_train))
identical(rownames(filtered_TEST_matrix), rownames(metadata_test))


######## B2: Percentile Normalization ----
# Add tiny jitter before percentile normalization
set.seed(42)

filtered_TRAIN_matrix <- filtered_TRAIN_matrix + matrix(
  runif(length(filtered_TRAIN_matrix), min = 0, max = 1e-9),
  nrow = nrow(filtered_TRAIN_matrix),
  ncol = ncol(filtered_TRAIN_matrix),
  dimnames = dimnames(filtered_TRAIN_matrix)
)

filtered_TEST_matrix <- filtered_TEST_matrix + matrix(
  runif(length(filtered_TEST_matrix), min = 0, max = 1e-9),
  nrow = nrow(filtered_TEST_matrix),
  ncol = ncol(filtered_TEST_matrix),
  dimnames = dimnames(filtered_TEST_matrix)
)

# make sure metadata matches filtered_matrix
set.seed(42)
identical(rownames(filtered_TRAIN_matrix), rownames(metadata_train))
identical(rownames(filtered_TEST_matrix), rownames(metadata_test))

# percentile normalization
percentile_norm <- function(values, controls) {
  if (length(controls) == 0) {
    return(rep(NA_real_, length(values)))
  }
  
  sapply(values, function(x) {
    mean(controls <= x) * 100
  })
}
# define control samples from train only (HC-Normalization)
control_idx <- which(metadata_train$disease_class == "HC")
# extract control matrix
control_matrix <- filtered_TRAIN_matrix[control_idx, , drop = FALSE]
# sanity check! 
length(control_idx)
# apply percentile normalization to train matrix
pnorm_train <- apply(filtered_TRAIN_matrix, 2, function(feature_values) {
  controls <- feature_values[control_idx]
  percentile_norm(feature_values, controls)
})
# check
dim(pnorm_train)
summary(as.vector(pnorm_train))
sum(is.na(pnorm_train))

# repeat with test data cautiously!(make sure using train controls)
pnorm_test <- sapply(seq_len(ncol(filtered_TEST_matrix)), function(j) {
  test_values <- filtered_TEST_matrix[, j]
  train_controls <- filtered_TRAIN_matrix[control_idx, j]
  percentile_norm(test_values, train_controls)
})

# fix names 
pnorm_test <- as.matrix(pnorm_test)
colnames(pnorm_test) <- colnames(filtered_TEST_matrix)
rownames(pnorm_test) <- rownames(filtered_TEST_matrix)

# check it 
dim(pnorm_test)
summary(as.vector(pnorm_test))
sum(is.na(pnorm_test))

#### B: Set Up ----
train_df <- data.frame(RF_Class = y_train, pnorm_train, check.names = FALSE)
test_df  <- data.frame(RF_Class = y_test, pnorm_test, check.names = FALSE)

# Class weights from training fold only
class_counts <- table(y_train)
class_weights <- sum(class_counts) / (length(class_counts) * class_counts) # Keep as table to keep them "named", ranger prefers this way
# Do not convert class_weights to numeric using "as.numeric", table is already seen as a numeric in the environment panel
names(class_weights) <- names(class_counts)

print(class_counts)
print(class_weights)

# Number of predictor features (exclude outcome column)
p <- ncol(train_df) - 1

# Set mtry to X% of total predictors, or 1, whichever is larger
mtry_val <- max(1, ceiling(params$mtry_fraction * p))

#### D: Model Generation ----
rf_model <- ranger(
  dependent.variable.name = "RF_Class", # outcome modeled by all predictors
  data = train_df,                      # training data only
  num.trees = params$num_trees,         
  max.depth = params$max_depth,
  class.weights = class_weights,
  mtry = mtry_val,                      
  min.node.size = params$min_node_size, 
  probability = TRUE,        # needed for multiclass probabilities / ROC
  importance = "impurity",   # feature importance
#  splitrule = "gini",        # closest standard classification impurity rule in ranger, entropy is not available in ranger
  num.threads = params$num_threads
)

#### E. Model Output ----
rf_model

#### F: Feature Importance (check which are most influential in splits, discussion/conclusion related?) ----

# Built In
# sort(importance(rf_model), decreasing = TRUE)[1:15]

# Count Splits
# split_counts <- unlist(
#   lapply(1:rf_model$num.trees, function(t) {
#     treeInfo(rf_model, tree = t)$splitvarName
#   })
# )
# head(sort(table(split_counts, useNA = "no"), decreasing = TRUE), 15)

# Step 7: Check against test set ##############################################
# explicitly drop the classifier so it doesn't get used in prediction as a safeguard
pred_probs <- predict(rf_model, data = test_df[, colnames(test_df) != "RF_Class"])$predictions 
pred_class <- colnames(pred_probs)[max.col(pred_probs)]
pred_class <- factor(pred_class, levels = levels(train_df$RF_Class))

#### A: Confusion Matrix ----
cm <- confusionMatrix(pred_class, y_test)

cm_df <- as.data.frame(cm$table)

# Add percentages (row-wise = % of actual class)
cm_df <- cm_df %>%
  group_by(Reference) %>%
  mutate(Percent = Freq / sum(Freq) * 100)

ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = paste0(Freq, "\n(", round(Percent, 1), "%)")), size = 4) +
  scale_fill_gradient(low = "white", high = "darkseagreen") +
  theme_bw() +
  labs(
    title = "Confusion Matrix: Ranger forest",
    x = "Actual",
    y = "Predicted", 
    caption = "Percentages sum to 100% vertically for the actual classifier."
  )
# Save summary metrics
byclass <- cm$byClass

results_summary <- data.frame(
  accuracy = as.numeric(cm$overall["Accuracy"]),
  kappa = as.numeric(cm$overall["Kappa"]),
  sensitivity = as.numeric(byclass["Sensitivity"]),
  specificity = as.numeric(byclass["Specificity"]),
  pos_pred_value = as.numeric(byclass["Pos Pred Value"]),
  neg_pred_value = as.numeric(byclass["Neg Pred Value"]),
  balanced_accuracy = as.numeric(byclass["Balanced Accuracy"])
)

# Save everything together
model_results <- list(
  confusion_matrix = cm,
  pred_probs = pred_probs,
  pred_class = pred_class,
  y_test = y_test,
  class_weights = class_weights,
  summary = results_summary
)

# print results plz 

print(cm)
print(results_summary)

#### B: F1 Score ----
if (is.matrix(cm$byClass)) {
  f1_df <- data.frame(
    Class = rownames(cm$byClass),
    F1 = cm$byClass[, "F1"],
    row.names = NULL
  )
} else {
  f1_df <- data.frame(
    Class = "Overall_or_Binary",
    F1 = unname(cm$byClass["F1"])
  )
}

f1_df

#### C: AUC-ROC Curve ----
# Set up 
roc_list <- list()       # for AUC
roc_df_list <- list()    # for plotting

for (class in colnames(pred_probs)) {
  
  binary_truth <- ifelse(test_df$RF_Class == class, 1, 0)
  
  # Skip if class not present
  if (sum(binary_truth) == 0 || sum(binary_truth) == length(binary_truth)) next
  
  roc_obj <- roc(binary_truth, pred_probs[, class])
  
  # Store ROC object
  roc_list[[class]] <- roc_obj
  
  # Store data for ggplot
  roc_df_list[[class]] <- data.frame(
    TPR = roc_obj$sensitivities,
    FPR = 1 - roc_obj$specificities,
    Class = class
  )
}

# Combine for ggplot
roc_df_all <- dplyr::bind_rows(roc_df_list)

# Plot the ROC Curves 
ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # random baseline
  theme_bw() +
  labs(
    title = "ROC Curves (One-vs-Rest)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) 

# Extract AUC Values
auc_values <- sapply(roc_list, function(x) auc(x))

auc_df <- data.frame(
  Class = names(auc_values),
  AUC = as.numeric(auc_values)
)

auc_df

#### D: Overfitting check: compare train vs test accuracy ----
# Get predicted probabilities  and convert to predicted class
train_pred_probs <- predict(rf_model, data = train_df)$predictions
train_pred_class <- colnames(train_pred_probs)[max.col(train_pred_probs)]

test_pred_probs <- predict(rf_model, data = test_df)$predictions
test_pred_class <- colnames(test_pred_probs)[max.col(test_pred_probs)]

# Compute accuracies
train_accuracy <- mean(train_pred_class == train_df$RF_Class)
test_accuracy  <- mean(test_pred_class == test_df$RF_Class)

# Build plotting data frame
overfit_df <- data.frame(
  Set = c("Train", "Test"),
  Accuracy = c(train_accuracy, test_accuracy)
)

# Plot
overfit_df
ggplot(overfit_df, aes(x = Set, y = Accuracy, fill = Set)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Overfitting Check") +
  theme_minimal()

# Step 6: Fit LODO RF model ####################################################
### A: Set up LODO storage ----
study_ids <- unique(metadata_filtered$study_name)

lodo_results <- list()
cm_plots <- list()
auc_plots <- list()

lodo_metrics <- data.frame(
  heldout_study = character(),
  n_test = integer(),
  accuracy = numeric(),
  kappa = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  pos_pred_value = numeric(),
  neg_pred_value = numeric(),
  balanced_accuracy = numeric(),
  stringsAsFactors = FALSE
)
### B: LODO loop ----
for (heldout_study in study_ids) {
  
  cat("\nHolding out study:", heldout_study, "\n")
  
  train_idx <- metadata_filtered$study_name != heldout_study
  test_idx  <- metadata_filtered$study_name == heldout_study
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  metadata_train <- metadata_filtered[train_idx, , drop = FALSE]
  metadata_test  <- metadata_filtered[test_idx, , drop = FALSE]
  
  # quick skip if held-out study has only one class
  if (length(unique(y_test)) < 2) {
    cat("Skipping", heldout_study, "- test set has only one class\n")
    next
  }
  # Skip if training fold has only one class
  if (length(unique(y_train)) < 2) {
    cat("Skipping", heldout_study, "- training set has only one class\n")
    next
  }
  #### I: Prep data for model ----
  ######## A2: Prevalence filtering ----
  keep_prev_features <- prevalence_filter(X_train, params$prev_filter)
  
  filtered_TRAIN_matrix <- X_train[, keep_prev_features, drop = FALSE]
  filtered_TEST_matrix  <- X_test[, keep_prev_features, drop = FALSE]
  
  dim(filtered_TRAIN_matrix)
  dim(filtered_TEST_matrix)
  
  # check metadata still aligned
  identical(rownames(filtered_TRAIN_matrix), rownames(metadata_train))
  identical(rownames(filtered_TEST_matrix), rownames(metadata_test))
  
  ######## B2: Percentile Normalization ----
  # Add tiny jitter before percentile normalization
  set.seed(42)
  
  filtered_TRAIN_matrix <- filtered_TRAIN_matrix + matrix(
    runif(length(filtered_TRAIN_matrix), min = 0, max = 1e-9),
    nrow = nrow(filtered_TRAIN_matrix),
    ncol = ncol(filtered_TRAIN_matrix),
    dimnames = dimnames(filtered_TRAIN_matrix)
  )
  
  filtered_TEST_matrix <- filtered_TEST_matrix + matrix(
    runif(length(filtered_TEST_matrix), min = 0, max = 1e-9),
    nrow = nrow(filtered_TEST_matrix),
    ncol = ncol(filtered_TEST_matrix),
    dimnames = dimnames(filtered_TEST_matrix)
  )
  
  # make sure metadata matches filtered_matrix
  set.seed(42)
  if (
    identical(rownames(filtered_TRAIN_matrix), rownames(metadata_train)) &&
    identical(rownames(filtered_TEST_matrix), rownames(metadata_test))
  ) {
    message("✅ Metadata predictors aligned with train/test samples.")
  } else {
    stop("❌ Metadata alignment issues.")
  }
  
  # define control samples from train only (HC-Normalization)
  control_idx <- which(metadata_train$disease_class == "HC")
  
  # extract control matrix
  control_matrix <- filtered_TRAIN_matrix[control_idx, , drop = FALSE]
  
  # apply percentile normalization to train matrix
  pnorm_train <- apply(filtered_TRAIN_matrix, 2, function(feature_values) {
    controls <- feature_values[control_idx]
    percentile_norm(feature_values, controls)
  })
  
  message("Training set percentile normalization complete.")
  
  # repeat with test data cautiously!(make sure using train controls)
  pnorm_test <- sapply(seq_len(ncol(filtered_TEST_matrix)), function(j) {
    test_values <- filtered_TEST_matrix[, j]
    train_controls <- filtered_TRAIN_matrix[control_idx, j]
    percentile_norm(test_values, train_controls)
  })
  
  message("Test set percentile normalization complete.")
  
  # fix names 
  pnorm_test <- as.matrix(pnorm_test)
  colnames(pnorm_test) <- colnames(filtered_TEST_matrix)
  rownames(pnorm_test) <- rownames(filtered_TEST_matrix)
  
  #### B: Run the ranger RF model ----
  train_df <- data.frame(RF_Class = y_train, pnorm_train, check.names = FALSE)
  test_df  <- data.frame(RF_Class = y_test, pnorm_test, check.names = FALSE)
  
  # Class weights from training fold only
  class_counts <- table(y_train)
  class_weights <- sum(class_counts) / (length(class_counts) * class_counts) # Keep as table to keep them "named", ranger prefers this way
  # Do not convert class_weights to numeric using "as.numeric", table is already seen as a numeric in the environment panel
  names(class_weights) <- names(class_counts)
  
  print(class_counts)
  print(class_weights)
  
  # Model
  rf_model_LODO <- ranger(
    dependent.variable.name = "RF_Class",
    data = train_df,
    num.trees = params$num_trees,
    max.depth = params$max_depth,
    mtry = max(1, floor(params$mtry_fraction * (ncol(train_df) - 1))), #number of predictor features (minus outcome) * %mtry_fraction
    min.node.size = params$min_node_size,
    probability = TRUE,
    importance = "impurity",
    class.weights = class_weights,
    num.threads = params$num_threads
  )
  
  message("Random Forest generation complete")
  
  #### C: Predict on test set ----
  pred_probs <- predict(rf_model_LODO, data = test_df[, -1])$predictions
  pred_class <- colnames(pred_probs)[max.col(pred_probs)]
  pred_class <- factor(pred_class, levels = levels(train_df$RF_Class))
  
  #### D: Evaluate ----
  ##### i: Confusion Matrix ----
  cm <- confusionMatrix(pred_class, y_test)
  
  byclass <- cm$byClass
  
  results_summary <- data.frame(
    accuracy = as.numeric(cm$overall["Accuracy"]),
    kappa = as.numeric(cm$overall["Kappa"]),
    sensitivity = as.numeric(byclass["Sensitivity"]),
    specificity = as.numeric(byclass["Specificity"]),
    pos_pred_value = as.numeric(byclass["Pos Pred Value"]),
    neg_pred_value = as.numeric(byclass["Neg Pred Value"]),
    balanced_accuracy = as.numeric(byclass["Balanced Accuracy"])
  )
  
  ##### ii. Per-Class F1 Score ----
  if (is.matrix(cm$byClass)) {
    f1_df <- data.frame(
      Class = rownames(cm$byClass),
      F1 = cm$byClass[, "F1"],
      row.names = NULL
    )
  } else {
    f1_df <- data.frame(
      Class = "Overall_or_Binary",
      F1 = unname(cm$byClass["F1"])
    )
  }
  
  ##### iii. AUC-ROC Curve ----
  roc_list <- list()
  roc_df_list <- list()
  
  for (class in colnames(pred_probs)) {
    binary_truth <- ifelse(test_df$RF_Class == class, 1, 0)
    
    if (sum(binary_truth) == 0 || sum(binary_truth) == length(binary_truth)) next
    
    roc_obj <- roc(binary_truth, pred_probs[, class])
    roc_list[[class]] <- roc_obj
    
    roc_df_list[[class]] <- data.frame(
      TPR = roc_obj$sensitivities,
      FPR = 1 - roc_obj$specificities,
      Class = class
    )
  }
  
  roc_df_all <- dplyr::bind_rows(roc_df_list)
  
  auc_values <- sapply(roc_list, function(x) auc(x))
  auc_df <- data.frame(
    Class = names(auc_values),
    AUC = as.numeric(auc_values)
  )
  
  message("Evaluation completed. Now saving results.")
  
  #### E: Save everything for this held-out study ----
  lodo_results[[heldout_study]] <- list(
    heldout_study = heldout_study,
    train_study_names = unique(metadata_train$study_name),
    class_counts = class_counts,
    class_weights = class_weights,
    confusion_matrix = cm,
    confusion_matrix_table = cm$table,
    pred_probs = pred_probs,
    pred_class = pred_class,
    y_test = y_test,
    f1_df = f1_df,
    roc_list = roc_list,
    roc_df_all = roc_df_all,
    auc_df = auc_df,
    summary = results_summary
  )
  
  lodo_metrics <- rbind(
    lodo_metrics,
    data.frame(
      heldout_study = heldout_study,
      n_test = length(y_test),
      accuracy = results_summary$accuracy,
      kappa = results_summary$kappa,
      sensitivity = results_summary$sensitivity,
      specificity = results_summary$specificity,
      pos_pred_value = results_summary$pos_pred_value,
      neg_pred_value = results_summary$neg_pred_value,
      balanced_accuracy = results_summary$balanced_accuracy,
      stringsAsFactors = FALSE
    )
  )
  
  ##### PLOTS: ----
  cm_df_plot <- as.data.frame(cm$table)
  
  cm_df_plot <- cm_df_plot %>%
    dplyr::group_by(Reference) %>%
    dplyr::mutate(Percent = Freq / sum(Freq) * 100) %>%
    dplyr::ungroup()
  
  cm_plots[[heldout_study]] <- ggplot(cm_df_plot, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = paste0(Freq, "\n(", round(Percent, 1), "%)")), size = 4) +
    scale_fill_gradient(low = "white", high = "darkseagreen") +
    theme_bw() +
    labs(
      #      title = "Confusion Matrix",
      subtitle = paste("Held-out study:", heldout_study),
      x = "Actual",
      y = "Predicted"
    )
  
  auc_plots[[heldout_study]] <- ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Class)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # random baseline
    theme_bw() +
    labs(
      #      title = "ROC Curve",
      subtitle = paste("Held-out study:", heldout_study),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )
}