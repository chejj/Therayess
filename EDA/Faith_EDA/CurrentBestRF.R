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
  mtry_fraction = 0.03, # % of features considered at each split
  min_node_size = 20,   # minimum samples per terminal node
  num_threads = 8       # Set to match HPC core
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

# Step 4: Feature Transformations #############################################
#### A: Prevalence filtering ----
prevalence_filter <- function(mat, threshold = 0.1) {
  prev <- rowSums(mat > 0) / ncol(mat)
  mat[prev >= threshold, , drop = FALSE]
}

filtered_matrix <- prevalence_filter(combined_matrix_filtered, 0.1)
dim(filtered_matrix)

#### B: Percentile Normalization ----
# Add tiny jitter to reduce ties for percentile normalization
set.seed(42)
filtered_matrix <- filtered_matrix + matrix(
  runif(length(filtered_matrix), min = 0, max = 1e-9),
  nrow = nrow(filtered_matrix),
  ncol = ncol(filtered_matrix),
  dimnames = dimnames(filtered_matrix)
)

# make sure metadata matches filtered_matrix
metadata_filtered <- metadata_filtered[colnames(filtered_matrix), , drop = FALSE]

percentile_norm <- function(feature_values, group, control_label = "HC") {
  controls <- feature_values[group == control_label]
  
  if (length(controls) == 0) {
    return(rep(NA_real_, length(feature_values)))
  }
  
  sapply(feature_values, function(x) {
    mean(controls <= x) * 100
  })
}

batch <- metadata_filtered$study_name
group <- metadata_filtered$disease_class
control_label <- "HC"

percentile_matrix <- matrix(
  NA_real_,
  nrow = nrow(filtered_matrix),
  ncol = ncol(filtered_matrix),
  dimnames = dimnames(filtered_matrix)
)

for (b in unique(batch)) {
  idx <- which(batch == b)
  
  sub_data <- filtered_matrix[, idx, drop = FALSE]
  sub_group <- group[idx]
  
  if (!any(sub_group == control_label)) next
  
  sub_percentile <- apply(sub_data, 1, function(feature) {
    percentile_norm(feature, sub_group, control_label)
  })
  
  percentile_matrix[, idx] <- t(sub_percentile)
}

# check result
dim(percentile_matrix)
summary(as.vector(percentile_matrix))
sum(is.na(percentile_matrix))

# Step 5: Create final grouped labels #########################################

metadata_filtered$RF_Class <- dplyr::case_when(
  metadata_filtered$disease_class %in% c("HC", "Other") ~ "NEGATIVE_POLYP",
  metadata_filtered$disease_class %in% c("CRC", "CRC+", "CRC-M", "PA", "PA+", "PA-M") ~ "POSITIVE_POLYP",
  TRUE ~ NA_character_
) %>% as.factor() %>% factor(levels = c("POSITIVE_POLYP", "NEGATIVE_POLYP")  # 1st = positive class, 2nd = reference
)


y2 <- droplevels(metadata_filtered$RF_Class)
table(y2)

# Step 6: Fit Random Forest ###################################################
#### A: Prep data for model ----
X <- t(percentile_matrix)
y <- y2

dim(X)
length(y)

#### B: Run test train split ----
set.seed(42)

trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]

prop.table(table(droplevels(y_train)))
prop.table(table(droplevels(y_test)))

## QC Check: Original CRC Progression Class Splits
disease_train <- metadata_filtered$disease_class[trainIndex]
disease_test  <- metadata_filtered$disease_class[-trainIndex]

# Check counts 
table(droplevels(disease_train))
table(droplevels(disease_test))

# Check Proportions (could make a ggplot for this to show balance?)
prop.table(table(droplevels(disease_train)))
prop.table(table(droplevels(disease_test)))

#### C: Convert to data.frame and bind y (RF_Class) (ranger prefers formula or full data) ----
train_df <- data.frame(RF_Class = y_train, X_train)
test_df  <- data.frame(RF_Class = y_test, X_test)

# Check and prep for class imbalance
class_counts <- table(train_df$RF_Class)
class_weights <- sum(class_counts) / (length(class_counts) * class_counts)

# Number of predictor features (exclude outcome column)
p <- ncol(train_df) - 1

# Set mtry to X% of total predictors, or 1, whichever is larger
mtry_val <- max(1, ceiling(params$mtry_fraction * p))

#### D: Model Generation ----
PNrf_model <- ranger(
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
PNrf_model

#### F: Feature Importance (check which are most influential in splits, discussion/conclusion related?) ----

# Built In
# sort(importance(PNrf_model), decreasing = TRUE)[1:15]

# Count Splits
# split_counts <- unlist(
#   lapply(1:PNrf_model$num.trees, function(t) {
#     treeInfo(PNrf_model, tree = t)$splitvarName
#   })
# )
# head(sort(table(split_counts, useNA = "no"), decreasing = TRUE), 15)

# Step 7: Check against test set ##############################################

# explicitly drop the classifier so it doesn't get used in prediction as a safeguard
pred_probs <- predict(PNrf_model, data = test_df[, colnames(test_df) != "RF_Class"])$predictions 

pred_class <- colnames(pred_probs)[max.col(pred_probs)]

#### A. Confusion Matrix ----

cm <- confusionMatrix(
  factor(pred_class, levels = levels(train_df$RF_Class)),
  factor(test_df$RF_Class, levels = levels(train_df$RF_Class))
)
cm

cm_df <- as.data.frame(cm$table)
ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() +
  labs(title = "Confusion Matrix: Random Forest")

#### B. Per-Class F1 Score ----
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

ggplot(f1_df, aes(x = Class, y = F1)) +
  geom_col() +
  theme_bw() +
  labs(title = "Per-Class F1 Score", y = "F1 Score")

#### C. AUC-ROC Curve (Multiclass) ----
# i. Set up ---
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

# ii. Plot the ROC Curves ---
ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # random baseline
  theme_bw() +
  labs(
    title = "ROC Curves (One-vs-Rest)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) 

# iii. Extract AUC Values ---
auc_values <- sapply(roc_list, function(x) auc(x))

auc_df <- data.frame(
  Class = names(auc_values),
  AUC = as.numeric(auc_values)
)

auc_df
ggplot(auc_df, aes(x = Class, y = AUC)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(title = "AUC by Class")

#### D. Overfitting check: compare train vs test accuracy ----
# i. Get predicted probabilities  and convert to predicted class
train_pred_probs <- predict(PNrf_model, data = train_df)$predictions
train_pred_class <- colnames(train_pred_probs)[max.col(train_pred_probs)]

test_pred_probs <- predict(PNrf_model, data = test_df)$predictions
test_pred_class <- colnames(test_pred_probs)[max.col(test_pred_probs)]

# ii. Compute accuracies
train_accuracy <- mean(train_pred_class == train_df$RF_Class)
test_accuracy  <- mean(test_pred_class == test_df$RF_Class)

# iii. Build plotting data frame
overfit_df <- data.frame(
  Set = c("Train", "Test"),
  Accuracy = c(train_accuracy, test_accuracy)
)

# iv. Plot
overfit_df
ggplot(overfit_df, aes(x = Set, y = Accuracy, fill = Set)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Overfitting Check") +
  theme_minimal()

# Step 8: LODO for Cross-Validation ############################################
studies <- unique(metadata_filtered$study_name)

set.seed(42)
lodo_preds <- lapply(studies, function(test_study) {
  
  test_idx  <- which(metadata_filtered$study_name == test_study)
  train_idx <- which(metadata_filtered$study_name != test_study)
  
  X_train_lodo <- X[train_idx, , drop = FALSE]
  X_test_lodo  <- X[test_idx, , drop = FALSE]
  
  y_train_lodo <- y[train_idx]
  y_test_lodo  <- y[test_idx]
  
  train_df <- data.frame(RF_Class = y_train_lodo, X_train_lodo)
  test_df  <- data.frame(X_test_lodo)
  
  ## OPTION A: DOWN SAMPLING ----
  # Skip if fewer than 2 classes remain in training
  class_sizes_lodo <- table(train_df$RF_Class)
  if (length(class_sizes_lodo) < 2) {
    return(NULL)
  }
  
  min_class_lodo <- min(class_sizes_lodo)
  
  balanced_train_df <- train_df %>%
    dplyr::group_by(RF_Class) %>%
    dplyr::slice_sample(n = min_class_lodo) %>%
    dplyr::ungroup()
  
  balanced_train_df$RF_Class <- droplevels(as.factor(balanced_train_df$RF_Class))
  
  # ## OPTION B: CLASS WEIGHTING (Like in 80/20 split) ----
  # # Check and prep for class imbalance
  # class_counts <- table(train_df$RF_Class)
  # class_weights <- sum(class_counts) / (length(class_counts) * class_counts)
  # 
  # # Number of predictor features (exclude outcome column)
  # p <- ncol(train_df) - 1
  # 
  # # Set mtry to X% of total predictors, or 1, whichever is larger
  # mtry_val <- max(1, ceiling(params$mtry_fraction * p))
  # balanced_train_df <- train_df
  
  rg_lodo <- ranger(
    dependent.variable.name = "RF_Class", # outcome modeled by all predictors
    data = balanced_train_df,                      # training data only
    num.trees = params$num_trees,         
    max.depth = params$max_depth,
    class.weights = class_weights,
#    mtry = mtry_val,          # See option B above            
    min.node.size = params$min_node_size, 
    probability = FALSE,        
    importance = "impurity",   # feature importance (Can see if country is still the most important?)
    #  splitrule = "gini",        # closest standard classification impurity rule in ranger, entropy is not available in ranger
    num.threads = params$num_threads
  )
  
  pred_lodo <- predict(rg_lodo, data = test_df)$predictions
  
  data.frame(
    Study = test_study,
    Truth = y_test_lodo,
    Pred = pred_lodo
  )
})

lodo_preds <- do.call(rbind, lodo_preds)

#### A: Create per study metrics ----
lodo_results2 <- lodo_preds %>%
  group_by(Study) %>%
  do({
    cm <- confusionMatrix(.$Pred, .$Truth)
    
    data.frame(
      N_test = nrow(.),
      Accuracy = as.numeric(cm$overall["Accuracy"]),
      P_Value = as.numeric(cm$overall["AccuracyPValue"]),
      Kappa = as.numeric(cm$overall["Kappa"]), 
      Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
      Specificity = as.numeric(cm$byClass["Specificity"]),
      Pos_Pred_Value = as.numeric(cm$byClass["Pos Pred Value"]),
      Neg_Pred_Value = as.numeric(cm$byClass["Neg Pred Value"])
    )
  }) %>%
  ungroup()

lodo_results2

#### B: Get the LODO means ----
### MEANS MIGHT NOT BE APPROPRIATE?? MIGHT WANT TO DO A WEIGHTED MEAN DO TO SAMPLE SIZE DISCREPANCIES?? Need to look into it more perhaps
### Would range be a better measure here to show the spread rather than the average?
mean(lodo_results2$Accuracy) 
mean(lodo_results2$Kappa)
mean(lodo_results2$Sensitivity)
mean(lodo_results2$Pos_Pred_Value)

classes <- levels(factor(lodo_preds$Truth))

kappa_by_class <- lapply(classes, function(cl) {
  
  truth_bin <- factor(ifelse(lodo_preds$Truth == cl, "Yes", "No"))
  pred_bin  <- factor(ifelse(lodo_preds$Pred == cl, "Yes", "No"))
  
  cm <- confusionMatrix(pred_bin, truth_bin)
  
  data.frame(
    Class = cl,
    Kappa = as.numeric(cm$overall["Kappa"])
  )
})

kappa_by_class <- do.call(rbind, kappa_by_class)
kappa_by_class

#### C: Plot Kappas for LODO analysis ----
ggplot(kappa_by_class, aes(x = Class, y = Kappa, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Kappa by Class (LODO + Percentile Normalization)") +
  theme_minimal() +
  guides(fill = "none")