# Model 1: ComBat + Binary (ie Neg_polyp, Pos_Polyp)

## SETUP #######################################################################
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TreeSummarizedExperiment)
  library(rsample) #study-aware train/test split
  library(ranger)
  library(caret)
  library(tidyverse)
  library(pROC)
  library(sva) #ComBat
})
###### PARAMETERS ##############################################################
params <- list(
  num_trees = 500,
  max_depth = 10,
  mtry_fraction = 0.03, # % of features considered at each split
  min_node_size = 15,   # minimum samples per terminal node
  num_threads = 8,      # Set to match HPC core
  splitrule = "gini",   # options: "gini", "extratrees" or "hellinger"
  prev_filter = 0.01,    # proportion filtered out in prevalence filtering step (i.e. 0.1 is 10%)
  train_split = 0.8,     # proportion of the split to delegate for training (i.e. 0.8 is 80/20 split)
  class_keep = c(
    "HC",
    "PA",
    "CRC",
    "Other",
    "PA+",
    "CRC+"
#    "CRC-H", 
#    "CRC-M", 
#    "PA-M"
  ), #options: "HC", "PA", "CRC", "Other", "PA+", "CRC+", "CRC-H", "CRC-M", "PA-M" 
  study_drop = c("HMP_2012"), 
  meta_drop = c(
    "RF_Class",
#    "country",
    "age",
    "disease_class",
    "keep_study",
    "study_name",
    "subject_id",
    "study_condition",
    "disease",
    "sequencing_platform",
    "DNA_extraction_kit",
    "PMID",
    "curator",
    "disease_stage",
    "disease_location"
  ),
  meta_completeness = 0.90
)

###--- OUTLINE / PLAN: ---######################################################
# Step 0: QC Check: Is data loaded and are sample IDs/order are identical across? 
# Step 1: Final Metadata Clean -> Final Class filters and grouping
# Step 2: Prep for Modeling
  # Step 2.1: Filtering Studies and feature tables
  # Step 2.2: Study-aware train/test split to avoid leakage
  # Step 2.3: Feature Transformations (ComBat)
    # NOTE: Transformations can only be applied to the training set, not the testing set, as ComBat does not have an easy "apply to new samples" method
    # Prevalence & Variance filtering
# Step 3: Fit Random Forest 
  # num.threads to same as HPC core
  # Class imbalance handling (class weights in ranger)
  # Cross validation for hyperparameter tuning
  # OOB Error for internal performance estimate
# Step 4: Check Random Forest final model on test set
  # Step 4.1: Confusion Matrix
  # Step 4.2: F1-Score
  # Step 4.3: AUC-ROC Curve
# Step 5: Inspect Feature Importance

## Step 0: QC Checks ----
if (!exists("merged_studies")) {
  stop("Error: `merged_studies` not loaded. Be sure to run 'Study-Merge+Clean.R' prior to running any of the ML scripts.")
}

if (identical(
  colnames(merged_studies[["relative_abundance"]]),  colnames(merged_studies[["pathway_abundance"]])
) && identical(
  colnames(merged_studies[["relative_abundance"]]),
  colnames(merged_studies[["pathway_coverage"]])
)) {
  message("✅ Sample IDs/order are identical across all three merged objects.")
} else {
  message("❌ Sample IDs/order are NOT identical across all three merged objects.")
}

## Step 1: Final Metadata Cleaning ----

# Create modeling columns to avoid removing original columns
studies_meta_df$RF_Class <- case_when(
  studies_meta_df$disease_class %in% params$class_keep ~ studies_meta_df$disease_class,
  TRUE ~ NA_character_
)

### This allows for class collapsing:
studies_meta_df$RF_Class <- dplyr::case_when(
  studies_meta_df$disease_class %in% c("HC", "Other") ~ "NEGATIVE_POLYP",
  studies_meta_df$disease_class %in% c("CRC", "CRC+", "CRC-M", "PA", "PA+", "PA-M") ~ "POSITIVE_POLYP",
  TRUE ~ NA_character_
) %>% as.factor() %>% factor(levels = c("POSITIVE_POLYP", "NEGATIVE_POLYP")  # 1st = positive class, 2nd = reference
)

studies_meta_df$keep_study <- !studies_meta_df$study_name %in% params$study_drop

## Step 2: Prep For Modeling ----
### 2.1: Filtering Datasets ----
# Get samples to keep
keep_samples <- rownames(studies_meta_df)[
  !is.na(studies_meta_df$RF_Class) &
    studies_meta_df$keep_study
]

#Filter all feature tables
studies_meta_df      <- studies_meta_df[keep_samples, , drop = FALSE]
studies_taxa         <- studies_taxa[keep_samples, , drop = FALSE]
studies_path_abund   <- studies_path_abund[keep_samples, , drop = FALSE]
studies_path_cov     <- studies_path_cov[keep_samples, , drop = FALSE]

# Make sure column names are safe for downstream processes
make_safe_names <- function(df, prefix) {
  colnames(df) <- paste0(prefix, make.names(colnames(df), unique = TRUE))
  return(df)
}

studies_taxa       <- make_safe_names(studies_taxa, "taxa__")
studies_path_abund <- make_safe_names(studies_path_abund, "pabund__")
studies_path_cov   <- make_safe_names(studies_path_cov, "pcov__")

### Step 2.2: Study-aware train/test split to avoid leakage ----
set.seed(42)

split_index  <- group_initial_split(
  studies_meta_df,
  group = study_name,  # one identifier per source study
  prop = params$train_split
)

train_meta <- training(split_index)
test_meta  <- testing(split_index)

# Sample IDs are the metadata row names
train_samples <- rownames(train_meta)
test_samples  <- rownames(test_meta)

# Optional: inspect which studies landed in each split
train_studies <- unique(train_meta$study_name)
test_studies  <- unique(test_meta$study_name)

## QC Check:
# Confirm no study is in both sets
stopifnot(length(intersect(train_studies, test_studies)) == 0)
# Confirm feature table frow names are sample IDs matching metadata row names
stopifnot(
  all(train_samples %in% rownames(studies_taxa)),
  all(test_samples  %in% rownames(studies_taxa)),
  all(train_samples %in% rownames(studies_path_abund)),
  all(test_samples  %in% rownames(studies_path_abund)),
  all(train_samples %in% rownames(studies_path_cov)),
  all(test_samples  %in% rownames(studies_path_cov))
)

# Check class coverage
prop.table(table(droplevels(train_meta$disease_class)))
prop.table(table(droplevels(test_meta$disease_class)))

# Apply the same sample split to each feature table
train_taxa <- studies_taxa[train_samples, , drop = FALSE]
test_taxa  <- studies_taxa[test_samples, , drop = FALSE]

train_path_abund <- studies_path_abund[train_samples, , drop = FALSE]
test_path_abund  <- studies_path_abund[test_samples, , drop = FALSE]

train_path_cov <- studies_path_cov[train_samples, , drop = FALSE]
test_path_cov  <- studies_path_cov[test_samples, , drop = FALSE]

### 2.3: Feature Transformations ----
#### Setup ----
batch <- studies_meta_df$study_name # expected batch effect
group <- studies_meta_df$RF_Class # groups we want to be able to detect
control_label <- "HC" # baseline group or control
run_combat <- function(control_mat,
                       metadata,
                       batch,
                       method_name = NULL) {
  
  # ComBat expects features x samples, so transpose in and out (currently in samples x features)
  combat_input <- t(control_mat)
  
  # preserve RF_Class so batch correction does not remove biology
  mod <- model.matrix(
    ~ RF_Class,
    data = metadata
  )
  
  combat_corrected <- ComBat(
    dat = combat_input,
    batch = factor(batch),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  message("[", method_name, "] COMBAT correction complete.")
  
  t(combat_corrected)
}
#### A: METADATA ----
# 1. Start from train/test metadata
train_meta_x <- train_meta
test_meta_x  <- test_meta

# 2. Drop unwanted metadata columns (from params)
train_meta_x <- train_meta_x[, !names(train_meta_x) %in% params$meta_drop, drop = FALSE]
test_meta_x  <- test_meta_x[, !names(test_meta_x) %in% params$meta_drop, drop = FALSE]

# 3. Apply metadata completeness filter (TRAIN ONLY)
# Keep columns where at least X% of values are non-missing in TRAIN
meta_keep <- colMeans(!is.na(train_meta_x)) >= params$meta_completeness

train_meta_x <- train_meta_x[, meta_keep, drop = FALSE]
test_meta_x  <- test_meta_x[, meta_keep, drop = FALSE]

# 4. Make metadata column names safe for modeling
train_meta_x <- make_safe_names(train_meta_x, "meta__")
test_meta_x  <- make_safe_names(test_meta_x,  "meta__")

# Ensure alignment preserved
if (
  identical(rownames(train_meta_x), rownames(train_meta)) &&
  identical(rownames(test_meta_x), rownames(test_meta))
) {
  message("✅ Metadata predictors aligned with train/test samples.")
} else {
  stop("❌ Metadata alignment issue after filtering.")
}

# Print summary
message("Metadata features kept: ", ncol(train_meta_x))

#### B: RELATIVE ABUNDANCE ----
## Keep biologically present features and statistically informative features, based only on training data.
# 1. Start from train/test relative abundance
train_taxa_X <- train_taxa
test_taxa_X  <- test_taxa

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
feature_prevalence <- colMeans(train_taxa_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= params$prev_filter])

# subset test and train using features selected from TRAIN to avoid leakage
train_X_prev <- train_taxa_X[, keep_prev_features, drop = FALSE]
test_X_prev  <- test_taxa_X[, keep_prev_features, drop = FALSE]

# 3. Variance filtering
# Remove near-zero variance features in TRAIN
nzv <- nearZeroVar(train_X_prev) # can tune with freqCut and uniqueCut

if (length(nzv) > 0) { # any near zero features found in TRAIN are removed from both TRAIN and TEST
  train_taxa_filt <- train_X_prev[, -nzv, drop = FALSE]
  test_taxa_filt  <- test_X_prev[, -nzv, drop = FALSE]
} else { # if nothing flagged, then data is unchanged
  train_taxa_filt <- train_X_prev
  test_taxa_filt  <- test_X_prev
}

# 4. Normalize with ComBat (Train Only)
train_taxa_filt <- run_combat(
  control_mat = train_taxa_filt,
  metadata = train_meta,
  batch = train_meta$study_name,
  method_name = "Taxa"
)

message("Relative Abundance Features Kept: ", ncol(train_taxa_filt))

#### C: PATHWAY ABUNDANCE ----
# 1. Start from train/test pathway abundance
train_pathab_X <- train_path_abund
test_pathab_X  <- test_path_abund

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
feature_prevalence <- colMeans(train_pathab_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= params$prev_filter])

# subset test and train using features selected from TRAIN to avoid leakage
train_X_prev <- train_pathab_X[, keep_prev_features, drop = FALSE]
test_X_prev  <- test_pathab_X[, keep_prev_features, drop = FALSE]

# 3. Variance filtering
# Remove near-zero variance features in TRAIN
nzv <- nearZeroVar(train_X_prev) # can tune with freqCut and uniqueCut

if (length(nzv) > 0) { # any near zero features found in TRAIN are removed from both TRAIN and TEST
  train_pathab_filt <- train_X_prev[, -nzv, drop = FALSE]
  test_pathab_filt  <- test_X_prev[, -nzv, drop = FALSE]
} else { # if nothing flagged, then data is unchanged
  train_pathab_filt <- train_X_prev
  test_pathab_filt  <- test_X_prev
}

# 4. Normalize with ComBat (Train Only)
train_pathab_filt <- run_combat(
  control_mat = train_pathab_filt,
  metadata = train_meta,
  batch = train_meta$study_name,
  method_name = "Pathway Abundance"
)

message("Pathway Abundance Features Kept: ", ncol(train_pathab_filt))

#### D: PATHWAY COVERAGE ----
# 1. Start from train/test pathway coverage
train_pathcov_X <- train_path_cov
test_pathcov_X  <- test_path_cov

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
feature_prevalence <- colMeans(train_pathcov_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= params$prev_filter])

# subset test and train using features selected from TRAIN to avoid leakage
train_X_prev <- train_pathcov_X[, keep_prev_features, drop = FALSE]
test_X_prev  <- test_pathcov_X[, keep_prev_features, drop = FALSE]

# 3. Variance filtering
# Remove near-zero variance features in TRAIN
nzv <- nearZeroVar(train_X_prev) # can tune with freqCut and uniqueCut

if (length(nzv) > 0) { # any near zero features found in TRAIN are removed from both TRAIN and TEST
  train_pathcov_filt <- train_X_prev[, -nzv, drop = FALSE]
  test_pathcov_filt  <- test_X_prev[, -nzv, drop = FALSE]
} else { # if nothing flagged, then data is unchanged
  train_pathcov_filt <- train_X_prev
  test_pathcov_filt  <- test_X_prev
}

# 4. Normalize with ComBat (Train Only)
train_pathcov_filt <- run_combat(
  control_mat = train_pathcov_filt,
  metadata = train_meta,
  batch = train_meta$study_name,
  method_name = "Pathway Coverage"
)

message("Pathway Coverage Features Kept: ", ncol(train_pathcov_filt))



# Step 3: Fit Random Forest --------------------------------------------------
## A. Model Input Preparation ----
# Potentially only keep those that have all three (FAITH CODE HAS THIS)

train_X <- cbind(
  train_meta_x, 
  train_taxa_filt, 
  train_pathab_filt, 
  train_pathcov_filt)
test_X  <- cbind(
  test_meta_x,  
  test_taxa_filt,  
  test_pathab_filt,  
  test_pathcov_filt)

train_df <- cbind(train_X, RF_Class = train_meta$RF_Class)
test_df  <- cbind(test_X,  RF_Class = test_meta$RF_Class)

# Ensure outcome is factor
train_df$RF_Class <- factor(train_df$RF_Class)
test_df$RF_Class  <- factor(test_df$RF_Class)

# Check and prep for class imbalance
class_counts <- table(train_df$RF_Class)
class_weights <- sum(class_counts) / (length(class_counts) * class_counts) # Keep as table to keep them "named", ranger prefers this way
# Do not convert class_weights to numeric using "as.numeric", table is already seen as a numeric in the environment panel
names(class_weights) <- names(class_counts)

print(class_counts)
print(class_weights)

# Number of predictor features (exclude outcome column)
p <- ncol(train_df) - 1

# Set mtry to X% of total predictors, or 1, whichever is larger
mtry_val <- max(1, ceiling(params$mtry_fraction * p))

## B. Model Generation ----
rf_fit <- ranger(
  dependent.variable.name = "RF_Class", # outcome modeled by all predictors
  data = train_df,                       # training data only
  num.trees = params$num_trees,         # start with 500 trees, approach 10000
  max.depth = params$max_depth,
  class.weights = class_weights,
  mtry = mtry_val,                      # % of features considered at each split
  min.node.size = params$min_node_size, # minimum samples per terminal node
  probability = TRUE,        # needed for multiclass probabilities / ROC
  importance = "impurity",   # feature importance
  splitrule = "gini"        # closest standard classification impurity rule in ranger, entropy is not available in ranger
#  num.threads = params$num_threads      # match HPC core request
)

## C. Model Output ----
rf_fit

# Step 4: Check Random Forest Model on Test Split ----
pred_probs <- predict(rf_fit, data = test_df[, colnames(test_df) != "RF_Class"])$predictions 
pred_class <- colnames(pred_probs)[max.col(pred_probs)]
pred_class <- factor(pred_class, levels = levels(train_df$RF_Class))

## Step 4.1: Confusion Matrix ----
cm <- confusionMatrix(
  factor(pred_class, levels = levels(train_df$RF_Class)),
  factor(test_df$RF_Class, levels = levels(train_df$RF_Class))
)
cm_df <- as.data.frame(cm$table)

# Add percentages (row-wise = % of actual class)
cm_df <- cm_df %>%
  group_by(Reference) %>%
  mutate(Percent = Freq / sum(Freq) * 100)

### A: Plots ----
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

#### B: Summary ----
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
  y_test = test_df$RF_Class,
  class_weights = class_weights,
  summary = results_summary
)

print(cm)
print(results_summary)

## Step 4.2: F1-Score ----
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

## Step 4.3: AUC-ROC Curve ----

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

## 4.4: Overfitting Check ----
#compare train vs test accuracy
# Get predicted probabilities  and convert to predicted class
train_pred_probs <- predict(rf_fit, data = train_df)$predictions
train_pred_class <- colnames(train_pred_probs)[max.col(train_pred_probs)]

test_pred_probs <- predict(rf_fit, data = test_df)$predictions
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

# Step 5: Inspect Feature Importance ----
# Built In
sort(importance(rf_fit), decreasing = TRUE)[1:15]

# Count Splits
split_counts <- unlist(
   lapply(1:rf_fit$num.trees, function(t) {
     treeInfo(rf_fit, tree = t)$splitvarName
   })
 )
head(sort(table(split_counts, useNA = "no"), decreasing = TRUE), 15)
