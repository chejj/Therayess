64
### --- Step 3: Transpose so rows = samples, cols = features####################

# This is needed because RF expects rows to be samples, but SEs are not set up this way
taxa_df <- as.data.frame(t(studies_taxa))
path_abund_df <- as.data.frame(t(studies_path_abund))
path_cov_df <- as.data.frame(t(studies_path_cov))

# QUALITY CHECK: Check sample alignment
if (
  identical(rownames(taxa_df), rownames(studies_meta_df)) &&
  identical(rownames(path_abund_df), rownames(studies_meta_df)) &&
  identical(rownames(path_cov_df), rownames(studies_meta_df))
) {
  message("✅ Sample alignment confirmed between metadata and all feature tables.")
} else {
  stop("❌ Sample alignment mismatch between metadata and one or more feature tables.")
}

# Step 3.1: Remove duplicate subjects to reduce pseudoreplication ---------

# Keep the first sample encountered for each subject_id
keep_samples <- rownames(studies_meta_df)[!duplicated(studies_meta_df$subject_id)]

# Apply to metadata and all feature tables
studies_meta_df <- studies_meta_df[keep_samples, , drop = FALSE]
taxa_df         <- taxa_df[keep_samples, , drop = FALSE]
path_abund_df   <- path_abund_df[keep_samples, , drop = FALSE]
path_cov_df     <- path_cov_df[keep_samples, , drop = FALSE]

# QUALITY CHECK
if (any(duplicated(studies_meta_df$subject_id))) {
  stop("❌ Duplicate subject IDs still remain after filtering.")
} else {
  message("✅ Duplicate subject IDs removed.")
}
if (is.null(rownames(studies_meta_df))) {
  stop("❌ Rownames were lost after subject filtering.")
} else {
  message("✅ Rownames preserved after subject filtering.")
}

### --- Step 4: Prepare for modeling ###########################################
# Step 4.1: Disease class and study filtering ----------------------------------
# choose disease classes to keep for ML model so to not influence the split (now done in params)

# Create modeling columns to avoid removing original columns
studies_meta_df$RF_Class <- case_when(
  studies_meta_df$disease_class %in% params$class_keep ~ studies_meta_df$disease_class,
  TRUE ~ NA_character_
)

### This allows for class collapsing, see below:
# studies_meta_df$RF_Class <- dplyr::case_when(
#   studies_meta_df$disease_class %in% c("CRC", "CRC+", "CRC-H") ~ "CRC",
#   studies_meta_df$disease_class %in% c("PA", "PA+") ~ "PA",
#   studies_meta_df$disease_class == "HC" ~ "HC",
#   TRUE ~ NA_character_
# )

studies_meta_df$keep_study <- !studies_meta_df$study_name %in% params$study_drop

# Get samples to keep
keep_samples <- rownames(studies_meta_df)[
  !is.na(studies_meta_df$RF_Class) &
    studies_meta_df$keep_study
]

# Step 4.2: Apply across all metadata and feature tables -----------------------
studies_meta_df <- studies_meta_df[keep_samples, , drop = FALSE]
taxa_df         <- taxa_df[keep_samples, , drop = FALSE]
path_abund_df   <- path_abund_df[keep_samples, , drop = FALSE]
path_cov_df     <- path_cov_df[keep_samples, , drop = FALSE]

# Make sure column names are safe for downstream processes
make_safe_names <- function(df, prefix) {
  colnames(df) <- paste0(prefix, make.names(colnames(df), unique = TRUE))
  return(df)
}

taxa_df       <- make_safe_names(taxa_df, "taxa__")
path_abund_df <- make_safe_names(path_abund_df, "pabund__")
path_cov_df   <- make_safe_names(path_cov_df, "pcov__")


### --- Step 5: 70/30 study-aware train/test split to avoid leakage#############
# Split studies, not samples
# we want to make sure that all classes chosen are represented in the test and train data
# find another way to evenly split with good representation
# createDataPartition() failed all attempts when using disease_class

### OPTION: do this manually
set.seed(123) #allows random operations to be reproducible

# Table for samples by study, ordered in descending order
study_sizes <- studies_meta_df %>%
  filter(keep_study & !is.na(RF_Class)) %>%
  group_by(study_name) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(desc(n))
# Shuffle order to prevent picking largest studies first for test set
study_sizes <- study_sizes[sample(nrow(study_sizes)), ]

total_n <- sum(study_sizes$n) # Total samples across studies
target_test_n <- params$target_test_pct * total_n # Target samples in TEST set

# Initialization
test_studies <- c()
running_total <- 0

# Loop in random order to add study to test set and stop once you hit target TEST
for (i in seq_len(nrow(study_sizes))) {
  test_studies <- c(test_studies, study_sizes$study_name[i])
  running_total <- running_total + study_sizes$n[i]
  if (running_total >= target_test_n) break
}
# Remaining samples go into training set
train_studies <- setdiff(study_sizes$study_name, test_studies)

# Get sample IDs belonging to those studies
train_samples <- rownames(studies_meta_df)[
  studies_meta_df$study_name %in% train_studies &
    !is.na(studies_meta_df$RF_Class) & # safety net that no NA class makes it through
    studies_meta_df$keep_study # safety net that no removed study makes it through
]

test_samples <- rownames(studies_meta_df)[
  studies_meta_df$study_name %in% test_studies &
    !is.na(studies_meta_df$RF_Class) &
    studies_meta_df$keep_study
]

# Metadata 
train_meta <- studies_meta_df[train_samples, , drop = FALSE]
test_meta  <- studies_meta_df[test_samples, , drop = FALSE]

# Taxa
train_taxa <- taxa_df[train_samples, , drop = FALSE]
test_taxa  <- taxa_df[test_samples, , drop = FALSE]

# Pathway abundance
train_path_abund <- path_abund_df[train_samples, , drop = FALSE]
test_path_abund  <- path_abund_df[test_samples, , drop = FALSE]

# Pathway coverage
train_path_cov <- path_cov_df[train_samples, , drop = FALSE]
test_path_cov  <- path_cov_df[test_samples, , drop = FALSE]

# QUALITY CHECK: Check sample alignment preservation
all_aligned <- all(
  identical(rownames(train_taxa), rownames(train_meta)),
  identical(rownames(test_taxa), rownames(test_meta)),
  identical(rownames(train_path_abund), rownames(train_meta)),
  identical(rownames(test_path_abund), rownames(test_meta)),
  identical(rownames(train_path_cov), rownames(train_meta)),
  identical(rownames(test_path_cov), rownames(test_meta))
)

if (all_aligned) {
  message("✅ Sample alignment preserved.")
} else {
  stop("❌ Sample alignment NOT preserved.")
}

# QUALITY CHECK: Train/Test Split
message("Kept classes: ", paste(unique(studies_meta_df$RF_Class), collapse = ", "))

train_prop <- nrow(train_meta) / nrow(studies_meta_df)
test_prop  <- nrow(test_meta)  / nrow(studies_meta_df)
train_classes <- table(train_meta$RF_Class)
test_classes  <- table(test_meta$RF_Class)

# Check 1: No study leakage
no_leakage <- length(intersect(
  unique(train_meta$study_name),
  unique(test_meta$study_name)
)) == 0
if (no_leakage) {
  message("✅ No study leakage between train and test sets.")
} else {
  stop("❌ Study leakage detected between train and test sets.")
}

# Check 2: Reasonable split (allow tolerance, e.g. 60–80%)
split_ok <- train_prop > 0.6 & train_prop < 0.8
if (split_ok) {
  message(sprintf("✅ Split looks reasonable: Train = %.2f | Test = %.2f", train_prop, test_prop))
  message(sprintf("Train samples: %d | Test samples: %d", nrow(train_meta), nrow(test_meta)))
} else {
  message(sprintf("⚠️ Split is skewed: Train = %.2f | Test = %.2f", train_prop, test_prop))
}

# Check 3: All classes present in test set
all_classes_present <- all(unique(train_meta$RF_Class) %in% unique(test_meta$RF_Class))
if (all_classes_present) {
  message("✅ All classes in training set are present in test set.")
} else {
  message("⚠️ Some classes are missing from the test set.")
}
train_dist <- prop.table(table(train_meta$RF_Class))
test_dist  <- prop.table(table(test_meta$RF_Class))
dist_df <- rbind(
  data.frame(Class = names(train_dist), Proportion = as.numeric(train_dist), Set = "Train"),
  data.frame(Class = names(test_dist),  Proportion = as.numeric(test_dist),  Set = "Test")
)
ggplot(dist_df, aes(x = Class, y = Proportion, fill = Set)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(title = "Class Distribution: Train vs Test")

# Step 5.1: Feature Transformations --------------------------------------------
# A: METADATA ------------------------------------------------------------------
# 1. Start from train/test metadata
train_meta_x <- train_meta
test_meta_x  <- test_meta

# 2. Drop unwanted metadata columns (from params)
train_meta_x <- train_meta_x[, !names(train_meta_x) %in% params$meta_drop, drop = FALSE]
test_meta_x  <- test_meta_x[, !names(test_meta_x) %in% params$meta_drop, drop = FALSE]

# 3. Apply metadata completeness filter (TRAIN ONLY)
# Keep columns where at least X% of values are non-missing in TRAIN
meta_completeness <- params$meta_completeness

meta_keep <- colMeans(!is.na(train_meta_x)) >= meta_completeness

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

# B: RELATIVE ABUNDANCE --------------------------------------------------------
## Keep biologically present features and statistically informative features, based only on training data.
# 1. Start from train/test relative abundance
train_taxa_X <- train_taxa
test_taxa_X  <- test_taxa

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
prev_threshold <- params$prev_threshold
feature_prevalence <- colMeans(train_taxa_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= prev_threshold])

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

# 4. Normalize + log transform (Percentile might be better)
train_rel <- sweep(train_taxa_filt, 2, colSums(train_taxa_filt), "/")
test_rel  <- sweep(test_taxa_filt, 2, colSums(test_taxa_filt), "/")

train_taxa_filt <- log10(train_rel + 1e-6)
test_taxa_filt  <- log10(test_rel + 1e-6)

message("Relative Abundance Features Kept: ", ncol(train_taxa_filt))

# C: PATHWAY ABUNDANCE ---------------------------------------------------------
# 1. Start from train/test pathway abundance
train_pathab_X <- train_path_abund
test_pathab_X  <- test_path_abund

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
prev_threshold <- params$prev_threshold
feature_prevalence <- colMeans(train_pathab_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= prev_threshold])

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

# 4. Normalize + log transform (Percentile might be better)
train_rel <- sweep(train_pathab_filt, 2, colSums(train_pathab_filt), "/")
test_rel  <- sweep(test_pathab_filt, 2, colSums(test_pathab_filt), "/")

train_pathab_filt <- log10(train_rel + 1e-6)
test_pathab_filt  <- log10(test_rel + 1e-6)

message("Pathway Abundance Features Kept: ", ncol(train_pathab_filt))

# D: PATHWAY COVERAGE ~~~~~~~~~~~~~~~~~~~~~
# 1. Start from train/test pathway coverage
train_pathcov_X <- train_path_cov
test_pathcov_X  <- test_path_cov

# 2. Prevalence filtering
# Keep features present (>0) in at least X% of samples
prev_threshold <- params$prev_threshold
feature_prevalence <- colMeans(train_pathcov_X > 0) #compute prevalence of each column in TRAIN

# keep features that surpass the threshold set previously
keep_prev_features <- names(feature_prevalence[feature_prevalence >= prev_threshold])

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

message("Pathway Coverage Features Kept: ", ncol(train_pathcov_filt))


# Step 5.2: Fit Random Forest --------------------------------------------------

# 1. Model Input Preparation
train_X <- cbind(train_meta_x, train_taxa_filt, train_pathab_filt, train_pathcov_filt)
test_X  <- cbind(test_meta_x,  test_taxa_filt,  test_pathab_filt,  test_pathcov_filt)

train_df <- cbind(train_X, RF_Class = train_meta$RF_Class)
test_df  <- cbind(test_X,  RF_Class = test_meta$RF_Class)

# Ensure outcome is factor
train_df$RF_Class <- factor(train_df$RF_Class)
test_df$RF_Class  <- factor(test_df$RF_Class)

# Check and prep for class imbalance
class_counts <- table(train_df$RF_Class)
class_weights <- sum(class_counts) / (length(class_counts) * class_counts)

# Number of predictor features (exclude outcome column)
p <- ncol(train_df) - 1

# Set mtry to X% of total predictors, or 1, whichever is larger
mtry_val <- max(1, ceiling(params$mtry_fraction * p))
  
# 2. Model Generation
rf_fit <- ranger(
  RF_Class ~ .,              # outcome modeled by all predictors
  data = train_df,           # training data only
  num.trees = params$num_trees,           # start with 500 trees, approach 10000
  class.weights = class_weights,
  mtry = mtry_val,           # % of features considered at each split
  min.node.size = params$min_node_size,         # minimum samples per terminal node
  probability = TRUE,        # needed for multiclass probabilities / ROC
  importance = "impurity",   # feature importance
  splitrule = "gini",        # closest standard classification impurity rule in ranger, entropy is not available in ranger
  num.threads = params$num_threads            # match HPC core request
)

# 3. Model Output
rf_fit

# Step 5.3: Check Random Forest final model on test set ------------------------
pred_probs <- predict(rf_fit, data = test_df)$predictions

pred_class <- colnames(pred_probs)[max.col(pred_probs)]

# 1. Confusion Matrix

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

# 2. Per-Class F1 Score
f1_scores <- cm$byClass[, "F1"]
f1_scores

f1_df <- data.frame(
  Class = names(f1_scores),
  F1 = f1_scores
)

ggplot(f1_df, aes(x = Class, y = F1)) +
  geom_col() +
  theme_bw() +
  labs(title = "Per-Class F1 Score", y = "F1 Score")


# 3. AUC-ROC Curve (Multiclass)
# a. Set up "One vs. Rest" due to inherent binary nature ---
roc_list <- list()       # for AUC
roc_df_list <- list()    # for plotting

for (class in colnames(pred_probs)) {
  
  binary_truth <- ifelse(test_df$RF_Class == class, 1, 0)
  
  # Skip if class not present
  if (sum(binary_truth) == 0) next
  
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

# b. Plot the ROC Curves ---
ggplot(roc_df_all, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # random baseline
  theme_bw() +
  labs(
    title = "ROC Curves (One-vs-Rest)",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) 

# c. Extract AUC Values ---
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

### Step 6: Leave-One-Study-Out as gold standard ###############################
#   For each study:
#   - hold that study out as the test set
#   - train on all remaining studies
#   - perform feature filtering using TRAIN ONLY
#   - fit random forest
#   - evaluate performance on the held-out study

# a: Setup
set.seed(123)  # ensures any random parts of RF are reproducible
LOSO_taxa <- taxa_df
LOSO_pathab <- path_abund_df
LOSO_pathcov <- path_cov_df

# Get study IDs to loop over
study_ids <- unique(studies_meta_df$study_name)

# Create storage objects
loso_accuracy_results <- list()
loso_f1_results <- list()
loso_auc_results <- list()
plots <- list()

# b: LOSO loop
for (test_study in study_ids) {
  
  message("Running LOSO for held-out study: ", test_study)
  
  # -------------------------
  # A. Define train/test samples by study
  # -------------------------
  # Train = all samples NOT in the held-out study
  # Test  = all samples IN the held-out study
  train_samples <- rownames(studies_meta_df)[
    studies_meta_df$study_name != test_study
  ]
  
  test_samples <- rownames(studies_meta_df)[
    studies_meta_df$study_name == test_study
  ]
  
  # -------------------------
  # B. Subset metadata and predictors
  # -------------------------
  train_meta <- studies_meta_df[train_samples, , drop = FALSE]
  test_meta  <- studies_meta_df[test_samples, , drop = FALSE]
  
  train_X <- X_df[train_samples, , drop = FALSE]
  test_X  <- X_df[test_samples, , drop = FALSE]
  
  # -------------------------
  # C. QC: confirm row alignment is preserved
  # -------------------------
  if (
    identical(rownames(train_X), rownames(train_meta)) &&
    identical(rownames(test_X), rownames(test_meta))
  ) {
    message("  ✅ Train/test sample alignment preserved.")
  } else {
    stop("  ❌ Train/test sample alignment NOT preserved.")
  }
  
  # -------------------------
  # D. Skip studies with too few classes in test set
  # -------------------------
  # For multiclass evaluation, a held-out study with only one class is not very informative.
  # This avoids failures in confusion matrix / ROC calculations.
  if (length(unique(test_meta$RF_Class)) < 2) {
    message("  ⚠️ Skipping ", test_study, ": test set contains fewer than 2 classes.")
    next
  }
  
  # -------------------------
  # E. Feature filtering (TRAIN ONLY - avoid leakage.)
  # -------------------------
  # Prevalence filtering
  prev_threshold <- params$prev_threshold
  # Keep features present (> 0) in at least 10% of training samples.
  feature_prevalence <- colMeans(train_X > 0)
  
  keep_prev_features <- names(feature_prevalence[feature_prevalence >= prev_threshold])
  
  train_X_prev <- train_X[, keep_prev_features, drop = FALSE]
  test_X_prev  <- test_X[, keep_prev_features, drop = FALSE]
  
  
  # E2. Near-zero variance filtering
  # Remove features with almost no variation in the training set.
  nzv <- nearZeroVar(train_X_prev)
  
  if (length(nzv) > 0) {
    train_X_filt <- train_X_prev[, -nzv, drop = FALSE]
    test_X_filt  <- test_X_prev[, -nzv, drop = FALSE]
  } else {
    train_X_filt <- train_X_prev
    test_X_filt  <- test_X_prev
  }
  
  # -------------------------
  # F. Build modeling data frames
  # -------------------------
  # Combine filtered predictors with the outcome column.
  train_df <- cbind(train_X_filt, RF_Class = train_meta$RF_Class)
  test_df  <- cbind(test_X_filt,  RF_Class = test_meta$RF_Class)
  
  # Make sure outcome is a factor for classification.
  train_df$RF_Class <- factor(train_df$RF_Class)
  test_df$RF_Class  <- factor(test_df$RF_Class, levels = levels(train_df$RF_Class))
  
  # -------------------------
  # G. QC printout
  # -------------------------
  message("  Train samples: ", nrow(train_df))
  message("  Test samples: ", nrow(test_df))
  message("  Train classes: ", paste(unique(train_df$RF_Class), collapse = ", "))
  message("  Test classes: ", paste(unique(test_df$RF_Class), collapse = ", "))
  message("Original features: ", ncol(train_X))
  message("After prevalence filter: ", length(keep_prev_features))
  message("Final features after variance filter: ", ncol(train_X_filt))
  
  # -------------------------
  # H. Fit random forest
  # -------------------------
  # Number of predictor features in this LOSO training fold
  p <- ncol(train_df) - 1
  
  # Set mtry to the chosen fraction of predictors, minimum of 1
  mtry_val <- max(1, ceiling(params$mtry_fraction * p))
  
  # Compute class weights for THIS LOSO training set
  class_counts <- table(train_df$RF_Class)
  class_weights <- sum(class_counts) / (length(class_counts) * class_counts)
  
  rf_fit <- ranger(
    RF_Class ~ .,
    data = train_df,
    num.trees = params$num_trees,          # start here; increase later if desired
    class.weights = class_weights,
    mtry = mtry_val,          # % of features per split
    min.node.size = params$min_node_size,        # minimum samples per terminal node
    probability = TRUE,       # needed for multiclass probabilities / ROC
    importance = "impurity",  # feature importance measure
    splitrule = "gini",       # standard classification split rule in ranger
    num.threads = params$num_threads           # set equal to requested HPC cores
  )
  
  # -------------------------
  # I. Predict on held-out study
  # -------------------------
  pred_probs <- predict(rf_fit, data = test_df)$predictions
  
  # Convert probability matrix into predicted class labels
  pred_class <- colnames(pred_probs)[max.col(pred_probs)]
  pred_class <- factor(pred_class, levels = levels(train_df$RF_Class))
  
  # -------------------------
  # J. Confusion matrix
  # -------------------------
  cm <- confusionMatrix(
    data = pred_class,
    reference = test_df$RF_Class
  )
  
  cm_df <- as.data.frame(cm$table)
  plots[[test_study]] <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_bw() +
      labs(title = paste("Confusion Matrix (LOSO) —", test_study)) 
  
  # -------------------------
  # K. Store overall accuracy
  # -------------------------
  loso_accuracy_results[[test_study]] <- data.frame(
    Study = test_study,
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Kappa = as.numeric(cm$overall["Kappa"]),
    Train_Samples = nrow(train_df),
    Test_Samples = nrow(test_df),
    Features = ncol(train_df) - 1
  )
  
  
  # -------------------------
  # L. Store per-class F1 scores
  # -------------------------
  # In multiclass settings, cm$byClass is usually a matrix with one row per class.
  if (is.matrix(cm$byClass)) {
    f1_df <- data.frame(
      Study = test_study,
      Class = rownames(cm$byClass),
      F1 = cm$byClass[, "F1"],
      row.names = NULL
    )
    
    loso_f1_results[[test_study]] <- f1_df
  }
  
  # -------------------------
  # M. Compute one-vs-rest AUC for each class
  # -------------------------
  roc_list <- list()
  auc_values <- c()
  
  for (class in colnames(pred_probs)) {
    
    # Create binary truth labels for this class vs all other classes
    binary_truth <- ifelse(test_df$RF_Class == class, 1, 0)
    
    # Skip if held-out test study has no positive examples for this class
    # or no negative examples for this class
    if (sum(binary_truth) == 0 || sum(binary_truth) == length(binary_truth)) {
      next
    }
    
    # Build ROC object
    roc_obj <- roc(binary_truth, pred_probs[, class], quiet = TRUE)
    roc_list[[class]] <- roc_obj
    
    # Store AUC value
    auc_values[class] <- as.numeric(auc(roc_obj))
  }
  
  # Save AUC results for this held-out study
  if (length(auc_values) > 0) {
    loso_auc_results[[test_study]] <- data.frame(
      Study = test_study,
      Class = names(auc_values),
      AUC = as.numeric(auc_values),
      row.names = NULL
    )
  }
  
  message("  ✅ Finished LOSO for ", test_study)
}

# --- Combine stored LOSO results into data frames ---

loso_accuracy_df <- bind_rows(loso_accuracy_results)
loso_f1_df <- bind_rows(loso_f1_results)
loso_auc_df <- bind_rows(loso_auc_results)

# --- QC: LOSO summary ---

message("✅ LOSO complete.")
message("Studies evaluated: ", nrow(loso_accuracy_df))
message("Mean LOSO accuracy: ", round(mean(loso_accuracy_df$Accuracy, na.rm = TRUE), 3))
message("Median LOSO accuracy: ", round(median(loso_accuracy_df$Accuracy, na.rm = TRUE), 3))

### PLOTTING

# Accuracy by Study
ggplot(loso_accuracy_df, aes(x = reorder(Study, Accuracy), y = Accuracy)) +
  geom_col() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  coord_flip() +
  labs(
    title = "LOSO Accuracy by Held-Out Study",
    x = "Held-Out Study",
    y = "Accuracy"
  )

# LOSO F1 by Class
ggplot(loso_f1_df, aes(x = Class, y = F1)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(
    title = "LOSO Per-Class F1 Distribution",
    x = "Class",
    y = "F1 Score"
  )

# LOSO AUC by Class
ggplot(loso_auc_df, aes(x = Class, y = AUC)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(
    title = "LOSO One-vs-Rest AUC Distribution",
    x = "Class",
    y = "AUC"
  )













# # Remove metadata columns you do not want as predictors
# rf_input <- test_df_small %>% 
#   filter(disease_class %in% c("HC", "PA", "CRC", "Other", "PA+", "CRC+")) %>% 
#   filter(country == "JPN") %>% 
#   select(-study_name, -study_condition, -non_westernized, -sequencing_platform, -PMID, -curator)
# 
# rf_input <- test_df_small %>%
#   filter(!is.na(disease_class)) %>% # Keep only rows with known outcome
#   mutate(disease_class = as.factor(disease_class)) %>% # Make outcome a factor
#   select(disease_class, starts_with("taxa__"), starts_with("pabund__"), starts_with("pcov__"))

