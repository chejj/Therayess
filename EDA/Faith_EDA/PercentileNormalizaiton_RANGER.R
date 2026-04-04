###################################
# Try ing to complete model building
####################################

#---------------------------------------------
# Begin with creating a matrix for each study 
#---------------------------------------------
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(ranger)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pROC)


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

#------------------------
#Now apply to all studies 
#------------------------
combined_study_list <- lapply(CRC_progression_studies, build_combined_features)

#---------------------------------
# find all the shared features
#---------------------------------
common_features <- Reduce(intersect, lapply(combined_study_list, rownames))

length(common_features)


combined_study_list_intersection <- lapply(combined_study_list, function(mat) {
  mat[common_features, , drop = FALSE]
})

combined_matrix <- do.call(cbind, combined_study_list_intersection)
dim(combined_matrix)

#---------------------
# Build metadata table
#---------------------
metadata_list <- lapply(CRC_progression_studies, function(study) {
  as.data.frame(colData(study$relative_abundance))
})

metadata <- bind_rows(metadata_list)
dim(metadata)

# optional alignment check
identical(colnames(combined_matrix), rownames(metadata))

#--------------------------------------------
# REMOVE HMP-2012 and unwanted classes EARLY
#--------------------------------------------
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

#--------------------
# Prevalence filtering
#--------------------
prevalence_filter <- function(mat, threshold = 0.1) {
  prev <- rowSums(mat > 0) / ncol(mat)
  mat[prev >= threshold, , drop = FALSE]
}

filtered_matrix <- prevalence_filter(combined_matrix_filtered, 0.1)
dim(filtered_matrix)

#----------------
# Relative abundance
#----------------
rel_abundance <- sweep(filtered_matrix, 2, colSums(filtered_matrix), "/")
dim(rel_abundance)

# Add tiny jitter to reduce ties for percentile normalization
set.seed(42)
rel_abundance <- rel_abundance + matrix(
  runif(length(rel_abundance), min = 0, max = 1e-9),
  nrow = nrow(rel_abundance),
  ncol = ncol(rel_abundance),
  dimnames = dimnames(rel_abundance)
)

# make sure metadata matches rel_abundance
metadata_filtered <- metadata_filtered[colnames(rel_abundance), , drop = FALSE]

#-------------------------------
# Percentile normalization
#-------------------------------
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
  nrow = nrow(rel_abundance),
  ncol = ncol(rel_abundance),
  dimnames = dimnames(rel_abundance)
)

for (b in unique(batch)) {
  idx <- which(batch == b)
  
  sub_data <- rel_abundance[, idx, drop = FALSE]
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

#--------------------------
# Create final grouped labels
#--------------------------
metadata_filtered$label2 <- dplyr::case_when(
  metadata_filtered$disease_class == "HC" ~ "HC",
  metadata_filtered$disease_class %in% c("PA", "PA+") ~ "PA",
  metadata_filtered$disease_class %in% c("CRC", "CRC+") ~ "CRC",
  metadata_filtered$disease_class == "Other" ~ "Other"
)

y2 <- droplevels(as.factor(metadata_filtered$label2))
table(y2)

# prep data for model
X <- t(percentile_matrix)
y <- y2

dim(X)
length(y)

#---------------------
# run test train split
#----------------------
set.seed(42)

trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]

# Convert to data.frame and bind y (ranger prefers formula or full data)
train_df <- data.frame(y = y_train, X_train)
test_df  <- data.frame(X_test)

PNrf_model <- ranger(
  dependent.variable.name = "y", # outcome modeled by all predictors
  data = train_df,
  num.trees = 200,
  importance = "none",     
  classification = TRUE    # explicit for clarity
)

# Kaitlyn's Pipeline Model. Use this to further tune the above with more parameters?
# rf_fit <- ranger(
#   dependent.variable.name = "RF_Class", # outcome modeled by all predictors
#   data = train_df,                      # training data only
#   num.trees = params$num_trees,         # start with 500 trees, approach 10000
#   max.depth = params$max_depth,         # How many splits the tree goes to. Reduces specificity but also reduces overfitting
#   class.weights = class_weights,        # Formula: sum(class_counts) / (length(class_counts) * class_counts)
#   mtry = mtry_val,                      # % of features considered at each split. Uses formula: max(1, ceiling(params$mtry_fraction * p))
#   min.node.size = params$min_node_size, # minimum samples per terminal node. Currently set to 20 in Kaitlyn's pipeline
#   probability = TRUE,                   # needed for multiclass probabilities / ROC plotting
#   importance = "impurity",              # feature importance, so we can interpret for discussion
#   splitrule = "gini",                   # closest standard classification impurity rule in ranger, entropy is not available in ranger
#   num.threads = params$num_threads      # match HPC core request, set to 8 in Kaitlyn's pipeline
# )

predPN <- predict(PNrf_model, data = test_df)$predictions
confusionMatrix(predPN, y_test)


########################################
# LODO on percentile-normalized data
########################################
studies <- unique(metadata_filtered$study_name)

set.seed(42)
lodo_preds <- lapply(studies, function(test_study) {
  
  test_idx  <- which(metadata_filtered$study_name == test_study)
  train_idx <- which(metadata_filtered$study_name != test_study)
  
  X_train_lodo <- X[train_idx, , drop = FALSE]
  X_test_lodo  <- X[test_idx, , drop = FALSE]
  
  y_train_lodo <- y[train_idx]
  y_test_lodo  <- y[test_idx]
  
  train_df <- data.frame(Truth = y_train_lodo, X_train_lodo)
  test_df  <- data.frame(X_test_lodo)
  
  # Skip if fewer than 2 classes remain in training
  class_sizes_lodo <- table(train_df$Truth)
  if (length(class_sizes_lodo) < 2) {
    return(NULL)
  }
  
  min_class_lodo <- min(class_sizes_lodo)
  
  balanced_train_df <- train_df %>%
    dplyr::group_by(Truth) %>%
    dplyr::slice_sample(n = min_class_lodo) %>%
    dplyr::ungroup()
  
  balanced_train_df$Truth <- droplevels(as.factor(balanced_train_df$Truth))
  
  rg_lodo <- ranger(
    dependent.variable.name = "Truth",
    data = balanced_train_df,
    num.trees = 200,
    classification = TRUE,
    probability = FALSE,
    importance = "none"
  )
  
  pred_lodo <- predict(rg_lodo, data = test_df)$predictions
  
  data.frame(
    Study = test_study,
    Truth = y_test_lodo,
    Pred = pred_lodo
  )
})

lodo_preds <- do.call(rbind, lodo_preds)

#Create per study metrics 
lodo_results2 <- lodo_preds %>%
  group_by(Study) %>%
  do({
    cm <- confusionMatrix(.$Pred, .$Truth)
    
    data.frame(
      N_test = nrow(.),
      Accuracy = as.numeric(cm$overall["Accuracy"]),
      Kappa = as.numeric(cm$overall["Kappa"])
    )
  }) %>%
  ungroup()

lodo_results2
#get the LODO means 
mean(lodo_results2$Accuracy)
mean(lodo_results2$Kappa)

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

#Plot Kappas for LODO analysis 
ggplot(kappa_by_class, aes(x = Class, y = Kappa, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Kappa by Class (LODO + Percentile Normalization)") +
  theme_minimal() +
  guides(fill = "none")

