###################################
# Try ing to complete model building
####################################

#---------------------------------------------
# Begin with creating a matrix for each study 
#---------------------------------------------
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(randomForest)
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

length(y_train)
length(y_test)
table(y_train)
table(y_test)

#Percentile normalization RF model
PNrf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  importance = FALSE
)

predPN <- predict(PNrf_model, X_test)
confusionMatrix(predPN, y_test)

