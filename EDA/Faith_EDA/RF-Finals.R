###################################
# Try ing to complete model building
####################################

#---------------------------------------------
# Begin with creating a matrix for each study 
#---------------------------------------------
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(caret)

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

#-------------------------
# Check if it worked
#-------------------------
length(combined_study_list)
dim(combined_study_list[[1]])
head(rownames(combined_study_list[[1]]))
head(colnames(combined_study_list[[1]]))

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

#Check unique sample IDs
any(duplicated(colnames(combined_matrix)))

#--------------------
#Prevalence filtering
#--------------------
#what percent of patients have this feature present
prevalence_filter <- function(mat, threshold = 0.1) {
  prev <- rowSums(mat > 0) / ncol(mat)
  #keep only features that appear in at least 10% of patients
  mat[prev >= threshold, , drop = FALSE]
}


#new dimensions after filtering 
filtered_matrix <- prevalence_filter(combined_matrix, 0.1)
dim(filtered_matrix)

#----------------
# Normalize data
#rel_abundance <- sweep(filtered_matrix, 2, colSums(filtered_matrix), "/")
#log_data <- log10(rel_abundance + 1e-6)
#dim(log_data)
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


#---------------------
#Build metadata table
#---------------------
metadata_list <- lapply(CRC_progression_studies, function(study) {
  as.data.frame(colData(study$relative_abundance))
})

metadata <- bind_rows(metadata_list)
dim(metadata)

## Make sure sample order matches combined matrix. 
identical(colnames(combined_matrix), rownames(metadata))

##Safeguard option to confirm alignment 
metadata_filtered <- metadata_filtered[colnames(log_data), , drop = FALSE]
#--------------------------
# remove CRC-M and PA-M
#--------------------------
metadata_filtered <- metadata[!metadata$disease_class %in% c("CRC-M", "PA-M"), ]

table(metadata_filtered$disease_class)

# Filter rel_abundance to match metadata_filtered
rel_abundance_filtered <- rel_abundance[, rownames(metadata_filtered), drop = FALSE]

# Drop HMP study due to HC only
metadata_filtered <- metadata_filtered[metadata_filtered$study_name != "HMP_2012", , drop = FALSE]
rel_abundance_filtered <- rel_abundance_filtered[, rownames(metadata_filtered), drop = FALSE]

# Final alignment check
identical(colnames(rel_abundance_filtered), rownames(metadata_filtered))
dim(rel_abundance_filtered)
dim(metadata_filtered)

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
  nrow = nrow(rel_abundance_filtered),
  ncol = ncol(rel_abundance_filtered),
  dimnames = dimnames(rel_abundance_filtered)
)

for (b in unique(batch)) {
  idx <- which(batch == b)
  
  sub_data <- rel_abundance_filtered[, idx, drop = FALSE]
  sub_group <- group[idx]
  
  if (!any(sub_group == control_label)) {
    next
  }
  
  sub_percentile <- apply(sub_data, 1, function(feature) {
    percentile_norm(feature, sub_group, control_label)
  })
  
  percentile_matrix[, idx] <- t(sub_percentile)
}

# Check result
dim(percentile_matrix)
summary(as.vector(percentile_matrix))
sum(is.na(percentile_matrix))

#--------------------------
# Create final label vector
#--------------------------
labels <- droplevels(as.factor(metadata_filtered$disease_class))
table(labels)

#prep data for model bc random forest expects
#rows = samples and columns = features
X <- t(percentile_matrix)
#Combine into only 4  groups, CRC, HC, PA, Other
y <- ifelse(labels %in% c("CRC", "CRC+", "CRC-H"), "CRC",
            ifelse(labels %in% c("PA", "PA+"), "PA",
                   ifelse(labels == "HC", "HC", "Other")))

y <- as.factor(y)
table(y)

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



















