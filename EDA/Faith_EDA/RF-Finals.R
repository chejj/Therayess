###################################
# Try ing to complete model building
####################################

#---------------------------------------------
# Begin with creating a matrix for each study 
#---------------------------------------------
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

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
#----------------
rel_abundance <- sweep(filtered_matrix, 2, colSums(filtered_matrix), "/")
log_data <- log10(rel_abundance + 1e-6)

dim(log_data)

#---------------------
#Build metadata table
#---------------------
library(dplyr)

metadata_list <- lapply(CRC_progression_studies, function(study) {
  as.data.frame(colData(study$relative_abundance))
})

metadata <- bind_rows(metadata_list)

dim(metadata)
## Make sure it all worked. 
identical(colnames(combined_matrix), rownames(metadata))
##Safeguard option to confirm allignment 
metadata_filtered <- metadata_filtered[colnames(log_data), , drop = FALSE]
#--------------------------
# remove CRC-M and PA-M
#--------------------------
metadata_filtered <- metadata[!metadata$disease_class %in% c("CRC-M", "PA-M"), ]

table(metadata_filtered$disease_class)

#Filter log data to match the removed/filtered metadata_filtered
log_data_filtered <- log_data[, rownames(metadata_filtered), drop = FALSE]
dim(log_data_filtered)


#Drop HMP study due to HC only
metadata_filtered <- metadata_filtered[metadata_filtered$study_name != "HMP_2012", ]
log_data_filtered <- log_data_filtered[, rownames(metadata_filtered)]

#Create final label vector
labels <- droplevels(as.factor(metadata_filtered$disease_class))
table(labels)


#prep data for model bc random forest expects
#rows = samples and columns = features
X <- t(log_data_filtered)
y <- labels

dim(X)
length(y)

#---------------------
# run test train split
#----------------------
library(caret)

set.seed(42)

trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]

length(y_train)
length(y_test)



















