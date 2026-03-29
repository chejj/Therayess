#This is my fave model all chopped up and what not but if we run
# this whole thing it works lol 

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
#Check unique sample IDs
any(duplicated(colnames(combined_matrix)))
#Check unique subject IDs
any(duplicated(metadata_filtered$subject_id))
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
#Prevalence filtering
#--------------------
#what percent of patients have this feature present
prevalence_filter <- function(mat, threshold = 0.1) {
  prev <- rowSums(mat > 0) / ncol(mat)
  #keep only features that appear in at least 10% of patients
  mat[prev >= threshold, , drop = FALSE]
}

#new dimensions after filtering 
filtered_matrix <- prevalence_filter(combined_matrix_filtered, 0.1)
dim(filtered_matrix)

#----------------
# Normalize data
#----------------
rel_abundance <- sweep(filtered_matrix, 2, colSums(filtered_matrix), "/")
log_data <- log10(rel_abundance + 1e-6)

dim(log_data)
# make sure metadata matches transformed data
metadata_filtered <- metadata_filtered[colnames(log_data), , drop = FALSE]
#---------------------
#Create labels
#---------------------
labels <- droplevels(as.factor(metadata_filtered$disease_class))
table(labels)

#prep data for model bc random forest expects
#rows = samples and columns = features
X <- t(log_data)
y <- labels

dim(X)
length(y)

##########################
metadata_filtered$label2 <- dplyr::case_when(
  metadata_filtered$disease_class == "HC" ~ "HC",
  metadata_filtered$disease_class %in% c("PA", "PA+") ~ "PA",
  metadata_filtered$disease_class %in% c("CRC", "CRC+") ~ "CRC",
  metadata_filtered$disease_class == "Other" ~ "Other"
)

table(metadata_filtered$label2)

#rebuild Y
y2 <- droplevels(as.factor(metadata_filtered$label2))

#---------------------
# run test train split
#----------------------
set.seed(42)

trainIndex <- createDataPartition(y2, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y2[trainIndex]
y_test  <- y2[-trainIndex]

length(y_train)
length(y_test)

#change focus of model to treat all class sizes equally instead of focusing on HC and CRC
class_sizes <- table(y_train)
class_sizes

min_class <- min(class_sizes[class_sizes > 0])


rf_model4 <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  sampsize = rep(min_class, length(class_sizes))
)
pred_final <- predict(rf_model4, X_test)
confusionMatrix(pred_final, y_test)
#--------------------------------------------
#create a clean table for feature importance 
#--------------------------------------------
importance_mat <- importance(rf_model4)
feature_importance <- data.frame(
  Feature = rownames(importance_mat),
  MeanDecreaseGini = importance_mat[, "MeanDecreaseGini"]
)

# sort by importance
feature_importance <- feature_importance[order(-feature_importance$MeanDecreaseGini), ]
head(feature_importance, 20)
top_features <- feature_importance[1:20, ]
top_features
#create a little plot 
library(ggplot2)

ggplot(top_features, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "pink") +
  coord_flip() +
  ggtitle("Top 20 Important Features") +
  theme_minimal()

#Lets see those results!!!!
cm <- confusionMatrix(pred_final, y_test)
df <- as.data.frame(cm$byClass)

df$Class <- gsub("Class: ", "", rownames(df))

df$F1 <- 2 * (df$`Pos Pred Value` * df$Sensitivity) /
  (df$`Pos Pred Value` + df$Sensitivity)

# Plot 1: Per class F1
ggplot(df, aes(x = Class, y = F1, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Per Class F1 Score") +
  theme_minimal()

# Plot 2: Per class sensitivity
ggplot(df, aes(x = Class, y = Sensitivity, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Per Class Sensitivity") +
  theme_minimal()

#Plot 3: Per class precision
ggplot(df, aes(x = Class, y = `Pos Pred Value`, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Per Class Precision") +
  theme_minimal()

#Plot 4: overall accuracy 
acc_df <- data.frame(
  Metric = "Accuracy",
  Value = as.numeric(cm$overall["Accuracy"])
)

ggplot(acc_df, aes(x = Metric, y = Value)) +
  geom_col(fill = "steelblue") +
  ylim(0, 1) +
  ggtitle("Overall Accuracy") +
  theme_minimal()

#Plot 5: confusion matrix heatmap 
cm_df <- as.data.frame(cm$table)

ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq)) +
  ggtitle("Confusion Matrix") +
  theme_minimal()

#Plot 6: Overfitting Check 
train_pred <- predict(rf_model4, X_train)

overfit_df <- data.frame(
  Set = c("Train", "Test"),
  Accuracy = c(
    mean(train_pred == y_train),
    mean(pred_final == y_test)
  )
)

ggplot(overfit_df, aes(x = Set, y = Accuracy, fill = Set)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Overfitting Check") +
  theme_minimal()

#ROC/AUC
prob <- predict(rf_model4, X_test, type = "prob")

auc_df <- data.frame(
  Class = colnames(prob),
  AUC = sapply(colnames(prob), function(cl) {
    roc(as.numeric(y_test == cl), prob[, cl])$auc
  })
)

ggplot(auc_df, aes(x = Class, y = AUC, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("ROC/AUC by Class") +
  theme_minimal()

########################################
#LODO Model
########################################
# make sure 4-class labels exist
metadata_filtered$label2 <- dplyr::case_when(
  metadata_filtered$disease_class == "HC" ~ "HC",
  metadata_filtered$disease_class %in% c("PA", "PA+") ~ "PA",
  metadata_filtered$disease_class %in% c("CRC", "CRC+") ~ "CRC",
  metadata_filtered$disease_class == "Other" ~ "Other"
)

y2 <- droplevels(as.factor(metadata_filtered$label2))

# check studies being used
table(metadata_filtered$study_name)

studies <- unique(metadata_filtered$study_name)

lodo_results <- lapply(studies, function(test_study) {
  
  test_idx  <- which(metadata_filtered$study_name == test_study)
  train_idx <- which(metadata_filtered$study_name != test_study)
  
  X_train_lodo <- X[train_idx, ]
  X_test_lodo  <- X[test_idx, ]
  y_train_lodo <- y2[train_idx]
  y_test_lodo  <- y2[test_idx]
  
  class_sizes_lodo <- table(y_train_lodo)
  min_class_lodo <- min(class_sizes_lodo)
  
  rf_lodo <- randomForest(
    x = X_train_lodo,
    y = y_train_lodo,
    ntree = 500,
    sampsize = rep(min_class_lodo, length(class_sizes_lodo))
  )
  
  pred_lodo <- predict(rf_lodo, X_test_lodo)
  cm_lodo <- confusionMatrix(pred_lodo, y_test_lodo)
  
  data.frame(
    Study = test_study,
    N_test = length(y_test_lodo),
    Accuracy = as.numeric(cm_lodo$overall["Accuracy"]),
    Kappa = as.numeric(cm_lodo$overall["Kappa"])
  )
})

lodo_results <- do.call(rbind, lodo_results)
lodo_results
mean(lodo_results$Accuracy)
mean(lodo_results$Kappa)