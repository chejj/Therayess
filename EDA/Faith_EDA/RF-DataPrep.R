#######################################################
# Start to combine studies into one matrix for analysis
#######################################################
#----------------
#Handle batch effects
#-----------------  
#batch effect correction
#library(sva)
#combat_data <- ComBat(dat = as.matrix(log_data), batch = metadata$study)



#---------------------------------------------
# Begin with creating a matrix for each study 
#---------------------------------------------
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(randomForest)
library(caret)
library(ggplot2)


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

#--------------------------
# remove CRC-M and PA-M and CRC-H
#--------------------------
metadata_filtered <- metadata[!metadata$disease_class %in% c("CRC-M", "PA-M", "CRC-H"), ]

table(metadata_filtered$disease_class)

#Filter log data to match the removed/filtered metadata_filtered
log_data_filtered <- log_data[, rownames(metadata_filtered), drop = FALSE]
dim(log_data_filtered)

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
set.seed(42)

trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test  <- y[-trainIndex]

length(y_train)
length(y_test)

####################################
#fit random forest model!!!
####################################
set.seed(42)

rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  importance = TRUE
)

print(rf_model)


#evaluate on this test set
pred <- predict(rf_model, X_test)

confusionMatrix(pred, y_test)


##########################
#RF test model 2
###########################
#change focus of model to treat all class sizes equally instead of focusing on HC and CRC
class_sizes <- table(y_train)
class_sizes

min_class <- min(class_sizes[class_sizes > 0])

rf_model_balanced <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 200,
  sampsize = rep(min_class, length(class_sizes)),
  importance = FALSE
)

#now evaluate
pred_bal <- predict(rf_model_balanced, X_test)
confusionMatrix(pred_bal, y_test)

##########################
#Third test model (middle-ground)
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

#re-split + model 
set.seed(42)
trainIndex <- createDataPartition(y2, p = 0.8, list = FALSE)

X_train <- X[trainIndex, ]
X_test  <- X[-trainIndex, ]
y_train <- y2[trainIndex]
y_test  <- y2[-trainIndex]

rf_model2 <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 200,
  importance = FALSE
)

pred2 <- predict(rf_model2, X_test)
confusionMatrix(pred2, y_test)
#This model = best overall accuracy. Great for HC and CRC 
#Weak for PA and Other

#-------------------------------
# building upon model 3 for improvement
#--------------------------------------
class_sizes <- table(y_train)
min_class <- min(class_sizes)

rf_model3 <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 200,
  sampsize = rep(min_class, length(class_sizes))
)

pred3 <- predict(rf_model3, X_test)
confusionMatrix(pred3, y_test)
#This model produces better early disease detection but lower accuracy at 56.7%
#its better for PA detection and other class detection (more fair across classes)


#Saving the current two best models! 
rf_model2  # best accuracy
rf_model3  # best balance

##############################################
# Now extract important features from model 2
##############################################
importance2<- importance(rf_model2)
varImpPlot(rf_model2)


##########################################
# Scale up  
##########################################
rf_model4 <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  sampsize = rep(min_class, length(class_sizes))
)
pred_final <- predict(rf_model4, X_test)
confusionMatrix(pred_final, y_test)

#Lets see those results!!!!
cm <- confusionMatrix(pred_final, y_test)
df <- as.data.frame(cm$byClass)

df$Class <- gsub("Class: ", "", rownames(df))

df$F1 <- 2 * (df$`Pos Pred Value` * df$Sensitivity) /
  (df$`Pos Pred Value` + df$Sensitivity)

# Plot 1
ggplot(df, aes(x = Class, y = F1, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Per-Class F1 Score") +
  theme_minimal()

# Plot 2
ggplot(df, aes(x = Class, y = Sensitivity, fill = Class)) +
  geom_col() +
  ylim(0, 1) +
  ggtitle("Per-Class Sensitivity") +
  theme_minimal()

