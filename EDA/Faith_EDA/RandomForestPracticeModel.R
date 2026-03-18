#Random forest modeling testing
#starting with a single study and data type to test if it works
#only using HC vs CRC to test 

#####################################################################
#Count samples per study to find what study is best for testing model
#####################################################################
study_counts <- sapply(CRC_progression_studies, function(st) {
  x <- st[["relative_abundance"]]
  if (is.null(x)) return(NA)
  
  table(colData(x)$disease_class)
})
#clean up this messy summary table to find what is our best test model study
study_summary <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  
  x <- CRC_progression_studies[[study]][["relative_abundance"]]
  if (is.null(x)) return(NULL)
  
  counts <- table(colData(x)$disease_class)
  
  data.frame(
    study = study,
    HC  = ifelse("HC" %in% names(counts), counts["HC"], 0),
    CRC = ifelse("CRC" %in% names(counts), counts["CRC"], 0)
  )
}))

study_summary
#filter for selected studies
good_studies <- study_summary %>%
  filter(HC >= 20, CRC >= 20)

good_studies
#now select the best study from the group
good_studies %>%
  mutate(total = HC + CRC) %>%
  arrange(desc(total))

#############################################################
# Visualize the best studies for practice training the model
#############################################################
library(ggplot2)

study_summary %>%
  tidyr::pivot_longer(cols = c(HC, CRC), names_to = "class", values_to = "count") %>%
  ggplot(aes(x = study, y = count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
## Random Forest starter model for CRC vs Healthy Controls
## Using one study + one product first
################################################################################

library(SummarizedExperiment)
library(dplyr)
library(randomForest)
library(caret)

# -----------------------------
# 1) Choose one study and one data type
# -----------------------------
study_name <- "YachidaS_2019"

# -----------------------------
# 2) Extract feature matrix + metadata
# -----------------------------
se <- CRC_progression_studies[[study_name]][["relative_abundance"]]

mat  <- assay(se)
meta <- as.data.frame(colData(se))

# transpose so rows = samples, columns = microbial features
model_df <- as.data.frame(t(mat))

# add disease labels from metadata
model_df$disease_class <- meta$disease_class

# keep only Healthy Controls and CRC
model_df <- model_df %>%
  filter(disease_class %in% c("HC", "CRC")) %>%
  mutate(disease_class = factor(disease_class))
#clean column names for modeling
colnames(model_df) <- make.names(colnames(model_df), unique = TRUE)
#check if clean worked 
head(colnames(model_df))

#Table
table(model_df$disease_class)
dim(model_df)

# -----------------------------
# 4) Optional cleanup
#    Remove features with zero variance
# -----------------------------
nzv <- nearZeroVar(model_df[, !names(model_df) %in% "disease_class"])
if (length(nzv) > 0) {
  model_df <- model_df[, -nzv]
}

# -----------------------------
# 5) Train/test split
# -----------------------------
set.seed(123)

train_indices <- createDataPartition(model_df$disease_class, p = 0.7, list = FALSE)
train_data <- model_df[train_indices, ]
test_data  <- model_df[-train_indices, ]

# -----------------------------
# 6) Train Random Forest
# -----------------------------
rf_model <- randomForest(
  disease_class ~ .,
  data = train_data,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data) - 1)),
  importance = TRUE,
  na.action = na.omit
)

print(rf_model)

# -----------------------------
# 7) Predict on test set
# -----------------------------
predictions <- predict(rf_model, newdata = test_data)

# -----------------------------
# 8) Evaluate performance
# -----------------------------
confusionMatrix(predictions, test_data$disease_class)

# -----------------------------
# 9) Variable importance
# -----------------------------
importance(rf_model)
varImpPlot(rf_model)

# count how many samples in each class in training data
table(train_data$disease_class)

# see if test set is balanced
table(test_data$disease_class)


############################################################
# creating a confusion matrix heatmap 
############################################################

library(ggplot2)

cm <- as.data.frame(confusionMatrix(predictions, test_data$disease_class)$table)

ggplot(cm, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() +
  labs(title = "Confusion Matrix: Random Forest (HC vs CRC)")

#what the confusion matrix tells us: 
# model demonstrates moderate classification performance ~68% accuracy 
# with higher specificity than sensistivity insufficient alone for CRC detection

############################################################
# variable importance plot
############################################################
imp <- importance(rf_model)
imp_df <- data.frame(
  feature = rownames(imp),
  importance = imp[,1]
)

imp_df %>%
  arrange(desc(importance)) %>%
  slice(1:20) %>%  # top 20 features
  ggplot(aes(x = reorder(feature, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Top Microbial Features Driving Classification",
       x = "Feature",
       y = "Importance")

#########################################################
# predicted vs actual distribution
#########################################################
results_df <- data.frame(
  Actual = test_data$disease_class,
  Predicted = predictions
)

ggplot(results_df, aes(x = Actual, fill = Predicted)) +
  geom_bar(position = "dodge") +
  theme_bw() +
  labs(title = "Prediction Distribution by True Class")

############################################################
# probability distribution plot
############################################################
probs <- predict(rf_model, newdata = test_data, type = "prob")

prob_df <- data.frame(
  Actual = test_data$disease_class,
  Prob_CRC = probs[, "CRC"]
)

ggplot(prob_df, aes(x = Prob_CRC, fill = Actual)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(title = "Model Confidence (Probability of CRC)")



