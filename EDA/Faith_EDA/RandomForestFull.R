##############################################
# Build one combined dataframe for all studies
##############################################
library(SummarizedExperiment)
library(dplyr)

product_type <- "relative_abundance"

# 1. Find common microbial features across all studies
common_features <- Reduce(
  intersect,
  lapply(CRC_progression_studies, function(st) {
    x <- st[[product_type]]
    if (is.null(x)) return(NULL)
    rownames(assay(x))
  })
)

length(common_features)

# 2. Build one combined dataframe
combined_list <- lapply(names(CRC_progression_studies), function(study) {
  se <- CRC_progression_studies[[study]][[product_type]]
  if (is.null(se)) return(NULL)
  
  # keep only shared features
  mat <- assay(se)[common_features, , drop = FALSE]
  meta <- as.data.frame(colData(se))
  
  df <- as.data.frame(t(mat))
  df$disease_class <- meta$disease_class
  df$study_name <- meta$study_name
  df$sequencing_platform <- meta$sequencing_platform
  df$DNA_extraction_kit <- meta$DNA_extraction_kit
  
  df
})

combined_df <- bind_rows(combined_list)

dim(combined_df)
table(combined_df$disease_class)
table(combined_df$study_name)

#####################################
# restrict to classes we want
#####################################

combined_df <- combined_df %>%
  filter(disease_class %in% c("HC", "CRC")) %>%
  mutate(
    disease_class = factor(disease_class),
    study_name = factor(study_name),
    sequencing_platform = factor(sequencing_platform),
    DNA_extraction_kit = factor(DNA_extraction_kit)
  )

# optional quick check
combined_df %>%
  count(study_name, disease_class)

#######################################
# make column names formula safe
#######################################
colnames(combined_df) <- make.names(colnames(combined_df), unique = TRUE)

#############################################
# model leaving out one study for testing
#############################################
library(randomForest)
library(caret)

study_ids <- unique(combined_df$study_name)

loso_results <- lapply(study_ids, function(test_study) {
  
  train_data <- combined_df %>% filter(study_name != test_study)
  test_data  <- combined_df %>% filter(study_name == test_study)
  
  # clean test and train data
  train_data_clean <- na.omit(train_data)
  test_data_clean  <- na.omit(test_data)
  
  #debugging step for yachida study 
  print(test_study)
  print(table(test_data$disease_class))
  print(table(test_data_clean$disease_class))
  
  # skip studies that don't contain both classes after cleaning
  if (length(unique(test_data_clean$disease_class)) < 2) return(NULL)
  
  # optional: remove study_name from predictors if you want pure biology + technical covariates
  rf_model <- randomForest(
    disease_class ~ . - study_name,
    data = train_data_clean,
    ntree = 500,
    mtry = floor(sqrt(ncol(train_data_clean) - 1)),
    importance = TRUE
  )
  
  preds <- predict(rf_model, newdata = test_data_clean)
  cm <- confusionMatrix(preds, test_data_clean$disease_class)
  
  data.frame(
    test_study = as.character(test_study),
    accuracy = cm$overall["Accuracy"],
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"]
  )
})

loso_results_df <- bind_rows(loso_results)
loso_results_df

#check if all studies ran or some were skipped
nrow(loso_results_df)
unique(loso_results_df$test_study)
#get a quick summary 
summary(loso_results_df)
mean(loso_results_df$accuracy, na.rm = TRUE)

#plot of study accuracy
library(ggplot2)

ggplot(loso_results_df, aes(x = reorder(test_study, accuracy), y = accuracy)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Leave-One-Study-Out Accuracy",
    x = "Held-out study",
    y = "Accuracy"
  )

#what studies were the best, what was the worst?
loso_results_df %>% arrange(desc(accuracy))
loso_results_df %>% arrange(accuracy)


combined_df %>%
  filter(study_name == "YachidaS_2019") %>%
  count(disease_class)

#sanity check after running 
unique(loso_results_df$test_study)



################################
# solving yachida study issue
###############################

yachida_test <- combined_df %>%
  filter(study_name == "YachidaS_2019")

nrow(yachida_test)
nrow(na.omit(yachida_test))

table(yachida_test$disease_class)
table(na.omit(yachida_test)$disease_class)







