##############################################
# Build one combined dataframe for all studies
##############################################
install.packages("curl")

library(SummarizedExperiment)
library(dplyr)
library(curl)
library(lime)
library(MMUPHin)
library(ggplot2)


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

# 2. Build one combined data frame at GENUS level
combined_list_genus <- lapply(names(CRC_progression_studies), function(study) {
  se <- CRC_progression_studies[[study]][[product_type]]
  if (is.null(se)) return(NULL)
  
  # keep only shared features
  mat <- assay(se)[common_features, , drop = FALSE]
  
  # extract genus from row names, then collapse duplicate genera
  genus_names <- sub(".*\\|g__([^|]+).*", "\\1", rownames(mat))
  mat_genus <- rowsum(mat, group = genus_names)
  
  meta <- as.data.frame(colData(se))
  
  df <- as.data.frame(t(mat_genus))
  df$disease_class <- meta$disease_class
  df$study_name <- meta$study_name
  df$studyID <- meta$studyID
  df$sequencing_platform <- meta$sequencing_platform
  df$DNA_extraction_kit <- meta$DNA_extraction_kit
  df$subject_id <- meta$subject_id
  df$body_site <- meta$body_site
  df$study_condition <- meta$study_conidition
  df$disease <- meta$disease
  df$age <- meta$age
  df$age_category <- meta$age_category
  df$gender <- meta$gender
  df$country <- meta$country
  df$sequencing_platform <- meta$sequencing_platform
  df$PMID <- meta$PMID
  df$number_reads <- meta$number_reads
  df$number_bases <- meta$number_bases 
  df$minimum_read_length <- meta$minimum_read_length
  df$median_read_length <- meta$median_read_length
  df$BMI <- meta$BMI
  df$age_decade <- meta$age_decade
  df$disease_class <- meta$disease_class
  
  df
})

combined_df_genus <- bind_rows(combined_list_genus)

dim(combined_df_genus)
table(combined_df_genus$disease_class)
table(combined_df_genus$study_name)


# Filter country for Japan
#combined_df_genus <- combined_df_genus %>%  filter(country == "JPN")

# Filter CRC and HC
#combined_df_genus <- combined_df_genus %>% filter(disease_class %in% c("HC", "CRC"))


# Response variable
y <- as.factor(combined_df_genus$disease)

# Pick the metadata columns you want to include as predictors
meta_vars <- c("age_decade", "gender", "BMI")
meta_x <- meta_df[, intersect(meta_vars, colnames(meta_df)), drop = FALSE]

# Convert types
meta_x[] <- lapply(meta_x, function(v) {
  if (is.character(v)) as.factor(v) else v
})

num_vars <- intersect(c("age", "BMI"), colnames(meta_x))
meta_x[num_vars] <- lapply(meta_x[num_vars], function(v) as.numeric(as.character(v)))

# Create factors for desired variables
meta_x$age_decade <- as.factor(meta_x$age_decade)
meta_x$gender <- as.factor(meta_x$gender)
meta_x$BMI <- as.numeric(meta_x$BMI)

str(meta_x)

# Actually running the model
keep <- complete.cases(x) & !is.na(y)
x <- x[keep, , drop = FALSE]
y <- y[keep]

is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

library(randomForest)
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 500, importance = TRUE)

print(rf_fit)


# Plot results
cm <- rf_fit$confusion[, 1:2]
cm_prop <- prop.table(cm, 1)

cm_df <- as.data.frame(as.table(cm_prop))
colnames(cm_df) <- c("Actual", "Predicted", "Proportion")

ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Proportion)) +
  geom_tile() +
  geom_text(aes(label = round(Proportion, 2)), size = 5) +
  labs(title = "Confusion Matrix (Row Proportions)") +
  theme_minimal()




# Chat GPT Testing
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

# 2. Build one combined dataframe at GENUS level
combined_list_genus <- lapply(names(CRC_progression_studies), function(study) {
  se <- CRC_progression_studies[[study]][[product_type]]
  if (is.null(se)) return(NULL)
  
  # keep only shared features
  mat <- assay(se)[common_features, , drop = FALSE]
  
  # collapse rows to genus
  genus_names <- sub(".*\\|g__([^|]+).*", "\\1", rownames(mat))
  mat_genus <- rowsum(mat, group = genus_names)
  
  meta <- as.data.frame(colData(se))
  
  # samples x genera
  df <- as.data.frame(t(mat_genus))
  
  # attach metadata
  df$subject_id <- meta$subject_id
  df$study_name <- meta$study_name
  df$body_site <- meta$body_site
  df$study_condition <- meta$study_condition
  df$disease <- meta$disease
  df$age <- meta$age
  df$age_category <- meta$age_category
  df$gender <- meta$gender
  df$country <- meta$country
  df$sequencing_platform <- meta$sequencing_platform
  df$DNA_extraction_kit <- meta$DNA_extraction_kit
  df$PMID <- meta$PMID
  df$number_reads <- meta$number_reads
  df$number_bases <- meta$number_bases
  df$minimum_read_length <- meta$minimum_read_length
  df$median_read_length <- meta$median_read_length
  df$BMI <- meta$BMI
  df$age_decade <- meta$age_decade
  df$disease_class <- meta$disease_class
  
  df
})

combined_list_genus <- Filter(Negate(is.null), combined_list_genus)
combined_df_genus <- bind_rows(combined_list_genus)

dim(combined_df_genus)
table(combined_df_genus$disease_class)
table(combined_df_genus$study_name)

# 3. Filter to Japan
combined_df_genus <- combined_df_genus %>%
  filter(country == "JPN")

# 4. Filter to HC and CRC only
combined_df_genus <- combined_df_genus %>%
  filter(disease_class %in% c("HC", "CRC")) %>%
  droplevels()

# 4.2 Filter age_decade to 50 and above
combined_df_genus <- combined_df_genus %>% 
  filter(age_decade %in% c("50-59", "60-69", "70-79"))
  
# 5. Response variable
y <- droplevels(as.factor(combined_df_genus$disease_class))

# 6. Metadata predictors
meta_vars <- c("age_decade", "gender", "BMI")
meta_x <- combined_df_genus[, intersect(meta_vars, colnames(combined_df_genus)), drop = FALSE]

meta_x[] <- lapply(meta_x, function(v) {
  if (is.character(v)) as.factor(v) else v
})

if ("BMI" %in% colnames(meta_x)) {
  meta_x$BMI <- as.numeric(as.character(meta_x$BMI))
}
if ("age_decade" %in% colnames(meta_x)) {
  meta_x$age_decade <- as.factor(meta_x$age_decade)
}
if ("gender" %in% colnames(meta_x)) {
  meta_x$gender <- as.factor(meta_x$gender)
}

str(meta_x)

# 7. Genus abundance predictors
meta_cols <- c(
  "subject_id", "study_name", "body_site", "study_condition", "disease",
  "age", "age_category", "gender", "country", "sequencing_platform",
  "DNA_extraction_kit", "PMID", "number_reads", "number_bases",
  "minimum_read_length", "median_read_length", "BMI", "age_decade",
  "disease_class"
)

genus_x <- combined_df_genus[, !colnames(combined_df_genus) %in% meta_cols, drop = FALSE]

# 8. Combine genus features + metadata predictors
x <- cbind(genus_x, meta_x)

# 9. Keep complete cases
keep <- complete.cases(x) & !is.na(y)
x <- x[keep, , drop = FALSE]
y <- droplevels(y[keep])

# 10. Remove constant columns
is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

# 11. Check classes before fitting
table(y)
levels(y)

# 12. Fit random forest
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 10000, importance = TRUE)

print(rf_fit)


# Plot results
cm <- rf_fit$confusion[, 1:2]
cm_prop <- prop.table(cm, 1)

cm_df <- as.data.frame(as.table(cm_prop))
colnames(cm_df) <- c("Actual", "Predicted", "Proportion")

ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Proportion)) +
  geom_tile() +
  geom_text(aes(label = round(Proportion, 2)), size = 5) +
  labs(title = "Confusion Matrix (Row Proportions)") +
  theme_minimal()






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



# Flatten Lists so they can be filtered
all_tse <- unlist(CRC_progression_studies, recursive = FALSE)

# Filter by variable
jpn_tse <- lapply(all_tse, function(x) {
  cd <- colData(x)
  
  if (!"country" %in% colnames(cd)) {
    return(NULL)
  }
  
  keep <- cd$country == "JPN"
  x[, keep, drop = FALSE]
})

jpn_tse <- Filter(Negate(is.null), jpn_tse)

jpn_tse <- Filter(function(x) ncol(x) > 0, jpn_tse)

# Double check to see if filter is correct
names(jpn_tse)
sapply(jpn_tse, ncol)

# This confirms proper country
unique(unlist(lapply(jpn_tse, function(x) as.character(colData(x)$country))))

# Filter down to relative abundance or any other set of data
jpn_ra <- jpn_tse[grep("relative_abundance$", names(jpn_tse))]

# Double check results
names(jpn_ra)
sapply(jpn_ra, dim)


# Combine into one feature table
common_features <- Reduce(intersect, lapply(jpn_ra, rownames))
length(common_features)

# Combine abundance data
common_features <- Reduce(intersect, lapply(jpn_ra, rownames))

jpn_ra2 <- lapply(jpn_ra, function(x) x[common_features, , drop = FALSE])

feature_mat <- do.call(cbind, lapply(jpn_ra2, assay))
feature_df <- as.data.frame(t(feature_mat))   # samples x taxa

# Combine metadata
meta_df <- bind_rows(Map(function(x, nm) {
  df <- as.data.frame(colData(x))
  df$dataset <- nm
  df
}, jpn_ra2, names(jpn_ra2)))

# Drop unwanted columns
meta_drop <- c("non_westernized","study_name", "curator", "ajcc", "dataset", 
               "disease_location", "brinkman_index", "alcohol_numeric")

meta_df2 <- meta_df[, !colnames(meta_df) %in% meta_drop, drop = FALSE]



# Combine feature_df and meta_df2
# Response variable
y <- as.factor(meta_df2$disease)

# Pick the metadata columns you want to include as predictors
meta_vars <- c("age_decade", "gender", "BMI", "study_condition")
meta_x <- meta_df[, intersect(meta_vars, colnames(meta_df)), drop = FALSE]

# Convert types
meta_x[] <- lapply(meta_x, function(v) {
  if (is.character(v)) as.factor(v) else v
})

num_vars <- intersect(c("age", "BMI"), colnames(meta_x))
meta_x[num_vars] <- lapply(meta_x[num_vars], function(v) as.numeric(as.character(v)))

# Create factors for desired variables
meta_x$age_decade <- as.factor(meta_x$age_decade)
meta_x$gender <- as.factor(meta_x$gender)
meta_x$BMI <- as.numeric(meta_x$BMI)

str(meta_x)

# This ensures study_condition is used as a predictor
y <- meta_x$study_condition

meta_x2 <- meta_x[, c("age_decade", "gender", "BMI"), drop = FALSE]

# Combines relative abundance data and metadata
x <- cbind(feature_df, meta_x2)

# Actually running the model
keep <- complete.cases(x) & !is.na(y)
x <- x[keep, , drop = FALSE]
y <- y[keep]

is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

library(randomForest)
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 500, importance = TRUE)

print(rf_fit)

# Check table counts
table(y)
prop.table(table(y))

# Drop small proportion of adenoma and carcinoma
keep <- y %in% c("CRC", "control")
x_bin <- x[keep, , drop = FALSE]
y_bin <- droplevels(y[keep])

table(y_bin)
prop.table(table(y_bin))

# Rerun and fit
set.seed(1)
rf_bin <- randomForest(
  x = x_bin,
  y = y_bin,
  ntree = 500,
  importance = TRUE
)

print(rf_bin)
varImpPlot(rf_bin)

# Check how the model performed
pred <- predict(rf_bin, x_bin)

correct_crc <- which(y_bin == "CRC" & pred == "CRC")
wrong_crc   <- which(y_bin == "CRC" & pred == "control")
wrong_ctrl  <- which(y_bin == "control" & pred == "CRC")

# First ensure lime understand the type of model
model_type.randomForest <- function(x, ...) {
  "classification"
}

predict_model.randomForest <- function(x, newdata, type, ...) {
  as.data.frame(predict(x, newdata = newdata, type = "prob"))
}

# Create explainer
explainer <- lime(x_bin, rf_bin)

# Create explanation
explanation_crc <- explain(
  x_bin[correct_crc[5], , drop = FALSE],
  explainer,
  n_features = 10,
  n_labels = 1
)

plot_features(explanation_crc)
sort(table(explanation_crc$feature), decreasing = TRUE)

# For global summary
library(iml)

pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  predict(model, newdata = newdata, type = "prob")[, "CRC"]
}

pred <- Predictor$new(
  model = rf_bin,
  data = as.data.frame(x_bin),
  y = y_bin,
  predict.function = pred_fun
)

imp <- FeatureImp$new(pred, loss = "ce")
plot(imp)

# Try going to genus level before 

# Testing safe approach and refitting 
orig_names <- colnames(x_bin)
safe_names <- make.names(orig_names, unique = TRUE)

name_key <- data.frame(
  safe = safe_names,
  original = orig_names,
  stringsAsFactors = FALSE
)

colnames(x_bin) <- safe_names

# Refitting
set.seed(1)
rf_bin_safe <- randomForest(
  x = x_bin,
  y = y_bin,
  ntree = 500,
  importance = TRUE
)

# Testing global explainer
pred_fun <- function(model, newdata) {
  newdata <- as.data.frame(newdata)
  predict(model, newdata = newdata, type = "prob")[, "CRC"]
}

pred <- Predictor$new(
  model = rf_bin_safe,
  data = as.data.frame(x_bin),
  y = y_bin,
  predict.function = pred_fun
)

imp <- FeatureImp$new(pred, loss = "ce")
plot(imp)

# Summarize results from iml
# feature-importance table
imp_df <- imp$results

# top features by median importance
imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]
head(imp_df, 20)


#############################################
# model leaving out one study for testing
#############################################
library(randomForest)
library(caret)

study_ids <- unique(jpn_tse$study_name)

loso_results <- lapply(jpn_tse, function(test_study) {
  
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



#Separate meta data from the rest of the metagenomic data

meta_df <- combined_df_genus %>%  select(71:last_col())
relative_df <- combined_df_genus %>%  select(1:70)


fit_adjust_batch <- adjust_batch(feature_abd = t(relative_df/100),
                                 batch = "study_name",
                                 covariates = "disease_class",
                                 data = meta_df,
                                 control = list(verbose = FALSE))

relative_df_adj <- fit_adjust_batch$feature_abd_adj

library(vegan, quietly = TRUE)

D_before <- vegdist(relative_df)
D_after <- vegdist(t(relative_df_adj))

set.seed(1)
fit_adonis_before <- adonis2(D_before ~ study_name, data = meta_df)
fit_adonis_after <- adonis2(D_after ~ study_name, data = meta_df)
print(fit_adonis_before)

print(fit_adonis_after)

# Combine metadata with adjusted relative abundance data
combined_adj <- cbind(t(relative_df_adj), meta_df)


# PCA before batch correction
pca_before <- prcomp(relative_df, scale. = TRUE)

# PCA after batch correction
pca_after <- prcomp(t(relative_df_adj), scale. = TRUE)

# Create a data frame with PCA results and metadata for plotting
pca_before_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  study_name = meta_df$study_name,
  disease_class = meta_df$disease_class,
  batch_status = "Before Adjustment"
)

pca_after_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  study_name = meta_df$study_name,
  disease_class = meta_df$disease_class,
  batch_status = "After Adjustment"
)

# Combine for faceted plotting
pca_combined <- rbind(pca_before_df, pca_after_df)

# Visualize
library(ggplot2)

pdf("pca_batch_correction.pdf", width = 12, height = 5)

ggplot(pca_combined, aes(x = PC1, y = PC2, color = study_name, shape = disease_class)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~batch_status) +
  theme_minimal() +
  labs(title = "PCA: MMUPHin Batch Correction Effect",
       color = "Study", shape = "Disease Class")

dev.off()


y <- as.factor(meta_df$disease)

# Pick the metadata columns you want to include as predictors
meta_vars <- c("age_decade", "gender", "BMI", "study_condition")
meta_x <- meta_df[, intersect(meta_vars, colnames(meta_df)), drop = FALSE]

# Convert types
meta_x[] <- lapply(meta_x, function(v) {
  if (is.character(v)) as.factor(v) else v
})

num_vars <- intersect(c("age", "BMI"), colnames(meta_x))
meta_x[num_vars] <- lapply(meta_x[num_vars], function(v) as.numeric(as.character(v)))

# Create factors for desired variables
meta_x$age_decade <- as.factor(meta_x$age_decade)
meta_x$gender <- as.factor(meta_x$gender)
meta_x$BMI <- as.numeric(meta_x$BMI)

str(meta_x)

# This ensures study_condition is used as a predictor
y <- meta_x$study_condition

meta_x2 <- meta_x[, c("age_decade", "gender", "BMI"), drop = FALSE]

# Combines relative abundance data and metadata
x <- cbind(t(relative_df_adj), meta_df)

# Actually running the model
keep <- complete.cases(x) & !is.na(y)
x <- x[keep, , drop = FALSE] 
y <- y[keep]

is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

library(randomForest)
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 10000, importance = TRUE)

print(rf_fit)
