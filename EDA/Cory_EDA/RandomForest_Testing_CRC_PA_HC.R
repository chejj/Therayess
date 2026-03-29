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

# filter out all 
combined_df_genus <- combined_df_genus %>%
  filter(disease_class %in% c("HC", "PA", "CRC")) %>%
  mutate(disease_class = droplevels(disease_class))  # or factor(disease_class)

# Now table will only show HC, PA, CRC columns
table(combined_df_genus$study_name, combined_df_genus$disease_class)


#Separate meta data from the rest of the metagenomic data

meta_df <- combined_df_genus %>%  select(71:last_col())
meta_df
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


# Use disease or disease_class (whichever has multiple classes you want)
y <- as.factor(meta_df$disease)  # or meta_df$disease_class

# Metadata predictors
meta_x2 <- meta_df[, c("age_decade", "gender", "BMI"), drop = FALSE]

# Convert types
meta_x2$age_decade <- as.factor(meta_x2$age_decade)
meta_x2$gender <- as.factor(meta_x2$gender)
meta_x2$BMI <- as.numeric(meta_x2$BMI)

# Combine 
x <- cbind(t(relative_df_adj), meta_x2)

# Check complete cases BEFORE filtering
cat("Before filtering - samples:", nrow(x), "y values:", length(y), "\n")
cat("y classes:", table(y), "\n")

# Filter
keep <- complete.cases(x) & complete.cases(y)
x <- x[keep, , drop = FALSE] 
y <- y[keep]

# Check AFTER filtering
cat("\nAfter filtering - samples:", nrow(x), "y values:", length(y), "\n")
cat("y classes:", table(y), "\n")

# Remove constant columns
is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

# Now try RF
library(randomForest)
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 10000, importance = TRUE)

print(rf_fit)


###################################################################
# HC and CRC
###################################################################
# filter out all 
combined_df_genus <- combined_df_genus %>%
  filter(disease_class %in% c("HC", "CRC")) %>%
  mutate(disease_class = droplevels(disease_class))  # or factor(disease_class)

# Now table will only show HC, PA, CRC columns
table(combined_df_genus$study_name, combined_df_genus$disease_class)


#Separate meta data from the rest of the metagenomic data

meta_df <- combined_df_genus %>%  select(71:last_col())
meta_df
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


# Use disease or disease_class (whichever has multiple classes you want)
y <- as.factor(meta_df$disease)  # or meta_df$disease_class

# Metadata predictors
meta_x2 <- meta_df[, c("age_decade", "gender", "BMI"), drop = FALSE]

# Convert types
meta_x2$age_decade <- as.factor(meta_x2$age_decade)
meta_x2$gender <- as.factor(meta_x2$gender)
meta_x2$BMI <- as.numeric(meta_x2$BMI)

# Combine. Use t(if relative abundance needs to be transposed)
x <- cbind(relative_df, meta_x2)

# Check complete cases BEFORE filtering
cat("Before filtering - samples:", nrow(x), "y values:", length(y), "\n")
cat("y classes:", table(y), "\n")

# Filter
keep <- complete.cases(x) & complete.cases(y)
x <- x[keep, , drop = FALSE] 
y <- y[keep]

# Check AFTER filtering
cat("\nAfter filtering - samples:", nrow(x), "y values:", length(y), "\n")
cat("y classes:", table(y), "\n")

# Remove constant columns
is_constant <- sapply(x, function(v) length(unique(v[!is.na(v)])) <= 1)
x <- x[, !is_constant, drop = FALSE]

# Now try RF
library(randomForest)
set.seed(1)
rf_fit <- randomForest(x = x, y = y, ntree = 10000, importance = TRUE)

print(rf_fit)
 