library("ranger")
library(SummarizedExperiment)
library(dplyr)

# Pick one study
study <- CRC_progression_studies[["ThomasAM_2018a"]]

# Pull assays
test_taxa <- assay(study[["relative_abundance"]])
test_path_abund <- assay(study[["pathway_abundance"]])
test_path_cov <- assay(study[["pathway_coverage"]])

# Transpose so rows = samples, cols = features
# This is needed because RF expects rows to be samples, but SEs are not set up this way
taxa_df <- as.data.frame(t(test_taxa))
path_abund_df <- as.data.frame(t(test_path_abund))
path_cov_df <- as.data.frame(t(test_path_cov))

# Add prefixes so feature names do not collide
colnames(taxa_df) <- paste0("taxa__", make.names(colnames(taxa_df), unique = TRUE))
colnames(path_abund_df) <- paste0("pabund__", make.names(colnames(path_abund_df), unique = TRUE))
colnames(path_cov_df) <- paste0("pcov__", make.names(colnames(path_cov_df), unique = TRUE))

# Metadata
test_meta_df <- as.data.frame(colData(study[["relative_abundance"]]))

# Check sample alignment
stopifnot(
  identical(rownames(taxa_df), rownames(test_meta_df)),
  identical(rownames(path_abund_df), rownames(test_meta_df)),
  identical(rownames(path_cov_df), rownames(test_meta_df))
)

keep_taxa <- colMeans(taxa_df > 0) > 0.05
keep_pabund <- colMeans(path_abund_df > 0) > 0.05
keep_pcov <- colMeans(path_cov_df > 0) > 0.05

test_df_small <- bind_cols(
  test_meta_df,
  taxa_df[, keep_taxa, drop = FALSE],
  path_abund_df[, keep_pabund, drop = FALSE],
  path_cov_df[, keep_pcov, drop = FALSE]
)

dim(test_df_small)

# Remove metadata columns you do not want as predictors
rf_input <- test_df_small %>% 
  filter(disease_class %in% c("HC", "PA", "CRC", "Other", "PA+", "CRC+")) %>% 
  filter(country == "JPN") %>% 
  select(-study_name, -study_condition, -non_westernized, -sequencing_platform, -PMID, -curator)

rf_input <- test_df_small %>%
  filter(!is.na(disease_class)) %>% # Keep only rows with known outcome
  mutate(disease_class = as.factor(disease_class)) %>% # Make outcome a factor
  select(disease_class, starts_with("taxa__"), starts_with("pabund__"), starts_with("pcov__"))

# Fit random forest
rf_fit <- ranger(
  disease_class ~ .,
  data = rf_input,
  num.trees = 500,
  importance = "impurity",
  probability = TRUE
)

rf_fit
