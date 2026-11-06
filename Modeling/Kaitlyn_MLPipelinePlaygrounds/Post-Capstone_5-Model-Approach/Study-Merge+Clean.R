# Merge and Clean Studies prior to Machine Learning Algorithms, variable Setup

## SETUP #######################################################################
suppressPackageStartupMessages({
  library(curatedMetagenomicData)  # for mergeData
  library(SummarizedExperiment)
  library(TreeSummarizedExperiment)
  library(tidyverse)
})

## Outline #####################################################################
# Step 1: Merge studies and create Objects
# Step 2: Separate out into Taxa, Pathway Abundance, Pathway Coverage, and Metadata, set as samples x features
# Step 3: Clean up metadata -> Create disease class and correct age classifiers

## Step 1: Merge Studies and create objects ----
types <- c("relative_abundance", "pathway_coverage", "pathway_abundance")

merged_studies <- lapply(types, function(type) {
  objs <- lapply(CRC_progression_studies, `[[`, type)
  objs <- objs[!vapply(objs, is.null, logical(1))]
  mergeData(objs)
})

# Name the list
names(merged_studies) <- types

## Step 2: Samples are rows, features are columns (samples x features) ----
studies_taxa <- as.data.frame(t(assay(merged_studies[["relative_abundance"]])))
studies_path_abund <- as.data.frame(t(assay(merged_studies[["pathway_abundance"]])))
studies_path_cov <- as.data.frame(t(assay(merged_studies[["pathway_coverage"]])))
studies_meta_df <- as.data.frame(colData(merged_studies[["relative_abundance"]]))

## Step 3: Metadata Cleanup ----
# Transform metadata to have disease class and correct age classifiers 
studies_meta_df <- studies_meta_df %>%
  mutate(
    disease_class = case_when(
      disease == "healthy" ~ "HC",
      str_detect(disease, "\\b(adenoma|polyp)\\b") & str_detect(disease, "metasta") ~ "PA-M",
      str_detect(disease, "\\bCRC\\b") & str_detect(disease, "metasta") ~ "CRC-M",
      str_detect(disease, "history") ~ "CRC-H",
      disease == "CRC" ~ "CRC",
      str_detect(disease, "\\bCRC\\b") ~ "CRC+",
      disease %in% c("adenoma", "few_polyps") ~ "PA",
      str_detect(disease, "\\b(adenoma|polyp)\\b") ~ "PA+",
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    disease_class = factor(
      disease_class,
      levels = c("Other", "HC", "PA", "PA+", "PA-M", "CRC", "CRC+", "CRC-M", "CRC-H")
    ), 
    age_decade = case_when(
      age < 10 ~ "0-9", 
      age < 20 ~ "10-19",
      age < 30 ~ "20-29",
      age < 40 ~ "30-39",
      age < 50 ~ "40-49",
      age < 60 ~ "50-59",
      age < 70 ~ "60-69",
      age < 80 ~ "70-79",
      TRUE  ~ "80+"
    ),
    age_decade = if_else(age_decade %in% c("10-19", "20-29"), "18-29", age_decade),
    age_decade = factor(
      age_decade,
      levels = c("18-29","30-39","40-49","50-59","60-69","70-79","80+"),
      ordered = TRUE
    )
  )