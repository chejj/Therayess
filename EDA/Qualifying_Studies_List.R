################################################################################
## ---- 1) Load packages
suppressPackageStartupMessages({
  library(BiocManager)
  library(curatedMetagenomicData)
  library(dplyr)
  library(stringr)
})
#View(sampleMetadata)
str(sampleMetadata)
################################################################################
## ---- 2) Quick environment sanity check
cat("\n=== Session Info ===\n")
cat("R:", R.version.string, "\n")
cat("Bioc:", as.character(BiocManager::version()), "\n")
cat("Library paths:\n")
print(.libPaths())
cat("\ncuratedMetagenomicData:", as.character(packageVersion("curatedMetagenomicData")), "\n\n")

sampleMetadata <- sampleMetadata %>% 
  mutate(
    age_decade = case_when(
      age < 10 ~ "0-9", 
      age > 9 & age < 20 ~ "10-19",
      age > 19 & age < 30 ~ "20-29",
      age > 29 & age < 40 ~ "30-39",
      age > 39 & age < 50 ~ "40-49",
      age > 49 & age < 60 ~ "50-59",
      age > 59 & age < 70 ~ "60-69",
      age > 69 & age < 80 ~ "70-79",
      age > 79 & age < 90 ~ "80-89",
      age > 89 ~ "90+"
    )
  ) 
################################################################################
qualifying_studies <- sampleMetadata %>%
  filter(body_site == "stool", age >= 18, !is.na(disease)) %>%
  group_by(study_name) %>%
  filter(
    any(str_detect(disease, "\\bCRC\\b") | str_detect(disease, "polyp|adenoma")) |
      any(study_name == "HMP_2012")
  ) %>%
  ungroup() %>%
  mutate(
    disease_class = case_when(
      disease == "healthy" ~ "HC",
      disease %in% c("adenoma", "few_polyps") ~ "PA",
      str_detect(disease, "\\b(adenoma|polyp)\\b") ~ "PA+",
      str_detect(disease, "\\b(adenoma|polyp)\\b") & str_detect(disease, "metasta") ~ "PA-M",
      disease == "CRC" ~ "CRC",
      str_detect(disease, "\\bCRC\\b") ~ "CRC+",
      str_detect(disease, "history") ~ "CRC-H",
      str_detect(disease, "\\bCRC\\b") & str_detect(disease, "metasta") ~ "CRC-M",
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    disease_class = factor(
      disease_class,
      levels = c(
        "Other", 
        "HC",
        "PA",
        "PA+",
        "PA-M",
        "CRC",
        "CRC+",
        "CRC-H",
        "CRC-M"
      )
    )
  ) %>% 
  mutate(
    age_decade = if_else(age_decade == "10-19" | age_decade == "20-29", "18-29", age_decade)
  ) %>% 
  select(
    study_name, subject_id, antibiotics_current_use,
    age, age_category, age_decade, gender, BMI, smoker, ever_smoker, alcohol, alcohol_numeric, diet, country, location, population,
    study_condition, disease, disease_class, disease_subtype, disease_stage, disease_location, days_from_first_collection, days_after_onset,
    sequencing_platform, DNA_extraction_kit, 
    cholesterol, wbc, rbc, stool_texture)
View(qualifying_studies)
cat("\n=== Qualifying Studies ===\n")
cat("Qualifying Studies samples (all):", nrow(qualifying_studies), "\n")
cat("Qualifying Studies samples (healthy):", qualifying_studies %>% filter(disease == "healthy") %>% nrow(), "\n")
cat("Qualifying Studies samples (disease):", qualifying_studies %>% filter(disease != "healthy") %>% nrow(), "\n")
cat("Qualifying Studies:", length(unique(qualifying_studies$study_name)), "\n")
cat("Study Sample Counts in Qualifying Studies:\n")
print(head(sort(table(qualifying_studies$study_name), decreasing = TRUE), 15))
cat("Top diseases in Qualifying Studies:\n")
print(head(sort(table(qualifying_studies$disease_class), decreasing = TRUE), 35))

