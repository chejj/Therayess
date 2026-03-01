################################################################################
## ---- 1) Load packages
suppressPackageStartupMessages({
  library(BiocManager)
  library(curatedMetagenomicData)
  library(dplyr)
  library(stringr)
  library(knitr)
  library(kableExtra)
})
#View(sampleMetadata)
#str(sampleMetadata)
################################################################################
## ---- 2) Quick environment sanity check
#cat("\n=== Session Info ===\n")
#cat("R:", R.version.string, "\n")
#cat("Bioc:", as.character(BiocManager::version()), "\n")
#cat("Library paths:\n")
#print(.libPaths())
#cat("\ncuratedMetagenomicData:", as.character(packageVersion("curatedMetagenomicData")), "\n\n")

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
      age > 79  ~ "80+"
    )
  ) %>% 
  mutate(
    age_decade = factor(
      age_decade,
      levels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59",
                 "60-69", "70-79", "80+"),
      ordered = TRUE
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
        "CRC-M",
        "CRC-H"
      )
    )
  ) %>% 
  mutate(
    age_decade = if_else(age_decade == "10-19" | age_decade == "20-29", "18-29", age_decade)
  ) %>% 
  mutate(
    age_decade = factor(
      age_decade,
      levels = c("18-29","30-39","40-49","50-59","60-69","70-79","80+"),
      ordered = TRUE
    )
  )
#View(qualifying_studies)
#cat("\n=== Qualifying Studies ===\n")
#cat("Qualifying Studies samples (all):", nrow(qualifying_studies), "\n")
#cat("Qualifying Studies samples (healthy):", qualifying_studies %>% filter(disease == "healthy") %>% nrow(), "\n")
#cat("Qualifying Studies samples (disease):", qualifying_studies %>% filter(disease != "healthy") %>% nrow(), "\n")
#cat("Qualifying Studies:", length(unique(qualifying_studies$study_name)), "\n")
#cat("Study Sample Counts in Qualifying Studies:\n")
#print(head(sort(table(qualifying_studies$study_name), decreasing = TRUE), 15))
#cat("Top diseases in Qualifying Studies:\n")
#print(head(sort(table(qualifying_studies$disease_class), decreasing = TRUE), 35))

