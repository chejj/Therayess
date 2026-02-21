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
      levels = c("HC", "PA", "PA+", "CRC", "CRC+", "CRC-M", "CRC-H")
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


#creating a metadata heatmap showing some variables that aren't used downstream
meta_vars <- c(
  # Demographics
  "age", "age_decade", "gender", "BMI",
  
  # Lifestyle / exposures
  "diet", "alcohol", "alcohol_numeric", "smoker", "ever_smoker",
  "antibiotics_current_use",
  
  # Clinical context
  "disease_class", "disease_stage", "disease_location",
  
  # Geography
  "country", "location", "population",
  
  # Technical
  "sequencing_platform", "DNA_extraction_kit",
  
  # Timing
  "days_from_first_collection", "days_after_onset"
)

# Confirm they exist
meta_vars[meta_vars %in% colnames(qualifying_studies)]

#compute metadata completeness (percent of non-mising data)
library(tidyr)

meta_completeness <- qualifying_studies %>%
  summarise(across(all_of(meta_vars), ~ mean(!is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "completeness"
  ) %>%
  mutate(completeness = completeness * 100)

print(meta_completeness)

#create a heatmap for metadata completeness 
library(ggplot2)
meta_completeness <- meta_completeness %>%
  arrange(completeness) %>%
  mutate(variable = factor(variable, levels = variable))

ggplot(meta_completeness,
       aes(x = "", y = variable, fill = completeness)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 100),
    name = "% samples\nwith data"
  ) +
  labs(
    title = "Metadata Completeness in CRC Related Metagenomic Samples",
    subtitle = "Percent of samples with non-missing values per variable",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )








