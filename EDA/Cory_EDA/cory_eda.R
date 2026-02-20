################################################################################
## ---- 1) Load packages
suppressPackageStartupMessages({
  library(BiocManager)
  library(curatedMetagenomicData)
  library(dplyr)
  library(stringr)
  library(phyloseq)
  library(mia)
  library(microViz)
  library(microbiome)
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
    study_name, sample_id, subject_id, antibiotics_current_use,
    age, age_category, age_decade, gender, BMI, smoker, ever_smoker, alcohol, alcohol_numeric, diet, country, location, population,
    study_condition, disease, disease_class, disease_subtype, disease_stage, disease_location, days_from_first_collection, days_after_onset,
    sequencing_platform, DNA_extraction_kit, 
    cholesterol, wbc, rbc, stool_texture)
#View(qualifying_studies)
cat("\n=== Qualifying Studies ===\n")
cat("Qualifying Studies samples (all):", nrow(qualifying_studies), "\n")
cat("Qualifying Studies samples (healthy):", qualifying_studies %>% filter(disease == "healthy") %>% nrow(), "\n")
cat("Qualifying Studies samples (disease):", qualifying_studies %>% filter(disease != "healthy") %>% nrow(), "\n")
cat("Qualifying Studies:", length(unique(qualifying_studies$study_name)), "\n")
cat("Study Sample Counts in Qualifying Studies:\n")
print(head(sort(table(qualifying_studies$study_name), decreasing = TRUE), 15))
cat("Top diseases in Qualifying Studies:\n")
print(head(sort(table(qualifying_studies$disease_class), decreasing = TRUE), 35))

# This is the command that curatedMetagenomicData understands to actually pull the relative abundance data
# Hence, dataType being "relative_abundance"
data <- returnSamples(
  qualifying_studies,
  dataType = "relative_abundance",
  counts = FALSE
)

# Extract components to put into Phyloseq object
# There are multiple microbiome
otu_mat <- assay(data, "relative_abundance")
sample_df <- as.data.frame(colData(data))
tax_df <- as.data.frame(rowData(data))

# 1. Fix the known bad row in tax_df BEFORE building phyloseq. This becomes an error when validating the Phyloseq object
# If we fix it here, it's not an issue later...
problem <- rownames(tax_df)[grepl("g__Enorma.*Collinsella_massiliensis", rownames(tax_df))]
tax_df  <- tax_df[!rownames(tax_df) %in% problem, ]
otu_mat <- otu_mat[!rownames(otu_mat) %in% problem, ]

# 2. Build the phyloseq object
gutBugs <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(sample_df),
  tax_table(as.matrix(tax_df))
)

# 3. Validate
gutBugs <- phyloseq_validate(gutBugs)

# 4. Now run tax_fix() — should be clean with no warnings
gutBugs <- gutBugs %>% tax_fix()

# 5. Verify results
tax_table(gutBugs) %>%
  as.data.frame() %>%
  dplyr::filter(is.na(family) | nchar(family) < 4)
# Should return an empty data frame if everything is resolved

# Create phyloseq object
gutBugs <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(sample_df),
  tax_table(as.matrix(tax_df))
)

# Verify
gutBugs
ntaxa(gutBugs)
nsamples(gutBugs)

# Check sample names 
sample_names(gutBugs) %>% head()

# Check taxa names
taxa_names %>% head()

# Check sample names
sample_variables(gutBugs)

# Retrieves tibble version of ps
samdat_tbl(gutBugs)

ps_validated <- phyloseq_validate(gutBugs)

tax_fix(gutBugs)

tax_fix_interactive(gutBugs)

ord_explore(gutBugs)

gutBugs %>% tax_fix()

# Creating bar plot

gutBugs <- gutBugs %>% phyloseq_validate() %>% tax_fix()

gutBugs %>%
  comp_barplot(
    "genus", n_taxa = 15,
    merge_other = FALSE,
    label = NULL,
    counts_warn = FALSE
  ) +
  facet_wrap(vars(study_condition), scales = "free") +
  coord_flip() +
  ggtitle(
    "Gut microbiota of different conditions",
    "Text"
  ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))




# Averages the mean across the condtition ggplot version
gutBugs %>%
  ps_filter(study_condition %in% c("control", "CRC"), .keep_all_taxa = TRUE) %>%
  tax_agg(rank = "genus") %>%
  tax_transform("compositional") %>%
  ps_melt() %>%
  dplyr::group_by(study_condition, genus) %>%
  dplyr::summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  dplyr::group_by(genus) %>%
  dplyr::mutate(mean_abund = mean(Abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genus = forcats::fct_lump_n(genus, n = 15, w = mean_abund)) %>%
  ggplot(aes(x = study_condition, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Gut microbiota: Control vs CRC", x = NULL, y = "Mean Relative Abundance") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# microViz version of bar plots
gutBugs %>%
  merge_samples(group = "study_condition") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  comp_barplot(
    "genus", n_taxa = 20,
    merge_other = FALSE,
    counts_warn = FALSE
  ) +
  coord_flip() +
  ggtitle("Gut microbiota by condition") +
  theme(axis.ticks.y = element_blank())

# Filter to CRC and control
gutBugs_crc_cont <- gutBugs %>%
  ps_filter(study_condition %in% c("control", "CRC"), .keep_all_taxa =  TRUE)

gutBugs_senior <- gutBugs_crc_cont %>%
  ps_filter(age_category %in% ("senior"), .keep_all_taxa = TRUE)

# Bar plot on filtered study condition
gutBugs_crc_cont %>%
  merge_samples(group = "study_condition") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  comp_barplot(
    "genus", n_taxa = 30,
    merge_other = FALSE,
    counts_warn = FALSE
  ) +
  coord_flip() +
  ggtitle("Gut microbiota of Seniors") +
  theme(axis.ticks.y = element_blank())

# Bar plot on filtered age_category
gutBugs_crc_cont %>%
  merge_samples(group = "age_category") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  comp_barplot(
    "genus", n_taxa = 30,
    merge_other = FALSE,
    counts_warn = FALSE
  ) +
  coord_flip() +
  ggtitle("Gut microbiota by age_category") +
  theme(axis.ticks.y = element_blank())

# Filtering down to just seniors
gutBugs_senior <- gutBugs_crc_cont %>%
  ps_filter(age_category %in% ("senior"), .keep_all_taxa = TRUE)

# Bar plot on filtered study condition
gutBugs_senior %>%
  merge_samples(group = "study_condition") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  comp_barplot(
    "genus", n_taxa = 30,
    merge_other = FALSE,
    counts_warn = FALSE
  ) +
  coord_flip() +
  ggtitle("Gut microbiota of Seniors by condition") +
  theme(axis.ticks.y = element_blank())

# Filtering down to just adults
gutBugs_adult <- gutBugs_crc_cont %>%
  ps_filter(age_category %in% ("adult"), .keep_all_taxa = TRUE)

# Bar plot on filtered study condition
gutBugs_adult %>%
  merge_samples(group = "study_condition") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  comp_barplot(
    "genus", n_taxa = 30,
    merge_other = FALSE,
    counts_warn = FALSE
  ) +
  coord_flip() +
  ggtitle("Gut microbiota of Adults by condition") +
  theme(axis.ticks.y = element_blank())


#####################################
# PCoA
#####################################
# PCoA BMI VS. study condition
gutBugs_crc_cont %>%
  tax_transform("compositional", rank = "genus") %>%
  dist_calc("bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    colour = "study_condition",
    size   = "BMI",
    alpha  = 0.5
  ) +
  stat_ellipse(aes(colour = study_condition), linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

# PCoA gender and study condition
gutBugs_crc_cont %>%
  tax_transform("compositional", rank = "genus") %>%
  dist_calc("bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    colour = "gender",
    shape   = "study_condition",
    size = "age",
    alpha  = 0.5
  ) +
  stat_ellipse(aes(colour = study_condition), linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

gutBugs_crc_cont %>%
  tax_transform("compositional", rank = "genus") %>%
  dist_calc("bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    colour = "gender",
    shape   = "study_condition",
    size = "age",
    alpha  = 0.5
  ) +
  stat_ellipse(aes(colour = study_condition), linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

gutBugs_senior %>%
  tax_transform("compositional", rank = "family") %>%
  dist_calc("bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    colour = "gender",
    shape   = "study_condition",
    size = "age",
    alpha  = 0.5
  ) +
  #stat_ellipse(aes(colour = study_condition), linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom")


# # View your current taxonomy to find the problem rows
# tax_table(ps) %>% as.data.frame() %>% View()
# 
# # Or filter to just the rows missing family
# tax_table(ps) %>%
#   as.data.frame() %>%
#   dplyr::filter(is.na(family) | family == "" | family == "f__")
# 
# tax_table(ps)["Gemmiger formicilis", "family"] <- "Oscillospiraceae"
# 
# # Search for any row containing "Gemmiger" (case-insensitive, partial match)
# tax_table(ps) %>%
#   as.data.frame() %>%
#   dplyr::filter(if_any(everything(), ~ grepl("Gemmiger", ., ignore.case = TRUE)))
# 
# # Find the exact rownames where genus is Gemmiger
# gemmiger_rows <- rownames(tax_table(ps))[tax_table(ps)[, "genus"] == "Gemmiger"]
# 
# gemmiger_rows
# 
# # Assign the family for all of them at once
# tax_table(ps)[gemmiger_rows, "family"] <- "Oscillospiraceae"


