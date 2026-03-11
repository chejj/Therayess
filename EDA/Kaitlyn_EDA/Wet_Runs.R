################################################################################
# Wet Runs → This can do dry runs as well, if `dryrun_mode = TRUE`
# --- Run Mode: Flip From TRUE to FALSE for the wet run ---
dryrun_mode <- FALSE
################################################################################

# --- Qualifying Studies List ---
source(here::here("EDA", "Qualifying_Studies_List.R"))

# --- Interested in: --- 
focus <- c("relative_abundance", "pathway_abundance", "pathway_coverage")

# --- Studies available --- 
studies <- qualifying_studies %>%
  distinct(study_name) %>%
  pull(study_name)
print(studies)

# --- Initialization --- 
CRC_progression_studies <- list()

# --- helpers --- # adds a backslash in front of any regex special character inside study names.
escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x) 
studies_rx <- escape_regex(studies)

################################################################################
# Per-study, per-type (not merged)
################################################################################
# For each study
# → For each data type (taxa, pathway abundance, pathway coverage)
# → Download exactly one matching dataset
# → Store it in a nested list

for (i in seq_along(studies)) { #seq_along() in R generates an integer sequence from 1 to the length of a given vector or list
  study <- studies[i] # current study as string
  study_rx <- studies_rx[i] # converts to regex safe version
  
  CRC_progression_studies[[study]] <- list() 
  
  for (type in focus) {
    pattern <- paste0(study_rx, "\\.", type, "$") # build strict regex with study name

    res <- curatedMetagenomicData(pattern, dryrun = dryrun_mode) # searches cMD3 for matching datasets
    
    if (!dryrun_mode) { # runs if downloading
      # enforce exactly one match before taking [[1]]
      if (length(res) != 1) {
        stop("Expected 1 dataset for pattern: ", pattern, " but got ", length(res))
      }
      CRC_progression_studies[[study]][[type]] <- res[[1]]
    } else {
      message("\n[DRYRUN] Pattern: ", pattern)
    }
  }
}

################################################################################
# --- Add Extra Metadata Layers & Filter: ---
# Do not run below if doing Dry Run
################################################################################
library(S4Vectors)
library(stringr)

filter_samples <- function(x) { #filters down to inclusion criteria
  
  cd <- as.data.frame(colData(x)) #convert to data frame
  
  required <- c("body_site", "age", "disease")
  if (!all(required %in% colnames(cd))) return(x)
  
  keep <- with(cd,
               body_site == "stool" & # stool samples only
                 age >= 18 & # adults only
                 !is.na(disease) #disease not blank
  )
  
  keep[is.na(keep)] <- FALSE
  
  x[, keep]
}

add_age_decade <- function(x) { # create a function to add age decade
  # If column already exists, do nothing
  if ("age_decade" %in% colnames(colData(x))) {
    return(x)
  }
  # If age column doesn't exist, do nothing
  if (!"age" %in% colnames(colData(x))) {
    return(x)
  }
  
  cd <- as.data.frame(colData(x)) # convert to a dataframe
  cd$age_decade <- ifelse( # create a new column and...
    is.na(cd$age), 
    NA_character_,
    paste0(floor(cd$age / 10) * 10, "-", floor(cd$age / 10) * 10 + 9) 
  ) # Example: 64 / 10 = 6.4 → floor(6.4) = 6 → 6 * 10 = 60, therefore paste0(60, "-", 69)
  colData(x) <- S4Vectors::DataFrame(cd) # convert back to SV4 Dataframe
  x #return object
}

add_disease_class <- function(x) {
  # If column already exists, do nothing
  if ("disease_class" %in% colnames(colData(x))) {
    return(x)
  }
  # If disease column doesn't exist, do nothing
  if (!"disease" %in% colnames(colData(x))) {
    return(x)
  }
  
  cd <- as.data.frame(colData(x))
  
  # normalize for matching
  d <- str_to_lower(as.character(cd$disease))
  
  cd$disease_class <- dplyr::case_when(
    is.na(d) ~ NA_character_,
    
    # Metastatic / history first (most specific)
    str_detect(d, "\\bcrc\\b") & str_detect(d, "metasta") ~ "CRC-M",
    str_detect(d, "history") ~ "CRC-H",
    str_detect(d, "\\b(adenoma|polyp)\\b") & str_detect(d, "metasta") ~ "PA-M",
    
    # Exact / common labels
    d == "healthy" ~ "HC",
    d %in% c("adenoma", "few_polyps") ~ "PA",
    d == "crc" ~ "CRC",
    
    # Broader catch-alls
    str_detect(d, "\\b(adenoma|polyp)\\b") ~ "PA+",
    str_detect(d, "\\bcrc\\b") ~ "CRC+",
    
    TRUE ~ "Other"
  )
  
  cd$disease_class <- factor(
    cd$disease_class,
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
  
  colData(x) <- S4Vectors::DataFrame(cd)
  x
}

# Filter inclusion Criteria
CRC_progression_studies <- lapply(CRC_progression_studies, function(st) lapply(st, filter_samples))
# Add Age_Decade 
CRC_progression_studies <- lapply(CRC_progression_studies, function(st) lapply(st, add_age_decade)) 
# Add Disease_Class 
CRC_progression_studies <- lapply(CRC_progression_studies, function(st) lapply(st, add_disease_class))

################################################################################
################################################################################
# Checks for Success-- Long version: cMD3 nested list (CRC_progression_studies)
# Purpose:
#   - Confirm every study has all data types
#   - Confirm each object has non-empty assays (features x samples)
#   - Confirm metadata aligns with assay columns (sample IDs match)
#   - Confirm our derived columns exist (age_decade, disease_class)
#   - Confirm inclusion criteria actually applied (stool, age >= 18, disease present)
################################################################################

library(SummarizedExperiment)

# Data types we expect in each study
focus <- c("relative_abundance", "pathway_abundance", "pathway_coverage")

# Metadata columns we expect to exist after our enrichment functions
required_meta <- c("study_name", "body_site", "age", "disease", "age_decade", "disease_class")

################################################################################
# 1) Coverage check: Do all studies have all types?
################################################################################

types_present <- lapply(CRC_progression_studies, names)

missing_types <- lapply(types_present, function(x) setdiff(focus, x))
missing_types <- missing_types[lengths(missing_types) > 0]

if (length(missing_types) == 0) {
  message("✅ All studies contain all expected types: ", paste(focus, collapse = ", "))
} else {
  message("❌ Some studies are missing expected types:")
  print(missing_types)
}

################################################################################
# 2) Dimension check: Are any objects empty or NA?
#    This catches: failed downloads, filtering leaving 0 samples, etc.
################################################################################

dim_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  st <- CRC_progression_studies[[study]]
  
  do.call(rbind, lapply(focus, function(tp) {
    x <- st[[tp]]
    if (is.null(x)) {
      return(data.frame(study = study, type = tp, n_features = NA, n_samples = NA))
    }
    d <- dim(assay(x))
    data.frame(study = study, type = tp, n_features = d[1], n_samples = d[2])
  }))
}))

# Show the report
print(dim_report)

# Flag issues (NA dims or 0 samples/features)
bad_dims <- subset(dim_report,
                   is.na(n_features) | is.na(n_samples) | n_features == 0 | n_samples == 0)

if (nrow(bad_dims) == 0) {
  message("✅ No empty/NA objects detected (features > 0 and samples > 0 for all).")
} else {
  message("❌ Found empty or missing objects:")
  print(bad_dims)
}

################################################################################
# 3) Alignment check: Do assay sample names match colData rownames?
#    This is critical. If FALSE, any metadata joins/plots may be wrong.
################################################################################

alignment_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  st <- CRC_progression_studies[[study]]
  
  do.call(rbind, lapply(focus, function(tp) {
    x <- st[[tp]]
    if (is.null(x)) {
      return(data.frame(study = study, type = tp, aligned = NA))
    }
    aligned <- identical(colnames(x), rownames(colData(x)))
    data.frame(study = study, type = tp, aligned = aligned)
  }))
}))

print(alignment_report)

bad_align <- subset(alignment_report, is.na(aligned) | aligned == FALSE)

if (nrow(bad_align) == 0) {
  message("✅ All objects have perfect alignment: colnames(x) == rownames(colData(x)).")
} else {
  message("❌ Alignment problems detected (metadata may not match samples):")
  print(bad_align)
}

################################################################################
# 4) Metadata presence check: Do we have our required columns everywhere?
#    This catches: functions not applied to some objects, columns missing in some studies.
################################################################################

meta_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  x <- CRC_progression_studies[[study]][["relative_abundance"]]  # pick one type for metadata checks
  if (is.null(x)) {
    return(data.frame(study = study, missing_cols = paste(required_meta, collapse = ", ")))
  }
  cols <- colnames(colData(x))
  missing <- setdiff(required_meta, cols)
  data.frame(study = study, missing_cols = if (length(missing) == 0) "" else paste(missing, collapse = ", "))
}))

print(meta_report)

bad_meta <- subset(meta_report, missing_cols != "")

if (nrow(bad_meta) == 0) {
  message("✅ All studies contain required metadata columns: ", paste(required_meta, collapse = ", "))
} else {
  message("❌ Some studies are missing required metadata columns:")
  print(bad_meta)
}

################################################################################
# 5) Inclusion criteria check: confirm filter was applied (spot-check per study)
#    We check on relative_abundance because filtering should have been applied to all types.
################################################################################

criteria_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  x <- CRC_progression_studies[[study]][["relative_abundance"]]
  if (is.null(x)) return(NULL)
  
  cd <- as.data.frame(colData(x))
  
  # logical checks for criteria
  stool_only <- all(tolower(cd$body_site) == "stool", na.rm = TRUE)
  adults_only <- all(cd$age >= 18, na.rm = TRUE)
  disease_present <- all(!is.na(cd$disease) & cd$disease != "", na.rm = TRUE)
  
  data.frame(
    study = study,
    n_samples = nrow(cd),
    stool_only = stool_only,
    adults_only = adults_only,
    disease_present = disease_present
  )
}))

print(criteria_report)

bad_criteria <- subset(criteria_report, stool_only == FALSE | adults_only == FALSE | disease_present == FALSE)

if (nrow(bad_criteria) == 0) {
  message("✅ Inclusion criteria appear to be enforced for all studies (stool, age>=18, disease present).")
} else {
  message("❌ Inclusion criteria violations detected (investigate these studies):")
  print(bad_criteria)
}

################################################################################
# 6) Optional: all-zero sample check (useful for pathway abundance especially)
#    This catches samples where the entire feature vector is 0.
################################################################################

zero_sample_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  x <- CRC_progression_studies[[study]][["pathway_abundance"]]
  if (is.null(x)) return(NULL)
  
  totals <- colSums(assay(x))
  n_zero <- sum(totals == 0)
  
  data.frame(study = study, n_zero_total_abundance = n_zero)
}))

print(zero_sample_report)
message("✅ Check complete.")

################################################################################
############################ --- CLEANUP --- ###################################
################################################################################
rm(
  missing_types,
  dim_report,
  bad_dims,
  alignment_report,
  bad_align,
  meta_report,
  bad_meta,
  criteria_report,
  bad_criteria,
  zero_sample_report
)

message("Temporary QC objects removed from environment.")