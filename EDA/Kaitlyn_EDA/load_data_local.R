################################################################################
## Rebuild CRC_progression_studies from manifest in course directory folder
## - Works off of Wet_Runs.R file & collect_copy.r
## - uses manifest to reproducibly reload the nested list originally made in Wet_Runs.R
QC_CHECK <- TRUE # TRUE = Full QC check, FALSE for no full QC check after load
################################################################################

# --- Library Loading ---
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# --- Manifest ---
out_dir <- "/data/courses/class_lsc585spring2026_ndixit5/Therayess/data"
manifest_path <- file.path(out_dir, "manifest.csv")

stopifnot(file.exists(manifest_path))

manifest_df <- read.csv(manifest_path, stringsAsFactors = FALSE)

# Only keep successfully saved objects
manifest_ok <- manifest_df %>%
  filter(status == "OK")

# --- Rebuild CRC_progression_studies nested list ---
CRC_progression_studies <- list()

for (i in seq_len(nrow(manifest_ok))) {
  
  study <- manifest_ok$study[i]
  type  <- manifest_ok$product[i]
  file  <- manifest_ok$se_rds[i]
  
  if (!file.exists(file)) {
    stop("Missing RDS file referenced in manifest: ", file)
  }
  
  if (!study %in% names(CRC_progression_studies)) {
    CRC_progression_studies[[study]] <- list()
  }
  
  CRC_progression_studies[[study]][[type]] <- readRDS(file)
}

message("Reload complete.")
length(CRC_progression_studies)
sapply(CRC_progression_studies, names)

################################################################################
################################################################################
# Checks for Success-- Post-load version: RDS -> cMD3 nested list (CRC_progression_studies)
# Purpose:
#   - Confirm every study has all data types
#   - Confirm each object has non-empty assays (features x samples)
#   - Confirm metadata aligns with assay columns (sample IDs match)
#   - Confirm our derived columns exist (age_decade, disease_class)
#   - Confirm inclusion criteria still hold (stool, age >= 18, disease present)
#   - Confirm we loaded full SE/TSE objects (not just assay matrices)
################################################################################
if (QC_CHECK) {
library(SummarizedExperiment)

# Data types we expect in each study
focus <- c("relative_abundance", "pathway_abundance", "pathway_coverage")

# Metadata columns we expect to exist after our enrichment functions
required_meta <- c("study_name", "body_site", "age", "disease", "age_decade", "disease_class")

################################################################################
# 0) Quick sanity check: does CRC_progression_studies exist and look right?
################################################################################

stopifnot(exists("CRC_progression_studies"))
stopifnot(is.list(CRC_progression_studies))
message("✅ CRC_progression_studies exists and is a list. N studies: ", length(CRC_progression_studies))

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
# 2) Class check: Did we load full SE/TSE objects?
#    This catches accidentally loading assay-only RDS files.
################################################################################

class_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  st <- CRC_progression_studies[[study]]
  do.call(rbind, lapply(focus, function(tp) {
    x <- st[[tp]]
    if (is.null(x)) {
      return(data.frame(study = study, type = tp, class = NA_character_))
    }
    data.frame(study = study, type = tp, class = paste(class(x), collapse = ";"))
  }))
}))

print(class_report)

bad_class <- subset(class_report, is.na(class) | !grepl("SummarizedExperiment", class))
if (nrow(bad_class) == 0) {
  message("✅ All loaded objects inherit from SummarizedExperiment (full objects loaded).")
} else {
  message("❌ Some objects do NOT look like SummarizedExperiment/TSE (check you loaded the correct RDS files):")
  print(bad_class)
}

################################################################################
# 3) Dimension check: Are any objects empty or NA?
#    This catches: missing files, failed loads, filtering leaving 0 samples, etc.
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

print(dim_report)

bad_dims <- subset(dim_report,
                   is.na(n_features) | is.na(n_samples) | n_features == 0 | n_samples == 0)

if (nrow(bad_dims) == 0) {
  message("✅ No empty/NA objects detected (features > 0 and samples > 0 for all).")
} else {
  message("❌ Found empty or missing objects:")
  print(bad_dims)
}

################################################################################
# 4) Alignment check: Do assay sample names match colData rownames?
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
# 5) Metadata presence check: Do we have our required columns everywhere?
#    We check on relative_abundance because metadata should be consistent across types.
################################################################################

meta_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  x <- CRC_progression_studies[[study]][["relative_abundance"]]
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
# 6) Inclusion criteria check: confirm filter still holds (spot-check per study)
#    We check on relative_abundance because filtering should have been applied to all types.
################################################################################

criteria_report <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
  x <- CRC_progression_studies[[study]][["relative_abundance"]]
  if (is.null(x)) return(NULL)
  
  cd <- as.data.frame(colData(x))
  
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
# 7) Optional: all-zero sample check (useful for pathway abundance especially)
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
# 8) Manifest comparison check:
#    Confirm that all successfully saved RDS files were loaded into memory.
################################################################################
out_dir <- "/data/courses/class_lsc585spring2026_ndixit5/Therayess/data"
manifest_path <- file.path(out_dir, "manifest.csv")

if (file.exists(manifest_path)) {
  
  manifest_df <- read.csv(manifest_path, stringsAsFactors = FALSE)
  
  # Only consider successfully saved objects
  manifest_ok <- subset(manifest_df, status == "OK")
  
  # Build in-memory index
  memory_index <- do.call(rbind, lapply(names(CRC_progression_studies), function(study) {
    st <- CRC_progression_studies[[study]]
    data.frame(
      study = study,
      product = names(st),
      stringsAsFactors = FALSE
    )
  }))
  
  # Expected combinations from manifest
  manifest_index <- manifest_ok[, c("study", "product")]
  
  # Expected combinations from manifest
  manifest_index <- manifest_ok[, c("study", "product")]
  manifest_index <- dplyr::distinct(manifest_index)
  
  # Find manifest entries that are missing in memory
  missing_in_memory <- dplyr::anti_join(
    manifest_index,
    memory_index,
    by = c("study", "product")
  )
  
  if (nrow(missing_in_memory) == 0) {
    message("✅ All manifest entries (status == OK) are loaded into CRC_progression_studies.")
  } else {
    message("❌ Some manifest entries were NOT loaded into memory:")
    print(missing_in_memory)
  }
  
} else {
  message("⚠️ Manifest file not found — skipping manifest comparison check.")
}

################################################################################
############################ --- CLEANUP --- ###################################
################################################################################

rm(
  types_present,
  missing_types,
  class_report,
  bad_class,
  dim_report,
  bad_dims,
  alignment_report,
  bad_align,
  meta_report,
  bad_meta,
  criteria_report,
  bad_criteria,
  zero_sample_report, 
  missing_in_memory, 
  manifest_index, 
  memory_index,
  manifest_ok, 
  out_dir, 
  file, 
  QC_CHECK, 
  required_meta, 
  study, 
  type, 
  manifest_path,
  i, 
  focus
)

message("Temporary QC objects removed from environment.")

} else {
  message("QC Check skipped.")
}