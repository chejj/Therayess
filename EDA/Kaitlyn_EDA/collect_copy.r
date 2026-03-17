################################################################################
## CMD3 download + filter + save to HPC Therayess Data Folder
## - Works off of Wet_Runs.R file
## - Copy of Jake's EDA folder collect.r file with modifications to work on Sol supercomputer in conjunction with Wet_Runs.R

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(stringr)
})

## ---- SETTINGS ---------------------------------------------------------------
out_dir <- "/data/courses/class_lsc585spring2026_ndixit5/Therayess/data"
RDS_COMPRESS <- "xz"     # "xz" smallest, "gzip" faster, "bzip2" in-between
SAVE_ASSAY_MATRIX <- FALSE     # TRUE = save assay(x) too (bigger, but convenient)
SAVE_METADATA <- TRUE          # TRUE = save colData(x) too
OVERWRITE <- FALSE             # TRUE = overwrite existing files

## Output folders
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(out_dir, "manifest.csv")

## Helper: safe filename ----------------------------------------------
safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    str_replace_all("_+", "_")
}

stopifnot(exists("CRC_progression_studies"))
stopifnot(is.list(CRC_progression_studies))

## Save loop  ---------------------------
manifest <- list()

for (study in names(CRC_progression_studies)) {
  st <- CRC_progression_studies[[study]]
  if (!is.list(st)) next
  
  for (type in names(st)) {
    x <- st[[type]]
    if (is.null(x)) next
    
    base <- safe_name(paste0(study, "__", type))
    
    se_path   <- file.path(out_dir, paste0(base, ".rds"))
    assay_path <- file.path(out_dir, paste0("assay__", base, ".rds"))
    meta_path  <- file.path(out_dir, paste0("colData__", base, ".rds"))
    
    if (!OVERWRITE && file.exists(se_path)) {
      manifest[[length(manifest) + 1]] <- data.frame(
        study = study,
        product = type,
        status = "SKIPPED_EXISTS",
        n_features = NA_integer_,
        n_samples = NA_integer_,
        se_rds = se_path,
        assay_rds = if (SAVE_ASSAY_MATRIX) assay_path else NA_character_,
        colData_rds = if (SAVE_METADATA) meta_path else NA_character_,
        se_bytes = file.info(se_path)$size,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Basic integrity / dimensions
    n_features <- nrow(assay(x))
    n_samples  <- ncol(assay(x))
    
    # Save the full object (recommended)
    saveRDS(x, file = se_path, compress = RDS_COMPRESS)
    
    # Optional saves
    if (SAVE_ASSAY_MATRIX) {
      saveRDS(assay(x), file = assay_path, compress = RDS_COMPRESS)
    }
    if (SAVE_METADATA) {
      saveRDS(as.data.frame(colData(x)), file = meta_path, compress = RDS_COMPRESS)
    }
    
    manifest[[length(manifest) + 1]] <- data.frame(
      study = study,
      product = type,
      status = "OK",
      n_features = n_features,
      n_samples = n_samples,
      se_rds = se_path,
      assay_rds = if (SAVE_ASSAY_MATRIX) assay_path else NA_character_,
      colData_rds = if (SAVE_METADATA) meta_path else NA_character_,
      se_bytes = file.info(se_path)$size,
      stringsAsFactors = FALSE
    )
  }
}

manifest_df <- bind_rows(manifest)
write.csv(manifest_df, file = manifest_path, row.names = FALSE)

session_file <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_file)

zipfile <- file.path(out_dir, "cmd3_stool_bundle.zip")
rds_files <- list.files(file.path(out_dir, "rds"), full.names = TRUE)
if (length(rds_files) > 0) {
  try(utils::zip(zipfile, files = rds_files), silent = TRUE)
}

message("\nDone.\nSaved to: ", normalizePath(out_dir, winslash = "/"))
message("Manifest: ", normalizePath(manifest_path, winslash = "/"))
message("Session info: ", normalizePath(session_file, winslash = "/"))
if (file.exists(zipfile)) {
  message("Zip: ", normalizePath(zipfile, winslash = "/"))
} else if (length(rds_files) == 0) {
  message("Zip not created: no files were saved under ", normalizePath(file.path(out_dir, "rds"), winslash = "/"))
} else {
  message("Zip not created: zip command failed. RDS files are still available under ", normalizePath(file.path(out_dir, "rds"), winslash = "/"))
}