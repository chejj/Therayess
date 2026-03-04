##########
# Wet Runs → This can do dry runs as well, if `dryrun_mode = TRUE`
# --- Run Mode: Flip From TRUE to FALSE for the wet run ---
dryrun_mode <- FALSE
##########

# --- Qualifying Studies List ---
source("./EDA/Qualifying_Studies_List.R")



# --- Interested in: --- 
focus <- c("relative_abundance", "pathway_abundance", "pathway_coverage")

# --- Studies available --- 
studies <- qualifying_studies %>%
  distinct(study_name) %>%
  pull(study_name)
print(studies)

# --- Initialization --- 
merged_by_focus <- list()
CRC_progression_studies <- list()

# --- helpers --- # adds a backslash in front of any regex special character inside study names.
escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x) 
studies_rx <- escape_regex(studies)

########################
# By Focus (All Studies) → Unlikely to use since we are doing LODO
########################
for (type in focus) {
  pattern <- paste0("(", paste(studies_rx, collapse = "|"), ").*\\.", type, "$") 
  # Studies, separated by "OR" ("|"), then ".*" [type] and the "$" meaning END LINE
  
  x_list <- curatedMetagenomicData(pattern, dryrun = dryrun_mode)
  
  if (!dryrun_mode) { # Merges the summarizedTreeExperiments into one, by focus
    merged_by_focus[[type]] <- mergeData(x_list)
  } else {
    message("\n[DRYRUN] Pattern: ", pattern)
  } 
}

################################
# Per-study, per-type (no merge)
################################
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

####################################
# --- Add Extra Metadata Layers: ---
# Do not run below if doing Dry Run
####################################
library(S4Vectors)
library(stringr)

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

CRC_progression_studies <- lapply(CRC_progression_studies, function(st) {
  st <- lapply(st, add_age_decade)
  st
}) # Apply Age_Decade to entire Tree Summarized Experiment

CRC_progression_studies <- lapply(CRC_progression_studies, function(st) {
  lapply(st, add_disease_class)
}) # Apply Disease_Class to entire Tree Summarized Experiment


################################################################################
# --- Is everything that we need loaded in? ---
focus <- c("relative_abundance", "pathway_abundance", "pathway_coverage")

# which studies loaded?
loaded_studies <- names(CRC_progression_studies)

# which types exist per study?
types_present <- lapply(CRC_progression_studies, names)

# studies missing any type?
missing_types <- lapply(types_present, function(x) setdiff(focus, x))
missing_types <- missing_types[lengths(missing_types) > 0]
missing_types # prints nothing / empty list → we have all types for all studies.

library(SummarizedExperiment)

dims <- lapply(CRC_progression_studies, function(st) {
  sapply(focus, function(tp) {
    x <- st[[tp]]
    if (is.null(x)) return(c(NA, NA))
    dim(assay(x))
  })
})
dims[1:12] # Peek and look for no NA, no 0 rows/cols.

x <- CRC_progression_studies[[1]][["relative_abundance"]]
stopifnot(identical(colnames(x), rownames(colData(x)))) # pass means metadata can be added by sample ID

####################
# Checks for Success
####################
names(CRC_progression_studies)
names(CRC_progression_studies[["ZellerG_2014"]])
CRC_progression_studies[["ZellerG_2014"]][["relative_abundance"]]