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

####################
# Checks for Success
####################
names(CRC_progression_studies)
names(CRC_progression_studies[["ZellerG_2014"]])
CRC_progression_studies[["ZellerG_2014"]][["relative_abundance"]]