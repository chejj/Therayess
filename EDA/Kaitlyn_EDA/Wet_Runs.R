##########
# Wet Runs
##########
# Run Qualifying Studies List
source("./EDA/Qualifying_Studies_List.R")

# Interested in:
################
focus <- c(
  "relative_abundance", # taxa
  "pathway_abundance",  # pathway abundance (metabolic)
  "pathway_coverage")   # pathway "presence" (metabolic)

# Studies available
###################
studies <- qualifying_studies %>%
  distinct(study_name) %>%
  pull(study_name)
print(studies)


# Single Study
##############
for (type in focus) {
  pattern <- paste0("(", paste(studies, collapse = "|"), ").*", type, "$") # Studies, separated by "OR" ("|"), then ".*" [type] and the "$" meaning END LINE
  x_list <- curatedMetagenomicData(pattern, dryrun = TRUE) # Change this to FALSE when time comes
  merged_by_focus[[type]] <- mergeData(x_list) # Merges the summarizedTreeExperiments into one, by focus
}

# All Studies
#############
# "For each data type I care about, 
# build a search pattern that matches all my studies, 
# download only that data type, 
# and merge all studies into one combined dataset."
CRC_progression_studies <- list()

for (study in studies) {
  for (type in focus) {
    
    pattern <- paste0(study, ".", type, "$") # Studies, separated by "OR" ("|"), then ".*" [type] and the "$" meaning END LINE
    
    CRC_progression_studies[[study]][[type]] <-
      curatedMetagenomicData(pattern, dryrun = TRUE)[[1]] # Change this to FALSE when time comes
  }
}