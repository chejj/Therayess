##########
# Dry Runs
##########
# Run Qualifying Studies List
source("./EDA/Qualifying_Studies_List.R")

# Studies available
###################
studies <- qualifying_studies %>%
  distinct(study_name) %>%
  pull(study_name)
print(studies)

# Information Available per Study
#################################
for (study in studies) {
  outfile <- paste0("./EDA/Kaitlyn_EDA/Info_Availability/", study, "_info_availability.txt")
  writeLines(curatedMetagenomicData(study, dryrun = TRUE), outfile)
}

# Option for Wet Run
####################
# Interested in:
focus <- c(
  "relative_abundance", # taxa
  "pathway_abundance",  # pathway abundance (metabolic)
  "pathway_coverage")   # pathway "presence" (metabolic)

# "For each data type I care about, build a search pattern that matches all my studies, download only that data type, and merge all studies into one combined dataset."
for (type in focus) {
  pattern <- paste0("(", paste(studies, collapse = "|"), ").*", type, "$") # Studies, separated by "OR" ("|"), then ".*" [type] and the "$" meaning END LINE
  x_list <- curatedMetagenomicData(pattern, dryrun = TRUE) # Change this to FALSE when time comes
  merged_by_focus[[type]] <- mergeData(x_list) # Merges thesummarizedTreeExperiments into one, by focus
}