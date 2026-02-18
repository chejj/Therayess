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
  curatedMetagenomicData(study, dryrun = TRUE)
}
