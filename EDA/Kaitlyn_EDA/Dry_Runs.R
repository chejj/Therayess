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

##############################################################
# Testing with Single Small Study (ThomasAM_2018b, 59 Samples)
##############################################################
Test <- curatedMetagenomicData("ThomasAM_2018b.+", 
                               dryrun = TRUE) # Doesn't download the data, dry runs it
class(Test)
length(Test)
head(Test)

# Exploring Taxa Abundance
##########################
ThomasAM2018b_taxa_test <- curatedMetagenomicData("ThomasAM_2018b.+", 
                                                  dryrun = FALSE) # Actually downloads the data
