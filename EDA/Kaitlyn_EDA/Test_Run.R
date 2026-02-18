##############################################################
# Testing with Single Small Study (ThomasAM_2018b, 59 Samples)
##############################################################
# Run Qualifying Studies List
source("./EDA/Qualifying_Studies_List.R")


Test <- curatedMetagenomicData("ThomasAM_2018b.+", 
                               dryrun = TRUE) # Doesn't download the data, dry runs it
class(Test)
length(Test)
head(Test)

# Interested in:
################
focus <- c(
  "relative_abundance", # taxa
  "pathway_abundance",  # pathway abundance (metabolic)
  "pathway_coverage")   # pathway "presence" (metabolic)


#############################
# Exploring Taxa and Pathways
#############################
ThomasAM2018b_list <- list()

for (type in focus) {   # Loop through what we want to focus on
  pattern <- paste0("ThomasAM_2018b.", type, "$") # What we are looking for in the "information available per study" generated above in .txt files
                                                  # "$" means "match this exactly"
  ThomasAM2018b_list[[type]] <- curatedMetagenomicData(  # create new list element by the "type" currently being looped through
    pattern,
    dryrun = FALSE
  )[[1]]
  
}

class(ThomasAM2018b_list)
length(ThomasAM2018b_list)
head(ThomasAM2018b_list)

