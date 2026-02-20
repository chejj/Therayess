##############################################################
# Testing with Single Small Study (ThomasAM_2018a, 80 Samples)
##############################################################
# Check Thomas2018a as it is small 
# w/ most categories (HC, PA, PA+, CRC, CRC+, Other)
# ----

# Run Qualifying Studies List
source("./EDA/Qualifying_Studies_List.R")

# Interested in:
################
focus <- c(
  "relative_abundance", # taxa
  "pathway_abundance",  # pathway abundance (metabolic)
  "pathway_coverage")   # pathway "presence" (metabolic)

#############################
# Exploring Taxa and Pathways
#############################
ThomasAM2018a_list <- list() 

for (type in focus) {   # Loop through what we want to focus on
  pattern <- paste0("ThomasAM_2018a.", type, "$") # What we are looking for in the "information available per study" generated above in .txt files
                                                  # "$" means "match this exactly"
  ThomasAM2018a_list[[type]] <- curatedMetagenomicData(  # create new list element by the "type" currently being looped through
    pattern,
    dryrun = FALSE
  )[[1]]
  
}

class(ThomasAM2018a_list)
length(ThomasAM2018a_list)
head(ThomasAM2018a_list)
