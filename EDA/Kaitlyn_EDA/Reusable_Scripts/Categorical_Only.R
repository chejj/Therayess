# Library and Setup
#########################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

source("EDA/Qualifying_Studies_List.R")
source("EDA/Kaitlyn_EDA/Reusable_Scripts/Utility_Functions.R")
################################################################################
# Variables and Label Set Up
##########################
df <- qualifying_studies

x_var <- trimws(readline(prompt = "Enter x (categorical) column name: "))
fill_var <- trim_or_null(readline(prompt = "Enter fill column name (optional; press enter to skip): "))

# Validate Inputs:
#~~~~~~~~~~~~~~~~~
assert_col(df, x_var) # check if exists or if provided
if (!is.null(fill_var)) assert_col(df, fill_var) 

x_levels <- dplyr::n_distinct(df[[x_var]], na.rm = TRUE)
if (x_levels > 40) {
  message(sprintf("Note: '%s' has %d levels; bar plot may be hard to read.", x_var, x_levels))
} # warn if a lot of categories



########################################

# Plotting
qualifying_studies %>%
  # Plotting
  # --------
  ggplot(aes(x = reorder(x, x, 
                       FUN = function(x) -length(x)),
           fill = y)) +
  geom_bar() +
  scale_fill_manual(values = progression_colors, 
                    name = "CRC Progression \nClassifiers") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, lineheight = 1.1)) +
  # Labeling
  #---------
  xlab("Author Last Name, First Initial, and Year") +
  ylab("Sample Count") +
  ggtitle("Qualifying Study Sample Count by CRC Progression Class")