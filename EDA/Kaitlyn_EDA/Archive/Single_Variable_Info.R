# Library and Setup
###################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(viridis)
})

source("EDA/Qualifying_Studies_List.R")
source("EDA/project_settings.R")
source("EDA/Kaitlyn_EDA/Reusable_Scripts/Utility_Functions.R")

############################
# Variables and Label Set Up
############################
df <- qualifying_studies

x_var <- trimws(readline(prompt = "Enter x (categorical) column name: "))

if (!x_var %in% names(df)) {
  stop("Variable not found.")
}

x_raw <- df[[x_var]]

#fill_var <- trim_or_null(readline(prompt = "Enter fill column name (optional; press enter to skip): "))

#~~~~~~~~~~~~~~~~~
# Validate Inputs:
#~~~~~~~~~~~~~~~~~
assert_col(df, x_var) # check if exists or if provided
if (!is.null(fill_var)) assert_col(df, fill_var) 

x_levels <- dplyr::n_distinct(df[[x_var]], na.rm = TRUE)
if (x_levels > 40) {
  message(sprintf("Note: '%s' has %d levels; bar plot may be hard to read.", x_var, x_levels))
} # warn if a lot of categories

#~~~~~~~~~~~~~~~~~~
# Filter NA Values:
#~~~~~~~~~~~~~~~~~~
plot_df <- df %>%
  filter(!is.na(.data[[x_var]])) %>% # remove rows where x is NA
  { if (is.null(fill_var)) . else filter(., !is.na(.data[[fill_var]])) } # if fill exists, filter NA from fill values

#################
# Quick Summaries
#################

summarize_cat(plot_df, x_var, top_n = 15) # uses function in "Utility functions" to summarize x variable

if (!is.null(fill_var)) { # summarizes by fill variable categories if it exists, by fill variable
  summarize_cat(plot_df, fill_var, top_n = 15)
  cat_vs_cat(plot_df, x_var, fill_var, max_cells = 400, run_test = FALSE) # Change to TRUE if you want chi square test to run
}

################################
# Count and reorder by frequency
################################

if (is.null(fill_var)) { # if there is no fill varibale
  
  count_df <- plot_df %>%
    count(.data[[x_var]], name = "n") %>% # new dataframe with the category and the frequency
    mutate(x_reordered = reorder(.data[[x_var]], n)) #reordered in descending order
  
  p <- ggplot(count_df, aes(x = x_reordered, y = n)) +
    geom_col() #tells ggplot to use the computed counts
  
} else { # if there is a fill variable
  
  count_df <- plot_df %>%
    count(.data[[x_var]], .data[[fill_var]], name = "n") %>% 
    group_by(.data[[x_var]]) %>% # groups to allow counting by category
    mutate(total = sum(n)) %>% # sum fill segments for bar height
    ungroup() %>%
    mutate(x_reordered = reorder(.data[[x_var]], total)) # reorder based on total count, not per fill level
  
  p <- ggplot(count_df, aes(x = x_reordered, y = n, fill = .data[[fill_var]])) + 
    geom_col()
}

###############
# Palette Logic
###############

if (!is.null(fill_var)) { # if the fill variable exists
  
  if (identical(fill_var, progression_col) && exists("progression_colors")) { # if the fill variable matches the progression
    
    # Warn if factor has levels without defined colors
    missing_cols <- setdiff(levels(plot_df[[fill_var]]), names(progression_colors))
    if (length(missing_cols) > 0) {
      warning("No colors defined for: ", paste(missing_cols, collapse = ", "))
    }
    
    p <- p + scale_fill_manual(
      values = progression_colors,
      name = progression_legend_title,
      drop = FALSE
    )
    
  } else {
    
    p <- p + scale_fill_viridis_d(    # Generic palette for arbitrary categorical variables
      option = "D",
      end = 0.95,
      name = fill_var
    )
  }
}

##############################
# Labels, Theme, & Print Plots
##############################
p <- p +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, lineheight = 1.1)) +
  labs(
    x = x_var,
    y = "Sample Count",
    title = if (is.null(fill_var)) {
      sprintf("Sample Count by %s", x_var)
    } else {
      sprintf("Sample Count by %s (filled by %s)", x_var, fill_var)
    }
  )

print(p) # Count Plot
