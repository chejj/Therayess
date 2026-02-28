################################################################################
# Functions / Helpers
#####################

# General
#########
trim_or_null <- function(s) { 
  s <- trimws(s) # Removes leading or trailing spaces
  if (identical(s, "")) NULL else s # If there's nothing in the string, return NULL, makes fill variable optional
}

assert_col <- function(df, col) { # Checks if column is in the data frame
  if (!col %in% names(df)) {  
    stop(sprintf(
      "Column '%s' not found.\nAvailable columns (first 30): %s",
      col, paste(head(names(df), 30), collapse = ", ")
    ), call. = FALSE) # Prevents the error message function call stack from printing
  }
}

# Categorical
#############
summarize_cat <- function(df, var, top_n = 15) {
  v <- df[[var]] 
  
  overview <- tibble( # Small Summary Table
    variable = var,
    n = length(v),
    n_missing = sum(is.na(v)),
    pct_missing = round(100 * mean(is.na(v)), 2),
    n_distinct = dplyr::n_distinct(v, na.rm = TRUE)
  )
  print(overview)
  
  top_levels <- df %>%
    count(.data[[var]], sort = TRUE, name = "n") %>% # ".data[[var]]" is tidyverse data frame column name inside a string
    slice_head(n = top_n)
  
  print(top_levels)
  invisible(list(overview = overview, top_levels = top_levels)) # returns results without printing the object
}

cat_vs_cat <- function(df, x_var, y_var, max_cells = 400, run_test = FALSE) { # Contingency table + optional chi-square (with guardrails)
  nx <- dplyr::n_distinct(df[[x_var]], na.rm = TRUE)
  ny <- dplyr::n_distinct(df[[y_var]], na.rm = TRUE)
  
  if (nx * ny > max_cells) {
    message(sprintf(
      "Skipping contingency table: %d x %d = %d cells exceeds max_cells=%d",
      nx, ny, nx * ny, max_cells
    ))
    return(invisible(NULL)) # exits and returns nothing
  }
  
  ct <- table(df[[x_var]], df[[y_var]]) # Contingency table
  print(addmargins(ct)) #adds totals then prints
  
  if (run_test) { # must set run test to TRUE in the function for it to run
    tst <- suppressWarnings(chisq.test(ct))
    if (any(tst$expected < 5)) {
      message("Warning: Some expected counts < 5; chi-square assumptions may be violated.")
    }
    print(tst)
  }
  
  invisible(ct)
}

# Numerical
###########


# Mixed
#######

