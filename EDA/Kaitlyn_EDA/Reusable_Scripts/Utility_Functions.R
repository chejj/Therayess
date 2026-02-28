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

eda_single_var <- function(df, var, top_n = 20) {
  
  # ---- Safety Checks ----
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.", call. = FALSE)
  }
  
  if (!is.character(var) || length(var) != 1) {
    stop("var must be a single string.", call. = FALSE)
  }
  
  if (!var %in% names(df)) {
    stop(sprintf("Variable '%s' not found in df.", var), call. = FALSE)
  }
  
  # ---- Pull Raw Vector (DO NOT FILTER) ----
  x <- df[[var]]
  
  # ---- Basic Counts ----
  n_total   <- length(x)
  n_missing <- sum(is.na(x))
  pct_miss  <- round(100 * mean(is.na(x)), 2)
  n_unique  <- length(unique(x[!is.na(x)]))
  
  # ---- Type Detection ----
  is_numeric  <- is.numeric(x) || is.integer(x)
  is_cat      <- is.factor(x) || is.character(x) || is.logical(x)
  is_ordered  <- is.factor(x) && is.ordered(x)
  
  # ---- Print Header ----
  cat("\n=====================================\n")
  cat("Variable:", var, "\n")
  cat("Class:", paste(class(x), collapse = ", "), "\n")
  cat("Total Observations:", n_total, "\n")
  cat("Missing:", n_missing, sprintf("(%.2f%%)", pct_miss), "\n")
  cat("Distinct (non-missing):", n_unique, "\n")
  
  if (is_numeric) {
    cat("Detected Type: Numeric\n")
  } else if (is_cat) {
    cat("Detected Type: Categorical\n")
    if (is_ordered) cat("Subtype: Ordered (Ordinal)\n")
  } else {
    cat("Detected Type: Other / Unclassified\n")
  }
  
  cat("=====================================\n\n")
  
  # ---- Numeric Summary ----
  if (is_numeric) {
    
    numeric_summary <- summary(x)
    
    print(numeric_summary)
    
    # Optional quick dispersion indicators
    cat("\nSD:", sd(x, na.rm = TRUE), "\n")
    cat("IQR:", IQR(x, na.rm = TRUE), "\n\n")
    
    return(invisible(list(
      var = var,
      type = "numeric",
      class = class(x),
      n_total = n_total,
      n_missing = n_missing,
      pct_missing = pct_miss,
      n_distinct = n_unique,
      summary = numeric_summary
    )))
  }
  
  # ---- Categorical Summary ----
  if (is_cat) {
    
    x_factor <- as.factor(x)
    
    if (is_ordered) {
      freq_table <- table(x_factor) # preserves factor level order
    } else {
      freq_table <- sort(table(x_factor), decreasing = TRUE)
    }
    
    prop_table <- round(prop.table(freq_table), 3)
    
    if (length(freq_table) > top_n) {
      cat(sprintf("Top %d Levels by Count:\n", top_n))
      print(head(freq_table, top_n))
      cat(sprintf("\n(Showing top %d of %d levels)\n\n", top_n, length(freq_table)))
    } else {
      print(freq_table)
      cat("\n")
    }
    
    cat("Proportions (Top Levels):\n")
    if (length(prop_table) > top_n) {
      print(head(prop_table, top_n))
    } else {
      print(prop_table)
    }
    
    cat("\n")
    
    return(invisible(list(
      var = var,
      type = "categorical",
      ordered = is_ordered,
      class = class(x),
      n_total = n_total,
      n_missing = n_missing,
      pct_missing = pct_miss,
      n_distinct = n_unique,
      freq = freq_table,
      prop = prop_table
    )))
  }
  
  # ---- Fallback ----
  return(invisible(list(
    var = var,
    type = "other",
    class = class(x),
    n_total = n_total,
    n_missing = n_missing,
    pct_missing = pct_miss,
    n_distinct = n_unique
  )))
}

# Categorical Variables
#######################
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

