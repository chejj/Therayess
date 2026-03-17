################################################################################
# Functions / Helpers
#####################
library(curatedMetagenomicData)  # for mergeData
library(SummarizedExperiment)

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

merge_studies <- function(nested, type) {
  objs <- lapply(nested, `[[`, type) 
  # extracts the objects for one type, objs becomes a list with one element per study
  objs <- objs[!vapply(objs, is.null, logical(1))] 
  # removes missing entries, returns logical TRUE/FALSE vector, that's flipped to drop nulls
  curatedMetagenomicData::mergeData(objs)
}

eda_overview_table <- function(df) {
  
  if (!is.data.frame(df)) {
    stop("Input must be a data frame.", call. = FALSE)
  }
  
  results <- lapply(names(df), function(var) {
    
    x <- df[[var]]
    
    n_total   <- length(x)
    n_missing <- sum(is.na(x))
    pct_miss  <- round(100 * mean(is.na(x)), 2)
    n_unique  <- length(unique(x[!is.na(x)]))
    
    # Variable type detection
    var_type <- if (is.numeric(x) || is.integer(x)) {
      "Numeric"
    } else if (is.factor(x) && is.ordered(x)) {
      "Categorical (Ordered)"
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      "Categorical (Nominal)"
    } else {
      paste("Other (", paste(class(x), collapse = ", "), ")", sep = "")
    }
    
    data.frame(
      Column_Name     = var,
      Total_Count     = n_total,
      Missing_Count   = n_missing,
      Percent_Missing = pct_miss,
      Unique_Values   = n_unique,
      Variable_Type   = var_type,
      stringsAsFactors = FALSE
    )
  })
  
  overview_table <- do.call(rbind, results)
  
  return(overview_table)
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
