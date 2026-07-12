# Load all needed ##############################################################
if (!exists("CRC_progression_studies")) {
  stop("CRC Progression Studies not loaded. If load_data.R or load_data_local.R is not functional, run the following: Qualifying_Studies_List.R & Wet_Runs.R")
}
## Libraries --------
suppressPackageStartupMessages({
  library(curatedMetagenomicData)  # for mergeData
  library(SummarizedExperiment)
  library(vegan)
  library(tidyverse)
  library(scater)
  library(sva)
  library(MMUPHin)
#  BiocManager::install("sva") # COMBAT
#  BiocManager::install("MMUPHin")
#  BiocManager::install("curatedMetagenomicData")
})

## Merge Studies and create objects --------
types <- c("relative_abundance", "pathway_coverage", "pathway_abundance")

merged_studies <- lapply(types, function(type) {
  objs <- lapply(CRC_progression_studies, `[[`, type)
  objs <- objs[!vapply(objs, is.null, logical(1))]
  mergeData(objs)
})
# Name the list
names(merged_studies) <- types

# Samples are rows, features are columns (samples x features)
studies_taxa <- as.data.frame(t(assay(merged_studies[["relative_abundance"]])))
studies_path_abund <- as.data.frame(t(assay(merged_studies[["pathway_abundance"]])))
studies_path_cov <- as.data.frame(t(assay(merged_studies[["pathway_coverage"]])))
studies_meta_df <- as.data.frame(colData(merged_studies[["relative_abundance"]]))

### Metadata Cleanup ----
# Transform metadata to have disease class and correct age classifiers 
studies_meta_df <- studies_meta_df %>%
  mutate(
    disease_class = case_when(
      disease == "healthy" ~ "HC",
      str_detect(disease, "\\b(adenoma|polyp)\\b") & str_detect(disease, "metasta") ~ "PA-M",
      str_detect(disease, "\\bCRC\\b") & str_detect(disease, "metasta") ~ "CRC-M",
      str_detect(disease, "history") ~ "CRC-H",
      disease == "CRC" ~ "CRC",
      str_detect(disease, "\\bCRC\\b") ~ "CRC+",
      disease %in% c("adenoma", "few_polyps") ~ "PA",
      str_detect(disease, "\\b(adenoma|polyp)\\b") ~ "PA+",
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    disease_class = factor(
      disease_class,
      levels = c("Other", "HC", "PA", "PA+", "PA-M", "CRC", "CRC+", "CRC-M", "CRC-H")
    ), 
    age_decade = case_when(
      age < 10 ~ "0-9", 
      age < 20 ~ "10-19",
      age < 30 ~ "20-29",
      age < 40 ~ "30-39",
      age < 50 ~ "40-49",
      age < 60 ~ "50-59",
      age < 70 ~ "60-69",
      age < 80 ~ "70-79",
      TRUE  ~ "80+"
    ),
    age_decade = if_else(age_decade %in% c("10-19", "20-29"), "18-29", age_decade),
    age_decade = factor(
      age_decade,
      levels = c("18-29","30-39","40-49","50-59","60-69","70-79","80+"),
      ordered = TRUE
    )
  )

# Numeric matrix version lists #################################################
# for control (no batch effect correction), log, % norm, COMBAT, MMUPHin
## Taxa ----
taxa <- as.matrix(studies_taxa)
storage.mode(taxa) <- "numeric"
taxa_versions <- list(
  control = taxa,
  log = NULL,
  percentile = NULL,
  combat = NULL,
  mmuphin = NULL
)

## Pathway Abundance ----
pab <- as.matrix(studies_path_abund)
storage.mode(pab) <- "numeric"
pab_versions <- list(
  control = pab,
  log = NULL,
  percentile = NULL,
  combat = NULL,
  mmuphin = NULL
)

## Pathway Coverage ----
pcov <- as.matrix(studies_path_cov)
storage.mode(pcov) <- "numeric"
pcov_versions <- list(
  control = pcov,
  log = NULL,
  percentile = NULL,
  combat = NULL,
  mmuphin = NULL
)

# Complete all transformations for each batch effect adjustment type

# TRANSFORMATIONS ##############################################################
### Setup ----
batch <- studies_meta_df$study_name # expected batch effect
group <- studies_meta_df$disease_class # groups we want to be able to detect
control_label <- "HC" # baseline group or control

## A: Log-transformation ----

### TAXA ----
# summary(rowSums(taxa)) # if values are all ~100, next lines can be skipped
## log_rel <- sweep(log_taxa, 2, colSums(log_taxa), "/") # Not required, as cMD3 already does this step (above line as proof)
taxa_versions$log <- log10(taxa_versions$control + 1e-6)
message("Taxa log transformation complete")

### PAB ----
summary(rowSums(pab)) # if values are all ~100, next line can be skipped
pab_versions$log <- pab_versions$control
log_rel <- sweep(pab_versions$log, 2, colSums(pab_versions$log), "/")
pab_versions$log <- log10(log_rel + 1e-6)
message("Pathway Abundance log transformation complete")

### PCOV ----
summary(rowSums(pcov)) # if values are all ~100, next line can be skipped
pcov_versions$log <- pcov_versions$control
log_rel <- sweep(pcov_versions$log, 2, colSums(pcov_versions$log), "/")
pcov_versions$log <- log10(log_rel + 1e-6)
message("Pathway Coverage log transformation complete")

## B: Percentile Normalization ----
percentile_norm <- function(feature_values, group, control_label = "HC") {
  controls <- feature_values[group == control_label]
  if (length(controls) == 0) {
    return(rep(NA_real_, length(feature_values)))
  }
  vapply(feature_values, function(x) {
    mean(controls <= x) * 100
  }, numeric(1))
}

percentile_normalize_by_study <- function(control_mat,
                                          batch,
                                          group,
                                          control_label = "HC",
                                          feature_name = NULL,
                                          verbose = TRUE) {
  mat <- as.matrix(control_mat)
  storage.mode(mat) <- "numeric"
  
  batch <- as.character(batch)
  group <- as.character(group)
  
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  missing_studies <- character(0)
  
  for (b in unique(batch)) {
    idx <- which(batch == b)
    sub_data <- mat[idx, , drop = FALSE]
    sub_group <- group[idx]
    
    if (!any(sub_group == control_label, na.rm = TRUE)) {
      missing_studies <- c(missing_studies, b)
      next
    }
    
    sub_percentile <- matrix(
      NA_real_,
      nrow = nrow(sub_data),
      ncol = ncol(sub_data),
      dimnames = dimnames(sub_data)
    )
    
    for (j in seq_len(ncol(sub_data))) {
      sub_percentile[, j] <- percentile_norm(
        feature_values = sub_data[, j],
        group = sub_group,
        control_label = control_label
      )
    }
    
    out[idx, ] <- sub_percentile
  }
  
  missing_studies <- unique(missing_studies)
  
  if (verbose) {
    message("[", feature_name, "] Percentile normalization complete.")
    
    if (length(missing_studies) > 0) {
      message(
        "[", feature_name, "] Percentile normalization could not be performed for: ",
        paste(missing_studies, collapse = ", ")
      )
    }
  }
  
  list(
    percentile = out,
    missing_studies = missing_studies
  )
}

### TAXA ----

taxa_pct <- percentile_normalize_by_study(
  control_mat = taxa_versions$control,
  batch = batch,
  group = group,
  control_label = control_label,
  feature_name = "Taxa"
)

taxa_versions$percentile <- taxa_pct$percentile

### PAB ----
pab_pct <- percentile_normalize_by_study(
  control_mat = pab_versions$control,
  batch = batch,
  group = group,
  control_label = control_label,
  method_name = "Pathway Abundance"
)

pab_versions$percentile <- pab_pct$percentile

### PCOV ----
pcov_pct <- percentile_normalize_by_study(
  control_mat = pcov_versions$control,
  batch = batch,
  group = group,
  control_label = control_label,
  method_name = "Pathway Coverage"
)

pcov_versions$percentile <- pcov_pct$percentile

# C: COMBAT ----
run_combat <- function(control_mat,
                       metadata,
                       batch,
                       method_name = NULL) {
  
  # ComBat expects features x samples, so transpose in and out (currently in samples x features)
  combat_input <- t(control_mat)
  
  # preserve disease_class so batch correction does not remove biology
  mod <- model.matrix(
    ~ disease_class,
    data = metadata
  )
  
  combat_corrected <- ComBat(
    dat = combat_input,
    batch = factor(batch),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  message("[", method_name, "] COMBAT correction complete.")
  
  t(combat_corrected)
}

### TAXA ----
taxa_versions$combat <- run_combat(
  control_mat = taxa_versions$control,
  metadata = studies_meta_df,
  batch = batch,
  method_name = "Taxa"
)

### PAB ----
pab_versions$combat <- run_combat(
  control_mat = pab_versions$control,
  metadata = studies_meta_df,
  batch = batch,
  method_name = "Pathway Abundance"
)


### PCOV ----
pcov_versions$combat <- run_combat(
  control_mat = pcov_versions$control,
  metadata = studies_meta_df,
  batch = batch,
  method_name = "Pathway Coverage"
)

# D: MMUPHin ----
run_mmuphin <- function(control_mat,
                        metadata,
                        method_name = NULL) {
  
  # MMUPHin expects features x samples
  mmuphin_input <- t(control_mat)
  
  # Convert percentages to proportions if needed
  max_val <- max(mmuphin_input, na.rm = TRUE)
  scale_factor <- 1
  
  if (max_val > 1 && max_val <= 100) {
    mmuphin_input <- mmuphin_input / 100
    scale_factor <- 100
  }
  
  mmuphin_fit <- adjust_batch(
    feature_abd = mmuphin_input,
    batch = "study_name",
    covariates = "disease_class",
    data = metadata,
    control = list(verbose = FALSE)
  )
  
  message("[", method_name, "] MMUPHin correction complete.")
  
  t(mmuphin_fit$feature_abd_adj) * scale_factor
}

### TAXA ----
taxa_versions$mmuphin <- run_mmuphin(
  control_mat = taxa_versions$control,
  metadata = studies_meta_df,
  method_name = "Taxa"
)


### PAB ----
pab_versions$mmuphin <- run_mmuphin(
  control_mat = pab_versions$control,
  metadata = studies_meta_df,
  method_name = "Pathway Abundance"
)

### PCOV ----
pcov_versions$mmuphin <- run_mmuphin(
  control_mat = pcov_versions$control,
  metadata = studies_meta_df,
  method_name = "Pathway Coverage"
)

# PCoA #########################################################################
dist_methods <- c(
  control = "bray",
  percentile = "bray",
  combat = "bray",
  mmuphin = "bray",
  log = "euclidean"
)

run_pcoa <- function(feature_matrix,
                     metadata,
                     distance = "bray", 
                     method_name = NULL) {
  
  # Ensure numeric matrix
  mat <- as.matrix(feature_matrix)
  storage.mode(mat) <- "numeric"
  
  ## Report missing values
  na_count <- sum(is.na(mat))
  if (na_count > 0) {
    message("[", method_name, "] Replacing ", na_count," NA values with 0.")
    mat[is.na(mat)] <- 0
  }
  
  ## Report empty samples
  empty_rows <- rowSums(mat) == 0
  n_empty <- sum(empty_rows)
  
  if (n_empty > 0) {
    message("[", method_name, "] Removing ",n_empty," empty samples before PCoA.")
    mat <- mat[!empty_rows, , drop = FALSE]
    metadata <- metadata[!empty_rows, , drop = FALSE]
  }
  
  # Compute Bray-Curtis distances between samples
  dist_mat <- vegdist(mat, method = distance)
  
  # Run PCoA
  pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  # Combine coordinates with metadata
  pcoa_df <- cbind(
    metadata,
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2]
  )
  
  if (!is.null(method_name)) {
    pcoa_df$method <- method_name
  }
  
  # Percent variance explained
  eig_vals <- pcoa$eig
  var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
  
  message(
    "[", method_name, "] PCoA complete (",
    nrow(mat), " samples, ",
    ncol(mat), " features)."
  )
  
  list(
    distance = dist_mat,
    pcoa = pcoa,
    pcoa_df = pcoa_df,
    variance_explained = var_explained
  )
}


taxa_pcoa_results <- list()
pab_pcoa_results <- list()
pcov_pcoa_results <- list()

for (method in names(taxa_versions)) {
  taxa_pcoa_results[[method]] <-
    run_pcoa(
      feature_matrix = taxa_versions[[method]],
      metadata = studies_meta_df,
      distance = dist_methods[[method]],
      method_name = method
    )
}

for (method in names(pab_versions)) {
  pab_pcoa_results[[method]] <-
    run_pcoa(
      feature_matrix = pab_versions[[method]],
      metadata = studies_meta_df,
      distance = dist_methods[[method]],
      method_name = method
    )
}

for (method in names(pcov_versions)) {
  pcov_pcoa_results[[method]] <-
    run_pcoa(
      feature_matrix = pcov_versions[[method]],
      metadata = studies_meta_df,
      distance = dist_methods[[method]],
      method_name = method
    )
}

## Joint Dataframe ----
taxa_pcoa_df <- purrr::imap_dfr(taxa_pcoa_results, function(res, method) {
  res$pcoa_df %>%
    mutate(
      method = method,
      PC1_lab = paste0("PCoA1 (", round(res$variance_explained[1] * 100, 1), "%)"),
      PC2_lab = paste0("PCoA2 (", round(res$variance_explained[2] * 100, 1), "%)")
    )
})

pab_pcoa_df <- purrr::imap_dfr(pab_pcoa_results, function(res, method) {
  res$pcoa_df %>%
    mutate(
      method = method,
      PC1_lab = paste0("PCoA1 (", round(res$variance_explained[1] * 100, 1), "%)"),
      PC2_lab = paste0("PCoA2 (", round(res$variance_explained[2] * 100, 1), "%)")
    )
}) 

pcov_pcoa_df <- purrr::imap_dfr(pcov_pcoa_results, function(res, method) {
  res$pcoa_df %>%
    mutate(
      method = method,
      PC1_lab = paste0("PCoA1 (", round(res$variance_explained[1] * 100, 1), "%)"),
      PC2_lab = paste0("PCoA2 (", round(res$variance_explained[2] * 100, 1), "%)")
    )
})

# Plotting #####################################################################
progression_colors <- c(
  "HC"    = "#539deb",
  "PA"    = "#e4e85b",
  "PA+"   = "#e8a25b",
  "PA-M"  = "#c97c2b",  
  "CRC"   = "#cf1919",
  "CRC+"  = "#8f0000",
  "CRC-M" = "#4d0404",
  "CRC-H" = "#520f76",
  "Other" = "#aaa3a3"
)

progression_legend_title <- "CRC Progression\nClassifiers"

## Individual -----
taxa_pcoa_results$combat$plot_df %>%
  ggplot(aes(PC1, PC2, color = disease_class)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = progression_colors, name = progression_legend_title) +
  theme_bw() +
  labs(
    title = "PCoA (ComBat), Colored by Disease",
    x = paste0("PCoA1 (", round(pcoa_results$combat$variance_explained[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_results$combat$variance_explained[2] * 100, 1), "%)")
  )

taxa_pcoa_results$combat$plot_df %>%
  ggplot(aes(PC1, PC2, color = study_name)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) +
  theme_bw() +
  labs(
    title = "PCoA (ComBat), Colored by Study",
    x = paste0("PCoA1 (", round(pcoa_results$combat$variance_explained[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_results$combat$variance_explained[2] * 100, 1), "%)")
  )

## Faceted Together ----
plot_pcoa <- function(pcoa_df,
                      color_by = c("disease_class", "study_name"),
                      dataset_name = "Relative Abundance",
                      disease_filter = c("HC", "PA", "PA+", "CRC", "CRC+", "Other"),
                      show_ellipses = TRUE) {
  
  color_by <- match.arg(color_by)
  
  plot_df <- pcoa_df %>%
    filter(disease_class %in% disease_filter)
  
  p <- ggplot(
    plot_df,
    aes(PC1, PC2, color = .data[[color_by]])
  ) +
    geom_point(alpha = 0.6)
  
  if (show_ellipses) {
    p <- p +
      stat_ellipse(level = 0.95)
  }
  
  p <- p +
    facet_wrap(~method, scales = "free") +
    theme_bw()
  
  if (color_by == "disease_class") {
    
    p <- p +
      scale_color_manual(
        values = progression_colors,
        name = progression_legend_title
      ) +
      labs(
        title = paste0(
          "PCoA of ",
          dataset_name,
          ", Colored by CRC Progression Classifier"
        ),
        x = "PCoA1",
        y = "PCoA2"
      )
    
  } else {
    
    p <- p +
      labs(
        title = paste0(
          "PCoA of ",
          dataset_name,
          ", Colored by Study"
        ),
        x = "PCoA1",
        y = "PCoA2"
      )
    
  }
  
  p
}

## TAXA ----
plot_pcoa(
  taxa_pcoa_df,
  color_by = "disease_class",
  dataset_name = "Relative Abundance"
)

plot_pcoa(
  taxa_pcoa_df,
  color_by = "study_name",
  dataset_name = "Relative Abundance"
)

## PAB ----
plot_pcoa(
  pab_pcoa_df,
  color_by = "disease_class",
  dataset_name = "Pathway Abundance"
)

plot_pcoa(
  pab_pcoa_df,
  color_by = "study_name",
  dataset_name = "Pathway Abundance"
)

## PCOV ----
plot_pcoa(
  pcov_pcoa_df,
  color_by = "disease_class",
  dataset_name = "Pathway Coverage"
)

plot_pcoa(
  pcov_pcoa_df,
  color_by = "study_name",
  dataset_name = "Pathway Coverage"
)

# Statistical Support ##########################################################
evaluate_batch_effect <- function(pcoa_result, metadata) {
  
  # Align metadata to the distance matrix
  sample_ids <- attr(pcoa_result$distance, "Labels")
  meta <- metadata[sample_ids, , drop = FALSE]
  
  disease_perm <- adonis2(
    pcoa_result$distance ~ disease_class,
    data = meta
  )
  
  study_perm <- adonis2(
    pcoa_result$distance ~ study_name,
    data = meta
  )
  
  list(
    summary = data.frame(
      disease_R2 = disease_perm$R2[1],
      disease_p = disease_perm$`Pr(>F)`[1],
      study_R2 = study_perm$R2[1],
      study_p = study_perm$`Pr(>F)`[1]
    ),
    disease = disease_perm,
    study = study_perm
  )
}
evaluate_method_set <- function(pcoa_results, metadata) {
  
  evaluation_list <- list()
  
  for (method in names(pcoa_results)) {
    evaluation_list[[method]] <- evaluate_batch_effect(
      pcoa_result = pcoa_results[[method]],
      metadata = metadata
    )
    
    message("PERMANOVA complete for ", method)
  }
  
  evaluation_df <- bind_rows(
    lapply(evaluation_list, `[[`, "summary"),
    .id = "method"
  )
  
  control_study_R2 <- evaluation_df$study_R2[evaluation_df$method == "control"]
  control_disease_R2 <- evaluation_df$disease_R2[evaluation_df$method == "control"]
  
  evaluation_df <- evaluation_df %>%
    mutate(
      study_reduction_pct = round(((control_study_R2 - study_R2) / control_study_R2) * 100, 2),
      disease_change_pct = round(((disease_R2 - control_disease_R2) / control_disease_R2) * 100, 2)
    )
  
  list(
    results = evaluation_list,
    summary = evaluation_df
  )
}

# Evaluations
taxa_eval <- evaluate_method_set(taxa_pcoa_results, studies_meta_df)
pab_eval  <- evaluate_method_set(pab_pcoa_results, studies_meta_df)
pcov_eval  <- evaluate_method_set(pcov_pcoa_results, studies_meta_df)

taxa_evaluation_df <- taxa_eval$summary
pab_evaluation_df  <- pab_eval$summary
pcov_evaluation_df  <- pcov_eval$summary

taxa_evaluation_df
pab_evaluation_df
pcov_evaluation_df