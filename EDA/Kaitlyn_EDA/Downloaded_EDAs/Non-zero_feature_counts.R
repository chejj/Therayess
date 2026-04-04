## Count how many Taxa/Pathways there are in each study

taxa_counts <- sapply(unique(studies_meta_df$study_name), function(study) {
  
  samples <- rownames(studies_meta_df)[studies_meta_df$study_name == study]
  
  sum(colSums(taxa_df[samples, , drop = FALSE] > 0) > 0)
})

path_ab_counts <- sapply(unique(studies_meta_df$study_name), function(study) {
  
  samples <- rownames(studies_meta_df)[studies_meta_df$study_name == study]
  
  sum(colSums(path_abund_df[samples, , drop = FALSE] > 0) > 0)
})

path_cov_counts <- sapply(unique(studies_meta_df$study_name), function(study) {
  
  samples <- rownames(studies_meta_df)[studies_meta_df$study_name == study]
  
  sum(colSums(path_cov_df[samples, , drop = FALSE] > 0) > 0)
})

study_feature_summary <- data.frame(
  Study = unique(studies_meta_df$study_name),
  Samples = as.numeric(table(studies_meta_df$study_name)),
  Taxa = taxa_counts,
  Pathway_Abundance = path_ab_counts,
  Pathway_Coverage = path_cov_counts
)

study_feature_summary %>% 
  arrange(-Samples)