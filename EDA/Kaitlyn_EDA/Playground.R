################################################################################
# --- How to access things inside the objects/nested list: ---
# Create a data frame from the nested list:
nested_list <- CRC_progression_studies[["FengQ_2015"]][["relative_abundance"]]

metadata_df <- as.data.frame(colData(nested_list))

# --- Now play with original things like in EDA: ---

# Summary Stats
table(metadata_df$study_condition, useNA = "ifany")
table(metadata_df$age_decade, useNA = "ifany")
table(colData(nested_list)$disease_class, useNA = "ifany")
table(colData(nested_list)$disease, colData(nested_list)$disease_class, useNA = "ifany")
summary(metadata_df$age)
summary(metadata_df$BMI)

# Plots
barplot(table(metadata_df$study_condition), las = 2, main = "Study condition counts")
hist(metadata_df$age, breaks = 20, main = "Age distribution", xlab = "Age")

mat <- assay(x)   # features x samples
dim(mat)
v <- mat[1, ]
hist(v, breaks = 50, main = rownames(mat)[1], xlab = "Value")
summary(colSums(mat))
range(mat, na.rm = TRUE)


################################################################################
# Merging nested data → this will be important for LODO!!

# --- To merge from the nested data  ---
library(mia)  # for mergeData
library(SummarizedExperiment)

merge_all_studies <- function(nested, type) {
  objs <- lapply(nested, `[[`, type) # extracts the objects for one type, objs becomes a list with one element per study
  objs <- objs[!vapply(objs, is.null, logical(1))] # removes missing entries, returns logical TRUE/FALSE vector, that's flipped to drop nulls
  mergeData(objs)
}

x_all_taxa <- merge_all_studies(CRC_progression_studies, "relative_abundance")
x_all_path_cov <- merge_all_studies(CRC_progression_studies, "pathway_coverage")
x_all_path_ab <- merge_all_studies(CRC_progression_studies, "pathway_abundance")

# Now inspect/plot as usual ####################################################
library(dplyr)
library(ggplot2)
library(forcats)
library(scater)
################################################################################
# --- Taxa: metadata + structure -----------------------------------------------
meta_taxa <- as.data.frame(colData(x_all_taxa))

length(unique(meta_taxa$study_name))
table(meta_taxa$study_name)
table(meta_taxa$disease_class, useNA = "ifany")
table(meta_taxa$disease, meta_taxa$disease_class, useNA = "ifany")

nrow(assay(x_all_taxa))   # features 
ncol(assay(x_all_taxa))   # samples
length(unique(meta_taxa$study_name))
summary(colSums(assay(x_all_taxa))) # sums around 100 (percent scale) or 1 (proportion scale)
dim(assay(x_all_taxa)) # dimensions of the feature matrix, meaning "x" taxa by "y" samples

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- PCoA or PCA depending on Transformation ---> 
# PCA on taxa requires compositional transformation (e.g., CLR)
# PCoA with Bray–Curtis is a standard ecological approach for relative abundance data.

library(vegan)
mat_taxa <- assay(x_all_taxa)

# Run PCoA and store results inside object
dist_taxa <- vegdist(t(mat_taxa), method = "bray")
pcoa_taxa <- cmdscale(dist_taxa, k = 2, eig = TRUE)

# Build plotting data frame
meta_taxa <- as.data.frame(colData(x_all_taxa))

df_pcoa_t <- cbind(
  meta_taxa,
  PC1 = pcoa_taxa$points[, 1],
  PC2 = pcoa_taxa$points[, 2]
)

# Approx % variance explained (Bray–Curtis is non-Euclidean, so we use positive eigenvalues only)
eig_vals <- pcoa_taxa$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])

# Plots ---
df_pcoa_t %>% 
#  filter(disease_class %in% c("HC", "CRC", "PA")) %>%  # Comment out to remove filter  
  ggplot(aes(PC1, PC2, color = disease_class)) +
  geom_point(alpha = 0.6) +
#  stat_ellipse(level = 0.95) + # comment out to remove ellipses
  facet_wrap(~disease_class) + #comment out to overlay
  scale_color_manual(values = progression_colors, 
                     name = progression_legend_title) +
  theme_bw() +
  labs(title = "PCoA (Bray–Curtis) of Taxonomic Data, Colored by CRC Progression Classifier",
       x = paste0("PCoA1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(var_explained[2] * 100, 1), "%)")
  )

ggplot(df_pcoa_t, aes(PC1, PC2, color = study_name)) +
  geom_point(alpha = 0.6) +
#  stat_ellipse(level = 0.95) + # comment out to remove ellipses
#  facet_wrap(~study_name) + #comment out to overlay
  theme_bw() +
  labs(title = "PCoA (Bray–Curtis) Colored by Study", 
       x = paste0("PCoA1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(var_explained[2] * 100, 1), "%)")
  )

# --- PERMANOVA (quantify separation)
adonis2(dist_taxa ~ disease_class, data = meta_taxa) # Gives R2 value for disease effect
adonis2(dist_taxa ~ study_name, data = meta_taxa) # Gives R2 value for study effect

################################################################################
# --- Pathway Coverage ---------------------------------------------------------
# Stats ---
meta_cov <- as.data.frame(colData(x_all_path_cov))

nrow(assay(x_all_path_cov)) # number of pathway features
ncol(assay(x_all_path_cov)) # number of samples
dim(assay(x_all_path_ab)) # pathways x samples
length(unique(meta_cov$study_name)) # number of studies
summary(colSums(assay(x_all_path_cov))) # not a percent scale; it’s a coverage/count-like metric

table(meta_cov$disease_class, useNA="ifany")
tab_cov <- table(meta_cov$study_name, meta_cov$disease_class)
tab_cov

richness_cov <- colSums(assay(x_all_path_cov) > 0)
summary(richness_cov)

# Plots ---
#options(device = "RStudioGD") # Sometimes required to fix rendering issues on HPC
df_counts <- meta_cov %>%
  count(study_name, disease_class)

ggplot(df_counts, # Barplot for study count with disease class fill
       aes(x = fct_reorder(study_name, n, .fun = sum, .desc = TRUE),
           y = n,
           fill = disease_class)) +
  geom_col() +
  labs(title = "Sample counts by study and disease class",
       x = "Study", y = "N samples") +
  scale_fill_manual(values = progression_colors, 
                    name = progression_legend_title) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_counts, # Heat map of study and disease class with count
       aes(x = disease_class, y = study_name, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 3) +
  labs(title = "Study × disease class counts", x = "Disease class", y = "Study") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(meta_cov, aes(x = richness_cov)) + 
  geom_histogram(bins=40) + 
  ggtitle("Distribution of Pathway Coverage Richness")

ggplot(meta_cov, # Box plot of pathway coverage by disease class
       aes(x = disease_class, y = richness_cov)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, color = "red") +
  labs(title = "Pathway Coverage Richness by Class", y = "Number of pathways present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- PCoA Analysis ---> Coverage presence/absence-like data, better analyzed via distance metrics like Jaccard
library(vegan)
mat_cov <- assay(x_all_path_cov)
mat_bin <- mat_cov > 0   # convert to presence/absence

# Run PCoA and store results inside object
dist_mat <- vegdist(t(mat_bin), method = "jaccard")
pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)

# Build plotting data frame
meta_cov <- as.data.frame(colData(x_all_path_cov))

df_pcoa <- cbind(
  meta_cov,
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

eig_vals <- pcoa$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])

# Plots ---
df_pcoa %>% 
  filter(disease_class %in% c("HC", "CRC", "PA")) %>%  # Comment out to remove filter  
  ggplot(aes(PC1, PC2, color = disease_class)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) + # comment out to remove ellipses
#  facet_wrap(~disease_class) + #comment out to overlay
  scale_color_manual(values = progression_colors, 
                     name = progression_legend_title) +
  theme_bw() +
  labs(title = "PCoA (Jaccard) of Pathway Coverage",
       x = paste0("PCoA1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(var_explained[2] * 100, 1), "%)")
      )

ggplot(df_pcoa, aes(PC1, PC2, color = study_name)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) + # comment out to remove ellipses
#  facet_wrap(~study_name) + #comment out to overlay
  theme_bw() +
  labs(title = "PCoA (Jaccard) Colored by Study", 
       x = paste0("PCoA1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(var_explained[2] * 100, 1), "%)")
       )

# --- PERMANOVA (quantify separation)
adonis2(dist_mat ~ disease_class, data = meta_cov) # Gives R2 value for disease effect
adonis2(dist_mat ~ study_name, data = meta_cov) # Gives R2 value for study effect

################################################################################
# --- Pathway Abundance --------------------------------------------------------
# Stats ---
meta_ab <- as.data.frame(colData(x_all_path_ab))
meta_ab$total_abundance <- colSums(assay(x_all_path_ab))

nrow(assay(x_all_path_ab)) # number of functional features
ncol(assay(x_all_path_ab)) # number of samples
dim(assay(x_all_path_ab)) # functional x samples
summary(colSums(assay(x_all_path_ab))) #Normalized pathway abundance values (possibly relative functional abundance)
# consistent with HUMAnN output
# if any are 0, that's suspicious (min value)
# which(total_abundance == 0) # uncomment to determine which are 0

# Plots ---
ggplot(meta_ab[meta_ab$total_abundance != 0, ], # Box plot of pathway abundance by disease class
       aes(x = disease_class, y = total_abundance)) + 
  # there is a single 0 value that needs to be filtered out, as it is an outlier (biologically impossible)
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, color = "red") +
  labs(title = "Total Pathway Abundance by Class", y = "Sum of pathway abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(meta_ab[meta_ab$total_abundance != 0, ],
       aes(x = total_abundance)) + 
  # there is a single 0 value that needs to be filtered out, as it is an outlier (biologically impossible)
  geom_histogram(bins=40) + 
  ggtitle("Distribution of total abundance")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- PCA Analysis ---> Pathway abundance is continuous, Euclidean-friendly data
mat <- log1p(assay(x_all_path_ab))   # features x samples

# Optional: variance filter to make PCA faster/cleaner (Important particularly for ML)
#vars <- apply(mat, 1, var)
#summary(vars)
#keep <- vars > quantile(vars, 0.75)   # top 25% most variable features
#mat <- mat[keep, ]

# Run PCA and store results inside object
pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)  # samples x features
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

# Build plotting data frame
meta_ab <- as.data.frame(colData(x_all_path_ab))
df_pca <- cbind(
  meta_ab, 
  PC1 = pca$x[,1], 
  PC2 = pca$x[,2]
  )

# Plots
df_pca %>% 
#  filter(disease_class %in% c("HC", "CRC", "PA")) %>%  # Comment out to remove filter
  ggplot(aes(PC1, PC2, color = disease_class)) +
  geom_point(alpha = 0.6) +
#  stat_ellipse(level = 0.95) + # comment out to remove ellipses
#  facet_wrap(~disease_class) + #comment out to overlay
  scale_color_manual(values = progression_colors, 
                    name = progression_legend_title) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
    title = "PCA of Pathway Abundance, Colored by CRC Progression Classifier"
  )

df_pca %>% 
  ggplot(aes(PC1, PC2, color = study_name)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) + # comment out to remove ellipses
  facet_wrap(~study_name) + #comment this line out to overlay
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
    title = "PCA of Pathway Abundance, Colored by Study")
