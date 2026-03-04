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

# Now inspect/plot as usual ----------------------------------------------------
library(dplyr)
library(ggplot2)
# --- Taxa: metadata + structure 
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

# --- Pathway Coverage 
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
df_counts <- meta_cov %>%
  count(study_name, disease_class)

ggplot(df_counts, aes(x = study_name, y = n, fill = disease_class)) +
  geom_col() +
  labs(title = "Sample counts by study and disease class",
       x = "Study", y = "N samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_counts, aes(x = disease_class, y = study_name, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 3) +
  labs(title = "Study × disease class counts", x = "Disease class", y = "Study") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(meta_cov, aes(x = disease_class, y = richness_cov)) +
  geom_boxplot() +
  labs(title = "Pathway Coverage Richness by Class", y = "Number of pathways present") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Pathway Abundance 
# Stats ---
meta_ab <- as.data.frame(colData(x_all_path_ab))

nrow(assay(x_all_path_ab)) # number of functional features
ncol(assay(x_all_path_ab)) # number of samples
dim(assay(x_all_path_ab)) # functional x samples
summary(colSums(assay(x_all_path_ab))) #Normalized pathway abundance values (possibly relative functional abundance)
# consistent with HUMAnN output
# if any are 0, that's suspicious (min value)
# which(total_abundance == 0) # uncomment to determine which are 0

# Plots ---
ggplot(x_all_path_ab, aes(x = disease_class, y = total_abundance)) +
  geom_boxplot() +
  labs(title = "Total Pathway Abundance by Class", y = "Sum of pathway abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(x_all_path_ab[x_all_path_ab$total_abundance != 0, ]) + 
  geom_histogram(bins=40) + 
  ggtitle("Distribution of total abundance")

library(scater)

x_pca <- log1p(x_all_path_ab)
x_pca <- runPCA(x_pca)

plotPCA(x_pca, colour_by = "disease_class")
plotPCA(x_pca, colour_by = "study_name")

# Variance filtering (important for ML)
vars <- apply(assay(x_all_path_ab), 1, var)
summary(vars)