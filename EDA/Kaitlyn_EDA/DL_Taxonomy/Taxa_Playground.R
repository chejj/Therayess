# Table
#################################################################################
# Single Study, nested list (Example)
Thomas2018a_taxa <- CRC_progression_studies[["ThomasAM_2018a"]][["relative_abundance"]]

# All studies combined to one list
taxa <- merge_studies(CRC_progression_studies, "relative_abundance")

# Metadata into dataframes from nested list
Thomas2018a_taxa_meta <- as.data.frame(colData(Thomas2018a_taxa))
taxa_meta <- as.data.frame(colData(taxa))

# Overview Tables
Thomas2018a_taxa_meta_overview <- eda_overview_table(Thomas2018a_taxa_meta)
taxa_meta_overview <- eda_overview_table(taxa_meta)

Thomas2018a_taxa_meta_overview %>%
  arrange(Percent_Missing) %>%
  mutate(Index = row_number()) %>%              
  select(Index, everything()) %>%               
  kbl(
    col.names = c("Index",
                  "Variable Name",
                  "Total Rows",
                  "Empty Rows",
                  "Percent Missing",
                  "Unique Values",
                  "Variable Type"),
    caption = "Structural Metadata Summary for ThomasAM_2018a Relative Abundance Data",
    format = "html"
  ) %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(2, bold = TRUE)

taxa_meta_overview %>%
  arrange(Percent_Missing) %>%
  mutate(Index = row_number()) %>%              
  select(Index, everything()) %>%               
  kbl(
    col.names = c("Index",
                  "Variable Name",
                  "Total Rows",
                  "Empty Rows",
                  "Percent Missing",
                  "Unique Values",
                  "Variable Type"),
    caption = "Structural Metadata Summary for All Relative Abundance Data",
    format = "html"
  ) %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(2, bold = TRUE)

################################################################################
# Heatmaps
################################################################################
# Check that Pathway col names and taxa col names are identical (metadata are the same)
identical(
  colnames(assay(CRC_progression_studies$FengQ_2015$relative_abundance)),
  colnames(assay(CRC_progression_studies$FengQ_2015$pathway_abundance))
)
#if TRUE, then the meta data analysis will be the same across all of them

# Heatmap of Metadata Completeness
taxa_meta_overview %>%
  arrange(desc(Percent_Missing)) %>%
  mutate(Column_Name = factor(Column_Name, levels = Column_Name)) %>%  # lock order
  ggplot(aes(x = "", y = Column_Name, fill = Percent_Missing)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "steelblue",
    high = "white",
    limits = c(0, 100)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Structural Metadata Completeness: All Studies",
    subtitle = "Percent rows with missing values",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7, lineheight = 2),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

Thomas2018a_taxa_meta_overview %>%
  arrange(desc(Percent_Missing)) %>%
  mutate(Column_Name = factor(Column_Name, levels = Column_Name)) %>%  # lock order
  ggplot(aes(x = "", y = Column_Name, fill = Percent_Missing)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "steelblue",
    high = "white",
    limits = c(0, 100)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Structural Metadata Completeness:\nThomasAM_2018a - Relative Abundance Data",
    subtitle = "Percent rows with missing values",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7, lineheight = 2),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5)
  )

################################################################################
# CRC Progression Classifier Counts
################################################################################
taxa_meta %>% 
  ggplot(aes(x=disease_class, stat = "count", fill = disease_class)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.21, size = 3.5, color = "black") +
  scale_fill_manual(values = progression_colors, 
                    name = progression_legend_title) +
  # Labeling
  #---------
labs(y = "Sample Count", x = "CRC Progression Classes") +
  ggtitle("Counts of CRC Progression Classes in Study")

################################################################################
# CRC Progression Classifier w/ Proportions
################################################################################
# 

######
taxa_hierarchy_table %>% 
  kbl(
    col.names = c("Kingdom", "Phylum", "Class", "Mean Abundance (%)"),
    caption = "Mean Relative Abundance by Kingdom, Phylum, and Class",
    format = "html",
    digits = 2,
    booktabs = TRUE
  ) %>% 
  kable_styling(full_width = FALSE) %>% 
  column_spec(1, bold = TRUE) %>% 
  column_spec(2, bold = TRUE) %>% 
  collapse_rows(columns = 1:2, valign = "top")

taxa_hierarchy_table_order %>%
  kbl(
    col.names = c("Kingdom", "Phylum", "Class", "Order", "Mean Abundance (%)"),
    caption = "Mean Relative Abundance by Kingdom → Phylum → Class → Order",
    format = "html",
    digits = 2
  ) %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(2, bold = TRUE) %>%
  column_spec(3, bold = TRUE) %>%
  collapse_rows(columns = 1:3, valign = "top")