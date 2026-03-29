##--- SETUP ---###############################
suppressPackageStartupMessages({
  library(tidyverse)
  library(knitr)
  library(kableExtra)
  library(patchwork)
  
  library(vegan)
  library(mia)
  library(scater)
  library(lefser)
  
  library(grid)
  library(gridExtra)
})

progression_order <- c(
  "HC",
  "PA",
  "PA+",
  "PA-M",
  "CRC",
  "CRC+",
  "CRC-H",
  "CRC-M", 
  "Other"
)

progression_colors <- c(
  "HC" = "#539deb",
  "PA" = "#e4e85b",
  "PA+" = "#e8a25b",
  "PA-M" = "#d4700d",
  "CRC" = "#cf1919",
  "CRC+" = "#8f0000",
  "CRC-M" = "#4d0404",
  "CRC-H" = "#520f76",
  "Other" = "#aaa3a3"
)

progression_legend_title <- "CRC Progression\nClassifiers"
# DATAFRAME CREATION ###########################################################
eda_meta_df <- studies_meta_df %>%  filter(study_name != "HMP_2012") # make a copy of filtered dataframe from step 3.1 in FullML_Playground.R

# METADATA COMPLETENESS ########################################################
meta_cols <- setdiff(colnames(eda_meta_df), params$meta_drop)

meta_completeness <- data.frame(
  variable = meta_cols,
  completeness = colMeans(!is.na(eda_meta_df[meta_cols])) * 100
)

meta_completeness %>%
  arrange(completeness) %>%
  mutate(variable = factor(variable, levels = variable)) %>%
  ggplot(aes(x = "", y = variable, fill = completeness)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Metadata Completeness",
    subtitle = "Percent non-missing values",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# AGE NORMALITY ####################################################################

# Histogram ---
normality1 <- eda_meta_df %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = age)) + 
  geom_histogram(bins=20) + 
  labs(title = "Distribution of Age: Overall", x = "Age", y = "Count") +
  theme(plot.title = element_text(hjust = 0, size = 10), 
        axis.title = element_text(size = 10))

# Box Plot: Normality by CRC Category ---
normality2 <- eda_meta_df %>%   
  filter(!is.na(age)) %>%   
  ggplot(aes(x=disease_class, y=age)) +    
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(alpha = 0.07, width = 0.15) +   
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, color = "red") +
  labs(title = "Distribution of Age: By CRC Progression Classifier", 
     x = progression_legend_title, y = "Age") +
  theme(plot.title = element_text(hjust = 0, size = 10), 
        axis.title = element_text(size = 10))

# --- CRC ACROSS AGE/GENDER ---
# Gendered Barplot: Proportions ---
decade1 <- eda_meta_df %>%
  filter(!is.na(age_decade), !is.na(gender)) %>%
  ggplot(aes(x = age_decade, fill = disease_class)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = progression_colors, 
                    name = progression_legend_title) +
  facet_wrap(~gender) +
  labs(title = "CRC Progression Proportion by Age Decade, Grouped by Gender", 
     x = "Age by Decade", y = "Proportion") +
  theme(plot.title = element_text(size = 12))

# Display Patchwork
layout1 <- c( # t, l , b = t, r = l
  area(1, 1, 2, 2), 
  area(1, 3, 2, 4),
  area(3, 1, 5, 4), 
  area(3, 5, 5, 5)
)

age_gender_caption <- str_wrap(
  "(A) Histogram demonstrating a bimodal age distribution with a secondary peak near 25 years and a primary peak near 60 years. 
   (B) Box plot of age across CRC progression classifiers, with jittered transluscent points to demonstrate spread of the sample. The red dot represents the mean, the line represents the median.
   (C) Bar plot displaying the gender faceted proportion of CRC progression classifiers across each decade. Notable differences include and visual increase in HC from 18-69 for females, and increased PA and CRC/CRC+ across most decades in males.", 
  width = 150
)

normality1 + normality2 + decade1 + guide_area() + 
  plot_layout(design = layout1, guides = "collect", widths = c(1,1,1,1,0.7)) +
  plot_annotation(
    # title = "Figure 1: Age Normality & CRC Variation across Age and Gender",
    tag_levels = "A",
    # caption = age_gender_caption
  ) &
  theme(
    plot.title = element_text(hjust = 0, size = 14),
    plot.caption = element_text(hjust = 0, size = 12), 
    legend.key.size = unit(0.4, "cm"), 
    legend.title = element_text(size = 12)
  )

# BMI NORMALITY / ASSOCIATION ##################################################

# Histogram ---
normality3 <- eda_meta_df %>%
  filter(!is.na(BMI)) %>%
  ggplot(aes(x = BMI)) + 
  geom_histogram(bins=20) + 
  labs(title = "Distribution of BMI: Overall", x = "BMI", y = "Count") +
  theme(plot.title = element_text(hjust = 0, size = 10), 
        axis.title = element_text(size = 10))

# Box Plot ---
normality4 <- eda_meta_df %>%   
  filter(!is.na(BMI)) %>%   
  ggplot(aes(x=disease_class, y=BMI)) +    
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(alpha = 0.15, width = 0.15) +   
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, color = "red") +
  labs(
    title = "Distribution of BMI: By CRC Progression Classifier", 
    x = "CRC Progression Classifier", 
    y = "BMI"
  ) +
  theme(plot.title = element_text(size = 12))

# Scatterplot by CRC Progression --
bmi1 <- eda_meta_df %>% 
  filter (!is.na(age), !is.na(BMI)) %>% 
  ggplot(aes(x=age, y=BMI)) +    
  geom_point(aes(color = disease_class), alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "grey40") +
  scale_color_manual(values = progression_colors, 
                     name = progression_legend_title) +
  facet_wrap(~disease_class) +
  labs(
    title = "Age by BMI, Colored by CRC Progression Classifier", 
    x = "Age", 
    y = "BMI"
  ) +
  theme(plot.title = element_text(size = 12))

BMI_caption <- str_wrap(
  "(A) Histogram of BMI demonstrating _______ distribution with ______________.
   (B) Box plot of BMI by CRC progression classifier. The red dot marks the mean, while the black line marks the median. All classifiers generally have a normal distribution.
   (C) Scatter plots of age by BMI, with a linear model regression line of y ~ x, faceted by CRC progression classifier. Notice that there is minimally notable clustering.", 
  width = 200
)

normality3 + normality4 + bmi1 + guide_area() + 
  plot_layout(design = layout1, guides = "collect", widths = c(1,1,1,1,0.7)) +
  plot_annotation(
    # title = "Figure 2: BMI Normality & CRC Variation Associations with Age",
    tag_levels = "A",
    # caption = BMI_caption
  ) &
  theme(
    plot.title = element_text(hjust = 0, size = 14),
    plot.caption = element_text(hjust = 0, size = 12), 
    legend.key.size = unit(0.4, "cm"), 
    legend.title = element_text(size = 12)
  )
# AGE ANOVA ####################################################################
