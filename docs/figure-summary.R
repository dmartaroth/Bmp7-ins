# ## ######################################## ## #
#                     FIGURES                    #
# ## ######################################## ## #

# Date: Wed May 29 10:17:01 2024 ------------------
# Updated by: Daniela M. Roth


# FIGURE 1 ----------------------------------------------------------------
# Development of midface hypoplasia and correlated midfacial suture anomalies in Bmp7 ncko mice

# Create directory for Figure 1
library(here)
dir.create(here(fig1 <-
                  here("Figure-1")), recursive = TRUE) 
# Load libraries and themes
source(here::here("docs/packages.R"))
source(here::here("docs/themes.R"))

# Load micro-computed tomography data collection spreadsheet
file_path <- here("raw-data/BMP7-CT-scans.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Convert relevant columns to factors
data <- data %>%
  mutate(wormian_bone = factor(wormian_bone, levels = c(0, 1), labels = c("absent", "present")),
         Genotype = factor(Genotype, levels = c("Bmp7 ctrl", "Bmp7 ncko")),
         Age = factor(Age, levels = c("P0", "P07", "P14", "P30")),
         Sex = factor(Sex,levels = c("Female","Male")))

## Frequency of interfrontal/Wormian bone ----------------------------------

data$wormian_bone2 <- relevel(data$wormian_bone, "present")
filtered_data <- data[!is.na(data$Genotype),]
filtered_data <- filtered_data[!is.na(filtered_data$Age),]
wbdata <- filtered_data[!is.na(filtered_data$wormian_bone),]
# Create the stacked bar plot for wormian bone frequency
plot_wormian_bone <- ggplot(wbdata, aes(x = Genotype, fill = wormian_bone2)) +
  geom_bar(position = "fill", aes(y = ..count../sum(..count..)), size = 0.5, color = "black") +
  facet_wrap(~ Age) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of Wormian Bone Presence by Age and Genotype",
       x = "Genotype",
       y = "Proportion",
       fill = "Wormian Bone") +
  custom_theme_bar() +
  scale_fill_manual(values = wormianbonecols)

# Save the plot for wormian bone frequency
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.pdf"), plot = plot_wormian_bone, width = 4.5, height = 3.5)
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.png"), plot = plot_wormian_bone, width = 4.5, height = 3.5)

## Sex distribution --------------------------------------------------------

# Create the plot for sex distribution without filtering NA in wormian_bone
data <- read_excel(file_path, sheet = "Sheet1")

sexdata <- data[!data$Sex =="x",]
plot_sex_distribution <- ggplot(sexdata, aes(x = Age, fill = Sex)) +
  geom_bar(position = "fill", aes(y = ..count../sum(..count..)), size = 0.5, color = "black") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Sex Distribution by Age",
       x = "Age",
       y = "Proportion",
       fill = "Sex") +
  custom_theme_bar() +
  scale_fill_manual(values = sexcols)

# Save the plot for sex distribution
ggsave(here("Figure-1/02_sex-distribution.pdf"), plot = plot_sex_distribution, width = 3)
ggsave(here("Figure-1/02_sex-distribution.png"), plot = plot_sex_distribution, width = 3)

## Interdigitation at P14 --------------------------------------------------
# Ensure required libraries are loaded
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

# Create directory for Figure 1
dir.create(here("Figure-1"), recursive = TRUE)

# Load micro-computed tomography data collection spreadsheet
file_path <- here("raw-data/BMP7-CT-scans.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")


# Ensure the columns 'fpmx_complexity_L', 'fpmx_complexity_R', and 'ID' exist and are correctly spelled
if ("fpmx_complexity_L" %in% colnames(data) & "fpmx_complexity_R" %in% colnames(data) & "ID" %in% colnames(data)) {
  
  # Reshape the data to have separate rows for L and R complexities
  reshaped_data <- data %>%
    pivot_longer(cols = c(fpmx_complexity_L, fpmx_complexity_R),
                 names_to = "Side",
                 values_to = "fpmx_complexity") %>%
    filter(!is.na(fpmx_complexity))
  
  # Filter out rows where Genotype or Age is NA
  reshaped_data <- reshaped_data %>%
    filter(!is.na(Genotype) & !is.na(Age))
  
  # Convert relevant columns to factors
  reshaped_data <- reshaped_data %>%
    mutate(Genotype = factor(Genotype, levels = c("Bmp7 ctrl", "Bmp7 ncko")),
           Age = factor(Age, levels = c("P0", "P07", "P14", "P30")),
           ID = as.factor(ID))
  
  # Calculate average and standard deviation for each genotype and side
  summary_stats <- reshaped_data %>%
    group_by(Genotype) %>%
    summarise(mean_complexity = mean(fpmx_complexity),
              sd_complexity = sd(fpmx_complexity),
              n = n())
  
  # Perform Wilcoxon rank-sum test to compare Bmp7 ctrl and Bmp7 ncko
  wilcox_test <- wilcox.test(fpmx_complexity ~ Genotype, data = reshaped_data)
  
  # Print Wilcoxon rank-sum test results
  cat("\nWilcoxon Rank-Sum Test Results:\n")
  print(wilcox_test)
  
  # Create box plot to show distribution of suture complexity by genotype
  # Create the box plot with modifications
  box_plot <- ggplot(reshaped_data, aes(x = Genotype, y = fpmx_complexity)) +
    geom_boxplot(color = "black",fill = "white",size = 1.7, lineend = "round") +  # No fill, black outline
    labs(title = "Distribution of Suture Complexity by Genotype",
         x = "Genotype",
         y = "Suture Complexity (fpmx_complexity)") +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "#E7FDA6",linewidth = 0.5),  
          panel.grid.minor = element_line(color = "#E7FDA6",linewidth = 0.3)) 
  
  # Calculate mean and standard deviation for annotations
  summary_stats <- reshaped_data %>%
    group_by(Genotype) %>%
    summarise(mean_complexity = mean(fpmx_complexity),
              sd_complexity = sd(fpmx_complexity))
  
  # Add annotations for mean and standard deviation
  box_plot <- box_plot +
    geom_text(data = summary_stats, aes(label = paste("Mean:", round(mean_complexity, 2))),
              x = c(1, 2), y = summary_stats$mean_complexity + 0.5, 
              hjust = 1.1, vjust = 0,size = 3) +
    geom_text(data = summary_stats, aes(label = paste("Stdev:", round(sd_complexity, 2))),
              x = c(1, 2), y = summary_stats$mean_complexity + 1.5, 
              hjust = 1.1, vjust = 2,size = 3)
  
  # Save the modified box plot as a PDF file
  ggsave(here("Figure-1", "03_suture-complexity_by-genotype_distribution.pdf"), plot = box_plot, width = 5, height = 6)
  
  # Save the modified box plot as a png file
  ggsave(here("Figure-1", "03_suture-complexity_by-genotype_distribution.png"), plot = box_plot, width = 5, height = 6)
  
  # Create color palettes
  blue_palette <- scales::seq_gradient_pal("lightblue", "darkblue")(seq(0, 1, length.out = length(unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ctrl"]))))
  red_palette <- scales::seq_gradient_pal("lightpink", "darkred")(seq(0, 1, length.out = length(unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ncko"]))))
  
  # Combine the color palettes
  color_mapping <- c(setNames(blue_palette, unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ctrl"])),
                     setNames(red_palette, unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ncko"])))
  
  # Create shapes
  shapes <- seq(0, length(unique(reshaped_data$ID)) - 1) %% 24 + 1  # 24 unique shapes
  
  # Create scatter plot with text annotations for standard deviation
  plot <- ggplot(reshaped_data, aes(x = Genotype, y = fpmx_complexity, color = ID, shape = ID)) +
    geom_point(position = position_jitter(width = 0.2, height = 0), size = 3, stroke = 1.5) +
    scale_color_manual(values = color_mapping) +
    scale_shape_manual(values = shapes) +
    labs(title = "Suture Complexity by Genotype",
         x = "Genotype",
         y = "Suture Complexity (fpmx_complexity)",
         color = "ID",
         shape = "ID") +
    theme_minimal()+
    theme(plot.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          text = element_text(size = 12))
  
  
  # Save the plot as a PDF file
  ggsave(here("Figure-1", "04_suture-complexity_by-genotype_scatterplot.pdf"), plot = plot, width = 4, height = 4)
  
  # Save the plot as a png file
  ggsave(here("Figure-1", "04_suture-complexity_by-genotype_scatterplot.png"), plot = plot, width = 4, height = 4)
  
  
} else {
  cat("The columns 'fpmx_complexity_L', 'fpmx_complexity_R', and/or 'ID' do not exist in the data frame.\n")
}


## Asymmetry violin plot ---------------------------------------------------

# Load the data from the spreadsheet
file_path <- "raw-data/BMP7-CT-scans.xlsx"
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where Genotype or fpmx_complexity_ratio is NA
filtered_data <- data[complete.cases(data$Genotype, data$fpmx_complexity_ratio), ]

# Convert Genotype to factor with desired levels
filtered_data$Genotype <- factor(filtered_data$Genotype, levels = c("Bmp7 ctrl", "Bmp7 ncko"))

# Create violin plot for fpmx_complexity_ratio by Genotype
violin_plot <- ggplot(filtered_data, aes(x = Genotype, y = fpmx_complexity_ratio, fill = Genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.3) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(fill = Genotype), size = 2.5, alpha = 1, shape = 21) +  # Add jittered points with dodge
  scale_fill_manual(values = c("red3", "dodgerblue3")) +  # Color palette for each genotype
  labs(
    title = "Degree of Asymmetry in Suture Interdigitation by Genotype",
    x = "Genotype",
    y = "fpmx_complexity_ratio",
    fill = "Genotype"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +  # Set y-axis limits from 0.4 to 1
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_line(color = "floralwhite"),
    axis.text = element_text(size = 11, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )


# Save the violin plot as a PDF file
ggsave("Figure-1/05_fpmx_complexity_ratio_violinplot.pdf", plot = violin_plot, width = 6, height = 6)
ggsave("Figure-1/05_fpmx_complexity_ratio_violinplot.png", plot = violin_plot, width = 6, height = 6)

# Perform Wilcoxon rank-sum test
wilcox_test <- wilcox.test(fpmx_complexity_ratio ~ Genotype, data = filtered_data)

# Print Wilcoxon rank-sum test results
cat("\nWilcoxon Rank-Sum Test Results:\n")
print(wilcox_test)

# FIGURE 2 ----------------------------------------------------------------

## Internasal suture/nasal septum deviation --------------------------------






# FIGURE 3 ----------------------------------------------------------------
library(here)
dir.create(here(fig3 <-
                  here("Figure-3")), recursive = TRUE) 
# Load libraries and themes
source(here::here("docs/packages.R"))
source(here::here("docs/themes.R"))

## Osteoblast asymmetry (P14) ----------------------------------------------------

# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# Set file path for the new spreadsheet
file_path <- here("raw-data/alpl-coverage_measurements.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype or alpl_ratio is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(alpl_ratio) & !is.na(ID))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype),
         ID = factor(ID))

# Create color palettes for IDs based on genotype
red_palette <- scales::seq_gradient_pal("lightsalmon", "darkred")(seq(0, 1, length.out = length(unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ctrl"]))))
blue_palette <- scales::seq_gradient_pal("powderblue", "darkblue")(seq(0, 1, length.out = length(unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ncko"]))))

# Combine the color palettes
color_mapping <- c(setNames(red_palette, unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ncko"])),
                   setNames(blue_palette, unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ctrl"])))

# Create violin plot for alpl_ratio by genotype
violin_plot <- ggplot(filtered_data, aes(x = genotype, y = alpl_ratio, fill = genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(color = ID), size = 2.5, alpha = 1) +  # Add jittered points with dodge
  scale_fill_manual(values = c("dodgerblue3", "red3")) +  # Adjust color palette for genotypes
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  labs(
    title = "Symmetry of Alpl+ by Genotype",
    x = "Genotype",
    y = "alpl_ratio",
    fill = "Genotype",
    color = "ID"
  ) +
  scale_y_continuous(limits = c(0, 1)) +  # Adjust y-axis limits as needed
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )

# Print the violin plot
print(violin_plot)

# Save the violin plot as a PDF file
ggsave(
  here("Figure-3/01_alpl_ratio_violinplot.pdf"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("Figure-3/01_alpl_ratio_violinplot.png"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)

# Statistics
# Perform Wilcoxon rank-sum test
wilcox_test <- wilcox.test(alpl_ratio ~ genotype, data = filtered_data)

# Print Wilcoxon rank-sum test results
cat("\nWilcoxon Rank-Sum Test Results:\n")
print(wilcox_test)
# Wilcoxon rank-sum test (aka Mann-Whitney U test) compares distribution of
# central tendencies of 2 independent samples without assuming normality
# If p>0.05: fail to reject null hypothesis, which indicates no statistically
# significant difference in central tendencies (medians) between groups.


# Fligner-Killeen test (non-parametric) compares variances between 2 or more groups
# Perform Fligner-Killeen test for homogeneity of variances
fligner_test <- fligner.test(alpl_ratio ~ genotype, data = filtered_data)

# Print Fligner-Killeen test results
cat("\nFligner-Killeen Test Results:\n")
print(fligner_test)
# p-value less than 0.05: reject null hypothesis; statistically significant
# difference in variances/spread of alpl_ratio between genotypes

# Point 955A appears to be an outlier. Re-analyze with exclusion
# Exclude the 955A datapoint
filtered_data_excluded <- filtered_data %>% filter(ID != "955A")

# Create violin plot for alpl_ratio by genotype excluding 955A
violin_plot_excluded <- ggplot(filtered_data_excluded, aes(x = genotype, y = alpl_ratio, fill = genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(color = ID), size = 2.5, alpha = 1) +  # Add jittered points with dodge
  scale_fill_manual(values = c("dodgerblue3", "red3")) +  # Adjust color palette for genotypes
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  labs(
    title = "Symmetry of Alpl+ by Genotype",
    x = "Genotype",
    y = "alpl_ratio",
    fill = "Genotype",
    color = "ID"
  ) +
  scale_y_continuous(limits = c(0, 1)) +  # Adjust y-axis limits as needed
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )

# Print the violin plot excluding 955A
print(violin_plot_excluded)

# Save the violin plot excluding 955A as a PDF file
ggsave(here("Figure-3/01b_alpl_ratio_violinplot_excluding_955A.pdf"), plot = violin_plot_excluded, width = 7.5, height = 6)
ggsave(here("Figure-3/01b_alpl_ratio_violinplot_excluding_955A.png"), plot = violin_plot_excluded, width = 7.5, height = 6)

# Statistics for the data excluding 955A
# Perform Wilcoxon rank-sum test excluding 955A
wilcox_test_excluded <- wilcox.test(alpl_ratio ~ genotype, data = filtered_data_excluded)

# Print Wilcoxon rank-sum test results excluding 955A
cat("\nWilcoxon Rank-Sum Test Results (excluding 955A):\n")
print(wilcox_test_excluded)
# Wilcoxon rank-sum test (aka Mann-Whitney U test) compares distribution of
# central tendencies of 2 independent samples without assuming normality
# If p>0.05: fail to reject null hypothesis, which indicates no statistically
# significant difference in central tendencies (medians) between groups.

# Perform Fligner-Killeen test for homogeneity of variances excluding 955A
fligner_test_excluded <- fligner.test(alpl_ratio ~ genotype, data = filtered_data_excluded)

# Print Fligner-Killeen test results excluding 955A
cat("\nFligner-Killeen Test Results (excluding 955A):\n")
print(fligner_test_excluded)
# p-value less than 0.05: reject null hypothesis; statistically significant
# difference in variances/spread of alpl_ratio between genotypes

## Osteoclast infiltration (P14) -------------------------------------------------


### TRAP stain quantification -----------------------------------------------


# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# Set file path for the new spreadsheet
file_path <- here("raw-data/TRAP-coverage_measurements.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype or TRAP_ratio is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(TRAP_ratio) & !is.na(ID))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype),
         ID = factor(ID))

# Create color palettes for IDs based on genotype
unique_ctrl_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ctrl"])
unique_ncko_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ncko"])

red_palette <- scales::seq_gradient_pal("lightsalmon", "darkred")(seq(0, 1, length.out = length(unique_ncko_ids)))
blue_palette <- scales::seq_gradient_pal("lightblue1", "navy")(seq(0, 1, length.out = length(unique_ctrl_ids)))

# Combine the color palettes
color_mapping <- c(setNames(red_palette, unique_ncko_ids),
                   setNames(blue_palette, unique_ctrl_ids))
# Create violin plot for TRAP_ratio by genotype
violin_plot <- ggplot(filtered_data, aes(x = genotype, y = TRAP_ratio, fill = genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(color = ID), size = 2.5, alpha = 1) +  # Add jittered points with dodge
  scale_fill_manual(values = c("dodgerblue3", "red3")) +  # Adjust color palette for genotypes
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  labs(
    title = "Symmetry of TRAP+ by Genotype",
    x = "Genotype",
    y = "TRAP_ratio",
    fill = "Genotype",
    color = "ID"
  ) +
  scale_y_continuous(limits = c(0, 1)) +  # Adjust y-axis limits as needed
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )

# Print the violin plot
print(violin_plot)

# Save the violin plot as a PDF file
ggsave(
  here("Figure-3/02_TRAP_ratio_violinplot.pdf"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("Figure-3/02_TRAP_ratio_violinplot.png"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)


### Ctsk IF quantification --------------------------------------------------
# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# Set file path for the new spreadsheet
file_path <- here("raw-data/Ctsk-IF-coverage_measurements.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype or Ctsk_ratio is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(Ctsk_ratio) & !is.na(ID))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype),
         ID = factor(ID))

# Create color palettes for IDs based on genotype
unique_ctrl_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ctrl"])
unique_ncko_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ncko"])

red_palette <- scales::seq_gradient_pal("lightsalmon", "darkred")(seq(0, 1, length.out = length(unique_ncko_ids)))
blue_palette <- scales::seq_gradient_pal("powderblue", "darkblue")(seq(0, 1, length.out = length(unique_ctrl_ids)))

# Combine the color palettes
color_mapping <- c(setNames(red_palette, unique_ncko_ids),
                   setNames(blue_palette, unique_ctrl_ids))

# Create violin plot for Ctsk_ratio by genotype
violin_plot <- ggplot(filtered_data, aes(x = genotype, y = Ctsk_ratio, fill = genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(color = ID), size = 2.5, alpha = 1) +  # Add jittered points with dodge
  scale_fill_manual(values = c("dodgerblue3", "red3")) +  # Adjust color palette for genotypes
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  labs(
    title = "Symmetry of Ctsk+ by Genotype",
    x = "Genotype",
    y = "Ctsk_ratio",
    fill = "Genotype",
    color = "ID"
  ) +
  scale_y_continuous(limits = c(0, 1)) +  # Adjust y-axis limits as needed
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )

# Print the violin plot
print(violin_plot)

# Save the violin plot as a PDF file
ggsave(
  here("Figure-3/03_Ctsk_ratio_violinplot.pdf"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("Figure-3/03_Ctsk_ratio_violinplot.png"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)




# FIGURE 4 ----------------------------------------------------------------

# Visium spatial transcriptomic definition of endo vs ectocranial domains

## Set up Visium objects ---------------------------------------------------
# Source import, processing, and segmentation script
source(here::here("docs/FIG4_01_import-visium-process.R"))


### Visualize Runx2, Cd200, Col2a1 expression -----------------------------------------

# Load packages and set directories
source(here::here("docs/packages.R"))

visium_dir <- here::here("raw-data/Visium")
dir.create(here(fig4 <-
                  here("Figure-4")), recursive = TRUE)



#crop regions 1, 3, and 4 for presentation in figures
insR1Par<-subset(WT1,subset = wt1_imagerow>120 & wt1_imagerow<220 &
                   wt1_imagecol>190 & wt1_imagecol<270)
insR2Par<-subset(WT1,subset = wt1_imagerow>130 & wt1_imagerow<190 &
                   wt1_imagecol>450 & wt1_imagecol<510)
insR3Par<-subset(WT1,subset = wt1_imagerow>355 & wt1_imagerow<415 &wt1_imagecol>445 & wt1_imagecol<505)

WT1allinsPar <- merge(insR1Par,c(insR2Par,insR3Par))

All <- readRDS(here::here("data-output/Visium/P3_All_endovsecto.Rds"))

p1 <- SpatialFeaturePlot(insR1Par, features = "Sost",image.alpha = 0.3,pt.size.factor = 6,stroke=NA)+
  theme(legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=8,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#C6E2FF", "#6CA6CD", "#436EEE", "#27408B"))

p2 <- SpatialFeaturePlot(insR1Par, features = "Cd200",image.alpha = 0.3,pt.size.factor = 6,stroke=NA)+
  theme(legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=8,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c( "white","lightgoldenrod2","gold1","red3","firebrick"))

p3 <- SpatialFeaturePlot(insR1Par, features = "Runx2",image.alpha = 0.3,pt.size.factor = 6,stroke=NA)+
  theme(legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=8,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#C1FFC1", "#4EEE94", "#3CB371", "#006400"))

p4 <- SpatialFeaturePlot(insR1Par, features = "Col2a1",image.alpha = 0.3,pt.size.factor = 6,stroke=NA)+
  theme(legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=8,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c( "white","#FFE1FF", "#DDA0DD", "#CD69C9", "#9932CC"))

plot <- plot_grid(p1,p2,p3,p4,ncol=2)

# Save as pdf
ggsave(here("Figure-4/07_visium_proof-of-concept_R1Par_WT.pdf"), plot = plot,width = 5, height = 5)

# Save as png
ggsave(here("Figure-4/07_visium_proof-of-concept_R1Par_WT.png"), plot = plot,width = 5, height = 5)


### Osteogenesis violin plots -----------------------------------------------
osteogenic <- c("Six2","Sparc","Runx2","Sp7","Dmp1","Sost")
plot <- Stacked_VlnPlot(All,features = osteogenic,colors_use = c("#E066FF", "#FFA54F"),pt.size = 0.5)&
  theme(text = element_text(size=10),
        axis.text.y =element_text(size=8),
        axis.text.x = element_text(size=10),
        plot.title = element_text(face = "bold.italic",size = 10),
        axis.title.y = element_text(angle=90,
                                    face = 'italic',
                                    vjust = -6,
                                    hjust = 0.4))

# Save as pdf
ggsave(here("Figure-4/08_visium_All_osteogenic-vln.pdf"), plot = plot,width = 1.75, height = 7)

# Save as png
ggsave(here("Figure-4/08_visium_All_osteogenic-vln.png"), plot = plot,width = 1.75, height = 7)

### Heatmap of marker genes (endo vs ecto) ----------------------------------

markers <- FindAllMarkers(All,only.pos=FALSE,log.fc.threshold=0.25,test.use = "poisson")
top10 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=10,order_by = avg_log2FC)


DoHeatmap(
  All,
  slot = "data",
  angle = 0,
  size = 3.75,
  group.bar.height = 0.08,
  hjust = 0.5,
  features = top10$gene,
  group.colors = c("#E066FF", "#FFA54F")
) + scale_fill_gradientn(
  colors = c(
    "dodgerblue",
    "lightsteelblue2",
    "oldlace",
    "mistyrose2",
    "tomato2",
    "red2"
  ),
  na.value = "white"
) + theme(
  axis.text.y = element_text(face = "italic"),
  text = element_text(size = 8),
  legend.position = "bottom"
)

### Visualization of selected genes (featplot and violin plots) -------------


## Left vs right segmentation ----------------------------------------------



# FIGURE 6 ----------------------------------------------------------------


## Load snRNA-seq data -----------------------------------------------------

## UMAP --------------------------------------------------------------------


## Marker gene dotplot -----------------------------------------------------

## Nebulosa density plot ---------------------------------------------------

## Bmp7+ cell subset -------------------------------------------------------

### Heatmap Bmp7+ cells -----------------------------------------------------


### ReviGO Molecular Function Bmp7+ cells -----------------------------------------------




# FIGURE 7 ----------------------------------------------------------------


## Bmp7 RNAscope quantification P0 -----------------------------------------
# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)
library(scales)  # For color palette
library(ggpubr)  # For statistical annotation
library(ggfun)

# Set file path for the new spreadsheet
file_path <- here("raw-data/P0_Bmp7_symmetry.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1") %>%
  select(ID, bmp7ratio)

# Convert ID to factor
data <- data %>%
  mutate(ID = factor(ID))

# Filter out rows where bmp7ratio is NA
data <- data %>%
  filter(!is.na(bmp7ratio))

# Perform Wilcoxon signed-rank test
wilcox_test <- wilcox.test(data$bmp7ratio, mu = 1, alternative = "two.sided")

# Create color palette for IDs (shades of blue)
blue_palette <- scales::seq_gradient_pal("powderblue", "darkblue")(seq(0, 1, length.out = length(unique(data$ID))))
color_mapping <- setNames(blue_palette, unique(data$ID))

# Create strip plot for Bmp7 ratio
bmp7_ratio_strip_plot <- ggplot(data, aes(x = ID, y = bmp7ratio, color = ID)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +  # Jittered points
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Horizontal line at ratio 1
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(
    title = "Bmp7+ L/R symmetry",
    x = NULL,
    y = "Bmp7+/DAPI+ symmetry",
    color = "ID",
    subtitle = paste("Wilcoxon signed-rank test p-value:", signif(wilcox_test$p.value, digits = 3))
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    plot.title = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8, face = "italic"),  # Subtitle size
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    text = element_text(size = 8),
    legend.background = element_blank(),
    legend.box.background = element_roundrect(colour = "black")
  )

# Print the strip plot
print(bmp7_ratio_strip_plot)

# Save the strip plot as a PDF file
ggsave(
  here("Figure-7/01_Bmp7_symmetry.pdf"),
  plot = bmp7_ratio_strip_plot,
  width = 3,
  height = 4
)

# Save the strip plot as a PNG file
ggsave(
  here("Figure-7/01_Bmp7_symmetry.png"),
  plot = bmp7_ratio_strip_plot,
  width = 3,
  height = 4
)

# Print the result of the Wilcoxon signed-rank test
print(wilcox_test)




## Osteocyte maturation quantification P0 ----------------------------------

# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

# Set file path for the new spreadsheet
file_path <- here("raw-data/osteocyte-maturation_measurements.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet2")

# View the first few rows of the dataframe
head(data)

# Filter out rows where genotype or total_L or total_R is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(total_L) & !is.na(total_R))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype),
         ID = factor(ID))

# Specify the column names for left and right sides
left_cols <- c("Dmp1_L", "Sost_L", "Dmp1_Sost_L", "Dmp1_Runx2_L", "Sost_Runx2_L", "Dmp1_Sost_Runx2_L")
right_cols <- c("Dmp1_R", "Sost_R", "Dmp1_Sost_R", "Dmp1_Runx2_R", "Sost_Runx2_R", "Dmp1_Sost_Runx2_R")

# Reshape the data to long format for left and right sides
long_data_L <- filtered_data %>%
  pivot_longer(cols = all_of(left_cols),
               names_to = "cell_type",
               values_to = "count",
               names_pattern = "(.*)_L") %>%
  mutate(total = total_L)

long_data_R <- filtered_data %>%
  pivot_longer(cols = all_of(right_cols),
               names_to = "cell_type",
               values_to = "count",
               names_pattern = "(.*)_R") %>%
  mutate(total = total_R)

# Add a side column to each long data
long_data_L <- long_data_L %>% mutate(side = "left")
long_data_R <- long_data_R %>% mutate(side = "right")

# Combine the left and right data
long_data <- bind_rows(long_data_L, long_data_R)

# Calculate proportions for each cell type
proportion_data <- long_data %>%
  mutate(proportion = count / total)

# Set the order of cell_type levels
proportion_data$cell_type <- factor(proportion_data$cell_type, 
                                    levels = c("Dmp1_Sost_Runx2", "Dmp1_Runx2", "Dmp1","Dmp1_Sost", "Sost", "Sost_Runx2"))

# View the proportion data
head(proportion_data)

# Define custom colors for cell types
custom_colors <- c("Dmp1" = "green",
                   "Dmp1_Sost" = "blue",
                   "Sost" = "violet",
                   "Sost_Runx2" = "lavender",
                   "Dmp1_Runx2" = "yellow",
                   "Dmp1_Sost_Runx2" = "red")

# Create the stacked bar plot with facets by ID
stacked_bar_plot <- ggplot(proportion_data, aes(x = interaction(genotype, side), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +  # Stacked bar plot with fill position
  scale_fill_manual(values = custom_colors) +  # Use custom colors
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  labs(
    title = "Proportions of Cell Types by Genotype and Side",
    x = "Genotype and Side",
    y = "Proportion of Cell Types",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold")
  ) +
  facet_wrap(~ID, ncol = 5)  # Create facets for each ID

# Print the stacked bar plot
print(stacked_bar_plot)

# Save the stacked bar plot as a PDF and PNG file
ggsave(here("Figure-7/02_osteocyte_maturation_stacked_bar_plot.pdf"), plot = stacked_bar_plot, width = 16, height = 12)
ggsave(here("Figure-7/02_osteocyte_maturation_stacked_bar_plot.png"), plot = stacked_bar_plot, width = 16, height = 12)


# Second visualization approach

# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# Set file path for the new spreadsheet
file_path <- here("raw-data/osteocyte-maturation_measurements.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype or asymm_oct_mat is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(asymm_oct_mat) & !is.na(ID))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype),
         ID = factor(ID))

# Create color palettes for IDs based on genotype
unique_ctrl_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ctrl"])  # Adjust "control" as needed
unique_ncko_ids <- unique(filtered_data$ID[filtered_data$genotype == "Bmp7 ncko"])  # Adjust "experimental" as needed

red_palette <- scales::seq_gradient_pal("lightsalmon", "darkred")(seq(0, 1, length.out = length(unique_ncko_ids)))
blue_palette <- scales::seq_gradient_pal("powderblue", "darkblue")(seq(0, 1, length.out = length(unique_ctrl_ids)))

# Combine the color palettes
color_mapping <- c(setNames(red_palette, unique_ncko_ids),
                   setNames(blue_palette, unique_ctrl_ids))

# Create violin plot for asymm_oct_mat by genotype
violin_plot <- ggplot(filtered_data, aes(x = genotype, y = asymm_oct_mat, fill = genotype)) +
  geom_violin(trim = FALSE, width = 0.8, alpha = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add box plot for clearer representation of median and quartiles
  geom_point(position = position_jitterdodge(), aes(color = ID), size = 2.5, alpha = 1) +  # Add jittered points with dodge
  scale_fill_manual(values = c("dodgerblue3", "red3")) +  # Adjust color palette for genotypes
  scale_color_manual(values = color_mapping) +  # Set the colors for each ID
  labs(
    title = "Symmetry of Osteocyte Maturation by Genotype",
    x = "Genotype",
    y = "asymm_oct_mat",
    fill = "Genotype",
    color = "ID"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),  # White grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "black"),  # Text size and color
    axis.title = element_text(size = 12, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 12, face = "bold")  # Plot title size and style
  )

# Print the violin plot
print(violin_plot)

# Save the violin plot as a PDF file
ggsave(
  here("Figure-7/03_osteocyte_maturation_asymmetry_violinplot.pdf"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("Figure-7/03_osteocyte_maturation_asymmetry_violinplot.png"),
  plot = violin_plot,
  width = 7.5,
  height = 6
)


# Perform Wilcoxon rank-sum test with continuity correction
wilcox_test <- wilcox.test(asymm_oct_mat ~ genotype, data = filtered_data, exact = FALSE)

# Print Wilcoxon rank-sum test results
cat("\nWilcoxon Rank-Sum Test Results:\n")
print(wilcox_test)
# Wilcoxon rank-sum test (aka Mann-Whitney U test) compares distribution of
# central tendencies of 2 independent samples without assuming normality
# If p > 0.05: fail to reject null hypothesis, which indicates no statistically
# significant difference in central tendencies (medians) between groups.

wilcox_result <- wilcox_test
# Extract U statistic and sample sizes
U <- wilcox_result$statistic
n1 <- sum(filtered_data$genotype == levels(filtered_data$genotype)[1])  # Sample size of genotype 1
n2 <- sum(filtered_data$genotype == levels(filtered_data$genotype)[2])  # Sample size of genotype 2

# Calculate common language effect size (CL)
CL <- U / (n1 * n2)

# Print the calculated CL
cat("Common Language Effect Size (CL):", CL, "\n")
# High effect size suggests notable difference between variables despite low
# p-value in Wilcoxon rank-sum test. For example, CL=0.9 indicates 90%
# probability that a randomly chosen measurement from Bmp7 ctrl will have a
# higher symmetry of asymm_oct_mat than one randomly selected from Bmp7 ncko


