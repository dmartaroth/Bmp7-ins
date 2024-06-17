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
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.pdf"), plot = plot_wormian_bone, width = 5)
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.png"), plot = plot_wormian_bone, width = 5)

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
  
  # Create color palettes
  blue_palette <- scales::seq_gradient_pal("lightblue", "darkblue")(seq(0, 1, length.out = length(unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ctrl"]))))
  red_palette <- scales::seq_gradient_pal("lightpink", "darkred")(seq(0, 1, length.out = length(unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ncko"]))))
  
  # Combine the color palettes
  color_mapping <- c(setNames(blue_palette, unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ctrl"])),
                     setNames(red_palette, unique(reshaped_data$ID[reshaped_data$Genotype == "Bmp7 ncko"])))
  
  # Create shapes
  shapes <- seq(0, length(unique(reshaped_data$ID)) - 1) %% 24 + 1  # 24 unique shapes
  
  # Create scatter plot with text annotations for standard deviation
  ggplot(reshaped_data, aes(x = Genotype, y = fpmx_complexity, color = ID, shape = ID)) +
    geom_point(position = position_jitter(width = 0.2, height = 0), size = 3, stroke = 1.5) +
    scale_color_manual(values = color_mapping) +
    scale_shape_manual(values = shapes) +
    labs(title = "Suture Complexity by Genotype",
         x = "Genotype",
         y = "Suture Complexity (fpmx_complexity)",
         color = "ID",
         shape = "ID") +
    theme_minimal()
  
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
    title = "Alpl Ratio by Genotype",
    x = "Genotype",
    y = "Alpl Ratio",
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
    axis.title = element_text(size = 14, face = "bold"),  # Axis title size and style
    plot.title = element_text(size = 16, face = "bold")  # Plot title size and style
  )

# Print the violin plot
print(violin_plot)

# Save the violin plot as a PDF file
ggsave(here("Figure-3/01_alpl_ratio_violinplot.pdf"), plot = violin_plot, width = 7.5, height = 6)
ggsave(here("Figure-3/01_alpl_ratio_violinplot.png"), plot = violin_plot, width = 7.5, height = 6)

## Osteoclast infiltration (P14) -------------------------------------------------

## Osteogenesis (P7) -------------------------------------------------------






# FIGURE 4 ----------------------------------------------------------------

# Visium spatial transcriptomic definition of endo vs ectocranial domains

## Set up Visium objects ---------------------------------------------------
# Source import, processing, and segmentation script
source(here::here("docs/FIG4_01_import-visium-process.R"))


### Visualize Runx2, Cd200, Col2a1 expression -----------------------------------------

## Visualize segmentation of endo vs ecto gems ---------------------------------------

### Osteogenesis violin plots -----------------------------------------------

### Heatmap of marker genes (endo vs ecto) ----------------------------------

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


## Osteocyte maturation quantification P0 ----------------------------------















