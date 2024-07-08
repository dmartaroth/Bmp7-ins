# ## ######################################## ## #
#                     FIGURES                    #
# ## ######################################## ## #

# Date: Wed May 29 10:17:01 2024 ------------------
# Updated by: Daniela M. Roth


# FIGURE 1 ----------------------------------------------------------------
# Development of midface hypoplasia and correlated midfacial suture anomalies in Bmp7 ncko mice

# Create directory for Figure 1
library(here)
dir.create(here("Figure-1"), recursive = TRUE) 

# Load libraries and themes
source(here::here("docs/packages.R"))
source(here::here("docs/themes.R"))


# Wormian bone observations -----------------------------------------------


# Load CT data
file_path <- here("raw-data/BMP7-CT-scans.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Adjust Age category for P08 entries and convert relevant columns to factors
data <- data %>%
  mutate(Age = ifelse(Age == "P08", "P07", Age),
         wormian_bone = factor(wormian_bone, levels = c(0, 1), labels = c("absent", "present")),
         Genotype = factor(Genotype, levels = c("Bmp7 ctrl", "Bmp7 ncko")),
         Age = factor(Age, levels = c("P0", "P07", "P14", "P30")),
         Sex = factor(Sex, levels = c("Female", "Male")))

# Filter out rows with NA in key columns
filtered_data <- data %>%
  filter(!is.na(Genotype) & !is.na(Age) & !is.na(wormian_bone))

# Relevel wormian_bone factor
filtered_data <- filtered_data %>%
  mutate(wormian_bone2 = relevel(wormian_bone, "present"))

# Function to check for normality and run appropriate tests
perform_stat_tests <- function(df) {
  results <- df %>%
    group_by(Age) %>%
    summarize(p_value = {
      # Create a contingency table
      tbl <- table(Genotype, wormian_bone)
      print(paste("Age:", unique(Age)))  # Debug: print current age group
      print(tbl)  # Debug: print contingency table
      
      # Check if the test should be Chi-squared or Fisher's exact test
      if (all(tbl >= 5)) {
        print("Using Chi-squared test")  # Indicate which test is used
        chisq_test_result <- chisq.test(tbl)
        chisq_test_result$p.value
      } else {
        print("Using Fisher's exact test")  # Indicate which test is used
        fisher_test_result <- fisher.test(tbl)
        fisher_test_result$p.value
      }
    }) %>%
    mutate(significant = ifelse(p_value < 0.05, "Yes", "No"))
  
  return(results)
}

# Apply the function to the data
stat_results <- perform_stat_tests(filtered_data)

# Report the number of measurements per genotype for each age
count_per_genotype <- filtered_data %>%
  group_by(Age, Genotype) %>%
  summarize(count = n(), .groups = 'drop')

# Print the statistical results and the counts
print(stat_results)
print(count_per_genotype)

# Create the stacked bar plot for wormian bone frequency
plot_wormian_bone <- ggplot(filtered_data, aes(x = Genotype, fill = wormian_bone2)) +
  geom_bar(position = "fill", aes(y = after_stat(count) / sum(after_stat(count))), linewidth = 0.5, color = "black") +
  facet_wrap(~ Age) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Wormian bone frequency",
       x = "Genotype",
       y = "Proportion",
       fill = "Observation") +
  custom_theme_bar() +
  scale_fill_manual(values = wormianbonecols)

# Save the plot for wormian bone frequency
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.pdf"), plot = plot_wormian_bone, width = 4.5, height = 3.5)
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.png"), plot = plot_wormian_bone, width = 4.5, height = 3.5)

## Sex distribution --------------------------------------------------------
library(ggplot2)
library(readxl)
library(here)
library(dplyr)

# Load the data
file_path <- here("raw-data/BMP7-CT-scans.xlsx")
data <- read_excel(file_path, sheet = "Sheet1")

# Adjust Age categories, filter out entries with 'Sex' marked as 'x' or blank
sexdata <- data %>%
  mutate(Age = case_when(
    Age == "P08" ~ "P07",
    Age == "P15" ~ "P14",
    TRUE ~ as.character(Age)
  )) %>%
  filter(Sex != "x", !is.na(Sex), Sex != "")  # Exclude invalid entries

# Group data by Age and Sex, calculate counts
plot_data <- sexdata %>%
  group_by(Age, Sex) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Total = sum(Count), Percentage = Count / Total) %>%
  arrange(Age, desc(Sex))

# Calculate position for labels
plot_data <- plot_data %>%
  group_by(Age) %>%
  mutate(ypos = cumsum(Percentage) - 0.5 * Percentage)

# Adjust the label positions differently for each age group
plot_data <- plot_data %>%
  mutate(ypos_adjusted = case_when(
    Age == "P07" ~ ypos * 10,  # Spread out the labels for P07 the most
    Age == "P14" ~ ypos * 4.5,  # Spread out the labels for P14 a little less
    Age == "P30" ~ ypos * 1.5      
  ))

# Perform binomial tests for each age group
binom_results <- plot_data %>%
  group_by(Age) %>%
  summarise(
    Male_Count = sum(Count[Sex == "M"]),
    Female_Count = sum(Count[Sex == "F"]),
    Total_Count = sum(Count),
    P_Value = binom.test(c(sum(Count[Sex == "M"]), sum(Count[Sex == "F"])), p = 0.5)$p.value,
    Significant = ifelse(P_Value < 0.05, "Yes", "No")
  )

# Print the binomial test results
print(binom_results)

# Create the plot for sex distribution using correct data structure
plot_sex_distribution <- ggplot(plot_data, aes(x = Age, y = Count, fill = Sex, label = Count)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_label(aes(y = ypos_adjusted),
             size = 3.5, fontface = "bold", color = "black", fill = "white", 
             show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Sex distribution",
       x = "Age",
       y = "Proportion",
       fill = "Sex") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        text = element_text(size = 8)) +
  scale_fill_manual(values = c("F" = "white", "M" = "black"))

# Save the plot for sex distribution
ggsave(here("Figure-1/02_sex-distribution.pdf"), plot = plot_sex_distribution, width = 3, height = 4)
ggsave(here("Figure-1/02_sex-distribution.png"), plot = plot_sex_distribution, width = 4, height = 3)


## Interdigitation at P14 --------------------------------------------------
# Ensure required libraries are loaded
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)

# Function to calculate Common Language Effect Size (CLES)
cles <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)
  U <- sum(outer(group1, group2, ">"))
  return(U / (n1 * n2))
}

# Create directory for Figure 1
dir.create(here("Figure-1"), recursive = TRUE)

# Load micro-computed tomography data collection spreadsheet
file_path <- here("raw-data/BMP7-CT-scans.xlsx")
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
  
  # Calculate Common Language Effect Size (CLES)
  bmp7_ctrl_values <- reshaped_data %>% filter(Genotype == "Bmp7 ctrl") %>% pull(fpmx_complexity)
  bmp7_ncko_values <- reshaped_data %>% filter(Genotype == "Bmp7 ncko") %>% pull(fpmx_complexity)
  cles_value <- cles(bmp7_ctrl_values, bmp7_ncko_values)
  
  # Print Wilcoxon rank-sum test results and CLES
  cat("\nWilcoxon Rank-Sum Test Results:\n")
  print(wilcox_test)
  cat("\nCommon Language Effect Size (CLES):", cles_value, "\n")
  
  # Create box plot to show distribution of suture complexity by genotype
  box_plot <- ggplot(reshaped_data, aes(x = Genotype, y = fpmx_complexity)) +
    geom_boxplot(color = "black", fill = "white", size = 1.7, lineend = "round") +  # No fill, black outline
    labs(title = "Distribution of Suture Complexity by Genotype",
         x = "Genotype",
         y = "Suture Complexity (fpmx_complexity)") +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "#E7FDA6", linewidth = 0.5),  
          panel.grid.minor = element_line(color = "#E7FDA6", linewidth = 0.3)) 
  
  # Calculate mean and standard deviation for annotations
  summary_stats <- reshaped_data %>%
    group_by(Genotype) %>%
    summarise(mean_complexity = mean(fpmx_complexity),
              sd_complexity = sd(fpmx_complexity))
  
  # Add annotations for mean and standard deviation
  box_plot <- box_plot +
    geom_text(data = summary_stats, aes(label = paste("Mean:", round(mean_complexity, 2))),
              x = c(1, 2), y = summary_stats$mean_complexity + 0.5, 
              hjust = 1.1, vjust = 0, size = 3) +
    geom_text(data = summary_stats, aes(label = paste("Stdev:", round(sd_complexity, 2))),
              x = c(1, 2), y = summary_stats$mean_complexity + 1.5, 
              hjust = 1.1, vjust = 2, size = 3)
  
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

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nortest)
library(irr)
library(here)


### Plotting ----------------------------------------------------------------


# Load the data
file_path <- here::here("raw-data/raw_scoring-deviation_KS-DMR.xlsx")
data <- read_excel(file_path)

# Remove unnecessary average columns
data <- data %>%
  select(-contains("average"))

# Automatically correct column names to valid R names
names(data) <- make.names(names(data), unique = TRUE)

# Rename columns for easier manipulation
data <- data %>%
  rename(
    ins_angle_KS1 = angle.of.internasal.suture.deviation_KS...5,
    ins_angle_KS2 = angle.of.internasal.suture.deviation_KS...6,
    ins_angle_KS3 = angle.of.internasal.suture.deviation_KS...7,
    ins_angle_DMR1 = angle.of.internasal.suture.deviation_DMR...9,
    ins_angle_DMR2 = angle.of.internasal.suture.deviation_DMR...10,
    ins_angle_DMR3 = angle.of.internasal.suture.deviation_DMR...11,
    suture_angle_KS1 = angle.of.suture.to.septum_KS...14,
    suture_angle_KS2 = angle.of.suture.to.septum_KS...15,
    suture_angle_KS3 = angle.of.suture.to.septum_KS...16,
    suture_angle_DMR1 = angle.of.suture.to.septum_DMR...18,
    suture_angle_DMR2 = angle.of.suture.to.septum_DMR...19,
    suture_angle_DMR3 = angle.of.suture.to.septum_DMR...20
  )

# Print the structure after renaming
print(str(data))

# Calculate the average for each technical replicate
data <- data %>%
  mutate(
    ins_angle_KS = rowMeans(select(., starts_with("ins_angle_KS")), na.rm = TRUE),
    ins_angle_DMR = rowMeans(select(., starts_with("ins_angle_DMR")), na.rm = TRUE),
    suture_angle_KS = rowMeans(select(., starts_with("suture_angle_KS")), na.rm = TRUE),
    suture_angle_DMR = rowMeans(select(., starts_with("suture_angle_DMR")), na.rm = TRUE)
  ) %>%
  select(block_id, section, age, genotype, ins_angle_KS, ins_angle_DMR, suture_angle_KS, suture_angle_DMR)

# Print the structure after calculating averages
print(str(data))

# Pivot longer to combine KS and DMR
data_long <- data %>%
  pivot_longer(
    cols = -c(block_id, section, age, genotype),
    names_to = c("measurement", "rater"),
    names_pattern = "(.+)_(.+)"
  ) %>%
  mutate(value = as.numeric(value))

# Print the structure after pivot longer
print(str(data_long))

# Clean the data
data_cleaned <- data_long %>%
  filter(complete.cases(.)) %>%
  mutate(
    genotype = recode(genotype, 'CTRL' = 'ctrl', 'NCKO' = 'ncko'),
    age = recode(age, 'P14' = '2-week', 'P15' = '2-week', 'P16' = '2-week', 'P7' = '1-week')
  )

# Calculate average value for each measurement and rater
data1 <- data_cleaned %>%
  group_by(block_id, age, genotype, measurement) %>%
  summarise(mean_value = mean(value), .groups = 'drop')

# Pivot wider to create a simpler table for plotting and analysis
data1 <- data1 %>%
  pivot_wider(names_from = measurement, values_from = mean_value)

# Ensure the data has the required columns and filter for complete cases
data1 <- data1 %>%
  filter(complete.cases(.))

# Print the structure after final pivot
print(str(data1))

# Function to create convex hull data
find_hull <- function(df) {
  df[chull(df$suture_angle, df$ins_angle), ]
}

# Calculate the convex hull for each combination of age and genotype
hulls <- data1 %>%
  group_by(age, genotype) %>%
  do(find_hull(.))

# Create the scatter plot
plot <- ggplot(data1, aes(x = suture_angle, y = ins_angle, color = genotype)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  geom_polygon(data = hulls, aes(fill = genotype, group = interaction(age, genotype)), alpha = 0.2) +
  scale_fill_manual(values = c("cornflowerblue","red1"))+
  geom_point(size=1)+ 
  facet_wrap(~ age,nrow = 2)+  
  scale_color_manual(values = c("cornflowerblue", "red1"))+
  theme_minimal() +
  theme(panel.background=element_rect(fill="seashell"),
        plot.background = element_rect(fill = "white", color = "white"),
        text = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10))+
  labs(
    x = "suture.septum.angle",
    y = "ins.angle",
    color = "Genotype",
    fill = "Genotype"
  ) +
  geom_rug()

# Print the plot to view it
print(plot)

# Optionally save the plot
ggsave(here("Figure-2", "01_genotype_measurement_comparison.png"), plot = plot, width = 2.5, height = 4)


### Sample size -------------------------------------------------------------


# Summarize the number of mice for each genotype and age group
mouse_summary <- data1 %>%
  group_by(age, genotype) %>%
  summarise(count = n_distinct(block_id), .groups = 'drop')

# Print the summary
print(mouse_summary)


### Normality ---------------------------------------------------------------



# Check normality for each group, switch tests based on sample size
normality_results <- data1 %>%
  pivot_longer(cols = c(suture_angle, ins_angle), names_to = "measurement", values_to = "value") %>%
  group_by(age, genotype, measurement) %>%
  summarise(
    count = n(),
    normality_test = if(count > 7) {
      ad.test(value)$p.value  # Anderson-Darling test for normality
    } else if(count >= 3) {
      shapiro.test(value)$p.value  # Shapiro-Wilk test for smaller samples
    } else {
      NA_real_  # Not enough data to perform any normality test
    },
    .groups = 'drop'
  ) %>%
  mutate(is_normal = ifelse(is.na(normality_test), FALSE, normality_test > 0.05))  # Common threshold for normality

print(normality_results)

# Based on normality, decide which test to apply
test_results <- data1 %>%
  pivot_longer(cols = c(suture_angle, ins_angle), names_to = "measurement", values_to = "value") %>%
  left_join(normality_results, by = c("age", "genotype", "measurement")) %>%
  group_by(age, measurement) %>%
  reframe(
    test_result = if_else(
      is_normal, 
      t.test(value ~ genotype, data = .)$p.value,  # Parametric test for normal data
      wilcox.test(value ~ genotype, data = ., exact = FALSE)$p.value  # Non-parametric test for non-normal data
    ),
    method = if_else(is_normal, "t-test", "Wilcoxon test")
  )

print(test_results, n = 38)


### Reliability (intra- and inter-rater) ------------------------------------


library(irr)

# Inter-rater reliability calculation for ins_angle
icc_ins_angle <- data %>%
  select(ins_angle_KS, ins_angle_DMR) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")

# Inter-rater reliability calculation for suture_angle
icc_suture_angle <- data %>%
  select(suture_angle_KS, suture_angle_DMR) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")

print(icc_ins_angle)
print(icc_suture_angle)

# Load the data again to get original replicates
data <- read_excel(file_path)

# Remove unnecessary average columns
data <- data %>%
  select(-contains("average"))

# Automatically correct column names to valid R names
names(data) <- make.names(names(data), unique = TRUE)

# Rename columns for easier manipulation
data <- data %>%
  rename(
    ins_angle_KS1 = angle.of.internasal.suture.deviation_KS...5,
    ins_angle_KS2 = angle.of.internasal.suture.deviation_KS...6,
    ins_angle_KS3 = angle.of.internasal.suture.deviation_KS...7,
    ins_angle_DMR1 = angle.of.internasal.suture.deviation_DMR...9,
    ins_angle_DMR2 = angle.of.internasal.suture.deviation_DMR...10,
    ins_angle_DMR3 = angle.of.internasal.suture.deviation_DMR...11,
    suture_angle_KS1 = angle.of.suture.to.septum_KS...14,
    suture_angle_KS2 = angle.of.suture.to.septum_KS...15,
    suture_angle_KS3 = angle.of.suture.to.septum_KS...16,
    suture_angle_DMR1 = angle.of.suture.to.septum_DMR...18,
    suture_angle_DMR2 = angle.of.suture.to.septum_DMR...19,
    suture_angle_DMR3 = angle.of.suture.to.septum_DMR...20
  )


#### Intra-rater reliability (KS) --------------------------------------------

icc_ins_angle_ks <- data %>%
  select(ins_angle_KS1, ins_angle_KS2, ins_angle_KS3) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")

icc_suture_angle_ks <- data %>%
  select(suture_angle_KS1, suture_angle_KS2, suture_angle_KS3) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")


#### Intra-rater reliability (DMR) -------------------------------------------

icc_ins_angle_dmr <- data %>%
  select(ins_angle_DMR1, ins_angle_DMR2, ins_angle_DMR3) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")

icc_suture_angle_dmr <- data %>%
  select(suture_angle_DMR1, suture_angle_DMR2, suture_angle_DMR3) %>%
  as.matrix() %>%
  icc(model = "twoway", type = "agreement", unit = "single")

print(icc_ins_angle_ks)
print(icc_suture_angle_ks)
print(icc_ins_angle_dmr)
print(icc_suture_angle_dmr)

# Perform ANOVA for ins_angle
anova_ins <- aov(ins_angle ~ genotype * age, data = data1)
summary(anova_ins)

# Perform ANOVA for suture_angle
anova_suture <- aov(suture_angle ~ genotype * age, data = data1)
summary(anova_suture)

# Post-hoc analysis if significant interaction is found
TukeyHSD(anova_ins)
TukeyHSD(anova_suture)


# Print the ANOVA results for ins_angle
print(summary(anova_ins))

# Print the ANOVA results for suture_angle
print(summary(anova_suture))

# Print post-hoc results if there are significant interactions
print(TukeyHSD(anova_ins))
print(TukeyHSD(anova_suture))

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

# Calculate the number of unique IDs per genotype
unique_ids_per_genotype <- filtered_data_excluded %>%
  group_by(genotype) %>%
  summarise(unique_ids = n_distinct(ID), .groups = 'drop')

# Print the summary
print(unique_ids_per_genotype)


# Function to calculate Common Language Effect Size (CLES)
cles <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)
  U <- sum(outer(group1, group2, ">"))
  return(U / (n1 * n2))
}

# Calculate the Common Language Effect Size (CLES)
cles_value <- cles(
  filtered_data_excluded$alpl_ratio[filtered_data_excluded$genotype == "Bmp7 ctrl"],
  filtered_data_excluded$alpl_ratio[filtered_data_excluded$genotype == "Bmp7 ncko"]
)

# Print the CLES
cat("\nCommon Language Effect Size (excluding 955A):\n")
print(cles_value)
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

# Function to calculate Common Language Effect Size (CLES)
cles <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)
  U <- sum(outer(group1, group2, ">"))
  return(U / (n1 * n2))
}

# Statistics for the Ctsk_ratio data
# Perform Wilcoxon rank-sum test
wilcox_test <- wilcox.test(Ctsk_ratio ~ genotype, data = filtered_data, exact = FALSE)

# Print Wilcoxon rank-sum test results
cat("\nWilcoxon Rank-Sum Test Results:\n")
print(wilcox_test)

# Perform Fligner-Killeen test for homogeneity of variances
fligner_test <- fligner.test(Ctsk_ratio ~ genotype, data = filtered_data)

# Print Fligner-Killeen test results
cat("\nFligner-Killeen Test Results:\n")
print(fligner_test)

# Calculate the number of unique IDs per genotype
unique_ids_per_genotype <- filtered_data %>%
  group_by(genotype) %>%
  summarise(unique_ids = n_distinct(ID), .groups = 'drop')

# Print the summary
print(unique_ids_per_genotype)

# Calculate the Common Language Effect Size (CLES)
cles_value <- cles(
  filtered_data$Ctsk_ratio[filtered_data$genotype == "Bmp7 ctrl"],
  filtered_data$Ctsk_ratio[filtered_data$genotype == "Bmp7 ncko"]
)

# Print the CLES
cat("\nCommon Language Effect Size:\n")
print(cles_value)


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
top30 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=30,order_by = avg_log2FC)

#top30.markers <- write.csv(top30,file = here(Figure-4/top30.endoectovisium.csv"))

# Using the top30 markers from FindAllMarkers, inputted genes for endo and ecto independently into <
#   https: /  / go.princeton.edu / cgi - bin / GOTermFinder > to predict Molecular Function from Gene Ontology.

# A plotting R script produced by the Revigo server at http://revigo.irb.hr/

library( ggplot2 )
library( scales )


revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0005488","binding",56.68847260979142,6.639871457866139,-0.025316426798752567,7.318248684861919,-5.379943301884653,1,-0),
                     c("GO:0005515","protein binding",5.073796515930947,-4.454026136073058,-3.7684425879131025,6.27008718737355,-6.409348626556189,0.9870947282679754,0.03236795),
                     c("GO:0016209","antioxidant activity",0.6670203825798291,2.3208400327934116,5.780463919658815,5.388894787172177,-3.9012394061847426,0.6671995787114399,0),
                     c("GO:0016684","oxidoreductase activity, acting on peroxide as acceptor",0.4512193500693483,2.4525682571111345,-5.9268634688881,5.21914424598055,-2.608700296857601,0.8647623242138107,-0),
                     c("GO:0031720","haptoglobin binding",0.0004985367537640342,-4.535435220802938,3.4672222472370393,2.2648178230095364,-2.2654784565911106,0.9870947282679754,-0));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );


library(ggrepel)

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(1) );
p1 <- p1 + scale_colour_gradientn( colours = c("#FFFFFF", "#FFDAB9", "#FF8247", "#FF7F24", "#EE7600"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("dodgerblue", 3) ));
p1 <- p1 + scale_size( range=c(3.5, 8)) + theme_bw();  #+ scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + ggtitle("Ectocranial ReviGO Molecular Function")+theme(legend.key = element_blank(),title = element_text(size=8)) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# output plot
p1

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).
#ggsave(here(Figure-4/P3insvisium_ecto_GOMF.pdf),width = 4, height = 4)

# Setting the max cut-off for spatial feature plots to be just above the maximum
## internasal suture value (from Violin Plots). Eg *Col9a1* goes up to 250 based
## on cartilage, but to interpret suture difference will need to decrease maximum
## for sensitivity (down to 40 in this case).

p1 <- SpatialFeaturePlot(insR1Par, 
                         features = "Col9a1",
                         image.alpha = 0.3,
                         pt.size.factor = 5,
                         stroke=NA,
                         max.cutoff = 40)+
  theme(legend.position="right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size=7,face = "italic"),
        legend.justification = "center",
        axis.text.y = element_text(face = "italic",
                                   size = 4,
                                   angle = 90),
        axis.text.x = element_text(face = "italic",
                                   size = 4,
                                   angle = 90))+
  ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#FFF0F5", "#FFAEB9", "#FF69B4", "#FF1493"))

p2 <- SpatialFeaturePlot(insR1Par, 
                         features = "Mgp",
                         image.alpha = 0.3,
                         pt.size.factor = 5,
                         stroke=NA,
                         max.cutoff = 30)+
  theme(legend.position="right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size=7,face = "italic"),
        legend.justification = "center",
        axis.text.y = element_text(face = "italic",
                                   size = 4,
                                   angle = 90),
        axis.text.x = element_text(face = "italic",
                                   size = 4,
                                   angle = 90))+
  ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#C1FFC1", "#4EEE94", "#3CB371", "#006400"))

p3 <- SpatialFeaturePlot(insR1Par, 
                         features = "Hmox1",
                         image.alpha = 0.3,
                         pt.size.factor = 5,
                         stroke=NA,
                         max.cutoff = NA)+
  theme(legend.position="right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size=7,face = "italic"),
        legend.justification = "center",
        axis.text.y = element_text(face = "italic",
                                   size = 4,
                                   angle = 90),
        axis.text.x = element_text(face = "italic",
                                   size = 4,
                                   angle = 90))+
  ggplot2::scale_fill_gradientn(colors=c( "white","lightgoldenrod2","gold1","red3","firebrick"))

p4 <- SpatialFeaturePlot(insR1Par, 
                         features = "Gpha2",
                         image.alpha = 0.3,
                         pt.size.factor = 5,
                         stroke=NA,
                         max.cutoff = NA)+
  theme(legend.position="right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size=7,face = "italic"),
        legend.justification = "center",
        axis.text.y = element_text(face = "italic",
                                   size = 4,
                                   angle = 90),
        axis.text.x = element_text(face = "italic",
                                   size = 4,
                                   angle = 90))+
  ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#C6E2FF", "#6CA6CD", "#436EEE", "#27408B"))

plot_grid(p1,p2,p3,p4,ncol=4)+theme(plot.margin = unit(c(0,0,0,0),"cm"))

ggsave(here(Figure-4/P3ins_visium_spatfeatureplots_GOMF.pdf),width = 10, height = 2.2)

my_genes <- c("Col9a1","Mgp","Hmox1","Gpha2")
Stacked_VlnPlot(All,
                features = my_genes,
                colors_use = c("#E066FF", "#FFA54F"),
                pt.size = 0.5)&
  theme(text = element_text(size=10),
        axis.text.y =element_text(size=8,),
        axis.text.x = element_text(size=10),
        plot.title = element_text(face = "bold.italic",size = 10),
        axis.title.y = element_text(angle=90,face = 'italic',vjust = -4))

#ggsave(here(Figure-4/P3ins_visium_vlnplot_GOMF.pdf),width = 2.5, height = 4.5)

## Left vs right segmentation ----------------------------------------------

# Future direction for further analysis

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
library(tidyr)
library(here)
library(scales)
library(ggpubr)
library(car)  # For Levene's test

# Function to calculate Common Language Effect Size (CLES)
cles <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)
  U <- sum(outer(group1, group2, function(x, y) x > y))
  return(U / (n1 * n2))
}

# Set file path for the new spreadsheet
file_path <- here("raw-data/P0_Bmp7_symmetry.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1") %>%
  select(ID, unnorm_ratio, DAPI_ratio)

# Convert ID to factor
data <- data %>%
  mutate(ID = factor(ID))

# Filter out rows where unnorm_ratio is NA
data <- data %>%
  filter(!is.na(unnorm_ratio))

# Descriptive statistics
summary_stats <- data %>%
  summarise(
    mean_bmp7 = mean(unnorm_ratio),
    sd_bmp7 = sd(unnorm_ratio),
    median_bmp7 = median(unnorm_ratio),
    mean_dapi = mean(DAPI_ratio),
    sd_dapi = sd(DAPI_ratio),
    median_dapi = median(DAPI_ratio)
  )
print(summary_stats)

# Normality tests
shapiro_test_bmp7 <- shapiro.test(data$unnorm_ratio)
shapiro_test_dapi <- shapiro.test(data$DAPI_ratio)

# Fligner-Killeen Test for homogeneity of variances
fligner_test <- fligner.test(list(data$unnorm_ratio, data$DAPI_ratio))

# Paired Wilcoxon signed-rank test for Bmp7 ratio vs. DAPI ratio
wilcox_test_paired <- wilcox.test(data$unnorm_ratio, data$DAPI_ratio, paired = TRUE, alternative = "two.sided")

# Wilcoxon rank-sum test for Bmp7 ratio vs. DAPI ratio
wilcox_test_comparison <- wilcox.test(data$unnorm_ratio, data$DAPI_ratio, alternative = "two.sided")


# Calculate CLES for paired data
cles_value_paired <- cles(data$unnorm_ratio, data$DAPI_ratio)

# Print the results
cat("\nShapiro-Wilk Normality Test Results (Bmp7 ratio):\n")
print(shapiro_test_bmp7)

cat("\nShapiro-Wilk Normality Test Results (DAPI ratio):\n")
print(shapiro_test_dapi)

cat("\nFligner-Killeen Test Results (Homogeneity of variances):\n")
print(fligner_test)

cat("\nWilcoxon Rank-Sum Test Results (Bmp7 ratio vs. DAPI ratio):\n")
print(wilcox_test_comparison)

cat("\nPaired Wilcoxon Signed-Rank Test Results (Bmp7 ratio vs. DAPI ratio):\n")
print(wilcox_test_paired)

cat("\nCommon Language Effect Size (CLES) for paired Bmp7 ratio vs. DAPI ratio:\n")
print(cles_value_paired)

# Rename unnorm_ratio to Bmp7_ratio
data <- data %>%
  rename(Bmp7_ratio = unnorm_ratio)

# Combine data into one dataframe with a grouping variable and reorder levels
long_data <- data %>%
  pivot_longer(cols = c(DAPI_ratio, Bmp7_ratio),  # Adjusted to use Bmp7_ratio
               names_to = "ratio_type",
               values_to = "ratio_value") %>%
  filter(!is.na(ratio_value)) %>%
  mutate(ratio_type = factor(ratio_type, levels = c("Bmp7_ratio", "DAPI_ratio")))  # Adjusted levels

# Debug: Print the structure and head of the long_data
print(str(long_data))
print(head(long_data))

# Debug: Print the unique values of ratio_type
print(unique(long_data$ratio_type))

# Create box plot for ratios with reordered levels
ratio_box_plot <- ggplot(long_data, aes(x = ratio_type, y = ratio_value, fill = ratio_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 1) +  # Jittered points
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Horizontal line at ratio 1
  scale_fill_manual(values = c("forestgreen", "violet")) +
  labs(
    y = "Symmetry",
    x = "",
    fill = "Ratio Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    text = element_text(size = 8),
    plot.background = element_rect(color = "white"),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

# Print the box plot
print(ratio_box_plot)

# Save the box plot as a PDF file
ggsave(
  here("Figure-7/01b_Ratio_Comparison.pdf"),
  plot = ratio_box_plot,
  width = 3.75,
  height = 2.5
)

# Save the box plot as a PNG file
ggsave(
  here("Figure-7/01b_Ratio_Comparison.png"),
  plot = ratio_box_plot,
  width = 3.75,
  height = 2.5
)


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


# Second visualization approach (main)

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

# Get the number of unique IDs per genotype
unique_ids_per_genotype <- filtered_data %>%
  group_by(genotype) %>%
  summarise(unique_ids = n_distinct(ID), .groups = 'drop')

# Print the summary of unique IDs per genotype
print(unique_ids_per_genotype)

# Get descriptive statistics for each genotype
descriptive_stats <- filtered_data %>%
  group_by(genotype) %>%
  summarise(
    mean_asymm = mean(asymm_oct_mat, na.rm = TRUE),
    sd_asymm = sd(asymm_oct_mat, na.rm = TRUE),
    median_asymm = median(asymm_oct_mat, na.rm = TRUE),
    min_asymm = min(asymm_oct_mat, na.rm = TRUE),
    max_asymm = max(asymm_oct_mat, na.rm = TRUE),
    .groups = 'drop'
  )

# Print the descriptive statistics
print(descriptive_stats)

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
  width = 7,
  height = 4
)
ggsave(
  here("Figure-7/03_osteocyte_maturation_asymmetry_violinplot.png"),
  plot = violin_plot,
  width = 7,
  height = 4
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


