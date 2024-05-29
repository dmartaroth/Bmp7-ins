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
file_path <- here("src/BMP7-CT-scans.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# View the first few rows of the dataframe
head(data)

# Filter out rows where wormian_bone is NA
filtered_data <- data[!is.na(data$wormian_bone),]
filtered_data <- filtered_data[!is.na(filtered_data$Genotype),]
filtered_data <- filtered_data[!is.na(filtered_data$Age),]

# Convert relevant columns to factors
filtered_data <- filtered_data %>%
  mutate(wormian_bone = factor(wormian_bone, levels = c(0, 1), labels = c("absent", "present")),
         Genotype = factor(Genotype, levels = c("Bmp7 ctrl", "Bmp7 ncko")),
         Age = factor(Age, levels = c("P0", "P07", "P14", "P30")))


# Check the structure of the filtered data
str(filtered_data)

filtered_data$wormian_bone2 <- relevel(filtered_data$wormian_bone,"present")
# Create the stacked bar plot
(plot <- ggplot(filtered_data, aes(x = Genotype, fill = wormian_bone2)) +
  geom_bar(position = "fill", aes(y = ..count../sum(..count..)),size = 0.5, color = "black") +
  facet_wrap(~ Age) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of Wormian Bone Presence by Age and Genotype",
       x = "Genotype",
       y = "Proportion",
       fill = "Wormian Bone") +
  custom_theme_bar()+scale_fill_manual(values = wormianbonecols))

# Save the plot using the here package for the file path
ggsave(here("Figure-1/01_wormian_bone_frequency_plot.pdf"), plot = plot,width = 5)


# Sex distribution --------------------------------------------------------
filtered_data <- filtered_data[!filtered_data$Sex =="x",]
(plot <- ggplot(filtered_data, aes(x = Age, fill = Sex)) +
   geom_bar(position = "fill", aes(y = ..count../sum(..count..)),size = 0.5, color = "black") +
   scale_y_continuous(labels = scales::percent) +
   labs(title = "Sex Distribution by Age",
        x = "Age",
        y = "Proportion",
        fill = "Sex") +
   custom_theme_bar()+scale_fill_manual(values = sexcols))

# Save the plot using the here package for the file path
ggsave(here("Figure-1/02_sex-distribution.pdf"), plot = plot,width = 3)


# Frequency of interfrontal/Wormian bone ----------------------------------


# Interdigitation at P14 --------------------------------------------------


