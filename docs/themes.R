# ## ######################################## ## #
#                      THEMES                    #
# ## ######################################## ## #

# Date: Wed May 29 11:18:33 2024 ------------------
# Updated by: Daniela M. Roth

# Colors
wormianbonecols <- c( "violet","floralwhite")
sexcols <- c("white","black")
my_colors <- c("goldenrod1", "red2", "dodgerblue", "green3", "orange")


error <- red $ bold
note1 <- red $ bold
note2 <- black $ bold
note3 <- blue $ bold

# Custom theme for barplot
custom_theme_bar <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(color = "white", fill = "white"),
      panel.grid.major = element_blank(),           # No major gridlines
      panel.grid.minor = element_blank(),           # No minor gridlines
      # legend.position = "none",                     # No legend
      axis.ticks = element_line(color = "black"),
      axis.ticks.x = element_line(color = NA),
      axis.title = element_text(color = "black",size = 8),
      axis.text = element_text(size = 8),           # Size of axis text
      axis.line.y = element_line(colour = NA), # Set y-axis line color
      axis.line.x = element_line(colour = NA),  
      axis.text.y = element_text(colour = "black"), # Set y-axis text color
      axis.ticks.length = unit(0.2, "cm"),          # Shorten tick length
      panel.border = element_blank(),   
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),         # No vertical gridlines
      text = element_text(colour = "black", size = 8),                        # Regular text
      plot.title = element_text(face = "bold", size = 8,hjust = 0),     # Bold plot title
      # plot.margin = margin(20, 20, 20, 20),         # Adjust plot margins
      plot.caption = element_text(color = "black")  # Caption color and position
    )
}


# Custom theme for density plot
custom_theme_density <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),                        # No vertical gridlines
      
      # Axis
      axis.ticks = element_line(color = "black"),                  # Color of axis ticks
      axis.line.y = element_line(colour = "black"),                # Color of y-axis line
      axis.line.x = element_line(colour = NA),                     # Remove x-axis line
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(size = 10),                         # Size of axis text
      axis.text.y = element_text(colour = "black"),                # Color of y-axis text
      axis.ticks.length = unit(0.2, "cm"),                         # Shorten tick length
      
      # Panel
      panel.border = element_rect(color = NA, fill = NA),     # Panel border
      
      # Text
      plot.title = element_text(face = "bold", vjust = 0.5),       # Bold plot title, centered vertically
      text = element_text(),                                       # Regular text
      
      # Legend and Caption
      legend.position = "right",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0)      # Caption color and position
    )
}

# Custom theme for scatter plots
custom_theme_scatter <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      
      # Axis
      axis.line = element_line(color = "black"),                   # Color of axis lines
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(color = "black", size = 8),         # Size and color of axis text
      
      # Legend and Caption
      legend.position = "none",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0.5),    # Caption color and position
      
      # Remove axis ticks
      axis.ticks = element_blank()                                 # No axis ticks
    )
}

# Custom violin plot theme
custom_theme_violin <- function() {
  theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend
      panel.background = element_rect(fill = "white"),  # White background for panels
      panel.grid.major.y = element_line(colour = "grey98"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),  # Color of axis lines
      axis.title = element_text(color = "black"),  # Color of axis titles
      axis.text = element_text(color = "black"),  # Color of axis text
      axis.ticks = element_line(color = "black"),  # Color of axis ticks
      text = element_text(color = "black"),  # Color of text
      plot.title = element_text(face = "bold", hjust = 0.5),  # Bold plot title, centered
      plot.caption = element_text(hjust = 0)  # Caption position
    )
}



# Define the custom theme for UMAP plots
umap_theme <- function() {
  theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),    # Background color
      axis.line = element_line(color = "black", linewidth = 0.5),          # Customize axis line
      axis.text = element_text(size = 10),                            # Customize axis text size
      axis.title = element_text(size = 10),                            # Customize axis title size
      axis.ticks = element_line(color = "black", linewidth = 0.5),         # Customize axis ticks
      panel.grid.major = element_blank(),                             # Remove major gridlines
      panel.grid.minor = element_blank(),                             # Remove minor gridlines
      plot.title = element_text( face = "bold",size = 10, hjust = 0.5),              # Title font size and alignment
      plot.caption = element_text(size = 10, hjust = 0.5, face = "bold"),            # Caption font size and alignment
      plot.margin = margin(10, 10, 10, 10),                           # Set plot margins
      plot.background = element_rect(fill = "white", color = NA),     # Background color
      plot.title.position = "plot",                                  # Title position
      plot.caption.position = "plot",                                # Caption position
      plot.tag.position = "plot",                                    # Tag position
      legend.position = "right",                                     # Move legend to the right
      legend.box.margin = margin(5, 5),                               # Adjust margin around legend box
      legend.spacing.y = unit(0.1, "cm"),                             # Adjust vertical spacing between legend items
      legend.text = element_text(size = 8),                           # Adjust legend text size
      legend.title = element_text(size = 10),                         # Adjust legend title size
      axis.title.y = element_text(angle = 90),          # Rotate y-axis label to be bottom-to-top
    )
}


# Custom feature plot theme
theme_feature_plot <- function(base_size = 12, base_family = "") {
  theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),    # Background color
      axis.line = element_line(color = "black", linewidth = 0.5),          # Customize axis line
      axis.text = element_text(size = rel(0.8)),                            # Customize axis text size
      axis.title = element_text(size = rel(0.9)),                            # Customize axis title size
      axis.ticks = element_line(color = "black", linewidth = 0.5),         # Customize axis ticks
      panel.grid.major = element_blank(),                             # Remove major gridlines
      panel.grid.minor = element_blank(),                             # Remove minor gridlines
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0.5),              # Title font size and alignment
      plot.caption = element_text(size = rel(0.8), hjust = 0.5, face = "bold"),            # Caption font size and alignment
      plot.margin = margin(10, 10, 10, 10),                           # Set plot margins
      plot.background = element_rect(fill = "white", color = NA),     # Background color
      plot.title.position = "plot",                                  # Title position
      plot.caption.position = "plot",                                # Caption position
      plot.tag.position = "plot",                                    # Tag position
      legend.position = "right",                                     # Move legend to the right
      legend.box.margin = margin(5, 5),                               # Adjust margin around legend box
      legend.spacing.y = unit(0.1, "cm"),                             # Adjust vertical spacing between legend items
      legend.text = element_text(size = rel(0.8)),                           # Adjust legend text size
      legend.title = element_text(size = rel(0.9)),                         # Adjust legend title size
      axis.title.y = element_text(angle = 90),          # Rotate y-axis label to be bottom-to-top
    )
}


# Custom dotplot theme
custom_dotplot_theme <- function() {
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(face = "italic", size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 8),  
        legend.title = element_text(size = 8)  
  )
}