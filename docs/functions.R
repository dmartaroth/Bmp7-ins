# ## ######################################## ## #
#                      FUNCTIONS                 #
# ## ######################################## ## #

# Date: Fri May 31 13:25:49 2024 ------------------
# Updated by: Daniela M. Roth

# Custom function to format labels in scientific notation if digits are more than 4
format_labels <- function(x) {
  ifelse(abs(x) >= 10000, format(x, scientific = TRUE), format(x))
}



# prepro.plots ------------------------------------------------------------


# Pre-processing of RNA sequencing data using Seurat v5
# Based on PopsicleR package (now incompatible)


prepro.plots <- function(data, output_dir) {
  plot_number <- 1  # Starting plot number
  
  # Check if 'data' is a Seurat object
  if (!inherits(data, "Seurat")) {
    stop("Input 'data' must be a Seurat object.")
  }
  
  # Check if default assay is set
  if (!"RNA" %in% DefaultAssay(data)) {
    stop("Default assay 'RNA' not found in the Seurat object.")
  }
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create subdirectory "01_Preprocessing"
  preprocessing_dir <- file.path(output_dir, "01_Preprocessing")
  if (!dir.exists(preprocessing_dir)) {
    dir.create(preprocessing_dir)
  }
  
  # Plot Cell Counts
  plot_cell_counts <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    cat(note1("Cell counts represent the number of unique cellular barcodes detected.\n"))
    
    plot <- ggplot(data, aes(x = orig.ident, fill = orig.ident)) +
      geom_bar(colour = "black", position = "dodge") +
      labs(title = "Cell Counts") +
      custom_theme_bar() +
      scale_fill_manual(values = my_colors)
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Cell_Counts.png", plot_number)), width = 3.5, height = 4.5, plot)
    
    # Increment plot number
    plot_number <<- plot_number + 1
  }
  
  # Plot UMI Counts
  plot_umi_counts <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    cat(note2("UMI counts per cell should be above 500...\n"))
    
    plot <- ggplot(data, aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "UMI Counts")
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_UMI_Counts.png", plot_number)), width = 5, height = 4, plot)
    
    # Increment plot number
    plot_number <<- plot_number + 1
  }
  
  # Plot Genes Detected per Cell
  plot_genes_detected_per_cell <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    p1 <- ggplot(data, aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "UMI Counts", fill = "orig.ident", color = "orig.ident")
    
    p2 <- ggplot(data, aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
      geom_boxplot() + 
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ggtitle("NCells vs NGenes")
    
    combined_plot <- plot_grid(p1, p2, rel_widths = c(2, 1.5), scale = c(1, 0.8), labels = "AUTO")
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Genes_Detected_Per_Cell.png", plot_number)), width = 7, height = 3, combined_plot)
    
    # Increment plot number
    plot_number <<- plot_number + 1
  }
  
  # Function to generate QC violin plots 
  generate_QC_violin_plots <- function(output_dir, data, plot_number) {
    cat((green("\nPlotting QC Violin plots \n")))
    
    mt_content <- PercentageFeatureSet(data, pattern = "^mt-")
    mt_ratio <- mt_content / 100
    
    data$percent.mt <- mt_ratio
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    plot_data <- as.data.frame(data@meta.data)
    
    generate_violin_plot <- function(data, feature, title) {
      violin_plot <- ggplot(data, aes(x = orig.ident, y = !!sym(feature), fill = orig.ident)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.4, aes(color = orig.ident)) +
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = my_colors) +
        scale_y_continuous(labels = format_labels) +
        custom_theme_violin() +
        ggtitle(title)
      
      density <- density(data[[feature]])
      max_density <- max(density$y)
      density_threshold <- 0.3 * max_density
      density_peaks <- which(density$y >= density_threshold)
      if (length(density_peaks) > 1) {
        peak_values <- density$x[density_peaks]
        upper_border <- max(peak_values)
      } else {
        upper_border <- density$x[density_peaks]
      }
      buffer <- 0.15 * (max(data[[feature]]) - min(data[[feature]]))
      zoomed_in_range <- c(0, upper_border + buffer)
      
      darker_colors <- alpha(my_colors, 0.8)
      
      zoomed_in_plot <- ggplot(data, aes(x = orig.ident, y = !!sym(feature), fill = orig.ident)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.6, aes(color = orig.ident)) +
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = darker_colors) +
        scale_y_continuous(labels = format_labels, limits = zoomed_in_range) +
        custom_theme_violin() +
        ggtitle(paste("zoom", title))
      
      rect_plot <- geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = min(zoomed_in_range), ymax = max(zoomed_in_range)), 
                             color = "black", fill = NA, linetype = "dashed", size = 0.1, 
                             inherit.aes = FALSE, 
                             lineend = "round", 
                             lwd = 0.4, lty = 2, 
                             dash = c(3, 5))
      
      violin_plot <- violin_plot + rect_plot
      
      return(list(violin_plot, zoomed_in_plot))
    }
    
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    
    pdf(file.path(output_dir, sprintf("%02d_b_QC_violin_plots.pdf", plot_number)), width=6, height=3, useDingbats=FALSE)
    for (i in 1:length(features)) {
      plots <- generate_violin_plot(plot_data, features[i], paste(features[i]))
      grid.arrange(plots[[1]], plots[[2]], nrow = 1, ncol = 2)
    }
    invisible(dev.off())
    
    # Increment plot number
    plot_number <<- plot_number + 1
  }
  
  # Call the plotting functions
  plot_cell_counts(data, plot_number)
  plot_umi_counts(data, plot_number)
  plot_genes_detected_per_cell(data, plot_number)
  generate_QC_violin_plots(preprocessing_dir, data, plot_number)
  
  # Print message
  cat("Plots saved in", preprocessing_dir, "\n")
}



# This function calculates the percentage of mitochondrial gene expression for each cell in a Seurat object and adds it to the metadata. It provides a convenient way to include mitochondrial content information in the Seurat object, which is useful for quality control and downstream analysis. The function takes a Seurat object as input and optionally allows specifying a pattern for mitochondrial genes. After execution, the `percent.mito` column is added to the metadata of the Seurat object, facilitating further analysis and visualization of mitochondrial content.
add_percent_mito <- function(seurat_object, pattern = "^mt-") {
  # Calculate mitochondrial content for the Seurat object
  mt_content <- PercentageFeatureSet(seurat_object, pattern = pattern)
  percent_mito <- mt_content / 100
  
  # Add percent.mito to the metadata of the Seurat object
  seurat_object$percent.mito <- percent_mito
  
  return(seurat_object)
}


# Gene-level filtering ----------------------------------------------------

# The presence of genes with zero counts can dramatically reduce avg expression for a cell; may remove from data with these transformations:
#   
# 1.  Remove genes with zero expression in **all** cells
# 2.  Perform filtering by prevalence, in which genes expressed in only a handful of cells are not particularly meaningful and still bring down averages for all cells not expressed in
# 
# In other words, keep only genes which are expressed in 10 or more cells. This should be adjusted depending on data.

genelvlfilt <- function(object) {
  
  # Output a logical vector for every gene on whether the more than zero counts per cell
  counts <- GetAssayData(object = object, assay = "RNA", layer = "counts") # extract counts
  nonzero <- counts > 0
  
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  
  # Only keeping those genes expressed in more than 10 cells
  filtered_counts <- counts[keep_genes, ]
  
  # Create a new Seurat object with filtered counts
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = object@meta.data)
  
  return(filtered_seurat)
}






# QC plots after filtering ------------------------------------------------------------


filt.plots <- function(data, output_dir) {
  plot_number <- 0  # Starting plot number
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create subdirectory "01_Preprocessing"
  preprocessing_dir <- file.path(output_dir, "02_QC-after-filtering")
  if (!dir.exists(preprocessing_dir)) {
    dir.create(preprocessing_dir)
  }
  
  # Plot Cell Counts
  plot_cell_counts <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    cat(note1("Cell counts represents the number of unique cellular barcodes detected.\n"))
    
    plot <- bind_rows(data = data, .id = "orig.ident") %>%
      ggplot(aes(x = orig.ident, fill = orig.ident)) +
      geom_bar(colour = "black", position = "dodge") +
      labs(title = "Cell Counts")+
      custom_theme_bar() +
      facet_wrap(~ orig.ident, scales = "free_x") +
      scale_fill_manual(values = my_colors)
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Cell_Counts_filtered.png", plot_number)), width = 3.5, height = 4.5, plot)
  }
  
  # Plot UMI Counts
  plot_umi_counts <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    cat(note2("UMI counts per cell should be above 500. If between 500 and 1000 counts, the cells may have benefited from deeper sequencing. The x intercept of this plot is set at 500.\n"))
    
    plot <- bind_rows(data = data, .id = "orig.ident") %>%
      ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "UMI Counts", fill = "orig.ident", color = "orig.ident")
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_UMI_Counts_filtered.png", plot_number)), width = 5, height = 4, plot)
  }
  
  
  # Plot Genes Detected per Cell
  plot_genes_detected_per_cell <- function(data, plot_number) {
    data <- as.data.frame(data@meta.data)
    
    p1 <- bind_rows(data = data, .id = "orig.ident") %>%
      ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "Filtered UMI Counts", fill = "orig.ident", color = "orig.ident")
    
    p2 <- bind_rows(data = data, .id = "orig.ident") %>% 
      ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
      geom_boxplot() + 
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ggtitle("Filtered NCells vs NGenes")
    
    combined_plot <- plot_grid(p1, p2, rel_widths = c(2, 1.5), scale = c(1, 0.8), labels = "AUTO")
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Genes_Detected_Per_Cell_Filtered.png", plot_number)),width = 7, height = 3, combined_plot)
  }
  
  # Plot Mitochondrial Content
  plot_mitochondrial_content <- function(data, plot_number) {
    
    # Calculate mitochondrial content for data
    mt_content_1 <- PercentageFeatureSet(data, pattern = "^mt-")
    mt_ratio_1 <- mt_content_1 / 100
    
    
    # Add mitochondrial content to data frames
    data$percent.mt <- mt_ratio_1
    
    # Add mitoRatio to meta.data
    data$mitoRatio <- mt_ratio_1
    
    # Plot mitochondrial content
    mitoratio_plot <- bind_rows(data = data@meta.data, .id = "orig.ident") %>%
      ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point(aes(color=mitoRatio)) +
      scale_colour_gradient(low = "gray90", high = "black") +
      geom_smooth(method=lm) +   
      scale_x_log10() +
      scale_y_log10() +
      theme_classic() +
      geom_vline(xintercept = 250) +
      geom_hline(yintercept = 250) +
      facet_wrap(~orig.ident)
    
    # Save plot
    ggsave(
      filename = file.path(preprocessing_dir, sprintf("%02d_Mitochondrial_Content_Filtered.png", plot_number)),
      mitoratio_plot,
      width = 6,
      height = 3
    )
  }
  
  # Function to generate QC violin plots 
  generate_QC_violin_plots <- function(output_dir, data, plot_number) {
    cat((green("\nPlotting QC Violin plots \n")))
    
    # Calculate mitochondrial content for data
    mt_content_1 <- PercentageFeatureSet(data, pattern = "^mt-")
    mt_ratio_1 <- mt_content_1 / 100
    
    
    # Add mitochondrial content to data frames
    data$percent.mt <- mt_ratio_1
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # Combine data frames for plotting
    plot_data <- bind_rows(
      data = data@meta.data,
      .id = "orig.ident"
    )
    
    # Function to generate violin plots with dynamically focused zoom
    generate_violin_plot <- function(data, feature, title) {
      # Main violin plot
      violin_plot <- ggplot(data, aes(x = orig.ident, y = !!sym(feature), fill = orig.ident)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +  # Adjust width of violins and position them closer
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.4, aes(color = orig.ident)) +  # Add dots with jitter, adjust transparency, and set color aesthetic
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = my_colors) +  
        scale_y_continuous(labels = format_labels) +  # Apply custom label formatting for y-axis
        custom_theme_violin() +
        ggtitle(title)
      
      # Dynamically determine zoom range based on data distribution
      density <- density(data[[feature]])
      max_density <- max(density$y)
      density_threshold <- 0.3 * max_density  # 50% of maximum density
      density_peaks <- which(density$y >= density_threshold)  # Find regions with density above threshold
      if (length(density_peaks) > 1) {
        peak_values <- density$x[density_peaks]
        upper_border <- max(peak_values)  # Focus on the highest peak
      } else {
        upper_border <- density$x[density_peaks]
      }
      buffer <- 0.15 * (max(data[[feature]]) - min(data[[feature]]))  # 15% buffer
      zoomed_in_range <- c(0, upper_border + buffer)  # Focus on upper border with a buffer
      
      # Darker colors for points
      darker_colors <- alpha(my_colors, 0.8)  # Adjust alpha for darker points
      
      # Zoomed-in violin plot
      zoomed_in_plot <- ggplot(data, aes(x = orig.ident, y = !!sym(feature), fill = orig.ident)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.6, aes(color = orig.ident)) +  # Adjust point size and alpha for darker points
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = darker_colors) +  # Use darker colors for points
        scale_y_continuous(labels = format_labels, limits = zoomed_in_range) +  # Adjust y-axis limits
        custom_theme_violin() +
        ggtitle(paste("zoom", title))
      
      # Red rounded rectangular outline to indicate zoom region on the original plot
      rect_plot <- geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = min(zoomed_in_range), ymax = max(zoomed_in_range)), 
                             color = "black", fill = NA, linetype = "dashed", size = 0.1, 
                             inherit.aes = FALSE, 
                             lineend = "round", 
                             lwd = 0.4, lty = 2, 
                             dash = c(3, 5))
      
      # Add the rectangle to the main violin plot
      violin_plot <- violin_plot + rect_plot
      
      return(list(violin_plot, zoomed_in_plot))
    }
    
    # Features to plot
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    
    # Create and save violin plots for each feature arranged in a grid
    pdf(file.path(output_dir, sprintf("%02d_b_QC_violin_plots_filtered.pdf", plot_number)), width=6, height=3, useDingbats=FALSE)  # Adjust width and height as needed
    for (i in 1:length(features)) {
      plots <- generate_violin_plot(plot_data, features[i], paste(features[i]))
      grid.arrange(plots[[1]], plots[[2]], nrow = 1, ncol = 2)
    }
    invisible(dev.off())  # Close the PDF device
  }
  
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_cell_counts(data, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_umi_counts(data, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_genes_detected_per_cell(data, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_mitochondrial_content(data, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  suppressWarnings(generate_QC_violin_plots(preprocessing_dir, data, plot_number))
  
  
  # Print message
  cat("Plots saved in", preprocessing_dir, "\n")
}



