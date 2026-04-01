#!/usr/bin/env Rscript

# Script to plot posterior distribution of trees using densiTree-style visualization
# with log-scaled node heights and coloring by sampling iteration

# Set working directory to the script's location
tryCatch({
  # Try to get script location from commandArgs first
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    script_dir <- dirname(script_path)
    setwd(script_dir)
    cat(sprintf("Working directory set to: %s\n", getwd()))
  }
}, error = function(e) {
  cat(sprintf("Using current directory: %s\n", getwd()))
})


# Load required libraries
library(ape)
library(phangorn)

# Function to log-transform node heights with correction for most recent sample
log_transform_tree_heights <- function(tree) {
  # Get node heights (distance from root to each node)
  node_heights <- node.depth.edgelength(tree)
  
  # Find the maximum height (tips/most recent samples)
  max_height <- max(node_heights)
  
  # Convert to ages (distance from present, with tips at 0)
  node_ages <- max_height - node_heights
  
  # Apply log transformation to ages (adding 1 to avoid log(0) for tips)
  log_ages <- log(node_ages + 1)
  
  # Convert back to heights (from root)
  max_log_age <- max(log_ages)
  log_heights <- max_log_age - log_ages
  
  # Update tree edge lengths based on log-transformed heights
  new_tree <- tree
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    new_tree$edge.length[i] <- abs(log_heights[parent] - log_heights[child])
  }
  
  return(new_tree)
}

# Custom densiTree function that allows coloring by iteration
plot_densitree_colored <- function(trees, colors, alpha = 0.1, ...) {
  # Get the first tree for structure
  tree1 <- trees[[1]]
  n_tips <- length(tree1$tip.label)
  
  # Plot the first tree to set up the plot area
  plot(tree1, show.tip.label = TRUE, edge.color = NA, ...)
  
  # Get plotting coordinates for all trees and overlay them
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    
    # Get coordinates
    coords <- plotPhyloCoor(tree, ...)
    
    # Draw edges
    for (j in 1:nrow(tree$edge)) {
      parent <- tree$edge[j, 1]
      child <- tree$edge[j, 2]
      
      x0 <- coords[parent, 1]
      y0 <- coords[parent, 2]
      x1 <- coords[child, 1]
      y1 <- coords[child, 2]
      
      # Draw the edge with specified color and transparency
      segments(x0, y0, x1, y1, 
               col = adjustcolor(colors[i], alpha.f = alpha), 
               lwd = 0.5)
    }
  }
}

# Read the trees file
cat("Reading trees file...\n")
trees_file <- "sub_rotavirus_genotype_25_rep5-VP1_rotavirusA_subsample_aa_1.trees"

# Create a temporary copy with END; added if missing
temp_file <- tempfile(fileext = ".trees")
file_content <- readLines(trees_file)

# Check if the last non-empty line is "End;" or "END;"
last_lines <- tail(file_content[nchar(trimws(file_content)) > 0], 5)
if (!any(grepl("^End;?$", last_lines, ignore.case = TRUE))) {
  cat("Adding END; to temporary file...\n")
  file_content <- c(file_content, "End;")
}

writeLines(file_content, temp_file)

# Read from the temporary file
trees <- read.nexus(temp_file)

# Clean up temporary file
unlink(temp_file)

# Check if trees is a list or single tree
if (class(trees) == "phylo") {
  trees <- list(trees)
}

n_trees <- length(trees)
cat(sprintf("Read %d trees\n", n_trees))

# Extract sampling iterations (assumes sequential sampling)
iterations <- seq_along(trees)

# Create color palette based on iterations
# Use a gradient from blue (early) to red (late)
color_palette <- colorRampPalette(c("blue", "red"))
colors <- color_palette(n_trees)

# Log-transform all trees
cat("Log-transforming tree heights...\n")
log_trees <- lapply(trees, log_transform_tree_heights)
class(log_trees) <- "multiPhylo"

# Create two plots side by side
cat("Creating densiTree plots...\n")

# Set up 2-panel plot
par(mfrow = c(1, 2), mar = c(2, 1, 4, 1))

# Plot 1: Original scale
densiTree(trees, 
          type = "cladogram",
          alpha = 0.02,
          consensus = NULL,
          scaleX = FALSE,
          col = colors,
          width = 1,
          cex = 0.001,
          underscore = FALSE)

title(main = "Original Scale",
      cex.main = 1.2)

# Add color legend
legend_colors <- color_palette(5)
legend_labels <- sprintf("Iteration %d", 
                        round(seq(1, n_trees, length.out = 5)))
legend("topleft", 
       legend = legend_labels,
       col = legend_colors,
       lwd = 2,
       bty = "n",
       cex = 0.8,
       title = "Sampling Iteration")

# Plot 2: Log-scaled
densiTree(log_trees, 
          type = "cladogram",
          alpha = 0.02,
          consensus = NULL,
          scaleX = FALSE,
          col = colors,
          width = 1,
          cex = 0.001,
          underscore = FALSE)

title(main = "Log-scaled Heights",
      cex.main = 1.2)

legend("topleft", 
       legend = legend_labels,
       col = legend_colors,
       lwd = 2,
       bty = "n",
       cex = 0.8,
       title = "Sampling Iteration")

cat("Done!\n")

