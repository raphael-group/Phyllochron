#if (!require("BiocManager", quietly = TRUE))
      #install.packages("BiocManager", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")

library(BiocManager, lib='../plotting/R/')

#BiocManager::install("timescape", lib="../../plotting/R/")
library(timescape, lib='../plotting/R/')
#install.packages("stringr", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(stringr, lib='../plotting/R/')
#install.packages("vctrs", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(vctrs, lib='../plotting/R/')
#install.packages("jsonlite", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(jsonlite, lib='../plotting/R/')
#install.packages("htmlwidgets", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(htmlwidgets, lib='../plotting/R/')
#install.packages("yaml", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(yaml, lib='../plotting/R/')
#install.packages("rmarkdown", lib='../../plotting/R/', repos = "http://cran.us.r-project.org")
library(rmarkdown, lib='../plotting/R/')

# Load required libraries

# Function to generate timescape plot
generate_timescape_plot <- function(s) {
  # Construct the file names using the input string s
  tree_edges_file <- paste0(s, "_tree_edges.csv")
  clonal_prev_file <- paste0(s, "_clone_prev.csv")
  timescape_file <- paste0(s, "_timescape_plot.html")
  timescape_png <- paste0(s, "_timescape_plot.png")

  # Read CSV files for tree edges and clonal prevalence
  tree_edges <- read.csv(tree_edges_file)
  clonal_prev <- read.csv(clonal_prev_file)
  
  unique_clones <- unique(c(tree_edges$parent, tree_edges$child))
  clone_ids <- as.character(unique_clones)
  num_clones <- length(clone_ids)
  color_palette <- rainbow(num_clones)
  
  # Define clone colours
  clone_colours <- data.frame(
    clone_id = clone_ids, 
    colour = color_palette  
  )

  # Generate timescape plot
  p <- timescape(
    clonal_prev = clonal_prev, 
    tree_edges = tree_edges, 
    clone_colours = clone_colours
  )
  htmlwidgets::saveWidget(p, timescape_file)
}

main <- function() {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if argument is provided
  #if (length(args) == 0) {
    #stop("Please provide the input string (file prefix) as a command line argument.")
  #}
  
  # Call the function with the input string from command line
  s <- "../data/AML/Phyllochron/AML-63"
  generate_timescape_plot(s)
}

main()
