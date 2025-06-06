##############################################################
## List of R packages requires in the 'minimal_output' mode ##
##############################################################

## CRAN packages:
required_packages_cran = c(
  "BiocManager",    # To install bioconductor packages
  "curl",           # Download files
  "dplyr",          # Data manipulation
  "data.table",     # Data manipulation
  "furrr",          # Run functions in parallel
  "optparse",       # Read command-line arguments
  "purrr",          # Iterations
  "rcartocolor",    # Cluster colors
  "reshape2",       # Dataframe manipulation
  "this.path",      # Create relative paths
  "tidyr",          # Data manipulation
  "dendsort",       # To draw heatmap
  "ggplot2",        ## For plotting
  "ggseqlogo",      # Draw logos
  "RColorBrewer",   # Heatmap cell colors
  "ape",            # Export hclust tree in newick format
  "RJSONIO",        # Export hclust tree in JSON format
  "circlize",       # Required to draw heatmaps
  "flexclust",      # Calculate adjusted rand index
  "htmlwidgets",    # Save plotly output as html
  "plotly",         # Interactive plots
  "jsonlite")       # To create the JSON file from the hclust outputs

message("; Installing these R packages from CRAN repository: ", required_packages_cran)
install.packages(
  required_packages_cran,
  dependencies = TRUE)

## Bioconductor packages:
required_packages_bioconductor <- c(
  "universalmotif",
  "ComplexHeatmap")

message("; Installing these R Bioconductor packages: ", required_packages_bioconductor)
BiocManager::install(
  required_packages_bioconductor)
