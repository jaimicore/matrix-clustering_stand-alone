##############################################################
## List of R packages requires in the 'minimal_output' mode ##
##############################################################
.libPaths( c( "/opt/software" , .libPaths() ) )

## CRAN packages:
required_packages_cran = c(
  # "devtools",       ## To install specific package versions
  "dplyr",          ## Data manipulation
  "data.table",     ## Read long matrices in a quick way
  "furrr",          ## Run functions in parallel
  "optparse",       ## Read command-line arguments
  "purrr",          ## Iterations
  "rcartocolor",    ## Cluster colors
  "reshape2",       ## Dataframe manipulation
  "this.path",      ## Create relative paths
  "tidyr",          ## Data manipulation
  "dendsort",       ## To draw heatmap
  "RColorBrewer",   ## Heatmap cell colors
  "ape",            ## Export hclust tree in newick format
  "RJSONIO",        ## Export hclust tree in JSON format
  "circlize",       ## Required to draw heatmaps
  "ggplot2",        ## Plotting
  "ggseqlogo",      ## Plotting
  "flexclust",      ## Calculate adjusted rand index
  "BiocManager",    ## Required for downloading Bio packages
  "jsonlite",
  "svglite")

message("; Installing these R packages from CRAN repository: ", required_packages_cran)
install.packages(required_packages_cran, repos="https://cran.uib.no/", lib="/opt/software")

## Intallation of specific package versions
# library(devtools)
# install_version("this.path", version = "1.2.0", repos = "https://cran.uib.no/")

## Bioconductor packages:
required_packages_bioconductor <- c(
  "universalmotif",
  "ComplexHeatmap")

message("; Installing these R Bioconductor packages: ", required_packages_bioconductor)
BiocManager::install(required_packages_bioconductor, lib="/opt/software")