# ----------------------- #
# Load required libraries #
# ----------------------- #

suppressPackageStartupMessages(library(dplyr))          # Data manipulation
suppressPackageStartupMessages(library(data.table))     # Read long matrices in a quick way
suppressPackageStartupMessages(library(furrr))          # Run functions in parallel
suppressPackageStartupMessages(library(ggplot2))        # Data visualization
suppressPackageStartupMessages(library(ggseqlogo))      # Sequence logo visualization
suppressPackageStartupMessages(library(htmlwidgets))    # Export plotly as HTML
suppressPackageStartupMessages(library(optparse))       # Read command-line arguments
suppressPackageStartupMessages(library(purrr))          # Iterations
suppressPackageStartupMessages(library(plotly))         # Interactive plots 
suppressPackageStartupMessages(library(reshape2))       # Dataframe operations
suppressPackageStartupMessages(library(rcartocolor))    # Nice color palettes
suppressPackageStartupMessages(library(svglite))        # SVG graphics
suppressPackageStartupMessages(library(this.path))      # Create relative paths
suppressPackageStartupMessages(library(tidyr))          # Data manipulation
suppressPackageStartupMessages(library(universalmotif)) # Motif analysis

# -------------- #
# Read arguments #
# -------------- #
option_list = list(
  
  make_option(c("-i", "--matrix_file_table"), type = "character", default = NULL, 
              help = "A text-delimited file where each line contain the next fields. 1: Motif file path; 2: Motif collection name; 3: Motif format (Mandatory). It does not expect a header, but it expects those columns in the indicated order.", metavar = "character"),
  
  make_option(c("-o", "--output_folder"), type = "character", default = NULL, 
              help = "Folder to save the results (Mandatory)", metavar = "character"),
  
  make_option(c("-m", "--comparison_metric"), type = "character", default = "Ncor", 
              help = "Comparison metric used to build the hierarchical tree. [Default: \"%default\" . Options: cor, Ncor]", metavar = "character"),  
  
  make_option(c("-l", "--linkage_method"), type = "character", default = "average", 
              help = "Linkage/agglomeration method to build the hierarchical tree. [Default: \"%default\" . Options: average, complete, single]", metavar = "character"),
  
  make_option(c("-w", "--number_of_workers"), type = "numeric", default = 2, 
              help = "number of cores to run in parallel. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("--export_newick"), type = "logical", default = FALSE, 
              help = "Export hierarchical tree in Newick format. [Default \"%default\"] ", metavar = "logical"), 
  
  make_option(c("--export_heatmap"), type = "logical", default = FALSE, 
              help = "Export heatmap with clusters in PDF. [Default \"%default\"] ", metavar = "logical"), 

  make_option(c("--radial_tree"), type = "logical", default = FALSE, 
              help = "Radial representation of all motifs aligned. In this alignment the thresholds are ignored, however the clusters are detected and highlighted in the radial tree [Default \"%default\"] ", metavar = "logical"), 
  
  make_option(c("--heatmap_color_palette"), type = "character", default = "RdGy", 
              help = "Cell colors in clusters heatmap. [Default: \"%default\"]. Options: any colorBrewer palette, see colorbrewer2.org for details", metavar = "character"),
  
  make_option(c("--color_palette_classes"), type = "numeric", default = 11, 
              help = "Number of classes to create color palette in clusters heatmap. [Default: \"%default\"]. Options: depends on the selected colorBrewer palette, see colorbrewer2.org for details", metavar = "number"),
  
  make_option(c("-c", "--cor_th"), type = "numeric", default = 0.75, 
              help = "Pearson correlation threshold. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-n", "--Ncor_th"), type = "numeric", default = 0.55, 
              help = "Normalized Pearson correlation threshold. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("-W", "--w_th"), type = "numeric", default = 5,
              help = "Minimum number of aligned positions between a pair of motifs. [Default \"%default\"] ", metavar = "number"),

  make_option(c("-M", "--minimal_output"), type = "logical", default = FALSE, 
              help = "When TRUE only returns the alignment, clusters and motif description tables. Comparison results, plots and trees are not exported. [Default \"%default\"] ", metavar = "logical"),

  make_option(c("-a", "--annotation_table"), type = "character", default = NULL,
              help = "Motif annotation table, when this file is provided with the '--radial_tree = TRUE' option add annotations to the radial tree. A tab-delimited file with at least the following columns : 1) motif_id, 2) class ID, 3) collection. If the input motifs are separated in many collections, concatenate all the annotations in a single file.", metavar = "character"),

  make_option(c("--title"), type = "character", default = "matrix-clustering",
              help = "Analysis title", metavar = "character"),
  
  make_option(c("--ARI"), type = "logical", default = FALSE, 
              help = "Calculate the Adjusted Rand Index (ARI) of the resulting clusters based in the provided annotation table (--annotation_table). [Default \"%default\"].", metavar = "logical")
  
);
message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


# Output file prefix
out.folder        <- opt$output_folder
matrix.file.table <- opt$matrix_file_table

# Mandatory input
if (!file.exists(matrix.file.table)) {
  stop("Mandatory input file not found: ", matrix.file.table)
}


# ARI flag
reference.clusters.flag <- opt$ARI

# In case of an annotation file is provided
motif.annotation.flag  <- FALSE
motif.annotation.file <- opt$annotation_table
if (!is.null(motif.annotation.file)) {
  if (!file.exists(motif.annotation.file)) {
    stop("Motif annotation file not found: ", motif.annotation.file)
  }
  motif.annotation.flag   <- TRUE
}

params.list <- list("export_newick"         = as.numeric(opt$export_newick),
                    "export_heatmap"        = as.numeric(opt$export_heatmap),
                    "heatmap_color_classes" = as.numeric(opt$color_palette_classes),
                    "heatmap_color_palette" = opt$heatmap_color_palette,
                    "comparison_metric"     = opt$comparison_metric,
                    "linkage_method"        = opt$linkage_method,
                    "w"                     = opt$w_th,
                    "cor"                   = as.numeric(opt$cor_th),
                    "Ncor"                  = as.numeric(opt$Ncor_th),
                    "nb_workers"            = as.numeric(opt$number_of_workers),
                    "min_output"            = opt$minimal_output,
                    "title"                 = opt$title,
                    "ARI"                   = reference.clusters.flag,
                    "annotation"            = motif.annotation.flag,
                    "radial_tree"           = opt$radial_tree)


# ------------------------- #
# Source here the functions #
# ------------------------- #

# 'here' is set to the path where matrix-clustering.R is located, then the libraries
# are sourced relative to where 'here' was set
source(this.path::here(.. = 0, "R", "General_utils.R"))
source(this.path::here(.. = 0, "R", "Hierarchical_clustering.R"))
source(this.path::here(.. = 0, "R", "Motif_alignment_utils.R"))
source(this.path::here(.. = 0, "R", "Motif_manipulation.R"))
source(this.path::here(.. = 0, "R", "Tree_partition_utils.R"))
source(this.path::here(.. = 0, "R", "D3_trees.R"))

# -------- #
# D3 files #
# -------- #
d3.radial.tree.template        <- this.path::here(.. = 0, "html", "templates", "Radial_tree_template.html")
html.interactive.tree.template <- this.path::here(.. = 0, "html", "templates", "Interactive_tree_template.html")
d3.interactive.tree.template   <- this.path::here(.. = 0, "html", "templates", "Interactive_tree_template.d3")
d3.min.lib                     <- this.path::here(.. = 0, "html", "js", "d3.v3.min.js")
jquery.lib                     <- this.path::here(.. = 0, "html", "js", "jquery.js")


# ----- #
# Debug #
# ----- #
## Example:

# matrix.file.table <- "Debug/Triana/motif_table_test.txt"
# out.folder        <- "results/Triana/Triana"


# Rscript matrix-clustering.R -i data/OCT4_datasets/OCT4_motif_table.txt -o results/OCT4_motifs_example/OCT4_motif_analysis --number_of_workers 8
# matrix.file.table <- "data/OCT4_datasets/OCT4_motif_table.txt"
# out.folder        <- "results/Oct4_debug_example/Oct4_debug_example"

# Rscript matrix-clustering.R -i data/JASPAR_2022/Jaspar_plants_motifs_tab.txt -o results/Jaspar_plants/Jaspar_plants --number_of_workers 8 --reference_cluster_annotation data/JASPAR_2022/Jaspar_2022_plants_TF_fam.tab
# matrix.file.table           <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/data/JASPAR_2022/Jaspar_plants_motifs_tab.txt"
# out.folder                  <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/results/Jaspar_plants/Jaspar_plants"
# reference.clusters.tab.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/data/JASPAR_2022/Jaspar_2022_plants_TF_fam.tab"

# matrix.file.table <- "data/JASPAR_2022/Jaspar_nematodes_motifs_tab.txt"
# out.folder        <- "results/JASPAR_nematodes_radial_debug/JASPAR_neamtodes_radial"

# params.list <- list("export_newick"         = 0,
#                     "export_heatmap"        = 0,
#                     "heatmap_color_classes" = NULL,
#                     "heatmap_color_palette" = NULL,
#                     "comparison_metric"     = "Ncor",
#                     "linkage_method"        = "average",
#                     "w"                     = 5,
#                     "cor"                   = 0.75,
#                     "Ncor"                  = 0.55,
#                     "nb_workers"            = 8,
#                     "min_output"            = FALSE,
#                     "ARI"                   = FALSE,
#                     "radial_tree"           = FALSE) #


# -------------------------------------------------------- #
# Initialize output + result lists + create output folders #
# -------------------------------------------------------- #

# Absolute path to result's main folder
# This path will be used to compute relative path required in the output html files
# 
# Example: this.path::as.rel.path(relative.to = results.main.dir, path = new.d3)
results.main.dir <- dirname(out.folder)

out.folder.list <- list(tables         = file.path(paste0(out.folder, "_tables")),
                        plots          = file.path(paste0(out.folder, "_plots")),
                        trees          = file.path(paste0(out.folder, "_trees")),
                        aligned_logos  = file.path(paste0(out.folder, "_aligned_logos")),
                        motifs         = file.path(paste0(out.folder, "_motifs")),
                        central_motifs = file.path(paste0(out.folder, "_motifs"), "central_motifs"),
                        root_motifs    = file.path(paste0(out.folder, "_motifs"), "root_motifs"),
                        indiv_motifs   = file.path(paste0(out.folder, "_motifs"), "individual_motifs_with_gaps"),
                        cluster_motifs = file.path(paste0(out.folder, "_motifs"), "motifs_sep_by_cluster"))

output.files.list <- list("Alignment_table"       = file.path(out.folder.list$tables, "alignment_table.tab"),
                          "Clusters_table"        = file.path(out.folder.list$tables, "clusters.tab"),
                          "Distance_table"        = file.path(out.folder.list$tables, "distance_table.tab"),
                          "Motif_compa"           = file.path(out.folder.list$tables, "pairwise_motif_comparison.tab"),
                          "Cluster_colors"        = file.path(out.folder.list$tables, "cluster_color_map.tab"),
                          "Summary_table"         = file.path(out.folder.list$tables, "summary_table.tab"),
                          "Motifs_transfac_tmp"   = file.path(out.folder.list$motifs, "input_motifs_parsed_id.tf.tmp"),
                          "Motifs_transfac"       = file.path(out.folder.list$motifs, "input_motifs_parsed_id.tf"),
                          "Root_motifs"           = file.path(out.folder.list$root_motifs, "Root_motifs.tf"),
                          "Heatmap_clusters"      = file.path(out.folder.list$plots, "Heatmap_clusters.pdf"),
                          "Clusters_vs_Reference" = file.path(out.folder.list$plots, "Clusters_vs_reference_contingency_table.pdf"),
                          "hclust_all"            = file.path(out.folder.list$trees, "tree.RData"),
                          "JSON_tree_all"         = file.path(out.folder.list$trees, "tree.json"),
                          "Newick_tree_all"       = file.path(out.folder.list$trees, "tree.newick"),
                          "JSON_radial_annotated" = file.path(out.folder.list$trees, "Annotated_tree_cluster_01.json"),
                          "D3_radial_tree"        = paste0(out.folder, "_D3_radial_tree.html"),
                          "D3_dynamic_tree"       = paste0(out.folder, "_D3_dynamic_tree.html"))

results.list <- list(Dist_table         = NULL,
                     Dist_matrix        = NULL,
                     Original_matrix    = NULL,
                     All_motifs_tree    = NULL,
                     Alignment_table    = NULL,
                     Cluster_color      = NULL,
                     Motif_info_tab     = NULL,
                     Motif_compa_tab    = NULL,
                     Reference_clusters = NULL,
                     Aligned_motifs_um  = NULL,
                     JSON_branch_nb     = NULL,
                     Clusters_files     = NULL)

## Create output folders
no.output <- sapply(out.folder.list, dir.create, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------- #
# Pre-process motif files + Generate motif description table #
# ---------------------------------------------------------- #

matrix.file.list <- check.status.motif.table(matrix.file.table = matrix.file.table)

## Returns the motif description table and the motifs (as a universalmotif object)
## with the modified IDs (unique identifiers with the motif number and motif collection name)
message("; Generating motif description table")
motif.info.and.motifs <- purrr::pmap(.l = matrix.file.list,
                                     .f = ~preprocess.one.motif.collection(motif.file      = ..1,
                                                                           motif.format    = ..3,
                                                                           collection.name = ..2))

## Concatenate the list into a single data.table
results.list$Motif_info_tab  <- rbindlist(purrr::map(motif.info.and.motifs, `[[`, "motif_info"))
all.motifs.um                <- unlist(purrr::map(motif.info.and.motifs, `[[`, "um_object")) ## All motifs as a universalmotif object


## Export transfac file with correct header to be read by compare-matrices-quick
write.transfac.parsed.header(old.tf.file = output.files.list$Motifs_transfac_tmp,
                             new.tf.file = output.files.list$Motifs_transfac,
                             um.object   = all.motifs.um) 

params.list[["Nb_motifs"]]      <- nrow(results.list$Motif_info_tab)
params.list[["Nb_collections"]] <- length(matrix.file.list$Motif_file)
message("; Analysis with ", params.list[["Nb_motifs"]], " motifs in ", params.list[["Nb_collections"]], " collections")


# --------------------------- #
# Read motif comparison table #
# --------------------------- #
results.list$Motif_compa_tab <- motif.comparison(transfac.file     = output.files.list$Motifs_transfac,
                                                 output.compa.file = output.files.list$Motif_compa)


# --------------------------------------------- #
# Convert distance table into a distance matrix #
# --------------------------------------------- #
message('; Converting correlation values to distances')
distances.objects <- build.distance.matrix(compa.table = results.list$Motif_compa_tab,
                                           metric      = params.list$comparison_metric)

results.list[["Dist_table"]]      <- data.table(as.data.frame.matrix(distances.objects$table))
results.list[["Dist_matrix"]]     <- distances.objects$matrix
results.list[["Original_matrix"]] <- data.table(as.data.frame.matrix(distances.objects$original))


# -------------------------------------- #
# Optional: read reference cluster table #
# -------------------------------------- #
if (params.list$ARI) {
  
  message('; Reading user-provided reference clusters table')
  reference.clusters.tab <- fread(motif.annotation.file, header = T) %>% 
                              dplyr::rename(ID         = motif_id,
                                            cluster    = class,
                                            Collection = collection)
  
  ## Parse the motif ID in the 'Motif Information Table' (remove the numeric suffix) to create Ids with simlar format
  motif.id.parsed <- gsub(x = results.list$Motif_info_tab$id, pattern = "_n\\d+$", replacement = "")
  
  ## Combine Id and collection columns 
  reference.clusters.tab <- reference.clusters.tab %>% 
                              mutate(id = paste0(Collection, "_", ID)) %>% 
                              within(rm(Collection, ID)) %>% 
                              dplyr::filter(id %in% motif.id.parsed) %>% 
                              select(id, cluster)

  results.list$Reference_clusters <- reference.clusters.tab
  
  ## Check that all IDs in the reference cluster table are found within the IDs obtained directly from the motifs
  if (!all(motif.id.parsed %in% reference.clusters.tab$id)) {
    stop("There are discrepancies in the id provided in the reference ID and the IDs obtained from the motif files. In addition, verify that the collection name provided with --matrix_file_table and --reference_cluster_annotation are identical.")
  }
}


# ---------------------------------------- #
# Compute the hierarchical clustering tree #
# ---------------------------------------- #

## Hierarchical clustering can only be applied with > 1 elements
if (params.list[["Nb_motifs"]] > 1) {
  
  results.list$All_motifs_tree <- hclust.motifs(results.list[["Dist_matrix"]],
                                                hclust.method = params.list$linkage_method)

  # --------------------- #
  # Find clusters section #
  # --------------------- #

  message("; Finding clusters")
  find.clusters.list <- find.motif.clusters(tree             = results.list$All_motifs_tree,
                                            comparison.table = results.list$Motif_compa_tab,
                                            parameters       = params.list)
  
  # Create a table with the path to the JSON file of each cluster
  results.list$Clusters_files <- data.table(Cluster     = names(find.clusters.list$clusters),
                                            JSON_folder = file.path(out.folder.list$trees, names(find.clusters.list$clusters)))
  
  results.list$Clusters_files <- results.list$Clusters_files |> 
    mutate(JSON_file           = file.path(JSON_folder, paste0(Cluster, "_tree.json")),
           JSON_annotated_file = file.path(JSON_folder, paste0(Cluster, "_tree_annotated.json")))
  
  purrr::walk(.x = results.list$Clusters_files$JSON_folder,
              .f = ~dir.create(path = .x, showWarnings = FALSE, recursive = TRUE))
 
  
  # When the --radial_tree option is indicated the find.motif.clusters function is launched two times.
  # 1st: Detect the clusters with the user-defined w, cor, and Ncor thresholds
  # 2nd: Set the w, cor, and Ncor thresholds to their minimum value, this will result in all motifs grouped and aligned in a single cluster
  if (params.list[["radial_tree"]]) {
    
    message("; Force the alignment of all motifs in a single cluster (--radial_tree option)")
    params.list.radial <- params.list
    params.list.radial$w    <- 0
    params.list.radial$cor  <- -1  
    params.list.radial$Ncor <- -1  
    
    find.clusters.list.radial <- find.motif.clusters(tree             = results.list$All_motifs_tree,
                                                     comparison.table = results.list$Motif_compa_tab,
                                                     parameters       = params.list.radial)
  }
  
  # ----------------------------- #
  # Calculate Adjusted Rand Index #
  # ----------------------------- #
  if (params.list$ARI) {
    
    message("; Calculating Adjusted Rand Index (ARI) based on the user-provided reference clusters")
    clustering.ari <- calculate.ARI(matrix.clustering.clusters = clusters.list.to.df(find.clusters.list$clusters, id.pttrn.rm = "_n\\d+$"),
                                    reference.clusters         = results.list$Reference_clusters)
    
    params.list[["ARI"]] <- clustering.ari$ARI
    params.list[["RI"]]  <- clustering.ari$RI
  }

  
  # ------------------------------------- #
  # Clusters and alignment in radial tree #
  # ------------------------------------- #
  if (params.list[["radial_tree"]]) {
    
    cluster.sizes                       <- purrr::map_dbl(find.clusters.list.radial$clusters, length)
    cl.singleton                        <- 0
    cl.many                             <- 1
    singleton.flag                      <- FALSE
    params.list[["Nb_clusters"]]        <- 1
    params.list.radial[["Nb_clusters"]] <- 1
    message("; Number of clusters: ", params.list[["Nb_clusters"]])
    
    r.cl.hclust.results <- hclust.cluster.ids(ids        = find.clusters.list.radial$clusters$cluster_1,
                                              compa      = results.list$Motif_compa_tab,
                                              parameters = params.list.radial)

    # ----------------------- #
    # Motif alignment section #
    # ----------------------- #
    
    ## Extract the hclust objects within the nested list
    r.cl.many.hclust <- r.cl.hclust.results$hclust
    
    ## Align motifs within each cluster (with 2 or more motifs, singletons are treated separately)
    ## This function is ran in parallel using furrr::future_map , see https://furrr.futureverse.org/
    plan(multisession, workers = params.list$nb_workers)
    options(future.globals.maxSize=400000000000000000)

    message("; Motif alignment step")
    r.aligment.clusters.tab <- align.motifs.in.cluster(tree       = r.cl.many.hclust,
                                                       compa      = results.list$Motif_compa_tab,
                                                       motif.info = results.list$Motif_info_tab,
                                                       parameters = params.list.radial)
    
    r.aligment.clusters.tab <- list(r.aligment.clusters.tab)

    # ----------------------------- #
    # Parse tables before exporting #
    # ----------------------------- #
    
    ## Prepare tables to export
    message("; Combining aligned clusters tables")
    
    ## Groups
    ## Add cluster name, rename columns and remove unnecessary columns
    r.alignment.clusters.tab.export <-  rbindlist(r.aligment.clusters.tab, idcol = "cluster") %>% 
                                            within(rm(N, Update_status)) %>% 
                                            dplyr::rename(strand            = Strand,
                                                          offset_up         = Offset_up,
                                                          offset_down       = Offset_down,
                                                          aligned_consensus = Oriented_consensus)

    results.list$Alignment_radial_table <- r.alignment.clusters.tab.export %>% 
                                            mutate(width = nchar(aligned_consensus)) 

    ## This is the column order of the original version
    results.list$Alignment_radial_table <- results.list$Alignment_radial_table %>% 
      select("id", "name", "cluster", "strand", "offset_up", "offset_down", "width", "aligned_consensus", "aligned_consensus_rc")

    # cat(results.list$Alignment_radial_table$aligned_consensus, sep = "\n")
    
    results.list$Clusters_table <- results.list$Alignment_radial_table %>% 
                                    group_by(cluster) %>% 
                                    summarise(id      = paste(id, collapse = ","),
                                              name    = paste(name, collapse = ","),
                                              .groups = "drop")    
  }

  
  
  
  # --------------------------------------- #
  # Clusters and alignments in motif forest #
  # --------------------------------------- #
  {
    cluster.sizes                <- purrr::map_dbl(find.clusters.list$clusters, length)
    cl.singleton                 <- names(cluster.sizes[cluster.sizes == 1])
    cl.many                      <- names(cluster.sizes[cluster.sizes > 1])
    singleton.flag               <- as.logical(length(cl.singleton))
    params.list[["Nb_clusters"]] <- length(find.clusters.list$clusters)
    all.singletons.flag          <- ifelse(params.list[["Nb_clusters"]] == params.list[["Nb_motifs"]], yes = 1, no = 0)
    message("; Number of clusters: ", params.list[["Nb_clusters"]])
    
  
    # ----------------------- #
    # Motif alignment section #
    # ----------------------- #
    
    ## Compute hierarchical clustering of each cluster
    cl.hclust.results <- purrr::map(.x = find.clusters.list$clusters,
                                    .f = ~hclust.cluster.ids(ids        = .x,
                                                             compa      = results.list$Motif_compa_tab,
                                                             parameters = params.list))
    ## Extract the hclust objects within the nested list
    cl.many.hclust <- purrr::map(cl.hclust.results[cl.many], `[[`, "hclust")
  
    ## Align motifs within each cluster (with 2 or more motifs, singletons are treated separately)
    ## This function is ran in parallel using furrr::future_map , see https://furrr.futureverse.org/
    if (!all.singletons.flag) {
      
      plan(multisession, workers = params.list$nb_workers)
      options(future.globals.maxSize=400000000000000000)

      message("; Motif alignment step") 
      aligment.clusters.tab <- furrr::future_map(.x = cl.many.hclust,
                                                 .f = ~align.motifs.in.cluster(tree       = .x,
                                                                               compa      = results.list$Motif_compa_tab,
                                                                               motif.info = results.list$Motif_info_tab,
                                                                               parameters = params.list))
    }

    
    ## Export the JSON file of the cluster
    json.trees <- furrr::future_map(.x = cl.many.hclust,
                                    .f = ~convert.hclust.to.JSON(hc = .x))
    message("; Exporting trees as a JSON files")
    
    # Subset to clusters with 2 or more motifs
    cl.many.files.df <- results.list$Clusters_files |>
                          dplyr::filter(Cluster %in% cl.many)
    
    furrr::future_walk2(.x = cl.many.files.df$JSON_file,
                        .y = json.trees,
                        .f = ~writeLines(.y, con = .x))

    # ----------------------------- #
    # Parse tables before exporting #
    # ----------------------------- #
    
    if (!all.singletons.flag) {
      ## Prepare tables to export
      message("; Combining aligned clusters tables")
      
      ## Groups
      ## Add cluster name, rename columns and remove unnecessary columns
      alignment.clusters.tab.export <-  rbindlist(aligment.clusters.tab, idcol = "cluster") %>% 
                                          within(rm(N, Update_status)) %>% 
                                          dplyr::rename(strand            = Strand,
                                                        offset_up         = Offset_up,
                                                        offset_down       = Offset_down,
                                                        aligned_consensus = Oriented_consensus)
    }
    
    
    ## Singleton section
    ## This if only applied when there is at least one singleton
    if (singleton.flag) {
      
      singleton.clusters.tab.export <- data.table(cluster = names(find.clusters.list$clusters[cl.singleton]),
                                                  id      = unlist(find.clusters.list$clusters[cl.singleton])) %>% 
        left_join(results.list$Motif_info_tab, by = "id") %>% 
        mutate(strand               = "D",
               offset_up            = 0,
               offset_down          = 0,
               aligned_consensus    = consensus,
               aligned_consensus_rc = rc_consensus) %>% 
        within(rm(n, width, IC, nb_sites, id_old)) 
      
      
      
      message("; Exporting trees as a JSON files")
      json.trees.singleton <- paste("{\n\"name\": \"\",\n\"children\":[\n{\n \"label\": \"", as.vector(unlist(find.clusters.list$clusters[cl.singleton])), "\"\n}\n]\n}", sep = "")
      
      # Subset to singletons
      cl.singleton.files.df <- results.list$Clusters_files |>
                                dplyr::filter(Cluster %in% cl.singleton)
      
      furrr::future_walk2(.x = cl.singleton.files.df$JSON_file,
                          .y = json.trees.singleton,
                          .f = ~writeLines(.y, con = .x))
      
      
      ## Combine groups + singleton clusters
      if (!all.singletons.flag) {
        results.list$Alignment_table <- rbind(alignment.clusters.tab.export, singleton.clusters.tab.export) %>% 
                                          mutate(width = nchar(aligned_consensus))
      } else {
        results.list$Alignment_table <- singleton.clusters.tab.export %>% 
                                          mutate(width = nchar(aligned_consensus))
      }
      
      
    } else {
      results.list$Alignment_table <- alignment.clusters.tab.export %>% 
                                        mutate(width = nchar(aligned_consensus))
    }
    
    ## This is the column order of the original version
    results.list$Alignment_table <- results.list$Alignment_table %>% 
                                      select("id", "name", "cluster", "strand", "offset_up", "offset_down", "width", "aligned_consensus", "aligned_consensus_rc")

    
    # ------------------- #
    # Clusters - ID table #
    # ------------------- #
    
    # At this point this variable contains 1 cluster when radial tree option is actiaved
    results.list$Clusters_table <- results.list$Alignment_table %>% 
      group_by(cluster) %>% 
      summarise(id      = paste(id, collapse = ","),
                name    = paste(name, collapse = ","),
                .groups = "drop")    
  }
  
} else {
  stop("The clustering analysis required at least 2 motifs.")
}



# ---------------------- #
# Export results section #
# ---------------------- #

# -------------------- #
# Export summary table #
# -------------------- #
summary.clustering.tab <- data.table(Nb_motifs         = params.list$Nb_motifs,
                                     Nb_collections    = params.list$Nb_collections,
                                     Nb_clusters       = params.list$Nb_clusters,
                                     Thresholds        = paste0("Ncor = ", params.list$Ncor, "; cor = ", params.list$cor, "; w = ", params.list$w),
                                     Linkage_method    = params.list$linkage_method,
                                     Similarity_metric = params.list$comparison_metric)

## Add columns when the Rand index is calculated
if (params.list$ARI) {
  
  ri.dt <- data.table(Rand_Index          = params.list$RI,
                      Adjusted_Rand_Index = params.list$ARI)
  
  summary.clustering.tab <- cbind(summary.clustering.tab, ri.dt)
}

message("; Exporting parameters + summary table: ", output.files.list$Summary_table)
fwrite(x         = summary.clustering.tab,
       file      = output.files.list$Summary_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")



# --------------------------------------------------------------------------------------- #
# Export the alignment + clusters table + motif files: these files are the minimal output #
# --------------------------------------------------------------------------------------- #
if (params.list[["radial_tree"]]) {
  align.tab <- results.list$Alignment_radial_table
} else {
  align.tab <- results.list$Alignment_table
}
# cat(align.tab$aligned_consensus, sep = "\n")


message("; Exporting alignment table: ", output.files.list$Alignment_table)
fwrite(x         = align.tab,
       file      = output.files.list$Alignment_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")


message("; Exporting clusters table: ", output.files.list$Clusters_table)
fwrite(x         = results.list$Clusters_table,
       file      = output.files.list$Clusters_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")


# ------------------------------- #
# Export motifs :                 #
# 1) Oriented motifs with gaps    #
# 2) Motifs separated by clusters #
# 3) Root motifs                  #
# 4) Central motifs               #
# ------------------------------- #

# -------------------------------------------------------------------------
# 1) Export motifs (without gaps) as transfac files in D and R orientation
message("; Exporting individual motif files in ", out.folder.list$indiv_motifs)


export.indiv.motif.files(un.motifs       = all.motifs.um,
                         alignment.table = align.tab,
                         outdir          = out.folder.list$indiv_motifs)

## universalmotif does not accepts rows containing only 0s, therefore, the insertion of gaps
## must be done directly to the transfac file, instead of using universalmotif functions
message("; Adding gaps to motif files")

# When the radial_tree option is active, all motifs are aligned in a single cluster
# The number of gaps was calculated in the 'Alignment_radial_table'

add.gaps.tab <- add.gaps.to.indiv.tf.files(motif.folder = out.folder.list$indiv_motifs,
                                           gap.info     = align.tab)

if (!all.singletons.flag) {
  
  add.gaps.tab <- add.gaps.tab |> na.omit()
  
  add.gaps.list <- list(File = as.vector(add.gaps.tab$file),
                        Up   = add.gaps.tab$offset_up,
                        Down = add.gaps.tab$offset_down)
  
  add.gaps.tab |> na.omit()
  
  plan(multisession, workers = params.list$nb_workers)
  options(future.globals.maxSize = 400000000000000000)
  
  furrr::future_pwalk(.l = add.gaps.list,
                      .f = ~add.gaps.transfac.motif(tf.file.in  = ..1,
                                                    gap.up      = ..2,
                                                    gap.down    = ..3,
                                                    tf.file.out = ..1))
}


# --------------------------------------- #
# Aligned motifs as UniversalMotif object #
# --------------------------------------- #
aligned.motif.files <- add.gaps.tab$file[grepl(add.gaps.tab$file, pattern = "oriented\\.tf")]
aligned.motifs.um <- purrr::map(.x = as.vector(aligned.motif.files),
                                .f = ~read.motif.file(motif.file   = .x,
                                                      motif.format = "tf"))

aligned.motif.rc.files <- add.gaps.tab$file[grepl(add.gaps.tab$file, pattern = "oriented_rc\\.tf")]
aligned.motifs.rc.um <- purrr::map(.x = as.vector(aligned.motif.rc.files),
                                   .f = ~read.motif.file(motif.file   = .x,
                                                         motif.format = "tf"))


#################################################################
## 2) Export motifs separated by clusters one file per cluster
message("; Exporting motifs separated by clusters: ", out.folder.list$cluster_motifs)
no.output <- sapply(file.path(out.folder.list$cluster_motifs, unique(results.list$Clusters_table$cluster)), dir.create, recursive = TRUE, showWarnings = FALSE)

motifs.files.per.cluster <- export.aligned.motifs.per.cluster(indiv.motis.folder    = out.folder.list$indiv_motifs,
                                                              cluster.motifs.folder = out.folder.list$cluster_motifs,
                                                              cluster.motif.id.tab  = find.clusters.list$clusters_df)


#################################################################
## 3) Export root motifs, one file per cluster
message("; Exporting root motifs: ", output.files.list$Root_motifs)
root.motifs.tf.vector <- purrr::map(.x = sort(motifs.files.per.cluster),
                                    .f = ~export.root.motif(cluster.tf.file = .x))
suppressMessages(root.motifs.tf.vector <- unlist(root.motifs.tf.vector))

## Print the vector as a text file
invisible(suppressWarnings(file.remove(output.files.list$Root_motifs, showWarnings = FALSE, recursive = TRUE)))
suppressMessages(writeLines(root.motifs.tf.vector, con = output.files.list$Root_motifs))



## When minimal output mode is not activated exports trees, heatmap, and cluster-color table
if (params.list$min_output == FALSE) {
  
  # -------------------------------- #
  # Export the tree as a newick file #
  # -------------------------------- #
  if (params.list$export_newick == 1) {
    message("; Exporting tree as a newick file", output.files.list$Newick_tree_all)
    newick.tree <- convert.hclust.to.newick(results.list$All_motifs_tree)
    write(newick.tree, file = output.files.list$Newick_tree)
  }
  
  
  # ------------------------------------------- #
  # Export the tree as a json file and as RData #
  # ------------------------------------------- #
  JSON.tree <- convert.hclust.to.JSON(hc = results.list$All_motifs_tree)
  message("; Exporting tree as a JSON file: ", output.files.list$JSON_tree_all)
  writeLines(JSON.tree, con = output.files.list$JSON_tree_all)
  
  
  message("; Exporting hclust object as a RData file: ", output.files.list$hclust_all)
  tree.rdata <- results.list$All_motifs_tree
  save(tree.rdata, file = output.files.list$hclust_all)
  
  
  # ------------------------- #
  # Export the distance table #
  # ------------------------- #
  message("; Exporting distance table: ", output.files.list$Distance_table)
  fwrite(x         = results.list[["Dist_table"]],
         file      = output.files.list$Distance_table,
         row.names = TRUE,
         col.names = TRUE,
         sep       = "\t")
  
  
  # -------------------------- #
  # Cluster - Hexa-Color table #
  # -------------------------- #
  
  results.list$Cluster_color <- cluster.color.map(cluster.tab = results.list$Clusters_table)
  cl.col                     <- results.list$Cluster_color
  message("; Exporting cluster-color table: ", output.files.list$Distance_table)
  fwrite(x         = cl.col,
         file      = output.files.list$Distance_table,
         row.names = FALSE,
         col.names = FALSE,
         sep       = "\t")
  
  
  # --------------------- #
  # Draw clusters heatmap #
  # --------------------- #
  if (params.list$export_heatmap == 1) {
    
    message("; Drawing cluster heatmap")
    heatmap.w.clusters <- draw.heatmap.motifs(dist.table    = results.list[["Original_matrix"]], 
                                              clusters      = results.list$Clusters_table, 
                                              metric        = "Ncor", 
                                              cluster.color = results.list$Cluster_color,
                                              color.palette = params.list$heatmap_color_palette,
                                              color.classes = params.list$heatmap_color_classes)
    
    message("; Exporting heatmap as PDF file: ", output.files.list$Heatmap_clusters)
    pdf(file = output.files.list$Heatmap_clusters, width = 15, height = 17)
    draw(heatmap.w.clusters)
    dev.off()
  }
  
  
  # ----------------------- #
  # Draw Rand Index heatmap #
  # ----------------------- #
  if (params.list$ARI) {
    
    message("; Drawing Adjusted Rand Index heatmap")
    heatmap.ari <- draw.heatmap.clusters.vs.ref(clusters.tab = clustering.ari$tab,
                                                comment      = paste0("ARI = ", round(clustering.ari$ARI, digits = 2)))
    
    export.heatmap.clusters.vs.ref(ht        = heatmap.ari,
                                   ht.matrix = clustering.ari$tab,
                                   pdf.file  = output.files.list$Clusters_vs_Reference)
    
    invisible(suppressWarnings(file.remove("Rplots.pdf", showWarnings = FALSE)))
  }
  
  
  # ------------ #
  # Export logos #
  # ------------ #
  message("; Exporting logos in ", out.folder.list["aligned_logos"])
  export.logos(um        = aligned.motifs.um,
               outfolder = out.folder.list["aligned_logos"],
               rev_tag   = FALSE)
  
  export.logos(um        = aligned.motifs.rc.um,
               outfolder = out.folder.list["aligned_logos"],
               rev_tag   = TRUE)
  
  # Add the columns with the path to the logos
  results.list$Motif_info_tab$Logo    <- file.path(out.folder.list["aligned_logos"], paste0(results.list$Motif_info_tab$id, ".svg"))
  results.list$Motif_info_tab$Logo_RC <- file.path(out.folder.list["aligned_logos"], paste0(results.list$Motif_info_tab$id, "_rc.svg"))

  
  # --------------------------------- #
  # Create radial tree HTML + D3 file #
  # --------------------------------- #
  if (params.list[["radial_tree"]]) {
    
    barplot.annotation.file <- NULL
    if (params.list[["annotation"]]) {
      
      motif.annotation.list <- create.color.annotation(motif.meta.file = motif.annotation.file,
                                                       ann.outdir      = out.folder.list$tables)
  
      barplot.annotation.file <- file.path(out.folder.list$plots, paste0(params.list.radial$title, "_barplot.svg"))
      
      # Create annotation barplot
      annotation.barplot(df            = motif.annotation.list$df,
                         class.colors  = motif.annotation.list$cpal,
                         barplot.title = params.list.radial$title,
                         barplot.file  = barplot.annotation.file)
      
      # Add url and color from the annotation table
      results.list$Motif_info_tab <- merge(results.list$Motif_info_tab, motif.annotation.list$df, by.x = "id_old", by.y = "motif_id") |> 
                                      within(rm(collection, class, class_nb))
      
    } else {
      results.list$Motif_info_tab$colour <- "#000000"
      results.list$Motif_info_tab$url    <- ""
    }

    # Update JSON file with annotations
    message("; Adding attributes to JSON file")
    Add_attributes_to_JSON_radial_tree(motif.description.tab = results.list$Motif_info_tab,
                                       clusters.list         = find.clusters.list,
                                       color.map             = cl.col,
                                       htree                 = results.list$All_motifs_tree,
                                       json.woa.file         = output.files.list$JSON_tree_all,
                                       json.wa.file          = output.files.list$JSON_radial_annotated,
                                       alignent.width        = max(results.list$Alignment_radial_table$width))
    
    # Fill the HTML + D3 template 
    create.html.radial.tree(json.file        = output.files.list$JSON_radial_annotated,
                            d3.template      = d3.radial.tree.template,
                            d3.outfile       = output.files.list$D3_radial_tree,
                            motif.info       = results.list$Motif_info_tab,
                            d3.lib           = d3.min.lib,
                            jq.lib           = jquery.lib,
                            outdir           = dirname(out.folder),
                            alignment.length = unique(results.list$Alignment_radial_table$width)[1],
                            html.legend      = motif.annotation.list$html,
                            barplot.ann      = barplot.annotation.file)
    
    if (params.list[["annotation"]]) {

      # Add color background to radial tree
      annotate.radial.tree(clusters         = find.clusters.list$clusters,
                           cluster2color    = results.list$Cluster_color,
                           tree             = results.list$All_motifs_tree,
                           motif.annotation = motif.annotation.list$df)
    }

  } else {

    # ------------------------ #
    # Create interactive trees #
    # ------------------------ #
    # save.image("Debug_radial.Rdata")
    plan(multisession, workers = params.list$nb_workers)
    options(future.globals.maxSize = 400000000000000000)

    message("; Adding attributes to JSON files (Interactive trees)")
    furrr::future_walk(.x = names(find.clusters.list$clusters),
                       .f = ~  Add_attributes_to_JSON_interactive_tree(cluster_name          = .x,
                                                                       motif.description.tab = results.list$Motif_info_tab,
                                                                       clusters.list         = find.clusters.list,
                                                                       color.map             = cl.col,
                                                                       htree                 = cl.hclust.results,
                                                                       json.file.df          = results.list$Clusters_files,
                                                                       alignment.df          = results.list$Alignment_table))

    
    create.html.interactive.tree(clusters      = find.clusters.list$clusters,
                                 clusters.df   = results.list$Clusters_files,
                                 cluster.color = cl.col,
                                 html.template = html.interactive.tree.template,
                                 hmtl.ready    = output.files.list$D3_dynamic_tree,
                                 d3.template   = d3.interactive.tree.template,
                                 d3.lib        = d3.min.lib,
                                 jq.lib        = jquery.lib,
                                 outdir        = dirname(out.folder))
    
  }
  
## Remove these folder when --minimal_output mode is activated
} else {
  
  # -------------------- #
  # Delete empty folders #
  # -------------------- #
  suppressWarnings(unlink(out.folder.list$plots, recursive = TRUE))
  suppressWarnings(unlink(out.folder.list$trees, recursive = TRUE))
  suppressWarnings(unlink(out.folder.list$central_motifs, recursive = TRUE))
  suppressWarnings(unlink(out.folder.list$root_motifs, recursive = TRUE))
  suppressWarnings(unlink(out.folder.list$cluster_motifs, recursive = TRUE))
  invisible(suppressWarnings(file.remove(output.files.list$Motifs_transfac, showWarnings = FALSE)))
  invisible(suppressWarnings(file.remove(output.files.list$Motif_compa, showWarnings = FALSE)))
  invisible(suppressWarnings(file.remove(output.files.list$Summary_table, showWarnings = FALSE)))

}
message("; End of program")