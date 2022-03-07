#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("dplyr",          ## Data manipulation
                      "data.table",     ## Read long matrices in a quick way
                      "furrr",          ## Run functions in parallel
                      "optparse",       ## Read command-line arguments
                      "purrr",          ## Iterations,
                      "this.path",      ## Create relative paths
                      "tidyr",          ## Data manipulation
                      "universalmotif") ## Motif manipulation (Bioconductor)

for (lib in required.packages) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}


####################
## Read arguments ##
####################
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
  
  make_option(c("--heatmap_color_palette"), type = "character", default = "RdGy", 
              help = "Cell colors in heatmap. [Default: \"%default\"]. Options: any colorBrewer palette, see colorbrewer2.org", metavar = "character"),
  
  make_option(c("--color_palette_classes"), type = "numeric", default = 11, 
              help = "Number of classes to create color palette. [Default: \"%default\"]. Options: any colorBrewer palette, see colorbrewer2.org", metavar = "number"),
  
  make_option(c("--cor_th"), type = "numeric", default = 0.75, 
              help = "Number of cores to run in parallel. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("--Ncor_th"), type = "numeric", default = 0.55, 
              help = "Number of cores to run in parallel. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("--minimal_output"), type = "logical", default = FALSE, 
              help = "When TRUE only returns the alignment, clusters and motif description tables. Comparison results, plots and trees are not exported. [Default \"%default\"] ", metavar = "logical")
  
);
message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Output file prefix
out.folder        <- opt$output_folder
matrix.file.table <- opt$matrix_file_table

## Mandatory input
if (!file.exists(matrix.file.table)) {
  stop("Mandatory input file not found: ", matrix.file.table)
}

params.list <- list("export_newick"         = as.numeric(opt$export_newick),
                    "export_heatmap"        = as.numeric(opt$export_heatmap),
                    "heatmap_color_classes" = as.numeric(opt$color_palette_classes),
                    "heatmap_color_palette" = opt$heatmap_color_palette,
                    "comparison_metric"     = opt$comparison_metric,
                    "linkage_method"        = opt$linkage_method,
                    # "w"                   = min.w,
                    "cor"                   = as.numeric(opt$cor_th),
                    "Ncor"                  = as.numeric(opt$Ncor_th),
                    "nb_workers"            = as.numeric(opt$number_of_workers),
                    "min_output"            = opt$minimal_output)


###############################
## Source here the functions ##
###############################

## 'here' is set to the path where matrix-clustering.R is located, then the libraries
## are sourced relative to where 'here' was set
source(this.path::here(.. = 0, "R", "General_utils.R"))
source(this.path::here(.. = 0, "R", "Hierarchical_clustering.R"))
source(this.path::here(.. = 0, "R", "Motif_alignment_utils.R"))
source(this.path::here(.. = 0, "R", "Motif_manipulation.R"))
source(this.path::here(.. = 0, "R", "Tree_partition_utils.R"))
# sourceCpp(file.path(params.list$clustering_lib_path, "Utils.cpp"))


## oriented instead of aligned + oriented_rc


###########
## Debug ##
###########
## Example:
## Rscript matrix-clustering.R -i data/OCT4_datasets/OCT4_motif_table.txt -o results/OCT4_motifs_example/OCT4_motif_analysis --number_of_workers 8 -q compare-matrices-quick/compare-matrices-quick --minimal_output FALSE -r ./R

# matrix.file.table <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/data/OCT4_datasets/OCT4_motif_table.txt"
# out.folder        <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/results/Oct4_example/Oct4_example"
# 
# params.list <- list("export_newick"         = 0,
#                     "export_heatmap"        = 0,
#                     "heatmap_color_classes" = NULL,
#                     "heatmap_color_palette" = NULL,
#                     "comparison_metric"     = "Ncor",
#                     "linkage_method"        = "average",
#                     # "w"                   = min.w,
#                     "cor"                   = 0.75,
#                     "Ncor"                  = 0.55,
#                     "nb_workers"            = 8,
#                     "min_output"            = TRUE)


##############################################################
## Initialize output + result lists + create output folders ##
##############################################################

out.folder.list <- list(tables = file.path(paste0(out.folder, "_tables")),
                        motifs = file.path(paste0(out.folder, "_motifs")),
                        plots  = file.path(paste0(out.folder, "_plots")),
                        trees  = file.path(paste0(out.folder, "_trees")))

output.files.list <- list("Alignment_table"     = file.path(out.folder.list$tables, "alignment_table.tab"),
                          "Clusters_table"      = file.path(out.folder.list$tables, "clusters.tab"),
                          "Distance_table"      = file.path(out.folder.list$tables, "distance_table.tab"),
                          "Motif_compa"         = file.path(out.folder.list$tables, "pairwise_motif_comparison.tab"),
                          "Cluster_colors"      = file.path(out.folder.list$tables, "cluster_color_map.tab"),
                          "Motifs_transfac_tmp" = file.path(out.folder.list$motifs, "input_motifs_parsed_id.tf.tmp"),
                          "Motifs_transfac"     = file.path(out.folder.list$motifs, "input_motifs_parsed_id.tf"),
                          "Heatmap"             = file.path(out.folder.list$plots, "Heatmap.pdf"),
                          "hclust_all"          = file.path(out.folder.list$trees, "tree.RData"),
                          "JSON_tree_all"       = file.path(out.folder.list$trees, "tree.json"),
                          "Newick_tree_all"     = file.path(out.folder.list$trees, "tree.newick"))


results.list <- list(Dist_table       = NULL,
                     Dist_matrix      = NULL,
                     Original_matrix  = NULL,
                     All_motifs_tree  = NULL,
                     Alignment_table  = NULL,
                     Cluster_color    = NULL,
                     Motif_info_tab   = NULL,
                     Motif_compa_tab  = NULL)

## Create output folders
no.output <- sapply(sapply(output.files.list, dirname), dir.create, showWarnings = FALSE, recursive = TRUE)


################################################################
## Pre-process motif files + Generate motif description table ##
################################################################

## Note: for the moment it only runs with one collection, but this should be adapted to run with many
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
export.tf.motifs.w.parsed.header(old.tf.file = output.files.list$Motifs_transfac_tmp,
                                 new.tf.file = output.files.list$Motifs_transfac,
                                 um.object   = all.motifs.um) 

params.list[["Nb_motifs"]]      <- nrow(results.list$Motif_info_tab)
params.list[["Nb_collections"]] <- length(matrix.file.list$Motif_file)
message("; Analysis with ", params.list[["Nb_motifs"]], " motifs in ", params.list[["Nb_collections"]], " collections")


#################################
## Read motif comparison table ##
#################################
results.list$Motif_compa_tab <- motif.comparison(transfac.file     = output.files.list$Motifs_transfac,
                                                 output.compa.file = output.files.list$Motif_compa)


#######################################################################
## Convert distance table into a distance matrix, required by hclust ##
#######################################################################
message('; Converting correlation values to distances')
distances.objects <- build.distance.matrix(compa.table = results.list$Motif_compa_tab,
                                           metric      = params.list$comparison_metric)

results.list[["Dist_table"]]      <- data.table(as.data.frame.matrix(distances.objects$table))
results.list[["Dist_matrix"]]     <- distances.objects$matrix
results.list[["Original_matrix"]] <- data.table(as.data.frame.matrix(distances.objects$original))


##############################################
## Compute the hierarchical clustering tree ##
##############################################

## Hierarchical clustering can only be applied with > 1 elements
if (params.list[["Nb_motifs"]] > 1) {
  
  results.list$All_motifs_tree <- hclust.motifs(results.list[["Dist_matrix"]],
                                                hclust.method = params.list$linkage_method)
  
  
  ###########################
  ## Find clusters section ##
  ###########################
  message("; Finding clusters")
  find.clusters.list <- find.motif.clusters(tree             = results.list$All_motifs_tree,
                                            comparison.table = results.list$Motif_compa_tab,
                                            parameters       = params.list)
  
  
  ## Identify singletons and clusters with many motifs
  cluster.sizes                <- purrr::map_dbl(find.clusters.list$clusters, length)
  cl.singleton                 <- names(cluster.sizes[cluster.sizes == 1])
  cl.many                      <- names(cluster.sizes[cluster.sizes > 1])
  singleton.flag               <- as.logical(length(cl.singleton))
  params.list[["Nb_clusters"]] <- length(find.clusters.list$clusters)
  message("; Number of clusters: ", params.list[["Nb_clusters"]])
  
  
  #############################
  ## Motif alignment section ##
  #############################
  
  ## Compute hierarchical clustering of each cluster
  cl.hclust.results <- lapply(find.clusters.list$clusters, hclust.cluster.ids, compa = results.list$Motif_compa_tab, parameters = params.list)
  
  ## Extract the hclust objects within the nested list
  cl.many.hclust <- purrr::map(cl.hclust.results[cl.many], `[[`, "hclust")
  
  ## Align motifs within each cluster (with 2 or more motifs, singletons are treated separately)
  ## This function is ran in parallel using furrr::future_map , see https://furrr.futureverse.org/
  plan(multisession, workers = params.list$nb_workers)
  message("; Motif alignment step")
  aligment.clusters.tab <- furrr::future_map(.x = cl.many.hclust,
                                             .f = ~align.motifs.in.cluster(tree       = .x,
                                                                           compa      = results.list$Motif_compa_tab,
                                                                           motif.info = results.list$Motif_info_tab,
                                                                           parameters = params.list))
  
  ###########################
  ## Debug : do not delete ##
  ###########################
  # tree.example <- cl.many.hclust$cluster_001
  # tree.example <- cl.many.hclust$cluster_02
  # tree.example <- cl.many.hclust$cluster_03
  # tree.example <- cl.many.hclust$cluster_04
  # tree.example <- cl.many.hclust$cluster_05
  # tree.example <- cl.many.hclust$cluster_06
  # tic()
  # aaaa <- align.motifs.in.cluster(tree = tree.example,
  #                         compa      = results.list$Motif_compa_tab,
  #                         motif.info = results.list$Motif_info_tab,
  #                         parameters = params.list)
  # toc()
  
  ###################################
  ## Parse tables before exporting ##
  ###################################
  
  ## Prepare tables to export
  message("; Combining aligned clusters tables")
  
  ## Groups
  ## Add cluster name, rename columns and remove unnecessary columns
  alignment.clusters.tab.export <-  rbindlist(aligment.clusters.tab, idcol = "cluster") %>% 
                                      within(rm(N, Update_status)) %>% 
                                      rename(strand            = Strand,
                                             offset_up         = Offset_up,
                                             offset_down       = Offset_down,
                                             aligned_consensus = Oriented_consensus)
  
  
  ## Singleton section
  ## This if only applied when there is at least one singleton
  if (singleton.flag) {
   
    singleton.clusters.tab.export <- data.table(cluster = names(find.clusters.list$clusters[cl.singleton]),
                                                id      = unlist(find.clusters.list$clusters[cl.singleton])) %>% 
                                        left_join(results.list$Motif_info_tab, by = "id") %>% 
                                        mutate(strand            = "D",
                                               offset_up         = 0,
                                               offset_down       = 0,
                                               aligned_consensus = consensus) %>% 
                                        within(rm(n, width, IC, nb_sites, id_old)) 
    
    ## Combine groups + singleton clusters
    results.list$Alignment_table <- rbind(alignment.clusters.tab.export, singleton.clusters.tab.export) %>% 
                                        mutate(alignment_width = nchar(aligned_consensus))
    
  } else {
    results.list$Alignment_table <- alignment.clusters.tab.export %>% 
                                      mutate(alignment_width = nchar(aligned_consensus))
  }

  
  #########################
  ## Clusters - ID table ##
  #########################
  results.list$Clusters_table <- results.list$Alignment_table %>% 
                                  group_by(cluster) %>% 
                                  summarise(id      = paste(id, collapse = ","),
                                            name    = paste(name, collapse = ","),
                                            .groups = "drop")
  
} else {
  stop("The clustering analysis required at least 2 motifs.")
}



############################
## Export results section ##
############################

#############################################################################################
## Export the alignment + clusters table + motif files: these files are the minimal output ##
#############################################################################################
message("; Exporting alignment table:", output.files.list$Alignment_table)
fwrite(x         = results.list$Alignment_table,
       file      = output.files.list$Alignment_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")


message("; Exporting clusters table:", output.files.list$Clusters_table)
fwrite(x         = results.list$Clusters_table,
       file      = output.files.list$Clusters_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")

save("Debug.Rdata")

## Export motifs as transfac files in D and R orientation
message("; Export individual motifs in transfac format")
export.indiv.motif.files(un.motifs = all.motifs.um,
                         outdir    = out.folder.list$motifs)


## When minimal output mode is not activated exports trees, heatmap, and cluster-color table
if (params.list$min_output == FALSE) {
  
  ######################################
  ## Export the tree as a newick file ##
  ######################################
  if (params.list$export_newick == 1) {
    message("; Exporting tree as a newick file", output.files.list$Newick_tree_all)
    newick.tree <- convert.hclust.to.newick(results.list$All_motifs_tree)
    write(newick.tree, file = output.files.list$Newick_tree)
  }
  
  
  #################################################
  ## Export the tree as a json file and as RData ##
  #################################################
  
  ## NOTE: JSON tree must be exported
  JSON.tree <- convert.hclust.to.JSON(results.list$All_motifs_tree)
  message("; Exporting tree as a JSON file: ", output.files.list$JSON_tree_all)
  writeLines(JSON.tree, con = output.files.list$JSON_tree_all)
  
  message("; Exporting hclust object as a RData file: ", output.files.list$hclust_all)
  tree.rdata <- results.list$All_motifs_tree
  save(tree.rdata, file = output.files.list$hclust_all)
  
  
  ###############################
  ## Export the distance table ##
  ###############################
  message("; Exporting distance table:", output.files.list$Distance_table)
  fwrite(x         = results.list[["Dist_table"]],
         file      = output.files.list$Distance_table,
         row.names = TRUE,
         col.names = TRUE,
         sep       = "\t")
  
  
  ###########################
  ## Cluster - Color table ##
  ###########################
  results.list$Cluster_color <- cluster.color.map(cluster.tab = results.list$Clusters_table)
  cl.col                     <- results.list$Cluster_color
  message("; Exporting cluster-color table:", output.files.list$Distance_table)
  fwrite(x         = cl.col,
         file      = output.files.list$Distance_table,
         row.names = FALSE,
         col.names = FALSE,
         sep       = "\t")
  
  ##################
  ## Draw heatmap ##
  ##################
  if (params.list$export_heatmap == 1) {
    
    message("; Drawing heatmap")
    heatmap.w.clusters <- draw.heatmap.motifs(dist.table    = results.list[["Original_matrix"]], 
                                              clusters      = results.list$Clusters_table, 
                                              metric        = "Ncor", 
                                              cluster.color = results.list$Cluster_color,
                                              color.palette = params.list$heatmap_color_palette,
                                              color.classes = params.list$heatmap_color_classes)
    
    message("; Exporting heatmap as PDF file: ", output.files.list$Heatmap)
    pdf(file = output.files.list$Heatmap, width = 15, height = 17)
    draw(heatmap.w.clusters)
    dev.off()
  }
  
## Remove these folder when --minimal_output mode is activated
} else {
  
  ##########################
  ## Delete empty folders ##
  ##########################
  unlink(out.folder.list$plots, recursive = TRUE)
  unlink(out.folder.list$trees, recursive = TRUE)
}
message("; End of program")