#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("dplyr",          ## Data manipulation
                      "data.table",     ## Read long matrices in a quick way
                      "furrr",          ## Run functions in parallel
                      "optparse",       ## Read command-line arguments
                      "purrr",          ## Iterations
                      "rcartocolor",    ## Cluster colors
                      "tidyr",          ## Data manipulation
                      "universalmotif") ## Motif manipulation


for (lib in required.packages) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


###################################
## Source here revised functions ##
###################################
source("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/General_utils.R")
source("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/hierarchical_clustering.R")
source("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/Motif_alignment_utils.R")
source("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/Motif_manipulation.R")
source("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/Tree_partition_utils.R")
# sourceCpp("/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/TFBMclust_revised/Utils.cpp")



####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-i", "--matrix_file_table"), type = "character", default = NULL, 
              help = "A text-delimited file where each line contain the next fields. 1: Motif file path; 2: Motif collection name; 3: Motif format (Mandatory). It does not expect a header, but it expects those columns in the indicated order.", metavar = "character"),
  
  make_option(c("-q", "--compare_matrices_quick_path"), type = "character", default = NULL, 
              help = "Path to the executable of compare-matrices-quick (Mandatory). From the main directory of the repository: compare-matrices/compare-matrices-quick", metavar = "character"),
  
  make_option(c("-o", "--output_folder"), type = "character", default = NULL, 
              help = "Folder to save the results (Mandatory)", metavar = "character"),
  
  # make_option(c("-t", "--title"), type = "character", default = "matrix-clustering_analysis", 
  #             help = "Analysis title, this is part of the ouput file names. [Default: \"%default\".", metavar = "character"),  
  
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
  
  make_option(c("--minimal_output"), type = "logical", default = TRUE, 
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



###########
## Debug ##
###########
## Example:
## Rscript /home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/matrix-clustering.R -i /home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/data/OCT4_datasets/OCT4_motif_table.txt -o /home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/results/OCT4_motifs_example --number_of_workers 8 -q /home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/compare-matrices/compare-matrices-quick --minimal_output FALSE

# matrix.file.tab <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/data/motif_table.txt"
# out.folder  <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/results/JASPAR_2022_example"

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
#                     "compare_matrices_path" = "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix_clustering/RSAT_matrix_clustering_tidyR/matrix-clustering/compare-matrices/compare-matrices-quick",
#                     "min_output"            = TRUE)



###########
## Flags ##
###########
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
                    "compare_matrices_path" = opt$compare_matrices_quick_path,
                    "min_output"            = opt$minimal_output)


##############################################################
## Initialize output + result lists + create output folders ##
##############################################################
output.files.list <- list("Alignment_table"     = file.path(paste0(out.folder, "_tables"), "alignment_table.tab"),
                          "Clusters_table"      = file.path(paste0(out.folder, "_tables"), "clusters.tab"),
                          "Distance_table"      = file.path(paste0(out.folder, "_tables"), "distance_table.tab"),
                          "Motif_compa"         = file.path(paste0(out.folder, "_tables"), "pairwise_motif_comparison.tab"),
                          "Cluster_colors"      = file.path(paste0(out.folder, "_tables"), "cluster_color_map.tab"),
                          "Motifs_transfac_tmp" = file.path(paste0(out.folder, "_motifs"), "input_motifs_parsed_id.tf.tmp"),
                          "Motifs_transfac"     = file.path(paste0(out.folder, "_motifs"), "input_motifs_parsed_id.tf"),
                          "Heatmap"             = file.path(paste0(out.folder, "_plots"), "Heatmap.pdf"),
                          "hclust_all"          = file.path(paste0(out.folder, "_trees"), "tree.RData"),
                          "JSON_tree_all"       = file.path(paste0(out.folder, "_trees"), "tree.json"),
                          "Newick_tree_all"     = file.path(paste0(out.folder, "_trees"), "tree.newick"))


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
                                                 output.compa.file = output.files.list$Motif_compa,
                                                 bin               = params.list$compare_matrices_path)


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
  aligment.clusters.tab <- furrr::future_map(.x = cl.many.hclust, .f = ~align.motifs.in.cluster(tree       = .x,
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
  aligment.clusters.tab.export <-  rbindlist(aligment.clusters.tab, idcol = "cluster") %>% 
                                      within(rm(N, Update_status)) %>% 
                                      rename(strand            = Strand,
                                             offset_up         = Offset_up,
                                             offset_down       = Offset_down,
                                             aligned_consensus = Oriented_consensus)
  
  
  ## Singletons
  singleton.clusters.tab.export <- data.table(cluster = names(find.clusters.list$clusters[cl.singleton]),
                                              id      = unlist(find.clusters.list$clusters[cl.singleton])) %>% 
                                      left_join(results.list$Motif_info_tab, by = "id") %>% 
                                      mutate(strand            = "D",
                                             offset_up         = 0,
                                             offset_down       = 0,
                                             aligned_consensus = consensus) %>% 
                                      within(rm(n, width, IC, nb_sites, id_old))
    
  
  ## Combine groups + singleton clusters
  clusters.tab.export <- rbind(aligment.clusters.tab.export, singleton.clusters.tab.export) %>% 
                            mutate(alignment_width = nchar(aligned_consensus))
  
  results.list$Alignment_table <- clusters.tab.export


  
  #########################
  ## Clusters - ID table ##
  #########################
  results.list$Clusters_table <- results.list$Alignment_table %>% 
                                    group_by(cluster) %>% 
                                    summarise(id   = paste(id, collapse = ","),
                                              name = paste(name, collapse = ","))
  
} else {
  stop("The clustering analysis required at least 2 motifs.")
}




############################
## Export results section ##
############################

################################################################################
## Export the alignment + clusters table:  these files are the minimal output ##
###############################################################################
message("; Exporting alignment table:", output.files.list$Alignment_table)
fwrite(x         = results.list$Alignment_table,
       file      = output.files.list$Alignment_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")


message("; Exporting clusters table:", output.files.list$Alignment_table)
fwrite(x         = results.list$Clusters_table,
       file      = output.files.list$Clusters_table,
       row.names = FALSE,
       col.names = TRUE,
       sep       = "\t")


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
}
message("; End of program")



# compatab <- results.list$Motif_compa_tab
# save(compatab, file = "Compa.Rdata")
# 
# load("Compa.Rdata")
# library(Rcpp)
# library(data.table)
# library(tictoc)
# sourceCpp("TFBMclust_revised/Utils.cpp")


# tic()
# aa <- MotifCompaEntry(compa_df = compatab,
#                 ids      = c("JASPAR_vertebrates_MA0004.1_n1", "JASPAR_vertebrates_MA0006.1_n2", "JASPAR_vertebrates_MA0019.1_n3"),
#                 # ids      = c("JASPAR_vertebrates_MA0004.1_n1", "JASPAR_vertebrates_MA0006.1_n2"),
#                 full     = FALSE,
#                 self     = FALSE)
# # toc()
# aa <- data.table(aa)
# aa
  
  ##############################################
  ## This block is only executed when the option
  ## -radial_tree_only is activated
  # if (radial.only == 1) { 
  # 
  # ################################################
  # ## Given a level of a hierarchical tree, find
  # ## the next level pointing the current level
  # ## NOTE: this function is only called within the function find.clusters
  # find.next.levels.in.tree <- function(x) {
  #   
  #   if (x == length(tree$merge)/2) {
  #     return(x)
  #     
  #   } else {
  #     
  #     ## Get the level
  #     level <- which(tree$merge == x)
  #     
  #     if (level > length(tree$merge)/2 & level <= length(tree$merge)) {
  #       return(level - length(tree$merge)/2)
  #     } else if (level <= length(tree$merge)/2) {
  #       return(level)
  #     }
  #   }
  # }
  # 
  # find.chained.levels <- function(x) {
  #   
  #   chained.levels <- NULL
  #   current.level <- NULL
  #   alignment.flag <- 1
  #   chained.levels <- append(chained.levels, x)
  #   
  #   while (alignment.flag == 1) {
  #     
  #     ## Search the next chained level
  #     current.level <- find.next.levels.in.tree(x)
  #     
  #     ## The case when the root tree is analyzed
  #     if (x == current.level) {
  #       return(x)
  #     }
  #     
  #     ## Get the alignment status at the level
  #     alignment.flag <- as.numeric(attributes.list[[paste("node_", current.level, sep = "")]][["alignment_flag"]])
  #     
  #     ## If the status at the current level is 0 then return the chained levels
  #     if (alignment.flag == 0) {
  #       return(chained.levels)
  #       
  #       ## Conversely, add the new level to the chain
  #     } else {
  #       chained.levels <- append(chained.levels, current.level)
  #       x <- current.level
  #     }
  #   }
  # }
  # 
  # 
  # chained.levels <- sapply(1:nrow(tree$merge), find.chained.levels)
  # alignment.flag <- as.numeric(sapply(alignment.attributes, function(n) {n[[1]]}))
  # 
  # node.to.cluster <- data.frame()
  # cluster.counter <- 0
  # 
  # ## Iterate over the clusters
  # x <- lapply(clusters, function(cl) {
  #   
  #   counter <- 0
  #   cluster.counter <<- cluster.counter + 1
  #   
  #   ## Iterate over the leaves per node
  #   xx <- lapply(leaves.per.node(tree), function(l) {
  #     
  #     counter <<- counter + 1
  #     
  #     ## Detect the nodes with the same number of motifs
  #     ## than the cluster
  #     if (length(l) >= length(cl)) {
  #       
  #       ## When all the elements of the cluster are found on the leave
  #       ## return the number of node (counter)
  #       if (all(cl %in% l)) {
  #         counter
  #       } 
  #     }
  #     
  #   })
  #   ## Select the first node where the cluster's members are contained
  #   selected.node <- as.numeric(unlist(xx))[1]
  #   
  #   ## Look for all the nodes pointing to the selected node
  #   all.nodes <- lapply(chained.levels, function(levels) { 
  #     if (is.element(selected.node, levels)) {
  #       levels
  #     }
  #   })
  #   
  #   all.nodes[sapply(all.nodes, function(x) length(all.nodes) == 0)] <- NA
  #   all.nodes <- unlist(all.nodes)
  #   
  #   ## Concatenate all the nodes, sort and remove repetitions
  #   all.nodes <- unique(sort(unlist(sapply(all.nodes, function(x) { 
  #     chained.levels[[x]]})
  #   )
  #   )
  #   )
  #   all.nodes <- paste(all.nodes, collapse = ",")
  #   all.leaves <- paste(leaves.per.node(tree)[[selected.node]], collapse = ",")
  #   node.alignment.status <- alignment.flag[selected.node]
  #   node.to.cluster <<- rbind(node.to.cluster, data.frame(selected.node, cluster.counter, all.nodes, all.leaves, node.alignment.status))
  #   
  # })
  # 
  # x1 <- sapply(as.vector(node.to.cluster$all.leaves), function(l) {
  #   as.numeric(unlist(strsplit(l, ",")))
  # })
  # names(x1) <- NULL
  # 
  # ## Find the motif number that appear in more than one cluster
  # ## This is to treat the node.to.cluster table
  # repeated.leaves <- as.vector(which(table(unlist(x1)) > 1))
  # 
  # th <-  sapply(repeated.leaves, function(r) {
  #   
  #   ## For each repeated number find the clusters where it is repeated
  #   co <- 0
  #   repetition.flag <- 0
  #   repeated.clusters <- sapply(x1, function(x) {
  #     
  #     co <<- co + 1
  #     if (r %in% x) {
  #       co
  #     }
  #     
  #   })
  #   repeated.clusters <- unlist(repeated.clusters)
  #   
  #   cou <- 0
  #   find.flag <- 0
  #   sapply(repeated.clusters, function(x) {
  #     
  #     status.alignment.clusters <- node.to.cluster[repeated.clusters,]$node.alignment.status
  #     
  #     cou <<- cou + 1
  #     
  #     if (status.alignment.clusters[cou] ==  1) {
  #       
  #     }
  #     
  #     if (cou == 1) {
  #       
  #       if (status.alignment.clusters[cou] ==  1) {
  #         NA
  #       } else {
  #         x1[[x]] <<- r
  #       }
  #       
  #     } else {
  #       x1[[x]][x1[[x]] == r] <<- NA
  #     }
  #   })
  #   
  # })
  # x1 <- lapply(x1, function(x) x[!is.na(x)])
  # x1 <- sapply(x1, paste, collapse = ",")
  # node.to.cluster$all.leaves <- x1
  # node.to.cluster$all.nodes <- as.vector(node.to.cluster$all.nodes)
  # node.to.cluster[which(node.to.cluster$node.alignment.status == 0),]$all.nodes <- 0
  # 
  # ## Export the table node -> cluster
  # ## Each node has its associated cluster
  # ## Nodes no associated to a cluster are indicated with 0
  # node.to.cluster.table <- data.frame()
  # nnodes <- nrow(tree$merge)
  # th <- sapply(node.to.cluster$all.nodes, function(n) {
  #   
  #   nodes <- as.numeric(unlist(strsplit(n, ",")))
  #   cluster <- node.to.cluster[which(node.to.cluster$all.nodes == n), "cluster.counter"]
  #   node.to.cluster.table <<- rbind(node.to.cluster.table, data.frame(nodes, cluster))
  #   
  # })
  # node.to.cluster.table <- unique(node.to.cluster.table)
  # 
  # all.nodes.in.tree <- 0:nnodes
  # nodes.in.clusters <- unique(sort(node.to.cluster.table$nodes))
  # missing.nodes <- setdiff(all.nodes.in.tree, nodes.in.clusters)
  # node.to.cluster.table <- rbind(node.to.cluster.table, data.frame(nodes = missing.nodes, cluster = 0))
  # 
  # node.to.cluster.table.file <- paste(sep = "", out.prefix, "_tables/node_to_cluster.tab")
  # write.table(node.to.cluster.table, file = node.to.cluster.table.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  # verbose(paste("Exporting node -> cluster table", node.to.cluster.table.file), 2)
  # 
  # ## Export the table node -> cluster
  # ## Each node has its associated cluster
  # ## Nodes no associated to a cluster are indicated with 0
  # leaf.to.cluster.table <- data.frame()
  # nleaves <- nrow(tree$merge) + 1
  # th <- sapply(node.to.cluster$all.leaves, function(n) {
  #   
  #   leaves <- as.numeric(unlist(strsplit(n, ",")))
  #   leaves <- lapply(leaves, function(x) {
  #     get.id(x)
  #   })
  #   leaves <- as.vector(unlist(leaves))
  #   cluster <- node.to.cluster[which(node.to.cluster$all.leaves == n), "cluster.counter"]
  #   leaf.to.cluster.table <<- rbind(leaf.to.cluster.table, data.frame(leaves, cluster))
  #   
  # })
  # 
  # leaf.to.cluster.table.file <- paste(sep = "", out.prefix, "_tables/leaf_to_cluster.tab")
  # write.table(leaf.to.cluster.table, file = leaf.to.cluster.table.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  # verbose(paste("Exporting leaf -> cluster table", leaf.to.cluster.table.file), 2)
  # 
  # clusters <<- list(cluster_1 = sort(unique(as.vector(unlist(clusters)))))
  # 
  # }
  
  # ############################################################################################  
  # clusters <<- lapply(clusters, function(x) {
  #   get.id(x)
  # })
 




  

#########################
## Produce the forests ##
#########################

# ## Creates a folder with where the separated information
# ## of each cluster will be stored
# clusters.info.folder <<- paste(out.prefix, "_clusters_information", sep = "")
# dir.create(clusters.info.folder, recursive=TRUE, showWarnings=FALSE)
# global.motifs.info <<- motifs.info
# forest.list <<- list()
# intermediate.levels.counter <- 0
# intermediate.levels <- vector()
# 
# ## Copy the global tables
# compa.table <- global.compare.matrices.table
# desc.tab <<- global.description.table
# 
# all.central.motifs <- vector()
# 
# central.motif.IDs <- vector()
# central.motif.names <- vector()
# central.motif.IDs.cluster <- vector()
# 
# # nb <- 1
# i <- sapply(1:length(clusters), function(nb) {
#   
#   alignment.cluster <<- list()
#   description.table <<- NULL
#   compare.matrices.table <<- NULL
#   tree <<- NULL
#   
#   ## Check the number of motifs within each cluster
#   ids <- clusters[[paste("cluster", nb, sep = "_")]]
#   if (length(ids) >= 2) {
#     case <- "case.2"
#   } else{
#     case <- "case.1"
#   }
#   
#   singleton.list <<- list()
#   
#   switch(case,
#          
#          ## If the cluster has only one element, create its JSON file and skip the
#          ## hierarchical clustering step
#          case.1 = {
#            
#            ids <- as.character(ids)
#            
#            ## New description table (with the ids of the current cluster)
#            global.description.table <<- desc.tab[desc.tab[,"id"] %in% ids, ]
#            
#            ## Get the central motif name
#            ## This will be used to rename the clusters
#            central.motif <- get.name(ids)[1]
#            central.motif <- gsub("\\(", "_", central.motif)
#            central.motif <- gsub("\\)", "_", central.motif)
#            central.motif <- gsub("\\.", "_", central.motif)
#            central.motif <- paste(central.motif, "cluster", nb, sep = "_")
#            central.motif <- paste("cluster", nb, sep = "_")
#            
#            ## Creates an individual folder for each cluster
#            cluster.folder <<- file.path(clusters.info.folder, central.motif)
#            dir.create(cluster.folder, recursive = TRUE, showWarnings = FALSE)
#            all.central.motifs <<- append(all.central.motifs, central.motif)
#            
#            ## Fill the cluster list with the data of the non-aligned motifs (singleton)
#            global.description.table <<- NULL
#            global.description.table <<- desc.tab
#            forest.list[[central.motif]][[ids]] <<- singleton.list
#            
#            ## Fill the elements of the 'central motif table'
#            central.motif.IDs.cluster[nb] <<- paste("cluster", nb, sep = "_")
#            central.motif.names[nb] <<- get.name(ids)[1]
#            central.motif.IDs[nb] <<- ids[1]
#            
#            ids <- as.character(ids)
#            
#            forest.list[[central.motif]][[ids]][["strand"]] <<- "D"
#            forest.list[[central.motif]][[ids]][["name"]] <<- get.name(ids)
#            forest.list[[central.motif]][[ids]][["consensus_d"]] <<- get.consensus(ids, RC = FALSE)
#            forest.list[[central.motif]][[ids]][["consensus_rc"]] <<- get.consensus(ids, RC = TRUE)
#            forest.list[[central.motif]][[ids]][["number"]] <<- as.numeric(1)
#            forest.list[[central.motif]][[ids]][["spacer.up"]] <<- as.numeric(0)
#            forest.list[[central.motif]][[ids]][["spacer.dw"]] <<- as.numeric(0)
#            
#            if (only.hclust == 0) {
#              
#              ## Create a JSON file for trees with a single node
#              ## In this situation this step is required because it is not possible to use the hclustToJson function
#              JSON.single.node <- paste("{\n\"name\": \"\",\n\"children\":[\n{\n \"label\": \"", ids, "\",\n}\n]\n}", sep = "")
#              json.file <- paste(out.prefix, "_trees/tree_", central.motif, ".json", sep = "")
#              verbose(paste("JSON tree file", json.file), 2)
#              writeLines(JSON.single.node, con=json.file)
#              
#              ## For consistency, print the empty file
#              ## It will be erased later
#              JSON.empty <- ";Empty_file\n"
#              JSON.clusters.table.file <- paste(sep = "", cluster.folder, "/levels_JSON_", central.motif, "_table.tab")
#              write.table(JSON.empty, file = JSON.clusters.table.file, sep = "\t", quote = FALSE, row.names = FALSE)
#            }
#            
#            intermediate.levels.counter <<- intermediate.levels.counter + 1
#            intermediate.levels[intermediate.levels.counter] <<- paste(intermediate.levels.counter, nb, 0, "NO_FILE", sep = "\t")
#          },
#          
#          ####################################################################
#          ## Align the internal cluster if they have more than a single node
#          ## NOTE: for each cluster found previously, its herarchical tree
#          ## is computed and each step of the alignment is exported in order
#          ## to create the branch-motifs (using perl)
#          case.2 = {
#            
#            ## New comparison table (with the ids of the current cluster)
#            ## Use dplyr::filter, this function is faster for large tables and cleaner code
#            ## Added by Jaime Castro: 25-02-2020
#            global.compare.matrices.table <<- compa.table[which( (compa.table[,"id1"] %in% ids) & (compa.table[,"id2"] %in% ids) ),]
#            # global.compare.matrices.table <<- compa.table %>% 
#            #                                    dplyr::filter(id1 %in% ids & id2 %in% ids)
#            
#            global.compare.matrices.table$id1 <<- as.vector(global.compare.matrices.table$id1)
#            global.compare.matrices.table$id2 <<- as.vector(global.compare.matrices.table$id2)
#            
#            ## New description table (with the ids of the current cluster)
#            global.description.table <<- desc.tab[desc.tab[,"id"] %in% ids, ]
#            global.description.table <<- global.description.table[order(as.vector(global.description.table$id)),]
#            global.description.table$n <<- 1:length(global.description.table$n)
#            
#            ## Convert distance table into a distance matrix, required by hclust
#            distances.objects <- build.distance.matrix(metric = metric)
#            dist.matrix <- distances.objects$matrix
#            
#            ## Calculate the central motif
#            ## This will be used to rename the clusters
#            mean.dist.per.motif <- apply(distances.objects[[1]], 1, mean)
#            central.motif <- names(which.min(mean.dist.per.motif)[1])
#            central.motif <- get.name(central.motif)
#            central.motif<- gsub("\\(", "_", central.motif)
#            central.motif<- gsub("\\)", "_", central.motif)
#            central.motif<- gsub("\\.", "_", central.motif)
#            central.motif <- paste(central.motif, "cluster", nb, sep = "_")
#            central.motif <- paste("cluster", nb, sep = "_")
#            
#            ## Fill the elements of the 'central motif table'
#            c.m <- names(which.min(mean.dist.per.motif)[1])
#            central.motif.IDs.cluster[nb] <<- paste("cluster", nb, sep = "_")
#            central.motif.names[nb] <<- get.name(c.m)
#            central.motif.IDs[nb] <<- c.m
#            
#            ## Creates an individual folder for each cluster
#            cluster.folder <<- file.path(clusters.info.folder, central.motif)
#            dir.create(cluster.folder, recursive=TRUE, showWarnings=FALSE)
#            
#            ## Build the tree by hierarchical clustering,
#            tree <<- hclust.motifs(dist.matrix, hclust.method = hclust.method)
#            
#            if (only.hclust == 0) {
#              
#              ## Creates and export the json file
#              JSON.tree <- convert.hclust.to.JSON(tree)
#              
#              json.file <- paste(out.prefix, "_trees/tree_", central.motif, ".json", sep = "")
#              verbose(paste("JSON tree file", json.file), 2)
#              writeLines(JSON.tree, con = json.file)
#              
#              tree.file.rdata <-  paste(out.prefix, "_trees/tree_", central.motif, ".Rdata", sep = "")
#              save(tree, file = tree.file.rdata)
#              
#              ## Creates a file indicating to which levels of the JSON tree correspond to the levels on the hclust tree
#              ## This step is required to assign a name to the branches in the JSON tree in order to create
#              ## the branch-motifs
#              JSON.clusters.table <- identify.JSON.tree.branches(tree)
#              JSON.clusters.table.file <- paste(sep = "", cluster.folder, "/levels_JSON_", central.motif,"_table.tab")
#              write.table(JSON.clusters.table, file = JSON.clusters.table.file, sep = "\t", quote = FALSE, row.names = FALSE)
#              
#              nodes <- as.vector(JSON.clusters.table$node)
#              nodes <- as.numeric(gsub("node_", "", nodes))
#              
#              ## Export the tree agglomeration order
#              tree.agg <- as.vector(tree[[1]])
#              
#              tree.agg.tab <- sapply(tree.agg, function(x) {
#                
#                if (x < 0) {
#                  new.x <- x*-1
#                  as.vector(global.description.table$id)[new.x]
#                } else {
#                  paste("node_", x, sep = "")
#                }
#              })
#              tree.agg.tab <- as.data.frame(matrix(tree.agg.tab, ncol = 2))
#              # tree.agg.tab$new <- rev(as.vector(JSON.clusters.table$node))
#              tree.agg.tab$order <- paste("node_", 1:nrow(tree.agg.tab), sep = "")
#              colnames(tree.agg.tab) <- c("child_1", "child_2", "merged_ID")
#              tree.agg.tab <- tree.agg.tab[,c("merged_ID", "child_1", "child_2")]
#              JSON.clusters.order.table.file <- paste(sep = "", cluster.folder, "/levels_JSON_", central.motif,"_table_linkage_order.tab")
#              write.table(tree.agg.tab, file = JSON.clusters.order.table.file, sep = "\t", quote = FALSE, row.names = FALSE)
#            }
#            
#            ## Align the motifs and retrieve the information of the intermediate alignments
#            if (radial.only == 1) {
#              thresholds <- list(Ncor = -1, cor = -1, w = 0)
#            }
#            
#            ##############################################
#            ## Forest: align the motif within each tree ##
#            ##############################################
#            verbose(paste("; Aligning each cluster individually"), v)
#            alignment.cluster <<- align.motifs(thresholds = thresholds,
#                                               method = hclust.method,
#                                               metric = metric,
#                                               align = TRUE,
#                                               nodes.attributes = TRUE,
#                                               intermediate.alignments = TRUE)
#            intern.alignment <- alignment.cluster$intermediate.alignments
#            
#            ## Export the table with the intermediates alignment information
#            sapply(names(intern.alignment), function(lev) {
#              
#              node.name <- gsub("merge_level", "node", lev, perl = TRUE)
#              level.info <- data.frame(t(intern.alignment[[node.name]]))
# 
#              # NOTE WSG. All columns in DF need to be converted from type
#              # list to character, as this leads to error in write.table in
#              # R 4.1.0
#              level.info <- apply(level.info, 2,unlist) %>% 
#                 as.data.frame()
# 
#              f <- paste(cluster.folder, "/levels_JSON_", central.motif, "_", node.name, "_dataframe.tab", sep = "")
#              write.table(level.info, file = f, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
#              create.dir.merge(node.name)
#            })
#            forest.list[[central.motif]] <<- alignment.cluster$motifs.alignment
#            all.central.motifs <<- append(all.central.motifs, central.motif)
#          }
#   )
# })
# 
# #gsub("^\\w+_(cluster_\\d+)", "\\1", all.central.motifs, perl = TRUE)
# 
# all.central.motifs <- as.vector(unlist(all.central.motifs))
# all.central.motifs.ids.df <-data.frame(all.central.motifs, 1:length(all.central.motifs))
# write.table(all.central.motifs.ids.df, file = paste(sep = "", out.prefix, "_cluster_IDs.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# ## Print a file with the Hexadecimals code for the colors of the clusters
# ## The color of the clusters showed in the heatmap will be the same
# ## colors in the D3 trees.
# ##
# ## Print as well a table with the central motif on each cluster
# if (only.hclust == 0) {
#   
#   n.colors <- NULL
#   if (radial.only == 1) {
#     n.colors <- original.number.clusters
#   } else {
#     n.colors <- length(clusters)
#   }
#   # colors <- rainbow(n.colors)
#   colors <- colorRampPalette(brewer.pal(8, "Dark2"), space="Lab")(n.colors)
#   
#   write.table(paste(paste("cluster_", 1:n.colors, sep = ""), " ", colors, sep = ""), file = paste(sep = "", out.prefix, "_hexa_colors.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
#   
# 
# }
# ## Create and export DF of central motifs
# central.motif.table <- data.frame(central.motif.IDs.cluster, central.motif.IDs, central.motif.names)
# colnames(central.motif.table) <- c("cluster", "ID", "name")
# write.table(central.motif.table, file = paste(sep = "", out.prefix, "_central_motifs_IDs.tab"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
# 
# 
# #################################
# ## Produce the alignment table
# if (forest.nb > 1) {
#   alignment.table <- sapply(forest.list, function(X) {
#     sapply(X, function(Y) {
#       return(c(Y[["number"]], Y[["strand"]], Y[["spacer.up"]], Y[["spacer.dw"]], Y[["consensus_d"]], Y[["consensus_rc"]], Y[["name"]]))
#     })
#   })
# } else{
#   alignment.table <- lapply(forest.list[[1]], function(X) {
#     return(c(X[["number"]], X[["strand"]], X[["spacer.up"]], X[["spacer.dw"]], X[["consensus_d"]], X[["consensus_rc"]], X[["name"]]))
#   })
# }
# alignment.table <- unlist(alignment.table)
# names(alignment.table) <- NULL
# alignment.table <- data.frame(matrix(alignment.table, ncol = 7, byrow = TRUE))
# colnames(alignment.table) <- c("number", "strand", "spacer.up", "spacer.dw", "consensus", "consensus_rc", "name")
# 
# ## Produce the column ID
# ids.names <- unlist(as.vector(sapply(forest.list, function(x) {names(x)})))
# names(ids.names) <- ids.names
# alignment.table$id <- ids.names
# 
# ## Produce the column Width
# width.tmp <- unlist(sapply(forest.list, function(X) {
#   sapply(X, function(Y) {
#     return( nchar(Y[["consensus_d"]]))
#   })
# }))
# width.tmp <- as.vector(width.tmp)
# alignment.table$width <- width.tmp
# 
# ## Produce the column Forest_ID
# forest.names <- names(forest.list)
# forest.id <- vector()
# for (name in forest.names) {
#   forest.id <- append(forest.id, rep(name, length(forest.list[[name]])))
# }
# alignment.table$cluster <- forest.id
# 
# ##  Re-order the table and export it
# alignment.table <- alignment.table[,c(8, 7, 10, 2:4, 9, 5:6)]
# colnames(alignment.table) <- c("#id", "name", "cluster", "strand", "offset_up", "offset_down", "width", "aligned_consensus", "aligned_consensus_rc")
# 
# alignment.file <- paste(sep = "", out.prefix, "_tables/alignment_table.tab")
# write.table(alignment.table, file = alignment.file, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# ## Print the table with the intermediate alignment (it wll be used in the perl code to create the branch-motifs)
# write.table(intermediate.levels, file = paste(out.prefix, "_tables/intermediate_alignments.tab", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
# 
