## Given the aggregation matrix (i.e., hclust_tree$merge object) and a level in the tree
## return the next (top) level  
find.next.levels.in.tree <- function(x          = NULL,
                                     tree.merge = NULL) {
  
  ## In case when the input value corresponds to the root level (nb nodes - 1)
  if (x == length(tree.merge)/2) {
    return(x)
    
  } else {
    
    ## Get the level
    level <- which(tree.merge == x)
    
    if (level > length(tree.merge)/2 & level <= length(tree.merge)) {
      return(level - length(tree.merge)/2)
    } else if (level <= length(tree.merge)/2) {
      return(level)
    }
  }
}



################################################################
## Collect the list of leaves (motifs) associated to each
## internal node of the tree.
##
## This method works for any hclust result tree, irrespective of the
## fact that it refers to motifs or anything else.
##
## Usage: generate a vector with the list of leaves (string-formatted)
## per internal node
##   treenodes <- leaves.per.node(tree)
#
## The input tree must be an hclust result
##
## Then to get a vector with the leaves associated to a given internal
## node (e.g. node 4):
##   leaves.for.node.4 <- as.numeric(unlist(strsplit(tree.nodes[4], split=" ")))
##
leaves.per.node <- function(tree = NULL) {
  
  merge.table <- tree$merge
  leave.lists <- list()
  
  for (i in 1:nrow(merge.table)) {
    branch1 <- merge.table[i, 1]
    branch2 <- merge.table[i, 2]
    
    ## Depending on whether the left branch points to a leave or an
    ## internal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch1 < 0) {
      nodes1 <- -branch1 ## branch one only contains one leave
    } else {
      nodes1 <- leave.lists[branch1]
    }
    
    ## Depending on whether the right branch points to a leave or an
    ## internal nodes, collect a single leave or the pre-defined list
    ## of leaves from this internal node
    if (branch2 < 0) {
      nodes2 <- -branch2 ## branch two only contains one leave
    } else {
      nodes2 <- leave.lists[branch2]
    }
    
    leave.lists[i] <- paste(nodes1, nodes2)
  }
  
  ## Transform the list to export it
  if (length(leave.lists) > 1) {
    leave.lists <- sapply(leave.lists, function(x){ as.numeric(unlist(strsplit(x, " "))) })
  } else{
    leave.lists <- list(as.numeric(unlist(strsplit(leave.lists[[1]], " "))))
  }
  return(leave.lists)
}



## Given a vector of motif IDs, return all entries in the comparison table (compa)
## with those motifs, removing repeated entries
##
## The argument 'hash.ind' must be a list where each element is a vector with the
## indexes of a motif in the comparison table
# motif.comparison.entries <- function(ids      = NULL,
#                                      compa    = NULL,
#                                      hash.ind = NULL,
#                                      full     = TRUE,
#                                      self     = FALSE) {
#   
#   ids               <- as.vector(unlist(ids))
#   comparison.subset <- compa
#   
#   
#   if (!all(ids %in% names(hash.ind))) {
#     
#     ids.indexes.all      <- sort(as.vector(unlist(hash.ind[ids])))
#     ids.indexes.selected <- ids.indexes.all[duplicated(ids.indexes.all)]
#     
#     comparison.subset <- comparison.subset[ids.indexes.selected] %>% 
#                             dplyr::mutate(Nb_compa = 1:n())
#     
#   } else {
#     comparison.subset <- comparison.subset %>% 
#                           dplyr::mutate(Nb_compa = 1:n())
#   }
#   
#   
#   ## Remove diagonal (motif compared to itself is always cor/Ncor = 1)
#   # if (self == FALSE) {
#   #   comparison.subset <- comparison.subset %>%
#   #     dplyr::filter(id2 != id1)
#   # }
#   
#   
#   ## Select relevant motifs + create comparison pair IDs
#   # comparison.subset <- comparison.subset %>%   
#   #   dplyr::filter(id1 %in% ids) %>% 
#   #   dplyr::filter(id2 %in% ids) %>% 
#   #   # dplyr::filter(id2 %in% ids & id1 %in% ids) %>%  ## Select table entries associated with the input motif IDs
#   #   dplyr::mutate(Nb_compa = 1:n())
#   
#   
#   ## Not full returns only the first instance of each comparison pair ID 
#   if (full == FALSE) {
#     comparison.subset <- comparison.subset  %>% 
#       group_by(Nb_compa) %>% 
#       # arrange(Compa_ID) %>% 
#       slice(1) %>%                                  ## To avoid lines with redundant information, take only the first time a comparison identifier appears
#       select(cor, Ncor, Nb_compa)
#   }
#   
#   
#   comparison.subset %>% 
#     ungroup()
# }
motif.comparison.entries <- function(ids      = NULL,
                                     compa    = NULL,
                                     hash.ind = NULL,
                                     full     = TRUE,
                                     self     = FALSE) {

  ids               <- as.vector(unlist(ids))
  comparison.subset <- compa



  ## Remove diagonal (motif compared to itself is always cor/Ncor = 1)
  if (self == FALSE) {
    comparison.subset <- comparison.subset %>%
                            dplyr::filter(id1 != id2)
  }


  ## Select relevant motifs + create comparison pair IDs
  comparison.subset <- comparison.subset %>%
                          ungroup() %>% 
                          dplyr::filter(id1 %in% ids) %>%
                          dplyr::filter(id2 %in% ids) %>%
                          # dplyr::filter(id2 %in% ids & id1 %in% ids) %>%  ## Select table entries associated with the input motif IDs
                          dplyr::mutate(Nb_compa = 1:n()) %>% 
                          group_by(Nb_compa) %>% 
                          mutate(Compa_ID = paste(sort(c(id1, id2)), collapse = ",")) %>% 
                          ungroup()

  ## Not full returns only the first instance of each comparison pair ID
  if (full == FALSE) {
    comparison.subset <- comparison.subset  %>%
                            group_by(Compa_ID) %>%
                            slice(1) %>%                ## To avoid lines with redundant information, take only the first time a comparison identifier appears
                            ungroup() %>% 
                            select(cor, Ncor, w)
                            
  }

  comparison.subset %>%
      ungroup() 

  # return(comparison.subset)
}


## Return the indexes of each motif ID in the comparison table
index.compa.table <- function(comparison.table = NULL) {
  
  ## Code taken from: https://stackoverflow.com/questions/22993637/efficient-r-code-for-finding-indices-associated-with-unique-values-in-vector
  id1.index <- as.data.table(comparison.table$id1)[, list(list(.I)), by = comparison.table$id1]
  setattr(id1.index$V1, 'names', id1.index$comparison.table)
  id1.index <- id1.index$V1
  
  id2.index <- as.data.table(comparison.table$id2)[, list(list(.I)), by = comparison.table$id2]
  setattr(id2.index$V1, 'names', id2.index$comparison.table)
  id2.index <- id2.index$V1
  

  index.hash <- Map(c, id1.index, id2.index)
  
  index.hash <- lapply(index.hash, function(x){
                  sort(unique(x))
                })
  
  return(index.hash)
}



#########################################################################
## Identify the "central"/"distal" leaf (motif) of a subtree (cluster)
## i.e. the pair of motif with the smallest/largest distance to
## among all the other motifs of these clusters.

# ids1 = "oct4_peak_motifs_m16_local_words_8nt_m1"
# ids2 = "oct4_peak_motifs_m17_local_words_8nt_m2"

## ids1, ids2 : motif IDs
## compa.tab  : motif comparison table
## metric     : comparison metric
## closest    : logical, should the function returns the closest comparison pair?
## compa.info : return the line corresponding to the highlighted comparison
closest.or.farthest.motifs.ids <- function(ids1       = NULL,
                                           ids2       = NULL,
                                           compa.tab  = NULL,
                                           metric     = "Ncor",
                                           closest    = TRUE,
                                           compa.info = FALSE){
  
  compa.subset <- compa.tab %>% 
                    dplyr::filter(id2 != id1) %>% 
                    dplyr::filter(id1 %in% ids1 & id2 %in% ids2)
  
  
  #########################################
  ## For correlations (Higher is better)
  cor.metrics <- c("Ncor", "cor")
  if (metric %in% cor.metrics) {
    if (closest) {
      compa.subset <-   compa.subset %>% 
                            ungroup() %>% 
                            dplyr::filter(!!sym(metric) == max(!!sym(metric))) %>% 
                            slice(1)
    } else {
      compa.subset <-   compa.subset %>% 
                            ungroup() %>% 
                            dplyr::filter(!!sym(metric) == min(!!sym(metric))) %>% 
                            slice(1)
    }
  }
  
  
  ## When the comparison information is not required, then returns the ids of the closest interaction
  if (!compa.info) {
    compa.subset <- compa.subset %>% 
                      select(id1, id2)
  }
  
  return(compa.subset)
}



## Calculate the mean Ncor and cor of each motif against the remaining ones
## Then calculate a single score multiplying the resulting mean Cor * mean Ncor
## The motif with the highest score is the one with overall highest similarity in
## the cluster
cluster.centroid.id <- function(compa.tab = cluster.compa.entries) {
  
  as.vector(unlist(compa.tab %>% 
                     group_by(id1) %>% 
                     summarise(Ncor    = median(Ncor),
                               cor     = median(cor),
                               .groups = "drop") %>% 
                     mutate(Score = Ncor * cor) %>% 
                     dplyr::filter(Score == max(Score)) %>% 
                     slice(1) %>% 
                     select(id1)))
  
}



## Given the cluster-id table
## Generate a dataframe with one color associated to each cluster
## Colors are generated using the 'Safe' palette from rcartcolors
cluster.color.map <- function(cluster.tab = NULL,
                              seed        = 130290) {

  suppressPackageStartupMessages(library("rcartocolor", character.only = TRUE, quietly = TRUE))
  message("; Generating cluster-color table")
  
  
  clusters.names <- unique(cluster.tab$cluster)
  nb.clusters    <- length(clusters.names)
  
  
  ## Cluster colors
  nb.classes.palette <- 12  ## We want to use the first 12 colors of the Safe palette
  nb.seed.colors     <- ifelse(nb.clusters < nb.classes.palette,
                               yes = nb.clusters,
                               no  = nb.classes.palette)
  
  ## Generate a carto palette and expand it to the number of clusters
  carto.pal.classes <- carto_pal(nb.seed.colors, "Safe")
  class.colors      <- colorRampPalette(carto.pal.classes, space = "Lab")(nb.clusters)
  
  ## Cluster-color map
  ## Randomize the color values
  set.seed(seed)
  cluster.color.map <- data.table(cluster = clusters.names,
                                  color   = sample(class.colors))
  
  return(cluster.color.map)
}



####################################################
## Draw a heatmap usig a particular agglomeration
## method.
draw.heatmap.motifs <- function(dist.table    = NULL, 
                                clusters      = NULL, 
                                metric        = "Ncor", 
                                cluster.color = NULL,
                                color.palette = "RdGy",
                                color.classes = 11) {
  
  ## List of packages to install from CRAN
  required.packages = c("circlize",      
                        "ComplexHeatmap",
                        "dendsort",
                        "RColorBrewer")
  
  for (lib in required.packages) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
  }
  
  
  ## Supported metrics
  supported.metrics <- c("cor",    ## cor    Pearson correlation (computed on residue occurrences in aligned columns)
                         "Ncor")   ## Ncor   Relative width-normalized Pearson correlation
  
  ## Define type of similarity metric
  ## This is required to prepare color palettes
  metric.definition <- NULL
  if (metric %in% supported.metrics) {
    metric.definition <- "correlation"
  } else {
    stop("; ", metric, " is not a supported metric. Supported metrics: ", paste(supported.metrics, collapse = ","))
  }
  
  
  
  ## Motif/clusters info
  id.per.cluster.tab <- clusters %>% 
    select(cluster, id) %>% 
    separate_rows(id, sep = ",")
  
  clusters.names <- unique(id.per.cluster.tab$cluster)
  nb.clusters    <- length(clusters.names)
  motif.ids      <- as.vector(unlist(id.per.cluster.tab %>%  select(id)))
  nb.motifs      <- length(motif.ids)
  
  ## Heatmap aesthetics
  col.row.title  <- paste0(nb.motifs, " motifs; ", nb.clusters, " clusters")
  row.col.order  <- dendsort(results.list$All_motifs_tree)
  
  
  ## Heatmap cell colors, depending in the similarity metric
  ## Correlation goes from -1 to +1
  ## Distance goes from 0 to infinite
  if (metric.definition == "correlation") {
    
    # col_fun = colorRamp2(c(-1, 0, 1), c("white", "#49006a","white"))
    palette <- rev(colorRampPalette(rev(brewer.pal(color.classes, color.palette)), space = "Lab")(color.classes))
    col_fun <- colorRamp2(seq(-1, 1, length.out = color.classes), palette)
    
  }
  ## TO DO: create color palette for distances (a sequential palette)
  
  
  
  ##############################
  ## Draw heatmap annotations ##
  ##############################
  row.col.order.w.cluster   <- data.table(id = results.list$All_motifs_tree$labels) %>% 
    left_join(id.per.cluster.tab, by = "id") %>% 
    left_join(cluster.color, by = "cluster")
  
  cluster.color.list        <- as.vector(cluster.color$color)
  names(cluster.color.list) <- cluster.color$cluster
  
  
  
  hac = HeatmapAnnotation(df               = row.col.order.w.cluster[,c("cluster")],
                          col              = list(cluster = cluster.color.list),
                          # gp             = gpar(col = "black"),
                          simple_anno_size = unit(1, "cm"),
                          which            = "column",
                          show_legend      = FALSE, show_annotation_name = FALSE)
  
  har = HeatmapAnnotation(df               = row.col.order.w.cluster[,c("cluster")],
                          col              = list(cluster = cluster.color.list),
                          # gp             = gpar(col = "black"),
                          simple_anno_size = unit(1, "cm"),
                          which            = "row",
                          show_legend      = FALSE, show_annotation_name = FALSE)
  
  
  
  ht1 = Heatmap(matrix            = dist.table,
                name              = metric,
                cluster_rows      = row.col.order,
                cluster_columns   = row.col.order,
                show_row_dend     = FALSE,
                # rect_gp           = gpar(col = "white", lwd = 0.75),
                column_title      = col.row.title, column_title_side = "bottom",
                row_title         = col.row.title,
                col               = col_fun,
                show_row_names    = FALSE,
                show_column_names = FALSE,
                top_annotation    = hac,
                left_annotation   = har)
  # dev.off()
  
  return(ht1)
}
