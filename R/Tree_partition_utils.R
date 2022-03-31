## Return the function that will be applied: mean (average linkage), min (single linkage),
## or max (complete linkage)
function.to.apply <- function(parameters = NULL){
  
  switch(params.list$linkage_method,
         average  = mean,
         single   = min,
         complete = max)
}



nodes.summary.stats <- function(comparison.list = NULL,
                                parameters      = NULL,
                                tree            = NULL) {
  
  ## Get the function that will by applied in the purrr::map function
  fun.apply <- function.to.apply(parameters = parameters)
  
  ## For each node, create a summary statistics
  stats.summary.per.node <- data.table(cor  = purrr::map(purrr::map(comparison.list, 1), fun.apply, na.rm = TRUE),
                                       Ncor = purrr::map(purrr::map(comparison.list, 2), fun.apply, na.rm = TRUE),
                                       w    = purrr::map(purrr::map(comparison.list, 3), fun.apply, na.rm = TRUE))
  
  
  ## Evaluate the correlation values in each node, return a TRUE/FALSE value depending
  ## whether the motifs in the node satisfy the threshold(s)
  stats.summary.per.node <- stats.summary.per.node %>% 
                              mutate(Merge = ifelse((cor  >= parameters$cor & Ncor >= parameters$Ncor) & w >= parameters$w,
                                                    yes = TRUE,
                                                    no  = FALSE))
  
  ## Add the tree linkage matrix to each node
  stats.summary.per.node <- cbind(stats.summary.per.node, data.frame(tree$merge)) %>% 
                                dplyr::rename(Merge1 = X1,
                                              Merge2 = X2)
  
  
  ## Update the merge status ('Merge' column) to set as FALSE all the parent nodes
  ## of those nodes that will not be merge
  ## comparisons.per.node$Merge[18] <- TRUE ## For debugging
  stats.summary.per.node.updated <- updated.parents.node.merge.status(tab  = stats.summary.per.node,
                                                                      tree = tree)
  
  return(stats.summary.per.node.updated)
}



## Given the aggregation matrix (i.e., hclust_tree$merge object) and a level in the tree
## return the next (top) level  
find.node.parents <- function(x          = NULL,
                              tree.merge = NULL) {
  
  node.parents <- NULL
  current.node <- NULL
  merge.flag   <- 1
  node.parents <- append(node.parents, x)
  
  while (merge.flag == 1) {
    
    ## Search the next chained level
    current.node <- find.next.levels.in.tree(x = x, tree.merge = tree.merge)
    
    ## The case when the root tree is analyzed
    if (x == current.node) {
      merge.flag <- 0
    }
    
    ## If the status at the current level is 0 then return the chained levels
    if (merge.flag == 0) {
      return(node.parents)
      
      ## Conversely, add the new level to the chain
    } else {
      node.parents <- append(node.parents, current.node)
      x <- current.node
    }
  }
}



## Given a logical vector indicating whether children in a tree node should be merged
## update such vector by setting as FALSE all parental nodes of the original FALSE
## values in the logical vector.
##
## This is required to avoid merging parent nodes that by chance satisfy the thresholds.
updated.parents.node.merge.status <- function(tab  = NULL,
                                              tree = NULL){
  
  ## Position in the vector that are not merged
  split.branches.ind <- which(!tab$Merge)
  
  ## These indexes should be set to FALSE, corresponding to the parrents of the 
  ## nodes that are not merged
  index.false <- sort(unique(unlist(purrr::map(split.branches.ind, find.node.parents, tree.merge = tree$merge))))
  
  tab$Merge[index.false] <- FALSE 
  
  tab
}



## Return the merge status (should this node be merged? Logical vector) of a particular
## node. The merge status vector must be provided
get.merge.status <- function(node.nb      = NULL,
                             merge.status = NULL) {
  
  ## 0 means that the node is a leaf
  if (node.nb == 0) {
    return(TRUE)
  } else {
    merge.status[node.nb]
  }
}


## Using a Merge status vector (provided in the merge.status.table),
## Return a vector with the elements corresponding to the clusters (the same notation as in the tree$merge object)
## Negative values indicate leaves, therefore these will be clusters of size 1 (a.k.a. singletons)
## Positive values correspond to nodes with two or more elements 
tree.partition <- function(tree               = NULL,
                           merge.status.table = NULL){
  
  nb.nodes  <- nrow(tree$merge)
  nodes.seq <- seq_len(nb.nodes)
  
  ## Return the node parents of each leaf/node in the tree
  node.parents <- purrr::map(seq_len(nb.nodes), find.node.parents, tree.merge = tree$merge)
  
  ## Identify the indexes (nodes) with a Merge status == FALSE (these are the nodes
  ## that will not be merged and therefore where the partition will occur)
  ## These indexes are selected to avoid iteration through the entire tree
  false.merge.ind <- which(merge.status.table$Merge == FALSE)
  

  ## This part of the code runs when there is at least one partition in the tree
  if (length(false.merge.ind) > 0) {
    
    ## A subset of the merge.status.table table only including those nodes with 
    ## FALSE in the 'Merge' column
    non.merged.nodes.tab        <- merge.status.table[false.merge.ind,]
    non.merged.nodes.tab$Merge1 <- ifelse(non.merged.nodes.tab$Merge1 < 0, yes = 0, no = non.merged.nodes.tab$Merge1)
    non.merged.nodes.tab$Merge2 <- ifelse(non.merged.nodes.tab$Merge2 < 0, yes = 0, no = non.merged.nodes.tab$Merge2)
    
    ## Obtain the Merge status (1 or 0) of the children of each node previously
    ## identified with a Merge == FALSE status
    ## Leaves were marked with 0s in the tree structure
    non.merged.nodes.tab <- non.merged.nodes.tab %>% 
                              cbind(data.table(Merge1_child_status = as.numeric(purrr::map_lgl(non.merged.nodes.tab$Merge1, get.merge.status, merge.status = merge.status.table$Merge)),
                                               Merge2_child_status = as.numeric(purrr::map_lgl(non.merged.nodes.tab$Merge2, get.merge.status, merge.status = merge.status.table$Merge)),
                                               Ind                 = false.merge.ind))
    
    ## Complete the table with the missing nodes (those were partitions do not occur)
    ## The children status of all missing values are set to FALSE
    missing.nodes <- nodes.seq[!nodes.seq %in% non.merged.nodes.tab$Ind]
    missing.ind   <- data.table(Merge1_child_status = rep(0, length(missing.nodes)),
                                Merge2_child_status = rep(0, length(missing.nodes)),
                                Ind                 = missing.nodes)
    
    ## This is the entire table with the Merging status of all nodes where partitions
    ## occur, the children of those nodes within a partition are all set to 0.
    ## The idea is to only focus in those positions in the tree where particion occurs
    ## assuming that all nodes above those partitions were not be merged (this is way be
    ## previously ran the updated.parents.node.merge.status function).
    non.merged.nodes.tab <- bind_rows(non.merged.nodes.tab, missing.ind) %>% 
                              arrange(Ind) %>% 
                              distinct() %>% 
                              select(Merge1_child_status, Merge2_child_status, Ind)
    
    
    merge.status.table <- cbind(merge.status.table, non.merged.nodes.tab) %>% 
      mutate(Partition1 = Merge1 * Merge1_child_status,
             Partition2 = Merge2 * Merge2_child_status)
    
    ## Returns an array with the entries (nodes or leaves) where the tree is partitioned
    ## Negative values indicate leaves, therefore these will be clusters of size 1 (a.k.a. singletons)
    ## Positive values correspond to nodes with two or more elements
    entries.with.clusters <- unique(c(merge.status.table$Partition1, merge.status.table$Partition2))
    entries.with.clusters <- sort(entries.with.clusters[entries.with.clusters != 0], decreasing = TRUE)
  
  ## When the tree is not partitioned return the number of the last node (i.e., the root)
  } else {
    entries.with.clusters <- nb.nodes
  }
  
  return(entries.with.clusters)
}



## This function returns a list with two elements:
## 1) clusters           : a list where each element represents a cluster, it contains the motif IDs
## 2) agglomeration.stats: the summary stats at each node in the tree and its merge status
find.motif.clusters <- function(tree             = NULL,
                                comparison.table = NULL,
                                parameters       = NULL){
  
  ## A list containing the leaves associated to each node of the tree
  motif.at.tree.level <- leaves.per.node(tree)
  
  ## A list containing the motif IDs associated to each label
  motif.IDs.per.tree.level <- lapply(motif.at.tree.level, function(x){
                                tree$labels[x]
                              })

  motif.ids.ind <- index.compa.table(comparison.table = comparison.table)

  ## Obtain the comparison entries among the motifs within each node in the tree
  ## Each node is stored as a dataframe nested in a list 
  ##
  ## NOTE: depending in the number of motifs, this function may be very slow
  plan(multisession, workers = parameters$nb_workers)
  comparisons.per.node <- furrr::future_map(.x = motif.IDs.per.tree.level, 
                                            .f = ~motif.comparison.entries(compa    = comparison.table,
                                                                           ids      = .x,
                                                                           hash.ind = motif.ids.ind,
                                                                           full     = FALSE,
                                                                           self     = FALSE))

  ## Return the Merge status of each node in the tree
  ## FALSE correspond to the nodes where the tree will be partitioned
  message(";     Calculating similarity among elements within each node")
  nodes.merge.status.tab <- nodes.summary.stats(comparison.list = comparisons.per.node,
                                                parameters      = parameters,
                                                tree            = tree)
  
  
  ## Return the values in the tree$merge object corresponding to the clusters
  message(";     Finding nodes where hierarchical tree will be partitioned")
  entries.with.clusters <- tree.partition(tree               = tree,
                                          merge.status.table = nodes.merge.status.tab)
  
  
  ## Get the IDs of the motifs in each cluster
  message(";     Associating tree nodes with motif IDs")
  clusters.list <- get.motifs.ids.in.clusters(tree                      = tree,
                                              clusters.leaves.nodes.ind = entries.with.clusters,
                                              motif.ids.per.node        = motif.IDs.per.tree.level)
  
  clusters.df <- clusters.list.to.df(clusters.list)
  
  
  return(list(clusters            = clusters.list,
              clusters_df         = clusters.df,
              agglomeration.stats = nodes.merge.status.tab,
              compa.per.node      = comparisons.per.node))
}



## Obtain the motifs IDs in each cluster
##
## tree                      : an hclust object
## clusters.leaves.nodes.ind : a vector indicating the leave numbers (negative ones) and node numbers (positives), each one corresponds to a cluster
## motif.ids.per.node        : a list containing the motifs IDs at each node in the hierarchical tree (output from leaves.per.node)
get.motifs.ids.in.clusters <- function(tree                       = NULL,
                                       clusters.leaves.nodes.ind  = NULL,
                                       motif.ids.per.node         = NULL) {
  
  ## Remove entries with 0s and sort the vector to ease the iteration
  clusters.leaves.nodes.ind <- sort(clusters.leaves.nodes.ind[clusters.leaves.nodes.ind != 0], decreasing = TRUE)
  
  ## Separate clusters (singletons: negative index; groups: positve index)
  clusters.leaves.nodes.ind.singleton <- clusters.leaves.nodes.ind[clusters.leaves.nodes.ind < 0]
  clusters.list.singleton             <- as.list(tree$labels[-clusters.leaves.nodes.ind.singleton])
  
  
  clusters.leaves.nodes.ind.group <- clusters.leaves.nodes.ind[clusters.leaves.nodes.ind > 0]
  clusters.list.group             <- motif.ids.per.node[clusters.leaves.nodes.ind.group]
  
  ## Combine singleton and group lists
  clusters.list        <- c(clusters.list.group,
                            clusters.list.singleton)
  
  ## Add a cluster name, add leading zeros in the cluster name
  nb.char.cluster.ids  <- nchar(as.character(length(clusters.list))) ## Number of digits
  sprint.string        <- paste0("%0", nb.char.cluster.ids, "d")     ## Sprintf options  
  names(clusters.list) <- paste0("cluster_", sprintf(sprint.string, seq_len(length(clusters.list))))
  
  clusters.list
}
