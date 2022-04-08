## Add upstrean/downstream gaps '-' to a given consensus (string) sequence
add.gaps.consensus <- function(consensus   = "wwATGCTAAtt",
                               offset.up   = 0,
                               offset.down = 0) {
  
  paste0(paste0(rep("-", times = offset.up), collapse = ""),
         consensus,
         paste0(rep("-", times = offset.down), collapse = ""), sep = "")
  
}



## Add gaps '-' in the downstream side of the alignment
## After this function is applied, all aligned consensuses have the same length
fill.alignment.downstream <- function(alignment.table = NULL,
                                      ids.subset      = NULL) {
  
  ## Select the consensus of the selected IDs
  consensus.ids.subset.tab <- alignment.table %>% 
                                dplyr::ungroup() %>% 
                                dplyr::filter(id %in% ids.subset) %>% 
                                mutate(Update_status = Update_status + 1)
  consensus.ids.subset <- as.vector(unlist(consensus.ids.subset.tab$Oriented_consensus))
  
  ## Calculate the missing gaps at the end of the alignment
  gaps.down.to.add <- abs(nchar(consensus.ids.subset) - max(nchar(consensus.ids.subset)))
  
  consensus.filled.downstream <- purrr::map2_chr(.x = consensus.ids.subset, .y = gaps.down.to.add, .f = ~add.gaps.consensus(consensus = .x, offset.down = .y))
  
  ## Count the gaps at the downstream part of the consensus
  gaps.down <- nchar(gsub(consensus.filled.downstream, pattern = "^-*[a-zA-Z]+", replacement = ""))
  
  consensus.ids.subset.tab$Oriented_consensus <- consensus.filled.downstream
  consensus.ids.subset.tab$Offset_down        <- gaps.down

  consensus.ids.subset.tab <- data.table(consensus.ids.subset.tab)

  return(consensus.ids.subset.tab)
  
}



## Return the hierarchical cluster of a set of motifs
hclust.cluster.ids <- function(ids        = NULL,
                               compa      = NULL,
                               parameters = NULL) {
  
  ## Check whether the input cluster is a singleton
  ## If it is not a singleton, then compute the hierarchical clustering with the input IDs
  if (length(ids) > 1) {
    
    ## Get the entries associated to the motifs in the cluster
    cluster.compa.entries <- motif.comparison.entries(ids   = ids,
                                                      compa = compa,
                                                      full  = TRUE,
                                                      self  = TRUE)
    
    ## Calculate the distance matrix
    distances.objects <- build.distance.matrix(compa.table = cluster.compa.entries,
                                               metric      = parameters$comparison_metric)
    
    ## Generate the hierarchical tree
    cluster.hclust <- hclust.motifs(distances.objects$matrix,
                                    hclust.method = parameters$linkage_method)
    
    return(list(hclust           = cluster.hclust,
                distance_objects = distances.objects))
    
    ## In case of a singleton, return the motif ID
  } else {
    return(ids)
  }
}



## Change the orientation and offset (gaps) of an input set of motif ids
invert.one.motif <- function(motif.line   = NULL,
                             alignment.df = NULL) {
  
  new.line    <- alignment.df[motif.line,]
  old.strand  <- new.line[, "Strand"]
  
  ## Invert offset values
  tmp.offset  <- new.line[, "Offset_up"]
  new.line[, "Offset_up"]   <- new.line[, "Offset_down"]
  new.line[, "Offset_down"] <- tmp.offset
  
  ## Invert strand values
  if (old.strand == "D") {
    new.line[, "Strand"]             <- "R"
    new.line[, "Oriented_consensus"] <- add.gaps.consensus(consensus   = new.line[, "rc_consensus"],
                                                           offset.up   = new.line[, "Offset_up"],
                                                           offset.down = new.line[, "Offset_down"])
  } else if (old.strand == "R") {
    new.line[, "Strand"]             <- "D"
    new.line[, "Oriented_consensus"] <- add.gaps.consensus(consensus   = new.line[, "consensus"],
                                                           offset.up   = new.line[, "Offset_up"],
                                                           offset.down = new.line[, "Offset_down"])
    # new.line[, "Oriented_consensus"] <- new.line[, "consensus"]
  }

  new.line[, "Update_status"] <- new.line[, "Update_status"] + 1
  
  new.line
}



## Use purrr::map to iterate across the motifs that need to be inverted
invert.alignment <- function(ids          = NULL,
                             alignment.df = NULL) {
  
  ## Row number in the table with the motifs of interest
  motif.lines <- which(alignment.df$id %in% ids)
  
  rbindlist(furrr::future_map(motif.lines, ~invert.one.motif(motif.line = .x, alignment.df = alignment.df)))
  
}



## Return the motif IDs according to the agglomeration matrix (tree$merge)
motifs.ids.agglomeration.tab <- function(node1         = NULL,
                                         node2         = NULL,
                                         tree          = NULL,
                                         ids.per.level = NULL) {
  
  node1.ids <- NULL
  node2.ids <- NULL
  
  ## If the number is positive (i.e., indicates a node), all the motifs
  ## in that node will be updated, motifs names are taken from the object motif.IDs.per.tree.level
  if (node1 > 0) {
    node1.ids <- ids.per.level[[node1]]
    
    ## If the number is negative (i.e., indicates a leaf), then the motif ID
    ## is taken directly from the tree$label 
  } else {
    node1.ids <- tree$labels[abs(node1)]
  }
  
  
  if (node2 > 0) {
    node2.ids <- ids.per.level[[node2]]
    
    ## If the number is negative (i.e., indicates a leaf), then the motif ID
    ## is taken directly from the tree$label 
  } else {
    node2.ids <- tree$labels[abs(node2)]
  }
  
  return(list(Node1 = node1.ids,
              Node2 = node2.ids))
}


## Given an x node in the tree, return the leaves within it
## This method uses the tree$merge structure
## Negative values: leaves
## Positive values: a group of leaves
leaves.at.x.node <- function(tree          = NULL,
                             ids.per.level = NULL) {
  
  tree.merge.df <- data.table(tree$merge)
  purrr::map2(.x = tree.merge.df$V1,
              .y = tree.merge.df$V2,
              .f = ~motifs.ids.agglomeration.tab(node1 = .x, node2 = .y, tree = tree, ids.per.level = ids.per.level))
  
}



## Returns a data.frame with the best comparison (most similar pair of motifs)
## of node1 vs node 2 at each level of the tree (provided in the tree$merge object)
best.comparison.per.node <- function(tree  = NULL,
                                     compa = NULL) {
  
  ## A list containing the leaves associated to each node of the tree
  motif.at.tree.level <- leaves.per.node(tree)
  
  ## A list containing the motif IDs associated to each label
  motif.IDs.per.tree.level <- lapply(motif.at.tree.level, function(x){
    tree$labels[x]
  })
  
  tree.merge        <- data.frame(tree$merge)
  motif.id.clusters <- tree$labels
  
  ## Return the motif IDs at each node in the tree, following the agglomeration
  ## order provided in the tree$merge object
  motifs.ids.per.tree.node <- apply(tree.merge, 1, function(tree.merge.row){
    
    ## Initialize a vector returned at each iteration
    evaluated.row <- rep(NA, 2)
    
    
    ## When the number in the agglomeration table is negative, indicates a leaf
    ## otherwise, a node with many leaves
    ## Nodes associated to negative values are obtained from the hclust object
    ## Nodes associated to positive vales are obtained from the motif.IDs.per.tree.level object
    
    ## Operations in column 1
    if (tree.merge.row[1] < 0) {
      evaluated.row[1] <- motif.id.clusters[abs(tree.merge.row[1])]
    } else {
      evaluated.row[1] <- paste(motif.IDs.per.tree.level[[tree.merge.row[1]]], collapse = ",")
    }
    
    
    ## Operations in column 2
    if (tree.merge.row[2] < 0) {
      evaluated.row[2] <- motif.id.clusters[abs(tree.merge.row[2])]
    } else {
      evaluated.row[2] <- paste(motif.IDs.per.tree.level[[tree.merge.row[2]]], collapse = ",") 
    }
    
    evaluated.row
  })
  
  motifs.ids.per.tree.node        <- data.table(t(motifs.ids.per.tree.node))
  motifs.ids.per.tree.node$Agg_ID <- seq_len(nrow(motifs.ids.per.tree.node))
  
  ## Separate the columns containing a list of comma-separated ids, this is required
  ## to ease the pairwise comparison
  motifs.ids.per.tree.node <- motifs.ids.per.tree.node %>%
                                separate_rows(V2, sep = ",") %>% 
                                separate_rows(V1, sep = ",")
  
  ## Separate the data.frame into a list, each element corresponds to a data.frame
  ## with all the comparisons associated to a node in the tree
  motifs.ids.per.tree.node.split <- split(motifs.ids.per.tree.node, f = motifs.ids.per.tree.node$Agg_ID)
  
  ## Each element in the list contains two elements corresponding to the motif IDs
  ## merged in each node in the tree$merge
  motifs.ids.per.tree.node.list <- lapply(motifs.ids.per.tree.node.split, function(l){
    list(N1 = as.vector(unique(l$V1)),
         N2 = as.vector(unique(l$V2)))
  })
  
  
  compa.ids.merge.list <- list(A = purrr::map(motifs.ids.per.tree.node.list, `[[`, "N1"),
                               B = purrr::map(motifs.ids.per.tree.node.list, `[[`, "N2"),
                               C = compa)
  # lapply(compa.ids.merge.list, length)
  
  
  ## Obtain the ids of the closest motif in each node, following the agglomeration order
  closest.ids.each.node <- furrr::future_pmap(.l = compa.ids.merge.list,
                                              .f = ~closest.or.farthest.motifs.ids(ids1       = ..1,
                                                                                   ids2       = ..2,
                                                                                   compa.tab  = ..3,
                                                                                   compa.info = TRUE))
  closest.ids.each.node <- rbindlist(closest.ids.each.node)

  return(list(closest_comparison = closest.ids.each.node,
              motif.id.per.node  = motif.IDs.per.tree.level))
}



#######################################
## Core function for motif alignment ##
#######################################
align.motifs.in.cluster <- function(tree       = NULL,
                                    compa      = NULL,
                                    motif.info = NULL,
                                    parameters = NULL) {
  
  
  ## A list containing the leaves associated to each node of the tree
  motif.at.tree.level <- leaves.per.node(tree)
  
  ## A list containing the motif IDs associated to each label
  motif.IDs.per.tree.level <- lapply(motif.at.tree.level, function(x){
    tree$labels[x]
  })
  

  ## A list where each element is a dataframe with the comparison entries among
  ## the nodes at each level of the input hierarchical tree
  plan(multisession, workers = parameters$nb_workers)
  comparisons.per.tree.node <- furrr::future_map(motif.IDs.per.tree.level, motif.comparison.entries, compa = compa, full = TRUE, self = FALSE)
  # compaaaa <- furrr::future_map(comparisons.per.tree.node, ~closest.or.farthest.motifs.ids(compa.tab = .x, metric = parameters$comparison_metric, closest = TRUE))
 
  
  ## A table with the paired motif comparison among the closest motifs in each
  ## tree node, following the agglomeration order
  ids.per.node.w.best.compa.list <- best.comparison.per.node(tree  = tree,
                                                             compa = comparisons.per.tree.node)
  best.comparison.pair.in.tree.nodes <- ids.per.node.w.best.compa.list$closest_comparison 
  ids.per.tree.node                  <- ids.per.node.w.best.compa.list$motif.id.per.node 
  
  
  ## Initialize alignment table
  ## All motif are in 'D' strand with 0 gaps at both sides
  alignment.tab <- motif.info %>% 
                    dplyr::filter(id %in% tree$labels) %>% 
                    select(id, name, consensus, rc_consensus) %>% 
                    dplyr::mutate(Strand             = "D",
                                  Offset_up          = 0,
                                  Offset_down        = 0,
                                  Oriented_consensus = consensus,
                                  N                  = 1:n(),
                                  Update_status      = 0)
  
  
  ## This list contains N entries (N = agglomeration steps in a tree)
  ## Each entry contains the motifs IDs at each merging step (node1 and node2)
  motif.ids.per.node.list <- leaves.at.x.node(tree          = tree,
                                              ids.per.level = ids.per.tree.node)

  
  # which(sapply(motif.at.tree.level, length) > 60)
  # tree$merge   View(tree$merge )  #380 190  96  48  24  12
  # 669 : 567 + 646
  # 646
  # nn <- 669
  # # l <- 28
  # for (l in 1:nn) {
  for (l in seq_len(nrow(best.comparison.pair.in.tree.nodes))) {
    
    #message("; Aligning cluster, node ", l)
    
    ## Paired Comparison information
    l.id1    <- as.vector(unlist(best.comparison.pair.in.tree.nodes[l, "id1"]))
    l.id2    <- as.vector(unlist(best.comparison.pair.in.tree.nodes[l, "id2"]))
    l.strand <- as.vector(unlist(best.comparison.pair.in.tree.nodes[l, "strand"]))
    l.offset <- as.vector(unlist(best.comparison.pair.in.tree.nodes[l, "offset"]))
    
    
    ## This is the current information about the IDs in the paired comparison
    id1.align.info <- subset(alignment.tab, id == l.id1)
    id2.align.info <- subset(alignment.tab, id == l.id2)
    
    
    ## l.id1, l.id2 are the IDs in the comparison table and may not necessarily
    ## correspond to the order in the agglomeration table
    ## With the next condition we ensure this correspondence exists
    ids.node1 <- motif.ids.per.node.list[[l]]$Node1
    ids.node2 <- motif.ids.per.node.list[[l]]$Node2
    
    if (!l.id1 %in% motif.ids.per.node.list[[l]]$Node1) {
      ids.node1 <- motif.ids.per.node.list[[l]]$Node2
      ids.node2 <- motif.ids.per.node.list[[l]]$Node1
    }
    
    
    
    ## Invert the alignment if needed
    ##
    ## Key   = Id1_orientation.Id2_orientation.Comparison_orientation
    ## D.D.D = ID1 orientations = D + ID2 orientation = D + Comparison orientation = D
    ##
    D.D.D <- id1.align.info$Strand == "D" & id2.align.info$Strand == "D" & l.strand == "D"  ## Nothing happens
    D.R.D <- id1.align.info$Strand == "D" & id2.align.info$Strand == "R" & l.strand == "D"  ## Invert ID2
    R.D.D <- id1.align.info$Strand == "R" & id2.align.info$Strand == "D" & l.strand == "D"  ## Invert ID1
    R.R.D <- id1.align.info$Strand == "R" & id2.align.info$Strand == "R" & l.strand == "D"  ## Invert ID1, ID2
    D.D.R <- id1.align.info$Strand == "D" & id2.align.info$Strand == "D" & l.strand == "R"  ## Invert ID2
    D.R.R <- id1.align.info$Strand == "D" & id2.align.info$Strand == "R" & l.strand == "R"  ## Nothing happens
    R.D.R <- id1.align.info$Strand == "R" & id2.align.info$Strand == "D" & l.strand == "R"  ## Invert ID1, ID2
    R.R.R <- id1.align.info$Strand == "R" & id2.align.info$Strand == "R" & l.strand == "R"  ## Invert ID1
    invert.these.ids <- NULL
    
    if (D.R.D | R.D.D | R.R.D | D.D.R | R.D.R | R.R.R) {
      
      ## Depending in the activated flag, select the set of motifs that will be inverted
      if (R.R.D | D.R.R) {
        # message("; Invert level: ", l, " - Case: 1")
        invert.these.ids <- c(ids.node1, ids.node2)
      } else if (R.D.D | R.R.R) {
        # message("; Invert level: ", l, " - Case: 2")
        invert.these.ids <- ids.node1
      } else if (D.R.D | D.D.R) {
        # message("; Invert level: ", l, " - Case: 3")
        invert.these.ids <- ids.node2
      }
      
      ## Invert the information of the selected motif ids in the alignment table
      inverted.entries <- invert.alignment(ids          = invert.these.ids,
                                           alignment.df = alignment.tab)
      
      ## Update the alignment table
      alignment.tab <- update.alignment.table(rbind(alignment.tab, inverted.entries))
      
      
      ## When the alignment is inverted, update the information of the reference
      ## motifs in the comparison table
      id1.align.info <- subset(alignment.tab, id == l.id1)
      id2.align.info <- subset(alignment.tab, id == l.id2)
    }

    
    ## Update offset: as this alignment is progressive (hierarchical) offset values
    ##                need to be updated across the iteration
    new.offset <- l.offset - id2.align.info$Offset_up + id1.align.info$Offset_up
    
    
    ## Get the consensus string after the conditional to update consensus orientation
    ## Negative offset: add gaps to id1 (in comparison table)
    ## Positive offset: add gaps to id2
    ids.to.update <- ids.node1
    if (new.offset >= 0) {
      ids.to.update <- ids.node2
    }
    
    
    ## The gaps are added AFTER the strand orientation
    ## Update the gaps
    subset.update.consensus <- alignment.tab %>% 
                                    dplyr::filter(id %in% ids.to.update) %>% 
                                    mutate(Update_status = Update_status + 1,
                                           Offset_up     = abs(new.offset) + Offset_up)

      
    ## Add the gaps
    subset.update.consensus.char <- gsub(subset.update.consensus$Oriented_consensus, pattern = "-", replacement = "")
    subset.update.consensus$Oriented_consensus <- purrr::map2_chr(.x = subset.update.consensus.char,
                                                                  .y = subset.update.consensus$Offset_up,
                                                                  ~add.gaps.consensus(consensus = .x, offset.up = .y))
      
    ## Update the alignment table, keep the most updated entries of each motif id
    alignment.tab <- update.alignment.table(rbind(alignment.tab, subset.update.consensus)) %>%
    		        data.table()
    
    
    ## Fill the gaps in the downstream side of all consensuses within this tree level      
    subset.update.consensus.down <- fill.alignment.downstream(alignment.table = alignment.tab,
                                                              ids.subset      = c(ids.node1, ids.node2))


    ## Update the alignment table, keep the most updated entries of each motif id
    alignment.tab <- update.alignment.table(tab = data.table(rbind(alignment.tab,
                                                                   subset.update.consensus.down))) %>%
    		        data.table()
  }
  

  ## Add the aligned consensus in RC
  consensus.rc.list <- list(consensus = alignment.tab$consensus,
                            gap_up    = alignment.tab$Offset_down,
                            gap_dw    = alignment.tab$Offset_up)
  
  alignment.tab$aligned_consensus_rc <- purrr::pmap_chr(.l = consensus.rc.list,
                                                        .f = ~add.gaps.consensus(consensus   = ..1,
                                                                                 offset.up   = ..2,
                                                                                 offset.down = ..3))
  
  alignment.tab <- data.frame(alignment.tab) %>%
                    slice(match(tree$labels[tree$order], id))
                    # %>%
                    # dplyr::filter(id %in% c(ids.node1, ids.node2))
  
  
  # data.frame(alignment.tab[,c("id", "Oriented_consensus", "Strand")]) %>%
  #   slice(match(tree$labels[tree$order], id)) %>%
  #   dplyr::filter(id %in% c(ids.node1, ids.node2))
  
  # alignment.tab %>%
  #   slice(match(tree$labels[tree$order], id)) %>%
  #   dplyr::filter(id %in% c(ids.node1, ids.node2))
  # 
  # alignment.tab <- alignment.tab %>%
  #                   slice(match(tree$labels[tree$order], id))
  
  return(alignment.tab)
}



## Group the entries by ID and returns the entry with the highest status
update.alignment.table <- function(tab = NULL) {
  
  tab %>% 
    dplyr::ungroup() %>% 
    group_by(id) %>%
    dplyr::filter(Update_status == max(Update_status)) %>%
    arrange(N)
}
