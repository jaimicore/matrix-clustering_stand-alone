###############################################################
## Build a distance matrix from the motif correlation values ##
###############################################################
build.distance.matrix <- function(compa.table = NULL,
                                  metric      = "Ncor") {
  
  ## Supported metrics
  supported.metrics <- c("cor",    ## cor    Pearson correlation (computed on residue occurrences in aligned columns)
                         "Ncor")   ## Ncor   Relative width-normalized Pearson correlation
  
  
  dist.table <- NULL
  distances.objects <- list()
  
  ## Extract metric values
  metric.values <- compa.table[, ..metric]
  
  if (metric %in% supported.metrics) {
    
    ## If required, convert similarities to distances
    ## Similarity sores bounded to 1
    metric.dist <- 1 - metric.values
    
  } else {
    
    stop("; ", metric, " is not a supported metric. Supported metrics: ", paste(supported.metrics, collapse = ","))
    
  }
  
  
  ## Add a column with metric column to the compare matrices table, will
  ## be used to generate a cross-table
  compa.table$distance <- as.vector(unlist(metric.dist))
  
  ## Build the distance table from the column metric
  dist.table <- xtabs(distance ~ id1+id2, compa.table)
  
  ## Ensure that symmetrical distances are defined
  dist.table.sym <- pmax(dist.table, t(dist.table))
  
  ## Cast the distance table into an object of class "dist"
  dist.matrix <- as.dist(dist.table.sym)
  
  original.values     <- xtabs(Ncor ~ id1+id2, compa.table)
  original.values.sym <- pmax(original.values, t(original.values))
    
  return(list(table    = dist.table.sym,
              matrix   = dist.matrix,
              original = original.values.sym))
}



#####################################
## Compute hierarchical clustering ##
#####################################
hclust.motifs <- function(dist.matrix, hclust.method = "average"){

   return(hclust(dist.matrix, method = hclust.method))
}



#####################################################
## Convert the hclust object to a character object ##
## with the lines ready to print a Newick file    ##
####################################################
convert.hclust.to.newick <- function(tree){
  
  ## Require ape (CRNA package) if it is required
  if (!require("ape")) {
    install.packages("ape")
  }
  suppressPackageStartupMessages(library("ape", character.only = TRUE, quietly = TRUE))
  message("; Converting hclust object to newick tree")
  
  ## Round the decimals of the branches distances
  newick.tree.obj <- as.phylo(tree) 
  
  ## Convert the hclust tree to Newick format
  ## If ‘flat=TRUE’ the result is a string can be written in a file      
  newick.tree.str <- write.tree(phy = newick.tree.obj)
  
  return(newick.tree.str)
}



#####################################################
## Convert the hclust object to a character object
## with the lines ready to print a JSON file
convert.hclust.to.JSON <- function(tree){
  
  if (!require("RJSONIO")) {
    install.packages("RJSONIO")
  }
  suppressPackageStartupMessages(library("RJSONIO", character.only = TRUE, quietly = TRUE))
  
  message("; Converting hclust object to a JSON tree")
  
  #######################################################
  ## Extract the tree from an hclust object
  createLeafNode <- function(hclust, i) {
    list(label = hclust$labels[[i]],
         order = hclust$order[[i]])
  }
  
  ########################################################
  ## Convert an hclust tree into a JSON format tree
  hclustToTree <- function(hclust) {
    if (length(hclust$merge) == 0)
      return(NULL)
    merges <- list()
    for (index in 1:nrow(hclust$merge)) {
      left <- hclust$merge[index, 1]
      right <- hclust$merge[index, 2]
      if (left < 0)
        left <- createLeafNode(hclust, -left)
      else
        left <- merges[[left]]
      if (right < 0)
        right <- createLeafNode(hclust, -right)
      else
        right <- merges[[right]]
      if (left$order > right$order) {
        tmp <- left
        left <- right
        right <- tmp
      }
      merges[[index]] <- list(
        children = list(
          left,
          right
        ),
        order = left$order
      )
    }
    return(merges[nrow(hclust$merge)])
  }
  
  ### Creates and parse the json string
  halfway.tree <- hclustToTree(tree)
  jsonTree <- toJSON(halfway.tree)
  
  
  ## Fix some little technical issues for JSON compatibility with the tree display javascript
  jsonTree <- gsub("\\],", "\\]", jsonTree, perl = TRUE)
  jsonTree <- paste("{\n\"name\": \"\",\n\"children\":", jsonTree, "}", sep = "")
  jsonTree <- gsub("\n\"order\":\\s+\\d+", "", jsonTree, perl = TRUE)
  
  return(jsonTree)
}


## Run compare-matrices-quick, this program returns the pairwise motif comparison
## among all motifs in the input collection and the offsets of the best alignment
## of each comparisons, this offsets will be used to align the motifs
motif.comparison <- function(transfac.file     = NULL,
                             output.compa.file = NULL) {
  
  bin <- here("compare-matrices-quick", "compare-matrices-quick")
  
  message("; Running compare-matrices-quick")
  ## Keep the lowest thresholds, this ensures all comparisons are shown (all comparisons satisfy the thresholds)
  compa.mat.quick.command <- paste0(bin," -v 2 -file1 ", transfac.file, " -file2 ", transfac.file, " -lth_ncor1 -1 -lth_ncor2 -1 -lth_ncor -1 -lth_cor -1 -lth_w 0 -mode matches -o ", output.compa.file)
  #message(compa.mat.quick.command)
  system(compa.mat.quick.command)
  
  if (file.exists(output.compa.file)) {
    
    message("; Pairwise motif comparison table: ", output.compa.file)
    motif.comparison.tab <- fread(output.compa.file)
    
    ## Rename and select relevant columns
    motif.comparison.tab <- motif.comparison.tab %>% 
      rename(id1 = '#id1') %>% 
      select(id1, id2, cor, Ncor, strand, offset)
    
    ## Check that the comparison table table contains the required metric and threshold columns
    required.columns <- c("cor", "Ncor")
    if (!all(required.columns %in% names(motif.comparison.tab))) {
      stop("; Missing columns in comparison table. Required columns: ", paste(required.columns, collapse = ", "))
    }
    
  } else {
    stop("Motif comparison table not found: ", output.compa.file)
  }
  
  return(motif.comparison.tab)
}
