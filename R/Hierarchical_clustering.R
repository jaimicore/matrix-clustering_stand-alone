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
  compa.table <- data.table(compa.table)
  
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
  
  bin <- this.path::here(.. = 0, "compare-matrices-quick", "compare-matrices-quick")
  
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
                              dplyr::rename(id1 = '#id1') %>% 
                              select(id1, id2, cor, Ncor, strand, offset, w) %>% 
                              data.table()
    
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



## Calculate the Adjusted-Rand-Index by comparing the RSAT matrix-clustering and 
## a user provided reference clusters
calculate.ARI <- function(matrix.clustering.clusters = NULL,
                          reference.clusters         = NULL) {
  
  if (!require("flexclust")) {
    install.packages("flexclust")
  }
  suppressPackageStartupMessages(library("flexclust", character.only = TRUE, quietly = TRUE))
  
  
  ## Combine matrix-clustering with reference-cluster tables
  mc.rf <- matrix.clustering.clusters %>% 
              left_join(reference.clusters, by = "id") %>% 
              dplyr::rename(matrix_clustering = cluster.x,
                            reference         = cluster.y)
  
  mc.rf.tab <- table(mc.rf$matrix_clustering, mc.rf$reference)
  mc.rf.ari <- randIndex(mc.rf.tab, correct = TRUE)
  mc.rf.ri  <- randIndex(mc.rf.tab, correct = FALSE)
  
  return(list(ARI = mc.rf.ari,
              RI  = mc.rf.ri,
              tab = mc.rf.tab))
  
}



## A small function to convert a list with vectors to a dataframe
clusters.list.to.df <- function(clusters.list = NULL,
                                id.pttrn.rm   = "") {
  
    # "_n\\d+$"
  
    reshape2::melt(clusters.list) %>% 
        dplyr::rename(id      = value,
                      cluster = L1) %>% 
        data.table() %>% 
        dplyr::mutate(id = gsub(x = id, pattern = id.pttrn.rm, replacement = ""))
}



draw.heatmap.clusters.vs.ref <- function(clusters.tab = NULL,
                                         comment      = "ARI = 0.5") {
  
  ## List of packages to install from CRAN
  required.packages = c("circlize",      
                        "ComplexHeatmap",
                        "RColorBrewer")
  
  for (lib in required.packages) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
  }
  
  
  ## Convert the actual values into frequencies
  ## This helps to better detect the classes in the heatmap
  clusters.tab.mt <- as.data.frame.matrix(clusters.tab) %>%
                        dplyr::rename(Unkwnon = V1) %>% 
                        as.matrix()
  clusters.tab.perc <- clusters.tab.mt

  
  palette <- rev(colorRampPalette(rev(brewer.pal(9, "PuRd")), space = "Lab")(max(clusters.tab.perc)))
  col_fun <- colorRamp2(seq(0, max(clusters.tab.perc), length.out = max(clusters.tab.perc)), palette)
  
  
  ## Sidebar annotations
  max.row.col <- max( max(colSums(clusters.tab.mt)), max(rowSums(clusters.tab.mt)))
  palette.ha  <- rev(colorRampPalette(rev(brewer.pal(9, "Reds")), space = "Lab")(max.row.col + 1))
  colfun.ha   <- colorRamp2(seq(0, max.row.col, length.out = max.row.col + 1 ), palette.ha)
  
  
  hac = HeatmapAnnotation(df                   = data.frame(Group_size = colSums(clusters.tab.mt)),
                          col                  = list(Group_size = colfun.ha),
                          simple_anno_size     = unit(0.5, "cm"),
                          which                = "column",
                          show_legend          = FALSE,
                          show_annotation_name = FALSE,
                          text                 = anno_text(colSums(clusters.tab.mt), rot = 0, location = 3, just = "center", gp = gpar(fontsize = 5)))
  
  har = HeatmapAnnotation(df                      = data.frame(Group_size = rowSums(clusters.tab.mt)),
                          col                     = list(Group_size = colfun.ha),
                          simple_anno_size        = unit(0.5, "cm"),
                          which                   = "row",
                          show_legend             = FALSE,
                          show_annotation_name    = FALSE,
                          annotation_legend_param = list(Group_size = list(direction      = "horizontal",
                                                                           legend_width   = unit(4, "cm"),
                                                                           title_position = "topcenter")),
                          text = anno_text(rowSums(clusters.tab.mt), location = -1.5, just = "center", gp = gpar(fontsize = 5)))
  
  ## Draw heatmap
  ht1 = Heatmap(matrix                  = clusters.tab.perc,
                name                    = "Fraction",
                row_title               = "RSAT matrix-clustering",
                column_title            = paste0("Reference groups - ", comment),
                column_title_side       = "top",
                col                     = col_fun,
                show_row_names          = TRUE,
                show_column_names       = TRUE,
                cluster_rows            = TRUE,
                show_row_dend           = FALSE,
                cluster_columns         = TRUE,
                show_column_dend        = FALSE,
                rect_gp                 = gpar(col = "white", lwd = 1),
                column_names_rot        = 65,
                bottom_annotation       = hac,
                right_annotation        = har, 
                column_names_max_height = max_text_width(colnames(clusters.tab.perc), gp = gpar(fontsize = 10)),
                row_names_gp            = gpar(fontsize = 13),
                heatmap_legend_param    = list(direction      = "horizontal", 
                                               legend_width   = unit(4, "cm"),
                                               title_position = "topcenter"),
                cell_fun          = function(j, i, x, y, width, height, fill) {
                  if (clusters.tab.perc[i, j] > 0) {
                    # grid.text(sprintf("%.2f", clusters.tab.perc[i, j]), x, y, gp = gpar(fontsize = 5))
                    grid.text(clusters.tab.perc[i, j], x, y, gp = gpar(fontsize = 5))
                  }
                })
  
  ht1
  
  # dev.off()
  
  return(ht1)
}


## Export a ComplexHeatmap object as a PDF. Precalculates the size to generate
## cells with same height and width
export.heatmap.clusters.vs.ref <- function(ht        = NULL,
                                           ht.matrix = NULL,
                                           pdf.file  = NULL) {
  
  nb.rows.contingecy.cl.vs.ref <- nrow(ht.matrix)
  nb.cols.contingecy.cl.vs.ref <- ncol(ht.matrix)
  y              <- NULL
  
  for (nr in c(nb.rows.contingecy.cl.vs.ref, nb.rows.contingecy.cl.vs.ref * 2)) {
    
    heatmap.ari.draw <- draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", height = unit(5, "mm") * nr, gap = unit(50, "mm"))
    ht_height        <- sum(component_height(heatmap.ari.draw)) + unit(4, "mm")
    ht_height        <- convertHeight(ht_height, "inch", valueOnly = TRUE)
    y                <- c(y, ht_height)
    dev.off()
  }
  lm.xy <- lm(y ~ c(nb.rows.contingecy.cl.vs.ref, nb.rows.contingecy.cl.vs.ref * 2))
  
  
  message("; Exporting cluster vs Reference contingency table as PDF file: ", pdf.file)
  ht_opt$message = FALSE
  pdf(file   = pdf.file,
      width  = ht_height/3.5,
      height = as.vector(lm.xy$coefficients[2]) * nb.rows.contingecy.cl.vs.ref + as.vector(lm.xy$coefficients[1]))
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", show_heatmap_legend = FALSE, height = unit(5, "mm") * nb.rows.contingecy.cl.vs.ref, gap = unit(50, "mm"))
  dev.off()
}
