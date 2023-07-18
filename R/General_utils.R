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
leaves.per.node <- function(tree   = NULL,
                            labels = FALSE) {
  
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
  
  if (labels) {
    leave.lists <- lapply(leave.lists, function(x){ tree$labels[x] })
  }
  
  return(leave.lists)
}



## Given a vector of motif IDs, return all entries in the comparison table (compa)
## with those motifs, removing repeated entries
##
## The argument 'hash.ind' must be a list where each element is a vector with the
## indexes of a motif in the comparison table
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



create.html.radial.tree <- function(json.file        = NULL,
                                    d3.template      = NULL,
                                    d3.outfile       = NULL,
                                    motif.info       = NULL,
                                    d3.lib           = NULL,
                                    jq.lib           = NULL,
                                    outdir           = NULL,
                                    alignment.length = 20,
                                    html.legend      = NULL) {
  
  # Read D3 template as an array an iterate over it
  d3.lines <- readLines(d3.template)
  d3.lines.updated <- d3.lines
  d3.line.counter  <- 0
  nb.motifs        <- nrow(motif.info)
  
  for (d3l in d3.lines) {
    
    d3.line.counter <- d3.line.counter + 1
    
    # Update tree width
    if (grepl(pattern = "--radial_w--", x = d3l)) {
      
      # Width is determined by number of motifs
      radial.w <- ifelse(nb.motifs > 200, yes = 2500, no = 2000)
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--radial_w--", x = d3l, replacement = radial.w)
    }
    
    if (!is.null(html.legend)) {
      if (grepl(pattern = '<!-- legend -->', x = d3l)) {
        d3.lines.updated[d3.line.counter] <- gsub(pattern = '<!-- legend -->', x = d3l, replacement = html.legend)
      }
    }

    
    # Update tree height
    if (grepl(pattern = "--radial_h--", x = d3l)) {
      
      # Width is determined by number of motifs
      radial.h <- ifelse(nb.motifs > 200, yes = 2500, no = 2000)
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--radial_h--", x = d3l, replacement = radial.h)
    }
    
    # Insert JSON updated file
    if (grepl(pattern = "--json_file--", x = d3l)) {
      
      json.rel.path <- this.path::relpath(relative.to = results.main.dir,
                                          path = json.file)
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--json_file--", x = d3l, replacement = json.rel.path)
    }
    
    
    # Set the tree radium according to the number of motifs
    # NOTE: the following radius were obtained empirically
    if (grepl(pattern = "--radium--", x = d3l)) {
      
      tree.radium <- 0
      
      if (nb.motifs <= 30) {
        tree.radium = 17
      } else if (nb.motifs > 30 & nb.motifs <= 50) {
        tree.radium = 200
      } else if (nb.motifs > 50 & nb.motifs <= 100) {
        tree.radium = 250
      } else if (nb.motifs > 100 & nb.motifs <= 150) {
        tree.radium = 275
      } else if (nb.motifs >= 150 & nb.motifs <= 200) {
        tree.radium = 300
      } else if (nb.motifs >= 200) {
        tree.radium = 450;
      }
      
      # Ring annotations tart/end coordinates depend on tree's radium
      start.annotation.pos <- tree.radium + 5
      end.annotation.pos   <- (max(nchar(motif.info$name)) * 7 + 5) + ((alignment.length + 2) * 5) + 0
      end.annotation.pos   <- start.annotation.pos + end.annotation.pos
      
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--radium--", x = d3l, replacement = tree.radium)
    }
    
    
    # Set the motif logo height according to the number of motifs
    if (grepl(pattern = "--h_motif--", x = d3l)) {
      
      logo.height = 10
      
      if (nb.motifs <= 30) {
        logo.height <- 30
      } else if (nb.motifs > 30 & nb.motifs < 100) {
        logo.height <- 20
      } else if (nb.motifs >= 100 & nb.motifs < 500) {
        logo.height <- 10
      } else if (nb.motifs >= 500) {
        logo.height <- 5
      }
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--h_motif--", x = d3l, replacement = logo.height);
    }
    
    
    # Set displacement (x axis) of the logos according to the number of motifs
    if (grepl(pattern = "--x_displ--", x = d3l)) {
      
      x.displacement <- max(nchar(motif.info$name)) * 7 + 5
      
      # x.displacement <- 75
      # if (nb.motifs <= 30) {
      #   x.displacement <- 75
      # } else if (nb.motifs > 30 & nb.motifs <= 100) {
      #   x.displacement <- 70
      # } else if (nb.motifs > 100 & nb.motifs < 500) {
      #   x.displacement <- 55
      # } else if (nb.motifs >= 500) {
      #   x.displacement <- 50
      # }
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--x_displ--", x = d3l, replacement = x.displacement);
    }
    
    
    # Set displacement (y axis) of the logos according to the number of motifs
    if (grepl(pattern = "--y-displ--", x = d3l)) {
      
      #y.displacement <- (logo.height/2) + 3;
      y.displacement <- 0
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--y-displ--", x = d3l, replacement = y.displacement);
    }
    
  

    # Set inner ring start/end variables
    if (grepl(pattern = "--innerRad_start--", x = d3l)) {
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--innerRad_start--", x = d3l, replacement = start.annotation.pos);
    }
    
    if (grepl(pattern = "--innerRad_end--", x = d3l)) {
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--innerRad_end--", x = d3l, replacement = end.annotation.pos);
    }
    
    
    # Add path to d3 library
    if (grepl(pattern = "--d3--", x = d3l)) {
      
      # Create a copy of the D3 library in the results folder
      d3.path <- cp.d3.lib(d3     = d3.lib,
                           folder = outdir)
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--d3--", x = d3l, replacement = d3.path);
    }
    
    # Add path to jquery library
      if (grepl(pattern = "--jquery--", x = d3l)) {
      
      # Create a copy of the D3 library in the results folder
      jquery.path <- cp.jquery.lib(jquery = jq.lib,
                                   folder = outdir)
      
      d3.lines.updated[d3.line.counter] <- gsub(pattern = "--jquery--", x = d3l, replacement = jquery.path);
    }
  }
  
  # Export JSON file with annotations
  message("; Exporting D3 Radial tree file: ", d3.outfile)
  writeLines(d3.lines.updated, con = d3.outfile)
}


cp.d3.lib <- function(d3     = NULL,
                      folder = NULL) {
  
  message("; Creating a copy of D3 library")
  file.copy(from = dirname(d3),
            to   = folder,
            recursive = TRUE)
  
  # new.d3           <- this.path::here(.. = 1, file.path(dirname(folder), "js", basename(d3)))
  # #results.main.dir <- this.path::here(.. = 1, folder)
  # 
  # d3.radial.js <- this.path::relpath(relative.to = folder,
  #                                    path        = new.d3)
  
  return("js/d3.v3.min.js")
}


cp.jquery.lib <- function(jquery = NULL,
                          folder = NULL) {
  
  message("; Creating a copy of Jquery library")
  file.copy(from = dirname(jquery),
            to   = folder,
            recursive = TRUE)
  
  return("js/jquery.js")
}


Add_attributes_to_JSON_radial_tree <- function(motif.description.tab = NULL,
                                               clusters.list         = NULL,
                                               color.map             = NULL,
                                               htree                 = NULL,
                                               json.woa.file         = NULL,
                                               json.wa.file          = NULL,
                                               alignent.width        = NULL)  {
  
  # Convert description table into a list
  motif.info.list <- motif.description.tab %>% purrr::transpose()
  names(motif.info.list) <- motif.description.tab$id
  
  # Same motif width for all logos
  motif.logo.size <- max(motif.description.tab$width)
  
  # Nodes (tree branches) to cluster association table
  node2cluster <- treenode2cluster(cluster_results = clusters.list,
                                   tree            = htree)
  
  # Motifs (tree leaves) to cluster association table
  leaf2cluster <- stack(clusters.list$clusters) %>% 
    rename("Motif"   = "values",
           "cluster" = "ind") %>% 
    data.table()
  leaf2cluster <- merge(leaf2cluster, color.map)
  
  # Read JSON file, it is stored as an array where each element corresponds to a line
  JSON.lines        <- readLines(json.woa.file)
  cluster.tree      <- "cluster_01"
  tree.branch       <- 0
  line.counter      <- 0
  tree.node         <- ""
  JSON.lines.parsed <- JSON.lines
  for (jl in JSON.lines) {
    
    # ll <- 7  # Node
    # ll <- 5  # Branch
    # ll <- 13 # LEave
    # jl <- JSON.lines[ll]
    # line.counter <- ll
    
    # Initialize
    json.flag <- 0
    line.counter <- line.counter + 1
    
    # ------------------------------------ #
    # Find the line indicating a tree leaf #
    # Update tree leaves                   #
    # ------------------------------------ #
    if (grepl(pattern = '"label":\\s*"(.+)"', jl)) {
      
      tree.label <- gsub(pattern = '"label":\\s*"(.+)"', replacement = "\\1", x = jl)
      tree.label <- gsub(pattern = "\\s+", replacement = "", x = tree.label)
      tree.label
      
      json.flag <- 1
      add.this  <- ""
      
      
      ## Define the URL of the logo files, relative to the location of the json file
      align.logo.link.relpath.F <- this.path::relpath(relative.to = results.main.dir,
                                                      path        = this.path::here(motif.info.list[[tree.label]]$Logo))
      align.logo.link.relpath.R <- this.path::relpath(relative.to = results.main.dir,
                                                      path        = this.path::here(motif.info.list[[tree.label]]$Logo_RC))
      
      ### Create the line that will be added to JSON file
      image.F.line      <- paste0(',\n "image" : "', align.logo.link.relpath.F, '"')
      image.R.line      <- paste0(',\n "image_rc" : "', align.logo.link.relpath.R, '"')
      url.line          <- paste0(',\n "url" : "', align.logo.link.relpath.F, '"')
      ic.line           <- paste0(',\n "ic" : ', round(motif.info.list[[tree.label]]$IC, digits = 3))
      size.line         <- paste0(',\n "size" : ', alignent.width)
      name.line         <- paste0(',\n "name" : "', motif.info.list[[tree.label]]$name, '"')
      link.ext.line     <- paste0(',\n "link_ext" : "', '', '"')
      branch.color.line <- paste0(',\n "branch_color" : "', as.vector(subset(leaf2cluster, Motif == tree.label)$color), '"')
      
      
      
      add.this <- paste0(image.F.line,
                         image.R.line,
                         url.line,
                         ic.line,
                         size.line,
                         name.line,
                         link.ext.line,
                         branch.color.line)
      
      # This needs to be completed
      # if (ID_link_flag) {
      #
      #   link.ext.line <- paste0(',\n "link_ext" : "', '', '"')
      #   color.line    <- paste0(',\n "color" : "', '', '"')
      #
      #   add.this <- paste0(add.this,
      #                      link.ext.line,
      #                      color.line)
      # }
      
      JSON.lines[line.counter] <- paste0(jl, add.this)
    }
    
    # -------------------------------------- #
    # Find the line indicating a tree branch #
    # Update tree branches                   #
    # -------------------------------------- #
    if (grepl(pattern = '"node"\\s*:', jl)) {
      tree.node <- gsub(pattern = '"node":"(node_\\d+)",', replacement = "\\1", x = jl)
      tree.node <- gsub(pattern = "\\s+", replacement = "", x = tree.node)
    }
    
    if (grepl(pattern = '"children"\\s*:', jl)) {
      
      # Update variables
      tree.branch     <- tree.branch + 1
      add.branch.line <- ""
      
      if (tree.branch > 1) {
        
        tree.label <- gsub(pattern = '"label":\\s*"(.+)"', replacement = "\\1", x = jl)
        tree.label <- gsub(pattern = "\\s+", replacement = "", x = tree.label)
        
        # Check this node_to_cluster_hash variable
        # Check that nodes without clusters have a different color
        branch.color <- '#ccc;'
        
        
        node.cluster <- as.vector(subset(node2cluster, node == tree.node)$cluster)
        #print(tree.node)
        # print(node.cluster)
        #print(line.counter)
        if (grepl(pattern = 'node_\\d+', x = tree.node)) {
          branch.color <- as.vector(subset(leaf2cluster, cluster == node.cluster)$color)
        }
        
        
        add.branch.line <- paste0('"branch_color" : "', branch.color, '",\n')
        
        #JSON.lines.parsed <- append(JSON.lines.parsed, add.branch.line, after = line.counter)
        JSON.lines[line.counter] <- paste0(add.branch.line, jl)
      }
    } # end of children grepl if
  } # End for loop
  
  # Export JSON file with annotations
  message("; Exporting JSON file with annotations: ", json.wa.file)
  writeLines(JSON.lines, con = json.wa.file)
  
}



annotate.radial.tree <- function(clusters         = NULL,
                                 cluster2color    = NULL,
                                 tree             = NULL,
                                 motif.annotation = NULL) {

  suppressMessages(library(jsonlite))

  # Motifs (tree leaves) to cluster association table
  treeleaf2cluster <- stack(clusters) %>%
                        rename("Motif"   = "values",
                               "cluster" = "ind") %>%
                        data.table()

  # Cluster::Motif::Cluster-color mapping table
  treeleaf2cluster <- merge(treeleaf2cluster, cluster2color)

  # Map motif IDs
  matrix_order.df <- tree$labels[tree$order]
  matrix_order.df <- treeleaf2cluster[match(matrix_order.df, treeleaf2cluster$Motif), ] |>
    select(Motif, color) |>
    rename(motif_id    = Motif,
           cluster_nb  = color) |>
    mutate(id_motif = 1:nrow(treeleaf2cluster))

  # Read and parse annotation file
  motif.annotation$class      <- gsub("'", "", as.character(motif.annotation$class))
  motif.annotation$collection <- gsub("\"", "", as.character(motif.annotation$collection))
  motif.annotation$motif_id   <- gsub("\"", "", as.character(motif.annotation$motif_id))


  # matrix_element <- annotation.df[1,]
  # rm(matrix_element)
  annotation.merge.df <- apply(motif.annotation, 1 ,function(matrix_element){

    # Create regex
    #regex_string    <- paste0("^", matrix_element["collection_name"], "_m\\d+_", matrix_element["matrix_name"],"_n+\\d+$")
    regex_string <- matrix_element["motif_id"]
    # message(regex_string)

    # NOTE WSG: spcial character from regex must be forbidden in
    # the file names. Need to add a sanity check
    # Safe check for special characters in regex
    #regex_string <- sub("\\+", "\\\\+", regex_string)

    # Find complete motif name
    match_matrix.df <- matrix_order.df[grepl(regex_string, as.character(matrix_order.df$motif_id)), ]
    # print(match_matrix.df)

    # Transpose match
    matrix_element <- as.data.frame(t(as.data.frame(matrix_element)))

    # Add annotation columns to matched matrix
    match_matrix.df <- cbind(match_matrix.df, matrix_element)

    return(match_matrix.df)
  }) %>% (plyr::rbind.fill)

  message("; Calculating degrees intervals for each annotation layer.")

  # Sort the data frame by id_motif
  annotation.merge.df       <- annotation.merge.df[order(annotation.merge.df$id_motif),]
  annotation.merge.df$start <- (0:(dim(annotation.merge.df)[1] - 1)) * (360/dim(annotation.merge.df)[1])
  annotation.merge.df$end   <- (1:(dim(annotation.merge.df)[1])) * (360/dim(annotation.merge.df)[1])

  # Convert annotation DF to JSON format
  annotation.json <- toJSON(annotation.merge.df)

  message("; Annotating JSON file")
  output.files.list$annotation_json_file  <- file.path(out.folder.list$trees, "annotation_matrix.json")
  write_json(annotation.json, output.files.list$annotation_json_file)

  prefix <- gsub(output.files.list$D3_radial_tree, pattern = "_D3_radial_tree.html", replacement = "")
  cmd_annotate_htmltree <- paste0("python3 annotate-html-radialtree.py ",
                                  "-i ", prefix)

  # Launch annotation script
  # Walter Santana wrote this script in python
  # It works for the moment but ideally we should convert the script to R code 
  #print(cmd_annotate_htmltree)
  system(cmd_annotate_htmltree)
}



create.color.annotation <- function(motif.meta.file = NULL,
                                    ann.outdir      = NULL) {
  
  # ---------------------------- #
  # Read and parse metadata file #
  # ---------------------------- #
  message("; Reading motif metadata table")
  motif.meta <- fread(motif.meta.file)
  
  ## For the dimers considers only the first TF class
  motif.meta$class <- gsub(motif.meta$class, pattern = ",.+$", replacement = "", perl = T)
  motif.meta$class <- gsub(motif.meta$class, pattern = "::.+$", replacement = "", perl = T)
  
  ## Add the 'Unknown' TF class to entries without a class
  motif.meta$class <- ifelse(motif.meta$class == "", yes = "Unknown", no = motif.meta$class)
  
  
  # --------------------- #
  # Create TF-class table #
  # --------------------- #
  
  
  ## Get the TF class names ordered by number of motifs
  TF.class.order <- motif.meta %>% 
    group_by(class) %>% 
    tally() %>% 
    arrange(desc(n)) %>% 
    dplyr::filter(class != "Unknown")
  TF.class.order <- TF.class.order$class  
  
  
  # ---------------------------- #
  # TF class - Color assignation #
  # ---------------------------- #
  message("; Assigning colours to TF classes")
  
  TF.known.classes    <- TF.class.order
  TF.known.classes.nb <- length(TF.known.classes)
  
  ## Grey color for Unknown class
  unknown.class.color <- "#888888"
  
  nb.classes.palette <- 12  ## We want to use the first 12 colors of the Safe palette
  nb.seed.colors     <- ifelse(TF.known.classes.nb < nb.classes.palette,
                               yes = TF.known.classes.nb,
                               no = nb.classes.palette)
  
  ## Generate a carto palette (remove gray value)
  carto.pal.classes  <- carto_pal(nb.seed.colors, "Safe")
  carto.pal.classes  <- carto.pal.classes[which(carto.pal.classes != "#888888")]
  
  ## Expand the color palette and add the gray color at the end
  class.colors        <- c(colorRampPalette(carto.pal.classes, space = "Lab")(TF.known.classes.nb),
                           unknown.class.color)
  names(class.colors) <- c(TF.known.classes, "Unknown")
  df.class.colour     <- data.frame(colour   = class.colors,
                                    class    = names(class.colors),
                                    class_nb = seq_len(TF.known.classes.nb + 1)) 
  
  
  
  # --------------------------- #
  # TF_class::Colour HTML table #
  # --------------------------- #
  
  ## Table header + tail
  head.tab <- "<div id='Color_class_tab' style='display: inline-block;float:left;position:relative;' class='color-legend' width='450px'><p style='font-size:12px;padding:0px;border:0px'><b></b></p><table id='Color_class_table' class='hover compact stripe' cellspacing='0' width='450px' style='padding:15px;align:center;'><thead><tr><th > Color </th> <th> TF Class </th> <th> Number</th> </tr></thead><tbody>"
  tab.lines <- paste("\n<tr><td class='color-box' style='background-color: --color--';></td> <td>--TFClass--</td> <td>--TFClass_ID--</td> </tr>", collapse = "")
  tail.tab <- "<tr><td class='non_validated'>*</td><td>Unvalidated</td></tr></tbody></table></div>"
  
  
  ## Table body
  table.body <- sapply(1:nrow(df.class.colour), function(r.nb){
    
    ## Set variables
    tab.lines.current.line <- tab.lines
    TF.class.Color <- df.class.colour[r.nb,1]
    TF.class.Name  <- df.class.colour[r.nb,2]
    TF.class.ID    <- df.class.colour[r.nb,3]
    
    tab.lines.current.line <- gsub("--TFClass--", TF.class.Name, tab.lines.current.line)
    tab.lines.current.line <- gsub("--color--", TF.class.Color, tab.lines.current.line)
    tab.lines.current.line <- gsub("--TFClass_ID--", TF.class.ID, tab.lines.current.line)
    
  })
  table.body <- paste(table.body, collapse = "")
  
  html.table <- paste(head.tab, table.body, tail.tab, collapse = "")
  html.annotation.table <- file.path(ann.outdir, "Radial_tree_legend.html")
  message("; Radial tree HTML legend : ", html.annotation.table)
  writeLines(html.table, html.annotation.table)
  
  
  # --------------------- #
  # Export metadata table #
  # --------------------- #
  
  ## Combine tables
  motif.meta.colour <- merge(motif.meta, df.class.colour, by = "class") |> 
    select(collection, motif_id, colour, class, class_nb)
  
  annotation.table.file <- file.path(ann.outdir, "annotation_table.txt")
  message("; Motif annotation table : ", annotation.table.file)
  fwrite(motif.meta.colour, file = annotation.table.file, sep = "\t")
  
  return(list(df   = motif.meta.colour,
              html = html.table))
}