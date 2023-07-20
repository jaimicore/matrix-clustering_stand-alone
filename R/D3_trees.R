
create.html.radial.tree <- function(json.file        = NULL,
                                    d3.template      = NULL,
                                    d3.outfile       = NULL,
                                    motif.info       = NULL,
                                    d3.lib           = NULL,
                                    jq.lib           = NULL,
                                    outdir           = NULL,
                                    alignment.length = 20,
                                    html.legend      = NULL,
                                    barplot.ann      = NULL) {
  
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
    
    if (!is.null(barplot.ann)) {
      
      barplot.ann.path <- relpath(path   = barplot.ann,
                                  relative.to = results.main.dir)
      
      barplot.ann.html <- paste0('<img class="barplot Color_class_table" src="', barplot.ann.path, '" >') 
      
      if (grepl(pattern = '<!-- barplot -->', x = d3l)) {
        d3.lines.updated[d3.line.counter] <- gsub(pattern = '<!-- barplot -->', x = d3l, replacement = barplot.ann.html)
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
  
  ID_link_flag <- FALSE
  if ("url" %in% colnames(motif.description.tab)) {
    ID_link_flag <- TRUE
  }
  
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
      branch.color.line <- paste0(',\n "branch_color" : "', as.vector(subset(leaf2cluster, Motif == tree.label)$color), '"')
      
      
      
      add.this <- paste0(image.F.line,
                         image.R.line,
                         url.line,
                         ic.line,
                         size.line,
                         name.line,
                         branch.color.line)
      
      # This needs to be completed
      if (ID_link_flag) {
        
        link.ext.line <- paste0(',\n "url" : "', motif.info.list[[tree.label]]$url, '"')
        color.line    <- paste0(',\n "color" : "', motif.info.list[[tree.label]]$colour, '"')
        
        add.this <- paste0(add.this,
                           link.ext.line,
                           color.line)
      }
      
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
  head.tab <- "<div id='Color_class_tab' style='display: inline-block;float:left;position:relative;' class='color-legend' width='450px'><p style='font-size:12px;padding:0px;border:0px'><b></b></p><table id='Color_class_table' class='Color_class_table hover compact stripe' cellspacing='0' width='450px' style='padding:15px;align:center;'><thead><tr><th > Color </th> <th> TF Class </th> <th> Number</th> </tr></thead><tbody>"
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
    select(collection, motif_id, colour, class, class_nb, url)
  
  annotation.table.file <- file.path(ann.outdir, "annotation_table.txt")
  message("; Motif annotation table : ", annotation.table.file)
  fwrite(motif.meta.colour, file = annotation.table.file, sep = "\t")
  
  return(list(df   = motif.meta.colour,
              html = html.table,
              cpal = class.colors))
}


which.collection <- function(id = NULL) {
  
  if (grepl(id, pattern = "MA\\d+\\.\\d+" )) {
    return("CORE")
  } else if (grepl(id, pattern = "UN\\d+\\.\\d+")) {
    return("UNVALIDATED")
  } else {
    return("None")
  }
}
wc <- Vectorize(which.collection)



# This function exports a barplot (jpeg) with the number of motifs in each class
annotation.barplot <- function(df            = NULL,
                               class.colors  = NULL,
                               barplot.title = NULL,
                               barplot.file  = NULL) {
  
  # Init
  collection.type <- NULL
  has.MA          <- NULL
  has.UM          <- NULL
  alpha.values    <- NULL
  y.axis.lab      <- NULL
  
  # Sort classes by number of TFs
  # Unkown classes at the end of the vector
  TF.class.sorted <- names(sort(table(df$class), decreasing = TRUE))
  TF.class.sorted <- TF.class.sorted[!TF.class.sorted %in% "Unknown"]
  TF.class.sorted <- c(TF.class.sorted, "Unknown")
  
  df$class <- factor(df$class, levels = rev(TF.class.sorted))
  
  # Check format of motif IDs (this is specific for JASPAR motifs)
  has.MA <- sum(as.vector(sapply(motif.annotation.list$df$motif_id, grepl, pattern = "MA\\d+\\.\\d+"))) > 0
  has.UM <- sum(as.vector(sapply(motif.annotation.list$df$motif_id, grepl, pattern = "UM\\d+\\.\\d+"))) > 0
  
  if (has.MA & !has.UM) {
    collection.type <- "CORE"
  } else if (has.MA & has.UM) {
    collection.type <- "UNVALIDATED"
  } else {
    collection.type <- "None"
  }
  
  ## This dataframe contains the labels to be displayed on each bar
  ## Adapt to collection type
  if (collection.type %in%  c("CORE", "UNVALIDATED")) {
    label.df <- df |> 
      mutate(type = wc(motif_id)) |> 
      group_by(class) %>% 
      mutate(Total_Class = n()) %>%
      group_by(class, type) |> 
      mutate(COREc = sum(type %in% "CORE"),
             UNVc  = sum(type %in% "UNVALIDATED")) %>% 
      select(class, COREc, UNVc) %>% 
      distinct() %>% 
      ungroup() %>% 
      group_by(class) %>% 
      mutate(CORE        = max(COREc),
             UNVALIDATED = max(UNVc)) %>% 
      select(class, CORE, UNVALIDATED) %>% 
      distinct()
    
  } else {
    label.df <- df |> 
      mutate(type = wc(motif_id)) |> 
      group_by(class) %>% 
      mutate(Total_Class = n()) %>%
      group_by(class, type) |> 
      mutate(COREc = sum(type %in% "None")) %>% 
      select(class, COREc) %>% 
      distinct() %>% 
      ungroup() %>% 
      group_by(class) %>% 
      mutate(CORE        = max(COREc)) %>% 
      select(class, CORE) %>% 
      distinct()
  }
  
  
  if (collection.type == "UNVALIDATED") {
    
    label.df <- label.df %>% 
      mutate(lab         = paste0(CORE, "/", UNVALIDATED),
             Total_Class = sum(c(CORE, UNVALIDATED))) %>% 
      data.table()
    
  } else if (collection.type %in% c("CORE", "None")) {
    
    label.df <- label.df %>% 
      mutate(lab         = CORE,
             Total_Class = sum(CORE)) %>% 
      data.table()
  }
  
  
  ## Use this variable to have enough space to insert the labels
  max.y <- nrow(df) + 25
  
  
  ## The alpha parameter and plot title change depending on the collection type
  if (collection.type == "UNVALIDATED") {
    alpha.values  <- c(0.4, 1)
    y.axis.lab    <- "Number of motifs (CORE/UNVALIDATED)"
  } else if (collection.type == "CORE") {
    alpha.values  <- 1
    y.axis.lab    <- "Number of motifs (CORE)"
  } else if (collection.type == "None") {
    alpha.values  <- 1
    y.axis.lab    <- "Number of motifs"
  }
  
  TF.class.barplot <- ggplot(df, aes(alpha = collection, x = class, fill = class)) +
    geom_bar() +
    coord_flip() +
    scale_fill_manual(values = class.colors) +
    theme_classic() +
    labs(x = "TF class", y = y.axis.lab, title = barplot.title) +
    scale_alpha_manual(values = alpha.values) +
    theme(text        = element_text(size = 15),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
          axis.text.x = element_text(hjust = 0.5, size = 15),
          legend.position = "none",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 9),
          legend.box      = "vertical",
          plot.title      = element_text(hjust = 0.5),
          panel.grid.major.x = element_line(color = "#969696",
                                            size = 0.25,
                                            linetype = 2)) +
    geom_text(data = label.df, aes(x = class, y = Total_Class, label = lab), vjust = 0.5, hjust = -0.2, inherit.aes = F, size = 3.5) +
    guides(fill = guide_legend(reverse = T))  +
    scale_y_continuous(limits = c(0, max.y), expand = c(0,2), breaks = seq(0, max.y, by = 50)[-1])
  
  ggsave(plot     = TF.class.barplot,
         filename = barplot.file,
         width    = 12,
         height   = 8)
  message("; JPEG barplot created: ", barplot.file)
}




# cl        <- "cluster_01"
# cl.motifs <- find.clusters.list$clusters[[cl]]
# motif.description.tab = results.list$Motif_info_tab |> dplyr::filter(id %in% cl.motifs)
# clusters.list         = find.clusters.list
# color.map             = cl.col
# htree                 = cl.hclust.results[[cl]]
#json.woa.file          = subset(results.list$Clusters_files, Cluster == cl)$JSON_file
#json.wa.file           = subset(results.list$Clusters_files, Cluster == cl)$JSON_annotated_file
# alignent.width        = max(subset(results.list$Alignment_table, cluster == cl)$width)
# 
# Add_attributes_to_JSON_interactive_tree <- function(motif.description.tab = NULL,  # Subset of the table including only the motifs in a cluster
#                                                     clusters.list         = NULL,
#                                                     color.map             = NULL,
#                                                     htree                 = NULL,
#                                                     json.woa.file         = NULL,
#                                                     json.wa.file          = NULL,
#                                                     alignent.width        = NULL)  {
#   
# 
#   # Convert description table into a list
#   motif.info.list <- motif.description.tab %>% purrr::transpose()
#   names(motif.info.list) <- motif.description.tab$id
#   
#   # Same motif width for all logos
#   motif.logo.size <- max(motif.description.tab$width)
#   
#   # Nodes (tree branches) to cluster association table
#   node2cluster <- treenode2cluster(cluster_results = clusters.list,
#                                    tree            = htree)
#   
#   # Motifs (tree leaves) to cluster association table
#   leaf2cluster <- stack(clusters.list$clusters) %>% 
#     rename("Motif"   = "values",
#            "cluster" = "ind") %>% 
#     data.table()
#   leaf2cluster <- merge(leaf2cluster, color.map)
#   
#   # Read JSON file, it is stored as an array where each element corresponds to a line
#   JSON.lines        <- readLines(json.woa.file)
#   cluster.tree      <- "cluster_01"
#   tree.branch       <- 0
#   line.counter      <- 0
#   tree.node         <- ""
#   JSON.lines.parsed <- JSON.lines
#   for (jl in JSON.lines) {
#     
#     # ll <- 7  # Node
#     # ll <- 5  # Branch
#     # ll <- 13 # LEave
#     # jl <- JSON.lines[ll]
#     # line.counter <- ll
#     
#     # Initialize
#     json.flag <- 0
#     line.counter <- line.counter + 1
#     
#     # ------------------------------------ #
#     # Find the line indicating a tree leaf #
#     # Update tree leaves                   #
#     # ------------------------------------ #
#     if (grepl(pattern = '"label":\\s*"(.+)"', jl)) {
#       
#       tree.label <- gsub(pattern = '"label":\\s*"(.+)"', replacement = "\\1", x = jl)
#       tree.label <- gsub(pattern = "\\s+", replacement = "", x = tree.label)
#       tree.label
#       
#       json.flag <- 1
#       add.this  <- ""
#       
#       
#       ## Define the URL of the logo files, relative to the location of the json file
#       align.logo.link.relpath.F <- this.path::relpath(relative.to = results.main.dir,
#                                                       path        = this.path::here(motif.info.list[[tree.label]]$Logo))
#       align.logo.link.relpath.R <- this.path::relpath(relative.to = results.main.dir,
#                                                       path        = this.path::here(motif.info.list[[tree.label]]$Logo_RC))
#       
#       ### Create the line that will be added to JSON file
#       image.F.line      <- paste0(',\n "image" : "', align.logo.link.relpath.F, '"')
#       image.R.line      <- paste0(',\n "image_rc" : "', align.logo.link.relpath.R, '"')
#       url.line          <- paste0(',\n "url" : "', align.logo.link.relpath.F, '"')
#       ic.line           <- paste0(',\n "ic" : ', round(motif.info.list[[tree.label]]$IC, digits = 3))
#       size.line         <- paste0(',\n "size" : ', alignent.width)
#       name.line         <- paste0(',\n "name" : "', motif.info.list[[tree.label]]$name, '"')
#       branch.color.line <- paste0(',\n "branch_color" : "', as.vector(subset(leaf2cluster, Motif == tree.label)$color), '"')
#       
#       
#       
#       add.this <- paste0(image.F.line,
#                          image.R.line,
#                          url.line,
#                          ic.line,
#                          size.line,
#                          name.line,
#                          branch.color.line)
#       
#       # This needs to be completed
#       if (ID_link_flag) {
#         
#         link.ext.line <- paste0(',\n "url" : "', motif.info.list[[tree.label]]$url, '"')
#         color.line    <- paste0(',\n "color" : "', motif.info.list[[tree.label]]$colour, '"')
#         
#         add.this <- paste0(add.this,
#                            link.ext.line,
#                            color.line)
#       }
#       
#       JSON.lines[line.counter] <- paste0(jl, add.this)
#     }
#     
#     # -------------------------------------- #
#     # Find the line indicating a tree branch #
#     # Update tree branches                   #
#     # -------------------------------------- #
#     if (grepl(pattern = '"node"\\s*:', jl)) {
#       tree.node <- gsub(pattern = '"node":"(node_\\d+)",', replacement = "\\1", x = jl)
#       tree.node <- gsub(pattern = "\\s+", replacement = "", x = tree.node)
#     }
#     
#     if (grepl(pattern = '"children"\\s*:', jl)) {
#       
#       # Update variables
#       tree.branch     <- tree.branch + 1
#       add.branch.line <- ""
#       
#       if (tree.branch > 1) {
#         
#         tree.label <- gsub(pattern = '"label":\\s*"(.+)"', replacement = "\\1", x = jl)
#         tree.label <- gsub(pattern = "\\s+", replacement = "", x = tree.label)
#         
#         # Check this node_to_cluster_hash variable
#         # Check that nodes without clusters have a different color
#         branch.color <- '#ccc;'
#         
#         
#         node.cluster <- as.vector(subset(node2cluster, node == tree.node)$cluster)
#         #print(tree.node)
#         # print(node.cluster)
#         #print(line.counter)
#         if (grepl(pattern = 'node_\\d+', x = tree.node)) {
#           branch.color <- as.vector(subset(leaf2cluster, cluster == node.cluster)$color)
#         }
#         
#         
#         add.branch.line <- paste0('"branch_color" : "', branch.color, '",\n')
#         
#         #JSON.lines.parsed <- append(JSON.lines.parsed, add.branch.line, after = line.counter)
#         JSON.lines[line.counter] <- paste0(add.branch.line, jl)
#       }
#     } # end of children grepl if
#   } # End for loop
#   
#   # Export JSON file with annotations
#   message("; Exporting JSON file with annotations: ", json.wa.file)
#   writeLines(JSON.lines, con = json.wa.file)
#   
# }