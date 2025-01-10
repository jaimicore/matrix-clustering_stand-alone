## Reset the ID (name) of universalmotif object
set.um.id.name <- function(um = NULL,
                           new.id = "New_ID") {
  um@name <- new.id
  return(um)
}


## Reset the 'Alternate name' field of universalmotif object
set.um.alt.name <- function(um       = NULL,
                            new.name = "New_alt_name") {
  um@altname <- new.name
  return(um)
}


## Reset the 'Alternate name' field of universalmotif object
set.um.nb.sites <- function(um       = NULL,
                            nb.sites = 100) {
  um@nsites <- nb.sites
  return(um)
}


## Reset the 'Motif' field of universalmotif object
set.um.motif <- function(um        = NULL,
                         new.count = NULL) {
  
  ## As the motif is modified, the universalmotif object needs to be updated
  new.motif <- create_motif(input    = new.count,
                            type     = um@type,
                            alphabet = "DNA")
  
  ## Rename it
  new.motif@name    <- um@name
  new.motif@altname <- um@altname

    return(new.motif)
}


## Subset the count matrix columns, this is applied when the motifs are trimmed
subset.matrix <- function(m   = NULL,
                          min = NULL,
                          max = NULL) {
  
  m <- m[, min:max, drop = FALSE]
  
  # In cases when subseting results in 1 column motifs reset the rownames
  rownames(m) <- NULL
  # rownames(m) <- c("A", "C", "G", "T")
  
  return(m)
}



## Count the number of sites in an input motif, in case there are columns with different
## numbers, the smallest number will be considered 
nb.sites.in.um <- function(um.motif = NULL){
  
  ## Calculate the number of sites in the motif by summin the columns in the matrix
  nbsites.vec       <- sort(unique(as.vector(colSums(um.motif))))
  nb.uniq.count.sum <- length(nbsites.vec)

  if (nb.uniq.count.sum > 1) {
    warning("; Motif colSum is different in ", nb.uniq.count.sum, " positions. Motif nbsites attibute will be set to ", min(nbsites.vec))
  }
  
  return(min(nbsites.vec))
}


# um <- cluster.buster.uo
calculate.nbsites.uo <- function(um = NULL) {
  
  um.motifs     <- furrr::future_map(um, `[`, "motif")
  motifs.colsum <- furrr::future_map_dbl(um.motifs, nb.sites.in.um)

  return(motifs.colsum)
}


# motif.file   = matrix.file.list$Motif_file

## This function reads motifs in cluster-buster format and returns a universalmotif object
read_cluster_buster <- function(motif.file = NULL) {
  
  ## Read motif file using universalmotif functions    
  ## The separator must be a blank line between motifs
  cluster.buster.uo <- universalmotif::read_matrix(file = motif.file, positions = "rows")
  # cluster.buster.uo <- cluster.buster.uo[[1]]
  
  ## Set nbsites attribute correctly
  ## By default it is set to 100 and this create a problem when generating RC because
  ## the counts are scaled to the number of sites
  cluster.buster.uo.nbsites <- calculate.nbsites.uo(cluster.buster.uo)
  
  ## Set nbsites attribute correctly
  cluster.buster.uo.new.id <- purrr::map2(.x = cluster.buster.uo,
                                          .y = cluster.buster.uo.nbsites,
                                          .f = ~set.um.nb.sites(um       = .x,
                                                                nb.sites = .y))
  
  ## Add name and alternate name fields
  cluster.buster.uo.id <-  purrr::map_chr(cluster.buster.uo, `[`, "name")
  cluster.buster.uo.id <- gsub(cluster.buster.uo.id, pattern = ">", replacement = "")
  
  
  cluster.buster.uo.new.id <- purrr::map2(.x = cluster.buster.uo,
                                          .y = cluster.buster.uo.id,
                                          .f = ~set.um.id.name(um     = .x,
                                                               new.id = .y))
  
  
  cluster.buster.uo.new.altname <- purrr::map2(.x = cluster.buster.uo.new.id,
                                               .y = cluster.buster.uo.id,
                                               .f = ~set.um.alt.name(um       = .x,
                                                                     new.name = .y))
  
  return(cluster.buster.uo.new.altname)
}



## Returns a data.frame with the description table of each motif in the input collection
preprocess.one.motif.collection <- function(motif.file      = NULL,
                                            motif.format    = NULL,
                                            collection.name = "Motif_collection_1") {
  
  ## Read motifs 
  message("; Reading motif collection: ", collection.name)
  motif.collection <- read.motif.file(motif.file   = motif.file,
                                      motif.format = motif.format, )
  
  # message("; Read CB")
  # read.motif.file(motif.file   = "Debug/Triana/matrix_clust.cb",
  #                 motif.format = "cluster-buster")
  
  
  ## This step is required when there are cases of motifs with empty columns, for some reason
  ## this generates an NA within the Universalmotif object and crashes the script
  #message("; Trimming motifs")
  # motif.collection <- universalmotif::trim_motifs(motif.collection, min.ic = 0.0001)
  motif.collection <- universalmotif::trim_motifs(motif.collection, min.ic = -1) ## To avoid omitting positions with equal IC distributions
  
  ## Generate the reverse-complement version of the input motif collection
  #message("; Reverse-complement motifs")
  motif.collection.rc <- universalmotif::motif_rc(motif.collection)
  
  ## In cases when there is only 1 motif in the motif collection, save the universalmotif
  ## object within a list
  if (length(motif.collection) == 1) {
    motif.collection    <- list(motif.collection)
    motif.collection.rc <- list(motif.collection.rc)
  }
  
  
  ## Returns a motif information data.table
  if (motif.format == "homer") {
    
    motif.info.dt <- data.table(id_old       = purrr::map_chr(motif.collection, `[`, "name"),
                                name         = purrr::map_chr(motif.collection, `[`, "name"),
                                rc_consensus = purrr::map_chr(motif.collection.rc, `[`, "consensus"),
                                consensus    = purrr::map_chr(motif.collection, `[`, "consensus"),
                                nb_sites     = round(purrr::map_dbl(motif.collection, `[`, "nsites")),
                                IC           = purrr::map_dbl(motif.collection, `[`, "icscore")) %>% 
      dplyr::mutate(width = nchar(consensus),
                    n     = 1:n()) %>% 
      dplyr::mutate(id = paste0(collection.name, "_", id_old, "_n", n))
    
  } else {
    
    motif.info.dt <- data.table(id_old       = purrr::map_chr(motif.collection, `[`, "name"),
                                name         = purrr::map_chr(motif.collection, `[`, "altname"),
                                rc_consensus = purrr::map_chr(motif.collection.rc, `[`, "consensus"),
                                consensus    = purrr::map_chr(motif.collection, `[`, "consensus"),
                                nb_sites     = round(purrr::map_dbl(motif.collection, `[`, "nsites")),
                                IC           = purrr::map_dbl(motif.collection, `[`, "icscore")) %>% 
                      dplyr::mutate(width = nchar(consensus),
                                    n     = 1:n()) %>% 
                      dplyr::mutate(id = paste0(collection.name, "_", id_old, "_n", n))
  }

  
  ## Change the motif IDs, this is required to map the collection of origin of each
  ## motif and to avoid repeated IDs
  # message("; Renaming motif IDs")
  motif.collection.new.id <- purrr::map2(.x = motif.collection,
                                         .y = motif.info.dt$id,
                                         .f = ~set.um.id.name(um     = .x,
                                                              new.id = .y))
  
  return(list(motif_info = motif.info.dt,
              um_object  = motif.collection.new.id))
}


## Export the motifs with updated ID as a transfac file
## The transfac file header exported by universalmotif::write_transfac is not correctly
## recognized by compare-matrices-quick (it only accepts transfac files) so the fields
## have to be renamed to be read by compare-matrices quick
write.transfac.parsed.header <- function(um.object   = NULL,
                                         old.tf.file = NULL,
                                         new.tf.file = NULL,
                                         verbose     = TRUE) {

  suppressWarnings(file.remove(new.tf.file, showWarnings = FALSE))
  suppressWarnings(dir.create(dirname(new.tf.file), recursive = TRUE, showWarnings = FALSE))
  
  if (verbose) {
    message("; Exporting motifs with updated ID: ", old.tf.file)
  }
  suppressMessages(universalmotif::write_transfac(motifs    = um.object,
                                                  file      = old.tf.file,
                                                  overwrite = TRUE))
  
  ## compare-matrices-quick expects the fields AC and ID (instead of ID and NA exported from universalmotif) 
  transfac.lines <- readLines(old.tf.file)
  transfac.lines <- gsub(transfac.lines, pattern = "^ID ", replacement = "AC ")
  transfac.lines <- gsub(transfac.lines, pattern = "^NA ", replacement = "ID ")
  suppressMessages(writeLines(text = transfac.lines,
                              con  = new.tf.file))
  invisible(suppressWarnings(file.remove(old.tf.file, showWarnings = FALSE)))
  
}


## Check that the input file, collection name and motif format are correct
check.input.motif.file <- function(motif.file       = NULL,
                                   motif.collection = NULL,
                                   motif.format     = NULL,
                                   file.line        = NULL) {

  
  ## Check the status of each element in the input table
  motif.file.flag       <- ifelse(file.exists(motif.file), yes = TRUE, no = FALSE)
  motif.collection.flag <- ifelse(motif.collection != "", yes = TRUE, no = FALSE)
  motif.format.flag     <- check.supported.formats(motif.format = motif.format)
  
  ## In case one of the element's flag is false, stop the script and report the line with the problem
  if (all(motif.file.flag, motif.collection.flag, motif.format.flag)) {
    return(TRUE)
  } else {
    
    if (!motif.file.flag) {
      stop("Line ", file.line, ": Input file not found: ", motif.file)
    }
    
    if (!motif.collection.flag) {
      stop("Line ", file.line, ": Invalid collection name")
    }
  }
}


## Verify the status of the files, collection names, and motif formats provided 
## in the --matrix_file_table command-line argument
## Also it reports when motif files or collection names are duplicated
check.status.motif.table <- function(matrix.file.table = NULL) {
  
  message("; Reading input motif file table: ", matrix.file.table)
  matrix.files <- fread(matrix.file.table,
                        header    = FALSE,
                        col.names = c("Motif_file", "Collection_name", "Format")) %>% 
                  dplyr::mutate(Line = 1:n())
  
  matrix.files.no.dup <- matrix.files %>% 
                          distinct()
  
  if (nrow(matrix.files) != nrow(matrix.files.no.dup)) {
    warning("There is a duplicated line in the input matrix file table. This line was removed and will not be considered in the analysis.")
  }
  
  
  matrix.files.list <- list(Motif_file   = matrix.files.no.dup$Motif_file,
                            Collection   = matrix.files.no.dup$Collection_name,
                            Motif_format = matrix.files.no.dup$Format,
                            Line         = matrix.files.no.dup$Line)
  
  ## Verify the input files exist and that motif formats are correct
  motif.lines.flags <- purrr::pmap_lgl(.l = matrix.files.list,
                                       .f = ~check.input.motif.file(motif.file       = ..1,
                                                                    motif.collection = ..2,
                                                                    motif.format     = ..3,
                                                                    file.line        = ..4))
  
  ## Check that input files and collection names are unique
  if (!all(!duplicated(matrix.files.no.dup$Motif_file))) {
    stop("There is a duplicated motif file name")
  }
  
  if (!all(!duplicated(matrix.files.no.dup$Collection_name))) {
    stop("There is a duplicated motif collection name")
  }
  
  return(matrix.files.list)
}



## This function returns two universalmotif objects, each object contains the motifs
## oriented according to the alignment, the second list is the reverse-complement
## 
## un.motifs         = Object of class 'universalmotif'
## orientation.table = a table containing the motif IDs and its final strand after alignment 
motifs.final.orientation <- function(un.motifs         = NULL,    
                                     orientation.table = NULL) {  
  
  ## Assign names to the entries of the universalmotif object
  ## By default the names are NULL, after these lines the object can be indexed through motif IDs
  
  motif.ids.um     <- purrr::map_chr(un.motifs, `[`, "name")
  names(un.motifs) <- motif.ids.um
  
  ## Obtain the RC
  all.motifs.um.D <- un.motifs
  
  ## Check why after this function the motifs are not the same
  all.motifs.um.R <- universalmotif::motif_rc(motifs = un.motifs)
  
  
  ## Get the IDs of motifs in D and R strand after alignment
  final.orientation.tab <- orientation.table %>% 
                            select(id, strand)
  
  motif.id.strandD <- subset(final.orientation.tab, strand == "D")$id
  motif.id.strandR <- subset(final.orientation.tab, strand == "R")$id
  
  # oriented.um    <- c(all.motifs.um.D[motif.id.strandD], all.motifs.um.R[motif.id.strandR])
  oriented.um    <- c(all.motifs.um.D[motif.ids.um %in% motif.id.strandD], all.motifs.um.R[motif.ids.um %in% motif.id.strandR])
  oriented.um.rc <- universalmotif::motif_rc(oriented.um)
  
  return(list(D = oriented.um,
              R = oriented.um.rc))
}



## Stop the script when the input format name is not in the list of supported formats
check.supported.formats <- function(motif.format = NULL) {
  
  supported.formats <- c("cluster-buster", "cisbp", "homer", "jaspar", "meme", "tf", "transfac", "uniprobe")
  
  if (motif.format %in% supported.formats) {
    return(TRUE)
  } else {
    stop(motif.format , " is not a supported format. Suported formats: ", paste(supported.formats, collapse = ", "))
  }
  
}



## Given the path of a motif file and its format load the content as a universalmotif object
read.motif.file <- function(motif.file  = NULL,
                            motif.format = NULL) {
  
  # message("; Reading input file in ", motif.format, " format: ", motif.file)
  
  ## Add here more supported motifs
  um.object <- switch(motif.format,
                      "cluster-buster" = read_cluster_buster(motif.file = motif.file),      ## This is a custom function to read cluster-buster format
                      "cisbp"          = universalmotif::read_cisbp(file = motif.file),
                      "homer"          = universalmotif::read_homer(file = motif.file),
                      "jaspar"         = universalmotif::read_jaspar(file = motif.file),
                      "meme"           = universalmotif::read_meme(file = motif.file),
                      "tf"             = universalmotif::read_transfac(file = motif.file),
                      "transfac"       = universalmotif::read_transfac(file = motif.file),
                      "uniprobe"       = universalmotif::read_uniprobe(file = motif.file))
  
  # Use this flag to avoid some bugs
  # When the input is one motif, R decompress the list and treat the object as a 
  # UniversalMotif object instead of a list, and purrr cannot iterate and crashes
  # Avoid this by re-generating a 1-length list wit the UniversalMotif object
  if (length(um.object) == 1) {
    um.object <- list(um.object)
  }
  
  # Check if the motif has only one site and prevent universalmotif from interpreting it as a frequency matrix
  um.object <- purrr::map(.x = um.object,
                          .f = ~check.one.site.motifs(.x))
  
  return(um.object)
}



# This function is used to check if the motif has only one site
# Universalmotif interprets this a frequency matrix, but it is a site matrix
check.one.site.motifs <- function(um) {
  
  motif.nb.sites     <- um@nsites
  motif.colsum.sites <- unique(universalmotif::colSums(um))
  
  if (motif.nb.sites == 1) {
    if (motif.nb.sites != motif.colsum.sites) {
      um@motif[um@motif > 0] <- 1
      warning("Motif has only one site, it will be interpreted as a count matrix")
    }
  }
  
  return(motifs.um)
}


## Export universalmotif object as cluster-buster motif file
write_cluster_buster <- function(motifs    = NULL,
                                 file      = NULL,
                                 overwrite = TRUE) {
  
  universalmotif::write_matrix(motifs    = motifs,
                               file      = file,
                               positions = "rows",
                               # sep       = "//",
                               overwrite = overwrite,
                               headers   = ">")
  
}



## Export a universalmotif file as a text-file motif
write.motif.file <- function(um.object    = NULL,
                             motif.format = NULL,
                             outfile.name = NULL) {
  
  dir.create(dirname(outfile.name), recursive = TRUE, showWarnings = FALSE)
  
  message("; Exporting motifs in ", motif.format, " format: ", outfile.name)
  
  if (motif.format == "homer") {
    universalmotif::write_homer(motifs    = um.object,
                                file      = outfile.name,
                                overwrite = TRUE)
    
  } else if (motif.format == "jaspar") {
    universalmotif::write_jaspar(motifs    = um.object,
                                 file      = outfile.name,
                                 overwrite = TRUE)
    
  } else if (motif.format == "meme") {
    universalmotif::write_meme(motifs    = um.object,
                               file      = outfile.name,
                               overwrite = TRUE)
    
  } else if (motif.format == "cluster-buster") {

    write_cluster_buster(motifs    = um.object,
                         file      = outfile.name,
                         overwrite = TRUE)
    
  } else if (motif.format %in% c("tf", "transfac")) {
    write.transfac.parsed.header(old.tf.file = paste0(outfile.name, ".tmp"),
                                 new.tf.file = outfile.name,
                                 um.object   = um.object,
                                 verbose     = FALSE) 
  } else {
    stop("; ", motif.format, " exporting function has not yet been implemented.")
  }
}



## Export one motif stored as a universalmotif object as a transfac file
## in a given directory
## NOTE: this function adds the suffix '_oriented' in the motif filename
export.one.motif.transfac <- function(un     = NULL,
                                      outdir = NULL,
                                      strand = "D") {

  ## Suffix based in strand
  strand.suffix <- switch(strand,
                          "D" = "",
                          "R" = "_rc")


  ind.motif.file.path.tmp <- file.path(outdir, paste0(un@name, "_oriented", strand.suffix,".tf.tmp"))
  ind.motif.file.path     <- file.path(outdir, paste0(un@name, "_oriented", strand.suffix,".tf"))

  ## Export transfac file with correct header to be read by compare-matrices-quick
  write.transfac.parsed.header(old.tf.file = ind.motif.file.path.tmp,
                               new.tf.file = ind.motif.file.path,
                               um.object   = un,
                               verbose     = FALSE)
}


## Export the motifs stored in the input universalmotif object as inidividual transfac files
## Each motif is exported in D and R orientations
export.indiv.motif.files <- function(un.motifs       = NULL,
                                     alignment.table = NULL,
                                     outdir          = NULL) {
  
  ## A list containing the oriented motifs in D and R orientations
  um.final.orientation <- motifs.final.orientation(un.motifs         = un.motifs,    
                                                   orientation.table = alignment.table)
  
  
  ## Iterate over the strands
  motifs.D <- um.final.orientation$D
  motifs.R <- um.final.orientation$R
  plan(multisession, workers = params.list$nb_workers)
  options(future.globals.maxSize = 400000000000000000)
  
  for (i in 1:2) {
    
    ## Select strand and universalmotif object
    if (i == 1) { 
      uo.motifs.oriented <- motifs.D
      strand.oriented    <- "D"
      message("; Exporting oriented individual motif files: D strand")
      
    } else if (i == 2) { 
      uo.motifs.oriented <- motifs.R
      strand.oriented    <- "R"
      message("; Exporting oriented individual motif files: R strand")
    }
    
    ## Iterate over each motif to export it as a transfac file
    furrr::future_walk2(.x = uo.motifs.oriented,
                        .y = rep(outdir, times = length(uo.motifs.oriented)),
                        .f = ~export.one.motif.transfac(un     = .x,
                                                        outdir = .y,
                                                        strand = strand.oriented))
  }
}



## Check whether the input file is indeed a transfac file
is.transfac.file <- function(AC  = NULL,
                             #ID  = NULL,
                             P0  = NULL,
                             SEP = NULL) {
  
  ## Minimal fields
  ## AC (beginning of a matrix), ID, P0 (Position 0, before the count matrix), SEP (end of matrix)
  AC  <- length(AC)
  #ID  <- length(ID)
  P0  <- length(P0)
  SEP <- length(SEP)
  
  
  #if (!all(AC, ID, P0, SEP)) {
  if (!all(AC, P0, SEP)) {
    
    stop("Incorrect transfac file format. Verify the following fields are present in your motif file: AC, ID, P0, //")
    
  } else {
    return(TRUE)
  }
  
}



## Return a 4 * nb matrix initialized with 0s
tf.empty.row <- function(nb = 1) {
  
  ## Matrix of 4 * nb initialized with 0s
  empty.rows <- matrix(nrow = nb, ncol = 4, data = 0)
  
  colnames(empty.rows) <- c("AA", "CC", "GG", "TT")
  rownames(empty.rows) <- NULL
  
  return(empty.rows)
  
}



## Add a given number of gaps (rows with 0s) to a transfac file containing a single motif
add.gaps.transfac.motif <- function(tf.file.in  = NULL,
                                    gap.up      = 0,
                                    gap.down    = 0,
                                    tf.file.out = NULL) {
  
  if (!file.exists(tf.file.in)) {
    stop("Transfac file not found: ", tf.file.in)
  }
  
  ## Read transfac file
  tf.lines <- readLines(tf.file.in)
  
  
  ## Save the important lines
  ac.line  <- which(grepl(x = tf.lines, pattern = "^\\s*AC"))
  id.line  <- which(grepl(x = tf.lines, pattern = "^\\s*ID"))
  p0.line  <- which(grepl(x = tf.lines, pattern = "^\\s*P0"))
  sep.line <- which(grepl(x = tf.lines, pattern = "^\\s*//"))
  
  ## Verify the input motif has the minimal transfac fields
  istf <- is.transfac.file(AC  = ac.line,
                           #ID  = id.line,
                           P0  = p0.line,
                           SEP = sep.line)
  
  ## Extract the motif from the transfac file
  xx.lines       <- which(grepl(x = tf.lines, pattern = "^\\s*XX"))
  end.motif.line <- xx.lines[min(which(xx.lines > p0.line))]
  motif.lines    <- tf.lines[(p0.line + 1):(end.motif.line - 1)]
  
  
  ## Convert the motif to a data.frame
  motif.df <- purrr::map(.x = motif.lines,
                         .f = ~strsplit(x     = .x,
                                        split = "\\s+"))
  
  motif.df <- matrix(unlist(purrr::map(motif.df, `[[`, 1)), byrow = TRUE, ncol = 6) %>% 
                  data.table() %>% 
                  dplyr::rename(AA = V2,
                                CC = V3,
                                GG = V4,
                                TT = V5) %>% 
                  select(AA, CC, GG, TT)

  
  ## Add upstream/downstream gaps
  gap.upstream   <- tf.empty.row(nb = gap.up)
  gap.downstream <- tf.empty.row(nb = gap.down)
  
  motif.w.gaps <- rbind(gap.upstream,
                        motif.df,
                        gap.downstream) %>% 
                    dplyr::mutate(LL = 1:n()) %>% 
                    select(LL, AA, CC, GG, TT)
  
  
  ## Reconstruct the transfac file, insert the matrix with gaps
  motif.w.gaps.text <- reconstruct.transfac.file.vector(AC  = tf.lines[ac.line],
                                                        ID  = tf.lines[id.line],
                                                        MAT = motif.w.gaps)
  
  ## Print the vector as a text file
  dir.create(dirname(tf.file.out), recursive = TRUE, showWarnings = FALSE)
  suppressMessages(file.remove(tf.file.out))
  writeLines(motif.w.gaps.text, con = tf.file.out)
  
}



## Return a dataframe where each row is a string
## These are the lines that will be printed in a text file
df.to.vector.text <- function(df = NULL) {
  
  ## Clean the column LL
  df$LL   <- as.numeric(df$LL)
  df.char <- apply(df, 2, as.character)
  
  # when the input df contains only one column, the previous operation returns a character vector
  # instead of a df (why R !?)
  # The following code reconstruct the DataFrame
  if (nrow(df) == 1) {
      df.char <- matrix(df.char) |> t() |> data.frame()
      colnames(df.char) <- c("LL", "AA", "CC",  "GG", "TT")
  } else {
    df.char <- apply(df.char, 2, function(x){
      gsub(x, pattern = "^ ", replacement = "")
    })
  }


  # df = data.frame(LL = 1, AA = 72, CC = 0, GG = 94, TT = 70)
  # df = data.frame(LL = 1:2, AA = 72:73, CC = 0:1, GG = 94:95, TT = 70:71)
  
  ## Convert each row in a text vector, then as a single string using paste
  nt.counts.lines <- apply(df.char, 1, function(x){
    nt.counts <- paste(x[2], x[3], x[4], x[5], sep = "   ")  ## One tab (/t) = 4 spaces
    paste(x[1], nt.counts, sep = "        ")                 ## Double tab
  })
  
  return(nt.counts.lines)
}



## Construct a transfac file as a character vector, each entry in the vector 
## corresponds to a line in the text file
reconstruct.transfac.file.vector <- function(AC  = NULL,
                                             ID  = NULL,
                                             MAT = NULL){
  
  ## Required fields
  tf.comment   <- "XX"
  tf.separator <- "//"
  p0           <- "P0       A     C     G     T"
  
  ## Concat the variables to create a vector
  new.tf.vec <- c(AC,
                  tf.comment,
                  ID,
                  tf.comment,
                  p0,
                  df.to.vector.text(MAT),
                  tf.comment,
                  tf.separator)
  
  return(new.tf.vec)
}




## Returns a table with the motifs that need to be modified by adding gaps (empty rows)
add.gaps.to.indiv.tf.files <- function(motif.folder = NULL,
                                       gap.info     = NULL){
  
  if (!dir.exists(motif.folder)) {
    stop("Motif folder not found: ", motif.folder)
  }
  
  ## Read transfac files in a directory
  transfac.files <- list.files(motif.folder)
  transfac.files <- transfac.files[transfac.files != "input_motifs_parsed_id.tf"]
  
  r.ind <- which(grepl(transfac.files, pattern = "_rc"))
  d.ind <- which(grepl(transfac.files, pattern = "oriented\\.tf"))
  
  ## Separate files by strand
  transfac.files.r <- file.path(motif.folder, transfac.files[r.ind])
  transfac.files.d <- file.path(motif.folder, transfac.files[d.ind])
  
  
  ## Get the motif ID (as part of the file name) to combine the alignment table 
  ## and identify the number of gaps that need to be added 
  gap.file.tab <- rbind(data.frame(file   = transfac.files.d,
                                   strand = "D",
                                   id     = gsub(basename(transfac.files.d), pattern = "_oriented.+$", replacement = "")),
                        data.frame(file   = transfac.files.r,
                                   strand = "R",
                                   id     = gsub(basename(transfac.files.r), pattern = "_oriented.+$", replacement = "")))
  
  gap.file.tab <- gap.file.tab %>% 
                    left_join(gap.info, by = "id") %>% 
                    rename(strand = strand.x) %>% 
                    select(file, strand, offset_up, offset_down) %>% 
                    data.table()
  
  
  ## Invert the offsets for the motifs in reverse complement
  gap.file.tab.list <- split(gap.file.tab, f = gap.file.tab$strand)
  tmp.off                         <- gap.file.tab.list$R$offset_up
  gap.file.tab.list$R$offset_up   <- gap.file.tab.list$R$offset_down
  gap.file.tab.list$R$offset_down <- tmp.off
  gap.file.tab <- rbindlist(gap.file.tab.list) %>% 
                      data.table()
  
  ## Restrict the table to those motifs that need gaps
  gap.file.tab <- gap.file.tab %>% 
                    mutate(sum_gaps = offset_up + offset_down)
                   # dplyr::filter(sum_gaps > 0) %>% 
                   # within(rm(sum_gaps))
  
  
  gap.file.tab$file_rc = gsub(gap.file.tab$file, pattern = "_oriented\\.tf", replacement = "_oriented_rc\\.tf")
  
  return(gap.file.tab)
}



## Read a transfac file and returns a list containing a count matrix in each element
concat.matrices <- function(um.object = NULL) {
  
  ## Calculate length of input motifs and check they have the same width
  motifs.length <- purrr::map_dbl(.x = um.object, `[`, .f = ~nchar(.x@consensus))
  
  if (length(unique(motifs.length)) != 1) {
    stop("The input motifs must have the same length before being merged.")
  }
  
  count.matrices <- purrr::map(.x = um.object, `[`, .f = ~.x@motif)
  
  return(count.matrices)
}




## Return a merged count matrix
## count.matrix.list is the output of the function 'concat.matrices'
reduce.count.matrix.list <- function(count.matrix.list = NULL,
                                     merge.method      = "sum",
                                     normalize.counts  = FALSE,
                                     normalize.factor  = 100) {
  
  ## Check that the combine method is supported
  supported.methods <- c("sum", "average")
  
  if (!merge.method %in% supported.methods) {
    stop("Combine method not recognized. Supported methods: ", paste(supported.methods, collapse = ", "))
    
  }
  
  
  ## Merge the count matrix according to the indicated method
  if (merge.method == "sum") {
    
    count.matrix.reduced <- Reduce('+', count.matrix.list)
    
  } else if (merge.method == "average") {
    
    count.matrix.reduced <- Reduce('+', count.matrix.list) / length(count.matrix.list)
    
  }
  
  colnames(count.matrix.reduced) <- NULL 
  count.matrix.reduced          <- round(count.matrix.reduced)
  
  
  ## Normalize the counts to a given factor (number)
  ## This may be needed because after merging the motifs some columns may not have the same sum
  if (normalize.counts) {
    
    ## Check the normalize factor is a numeric value
    if (!is.numeric(normalize.factor) |
        is.null(normalize.factor)    |
        is.na(normalize.factor)      |
        normalize.factor < 1 ) {
      stop("The normalize factor must be a positive round number")
    }
    normalize.factor <- round(normalize.factor)
    
    
    ## Normalize the matrix
    count.matrix.reduced.norm <- round(scale(count.matrix.reduced,
                                             center = FALSE,
                                             scale  = colSums(count.matrix.reduced)) * normalize.factor)
    count.matrix.reduced      <- count.matrix.reduced.norm
  }
  
  count.matrix.reduced <- data.frame(count.matrix.reduced)
  rownames(count.matrix.reduced) <- c("AA", "CC", "GG", "TT")
  
  if (ncol(count.matrix.reduced) == 1) {
    if (colnames(count.matrix.reduced) == "count.matrix.reduced") {
      colnames(count.matrix.reduced) <- "X1"
    }
  }


  
  return(count.matrix.reduced)
}



## Reduce all the motifs (count matrices) within a transfac file to a single count matrix
## applying a method (sum, average)
create.root.motif <- function(aligned.motif.cluster.file = NULL) {
  
  ## If there is only one motif in the transfac file, the universalmotif object
  ## needs to be converted into a list
  aligned.motifs.um <- universalmotif::read_transfac(file = aligned.motif.cluster.file)

  if (length(aligned.motifs.um) == 1) {
    aligned.motifs.um <- list(aligned.motifs.um)
  }
  
  ## Return the count matrices as a list (one element per cluster member)
  aligned.count.matrices <- concat.matrices(um.object = aligned.motifs.um)
  
  
  ## Reduce the aligned count matrices into a single matrix
  reduced.count.matrix <- reduce.count.matrix.list(count.matrix.list = aligned.count.matrices,
                                                   merge.method      = "sum")
  
  return(reduced.count.matrix)  
}



## Calculate the IC of each positon of a count matrix
calculate.col.IC.count.matrix <- function(count.matrix = matrix(rep(c(10,0,0,0), times = 6), nrow = 6)) {
  
  ## Convert count -> frequency matrix
  freq.matrix <- scale(count.matrix,
                       center = FALSE,
                       scale  = colSums(count.matrix)) 
  
  freq.matrix.log2 <- -log2(freq.matrix)
  freq.matrix.log2[is.infinite(freq.matrix.log2)] <- 0
  
  ## total IC - Sum of p(x) * log2(p(x))
  ## DNA alphabet = IC 2
  2 - colSums(Reduce("*", list(freq.matrix, freq.matrix.log2)))
  
}


## Given a numeric vector representing the IC at each position of a TF binding motif
## return the start and end of the motif after trimming the flanks with low IC (th argument)
## in a window of +/- k positions (this allows to remove IC spikes)
## IC spikes have a separate IC th (sp.th), but by default it is the same as th
motif.trimming <- function(IC.vector = NULL, 
                           th        = 0.25,
                           sp.th     = 0.25,
                           k         = 1) {
  options(stringsAsFactors = FALSE)
  flank_size        <- length(IC.vector)
  ## Positions above general threshold:
  pos.high.IC       <- which(IC.vector >= th)
  pos.binarized.IC  <- IC.vector >= th

  checked_positions <- data.frame()
  
  for (i in seq_along(pos.high.IC)) {
    
    ## Exact IC value for a position:
    pos.idx    <- pos.high.IC[i]
    IC.for.pos <- IC.vector[pos.idx]

    ## If a position is above spike thr, it is immediately kept, otherwise we
    ## run spike detection:
    if (IC.for.pos < sp.th | sp.th == 2) {

      pos.to.check <- pos.idx + seq(-k,k)
      pos.to.check <- pos.to.check[pos.to.check >= 1]
      pos.to.check <- pos.to.check[pos.to.check <= flank_size] ## Selecting the positions in our range only
      pos.to.check <- pos.to.check[pos.to.check != pos.idx]
    
      decision.sum <- sum(pos.binarized.IC[pos.to.check])

    } else {
      decision.sum <- 2*k
    }

    desicion.sum.fraction <- decision.sum/(2*k)

    position_sum_fraction <- cbind(pos.idx, decision.sum, desicion.sum.fraction)
    checked_positions     <- rbind(checked_positions, position_sum_fraction)
    
  }
  
  passed_thr_pos <- checked_positions[which(checked_positions$desicion.sum.fraction >= 0.5), ]
  comment.trim   <- ""
  
  # If no position passed the threshold, raise a Warning and return the whole motif
  if (nrow(passed_thr_pos) == 0) {
    
    min.trim.idx <- 1
    max.trim.idx <- flank_size
    comment.trim <- "Warning: All positions trimmed"
    warning("No position left after trimming. The original motif will be kept.")
    
  # Otherwise trim the motif
  } else {
    
    min.trim.idx <- min(passed_thr_pos$pos.idx)
    max.trim.idx <- max(passed_thr_pos$pos.idx)
    
    ## Trimming:
    # IC_trimmed   <- IC.vector[min.trim.idx:max.trim.idx]
  }
  
  return(data.table(min         = min.trim.idx,
                    max         = max.trim.idx,
                    matrix_size = flank_size,
                    comment     = comment.trim))
}


export.aligned.motifs.per.cluster <- function(indiv.motis.folder    = NULL,
                                              cluster.motifs.folder = NULL,
                                              cluster.motif.id.tab  = NULL) {

					    
  dir.create(indiv.motis.folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(cluster.motifs.folder, showWarnings = FALSE, recursive = TRUE)
  
  ## This table contains the motif files (D and R) with their respective clusters
  motif.id.file.tab <- data.table(id     = unique(gsub(list.files(indiv.motis.folder), pattern = "_oriented\\w{,3}.tf$", replacement = "")),
                                  File_D = file.path(indiv.motis.folder, list.files(indiv.motis.folder, pattern = "_oriented.tf"))) %>% 
                          left_join(cluster.motif.id.tab, by = "id") %>% 
                          mutate(New_File_D = file.path(cluster.motifs.folder, cluster, basename(File_D)))
  
  ## Creates a copy of each motif in the cluster in its corresponding folder
  furrr::future_walk2(.x = motif.id.file.tab$File_D,
                      .y = motif.id.file.tab$New_File_D,
                      .f = ~file.copy(from      = .x,
                                      to        = .y,
                                      overwrite = TRUE))
  
  ## For each cluster, read its motifs and export all together in a single transfac file
  clusters.motif.files.dir <- sort(unique(dirname(motif.id.file.tab$New_File_D)))
  motif.files <- sapply(clusters.motif.files.dir, function(d) {
    
    ## Read all the motif files (transfac format associated to a cluster)
    clusters.transfac.files.um <- unlist(purrr::map(.x = file.path(d, list.files(d, pattern = "_oriented\\w{,3}.tf")),
                                                    .f = ~read_transfac(file = .x)))
    
    cluster.motifs.file <- file.path(d, paste0(basename(d), ".tf"))
    
    ## Export transfac file with correct header to be read by compare-matrices-quick
    # message("Exporting motifs for ", basename(d), " : ", cluster.motifs.file)
    write.transfac.parsed.header(old.tf.file = paste0(cluster.motifs.file, ".tmp"),
                                 new.tf.file = cluster.motifs.file,
                                 um.object   = clusters.transfac.files.um,
                                 verbose     = FALSE)
    
    cluster.motifs.file
    
  })
  
  names(motif.files) <- NULL
  return(motif.files)
}



## The input of this function is a transfac file with all motifs are aligned (may
## contain gaps) and all of them have the same width
## and are aligned
export.root.motif <- function(cluster.tf.file = NULL) {
  
  ## Combine the count matrices to generate a reduced (root) motif
  root.motif <- create.root.motif(aligned.motif.cluster.file = cluster.tf.file)
  
  ## Assign AC and ID fields as cluster names
  root.motif.name <- gsub(basename(cluster.tf.file), pattern = "\\.tf$", replacement = "")
  
  ## Convert the root motif (a matrix) to the format required in the transfac format
  root.motif.tf <- as.matrix.data.frame(root.motif) %>%
                    t() %>% 
                    data.table() %>% 
                    dplyr::mutate(LL = 1:n()) %>% 
                    dplyr::select(LL, AA, CC, GG, TT)
  
  ## Reconstruct a transfac file with the minimal fields
  root.motif.file <- reconstruct.transfac.file.vector(AC  = paste0("AC ", root.motif.name),
                                                      ID  = paste0("ID ", root.motif.name),
                                                      MAT = root.motif.tf)
  return(root.motif.file)
}


## Trimm the input universalmotif object using as reference an IC threshold
## This method also considers the window (+/- k) surrounding a given position to remove 
## spikes in the IC.
trim.motifs.window <- function(um                 = NULL,
                               ic.threshold       = 0.25,
                               ic.spike.threshold = 0.25,
                               window.k           = 1) {

  message("; Trimming motifs. General IC threshold: ", ic.threshold, "; IC threshold for spikes: ", ic.spike.threshold, " in +/- ", window.k, " positions") 
  options(stringsAsFactors = FALSE)
  
  ## Get names of all motifs:
  list_of_names <- purrr::map(um, `[`, "name")
  
  ## Get the count matrices within the Universalmotif object
  ## Then obtain the IC content per column in each count matrix
  count.matrices.list    <- purrr::map(um, `[`, "motif")
  count.matrices.list.ic <- purrr::map(count.matrices.list, calculate.col.IC.count.matrix)
  count.matrices.list.ic <- lapply(count.matrices.list.ic, as.vector)


  ## Calculate positions that should be kept after trimming motifs
  trim.positions <- purrr::map_df(.x = count.matrices.list.ic,
                                  .f = ~motif.trimming(IC.vector = .x,
                                                       th        = ic.threshold,
                                                       sp.th     = ic.spike.threshold,
                                                       k         = 1))
  
  
  positions.trim.list <- list(count_matrices = count.matrices.list,
                              from           = trim.positions$min,
                              to             = trim.positions$max)
  
  ## Subset the motif using the previously calculated positions
  count.matrices.trimmed.list <- purrr::pmap(.l = positions.trim.list,
                                             .f = ~subset.matrix(m   = ..1,
                                                                 min = ..2,
                                                                 max = ..3))
  suppressWarnings(
    count.matrices.trimmed.um <- purrr::map2(.x = um,
                                             .y = count.matrices.trimmed.list,
                                             .f = ~set.um.motif(um        = .x,
                                                                new.count = .y))

  )
  

  trim_values_df <- cbind(unlist(list_of_names),
                          trim.positions$min,
                          trim.positions$max,
                          trim.positions$matrix_size,
                          trim.positions$comment)
  colnames(trim_values_df) <- c("name", "min_position", "max_position", "motif_size", "comment")
  
  # print(trim_values_df)
  trim_values_df <- data.frame(trim_values_df) %>%
                      dplyr::mutate(trim_left  = as.numeric(min_position) -1 ,
                                    trim_right = as.numeric(motif_size) - as.numeric(max_position))
  # print(trim_values_df)
  
  # return(count.matrices.trimmed.um)
  return(list(trimmed_matrix = count.matrices.trimmed.um,
              trim_values    = trim_values_df))
}


export.one.logo <- function(um.motif = NULL,
                            logofile = NULL) {
  
  # Format conversion
  motif <- convert_type(um.motif, "PPM")
  motif <- motif["motif"]
  
  invisible(suppressMessages(suppressWarnings(
    motif.gg <- ggplot() +
                  geom_logo(motif, method = 'bits', stack_width = 1) +
                  theme_logo() +
                  scale_x_continuous(expand = c(0, 0)) +
                  scale_y_continuous(limits = c(0,2), expand = c(0, 0)) +
                  guides(fill = "none")
  )))
  
  ggsave(filename = logofile, plot = motif.gg, bg = "white", width = 10, height = 6, dpi = 400)
}


export.logos <- function(um        = NULL,
                         outfolder = NULL,
                         rev_tag   = FALSE) {
  
  rtag = ifelse(rev_tag, yes = "_rc", no = "")
  
  # The logos must be exported in both orientations
  # The logo name is: motif_ID.jpeg and motif_ID_RC.jpeg
  
  # ------------------------- #
  # Logos Forward orientation #
  # ------------------------- #
  logos.F.name <- file.path(outfolder, paste0(purrr::map_chr(um, `[`, "name"), rtag, ".svg"))
  
  # plan(multisession, workers = params.list$nb_workers)
  purrr::walk2(.x = logos.F.name,
               .y = um,
               .f = ~export.one.logo(um.motif = .y,
                                     logofile = .x))
}
