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


## This function reads motifs in cluster-buster format and returns a universalmotif object
read_cluster_buster <- function(motif.file = NULL) {
  
  ## Read motif file using universalmotif functions    
  cluster.buster.uo <- list(universalmotif::read_matrix(file = motif.file, positions = "rows", sep = "//"))
  
  ## Add name and altternate name fields
  cluster.buster.uo.id <-  purrr::map_chr(cluster.buster.uo, `[`, "name")
  
  cluster.buster.uo.new.altname <- purrr::map2(.x = cluster.buster.uo,
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
  motif.collection <- switch(motif.format,
                             "cluster-buster" = read_cluster_buster(motif.file = motif.file),
                             "cisbp"          = universalmotif::read_cisbp(file = motif.file),
                             "homer"          = universalmotif::read_homer(file = motif.file),
                             "jaspar"         = universalmotif::read_jaspar(file = motif.file),
                             "meme"           = universalmotif::read_meme(file = motif.file),
                             "tf"             = universalmotif::read_transfac(file = motif.file),
                             "transfac"       = universalmotif::read_transfac(file = motif.file),
                             "uniprobe"       = universalmotif::read_uniprobe(file = motif.file))
  
  
  ## This step is required when there are cases of motifs with empty columns, for some reason
  ## this generates an NA within the Universalmotif object and crashes the script
  motif.collection <- universalmotif::trim_motifs(motif.collection, min.ic = 0.0001)
  
  ## Generate the reverse-complement version of the input motif collection
  motif.collection.rc <- universalmotif::motif_rc(motif.collection)
  
  ## In cases when there is only 1 motif in the motif collection, save the universalmotif
  ## object within a list
  if (length(motif.collection) == 1) {
    motif.collection <- list(motif.collection)
  }
  
  ## Returns a motif information data.table
  motif.info.dt <- data.table(id_old       = purrr::map_chr(motif.collection, `[`, "name"),
                              name         = purrr::map_chr(motif.collection, `[`, "altname"),
                              rc_consensus = purrr::map_chr(motif.collection.rc, `[`, "consensus"),
                              consensus    = purrr::map_chr(motif.collection, `[`, "consensus"),
                              nb_sites     = round(purrr::map_dbl(motif.collection, `[`, "nsites")),
                              IC           = purrr::map_dbl(motif.collection, `[`, "icscore")) %>% 
                    mutate(width = nchar(consensus),
                           n     = 1:n()) %>% 
                    mutate(id = paste0(collection.name, "_", id_old, "_n", n))
  
  
  ## Change the motif IDs, this is required to map the collection of origin of each
  ## motif and to avoid repeated IDs
  message("; Renaming motif IDs")
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
export.tf.motifs.w.parsed.header <- function(um.object   = NULL,
                                             old.tf.file = NULL,
                                             new.tf.file = NULL) {
  
  message("; Exporting motifs with updated ID: ", old.tf.file)
  universalmotif::write_transfac(motifs    = um.object,
                                 file      = old.tf.file,
                                 overwrite = TRUE)
  
  ## compare-matrices-quick expects the fields AC and ID (instead of ID and NA exported from universalmotif) 
  transfac.lines <- readLines(old.tf.file)
  transfac.lines <- gsub(transfac.lines, pattern = "^ID ", replacement = "AC ")
  transfac.lines <- gsub(transfac.lines, pattern = "^NA ", replacement = "ID ")
  writeLines(text = transfac.lines,
             con  = new.tf.file)
  file.remove(old.tf.file)
  
}


## Check that the input file, collection name and motif format are correct
check.input.motif.file <- function(motif.file       = NULL,
                                   motif.collection = NULL,
                                   motif.format     = NULL,
                                   file.line        = NULL) {
  
  supported.motif.formats <- c("cluster-buster",
                               "cisbp",
                               "homer",
                               "jaspar",
                               "meme",
                               "tf",
                               "transfac",
                               "uniprobe")
  
  ## Check the status of each element in the input table
  motif.file.flag       <- ifelse(file.exists(motif.file), yes = TRUE, no = FALSE)
  motif.collection.flag <- ifelse(motif.collection != "", yes = TRUE, no = FALSE)
  motif.format.flag     <- ifelse(motif.format %in% supported.motif.formats, yes = TRUE, no = FALSE)
  
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
    
    if (!motif.format.flag) {
      stop("Line ", file.line, ": Invalid motif format. upported formats: ", paste(supported.motif.formats, collapse = ", "))
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
    mutate(Line = 1:n())
  
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
