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



calculate.nbsites.uo <- function(um = NULL) {
  
  um.motifs     <- furrr::future_map(um, `[`, "motif")
  motifs.colsum <- furrr::future_map_dbl(um.motifs, nb.sites.in.um)

  return(motifs.colsum)
}



## This function reads motifs in cluster-buster format and returns a universalmotif object
read_cluster_buster <- function(motif.file = NULL) {
  
  ## Read motif file using universalmotif functions    
  cluster.buster.uo <- list(universalmotif::read_matrix(file = motif.file, positions = "rows", sep = "//"))
  
  ## Set nbsites attribute correctly
  ## By default it is set to 100 and this create a problem when genereting RC because
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
                                      motif.format = motif.format)
  
  
  ## This step is required when there are cases of motifs with empty columns, for some reason
  ## this generates an NA within the Universalmotif object and crashes the script
  motif.collection <- universalmotif::trim_motifs(motif.collection, min.ic = 0.0001)
  
  ## Generate the reverse-complement version of the input motif collection
  motif.collection.rc <- universalmotif::motif_rc(motif.collection)
  
  ## In cases when there is only 1 motif in the motif collection, save the universalmotif
  ## object within a list
  if (length(motif.collection) == 1) {
    motif.collection    <- list(motif.collection)
    motif.collection.rc <- list(motif.collection.rc)
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
write.transfac.pased.header <- function(um.object   = NULL,
                                        old.tf.file = NULL,
                                        new.tf.file = NULL,
                                        verbose     = TRUE) {
  
  if (verbose) {
    message("; Exporting motifs with updated ID: ", old.tf.file)
  }
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
  
  oriented.um    <- c(all.motifs.um.D[motif.id.strandD], all.motifs.um.R[motif.id.strandR])
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
  
  message("; Reading input file in ", motif.format, " format: ", motif.file)
  
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
  
  return(um.object)
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
    write.transfac.pased.header(old.tf.file = paste0(outfile.name, ".tmp"),
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
  write.transfac.pased.header(old.tf.file = ind.motif.file.path.tmp,
                              new.tf.file = ind.motif.file.path,
                              um.object   = un,
                              verbose     = FALSE) 
}



## Export the motifs stored in the input universalmotif object as inidividual transfac files
## Each motif is exported in D and R orientations
export.indiv.motif.files <- function(un.motifs = NULL,
                                     outdir    = NULL) {
  
  ## A list containing the oriented motifs in D and R orientations
  um.final.orientation <- motifs.final.orientation(un.motifs         = un.motifs,    
                                                   orientation.table = results.list$Alignment_table)
  
  
  ## Iterate over the strands
  motifs.D <- um.final.orientation$D
  motifs.R <- um.final.orientation$R
  plan(multisession, workers = params.list$nb_workers)
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
                        .y = rep(outdir, times = length(motifs.D)),
                        .f = ~export.one.motif.transfac(un     = .x,
                                                        outdir = .y,
                                                        strand = strand.oriented))
  }
}