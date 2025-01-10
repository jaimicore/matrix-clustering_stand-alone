#!/usr/bin/env Rscript

#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("data.table",
                      "dplyr",          ## For data manipulation
                      "furrr",          ## Run functions in parallel
                      "ggplot2",
                      "ggseqlogo",      ## Draw logos
                      "optparse",       ## Read arguments from commmand-line
                      "this.path",      ## Create relative paths
                      "universalmotif") ## Motif manipulation (Bioconductor)

for (lib in required.packages) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}


## 'here' is set to the path where matrix-clustering.R is located, then the libraries
## are sourced relative to where 'here' was set
source(this.path::here(.. = 0, "R", "Motif_manipulation.R"))


####################
## Read arguments ##
####################
option_list = list(
  
  make_option(c("-i", "--matrix_file"), type = "character", default = NULL, 
              help = "A file containing TF binding motifs (Mandatory). ", metavar = "character"),
  
  make_option(c("-o", "--output_file"), type = "character", default = NULL, 
              help = "Folder to save the results (Mandatory)", metavar = "character"),
  
  make_option(c("-l", "--logos"), type = "logical", default = FALSE, 
              help = "Indicates whether the logos are exported or not", metavar = "logical"),
  
  make_option(c("--from"), type = "character", default = NULL, 
              help = "Format of matrices in --matrix_file (Mandatory). Options: cluster-buster, cisbp, homer, jaspar, meme, tf, transfac, uniprobe].", metavar = "character"),
  
  make_option(c("--to"), type = "character", default = "transfac", 
              help = "Output format of motifs. [Default: \"%default\" . Options: cluster-buster, cisbp, homer, jaspar, meme, tf, transfac, uniprobe]", metavar = "character"),
  
  make_option(c("--rc"), type = "logical", default = FALSE, 
              help = "Return the reverse-complement of the input motifs. [Default \"%default\"] ", metavar = "logical"),
  
  make_option(c("-t", "--trim"), type = "logical", default = FALSE,
              help = "Trim the motif flanks. [Default \"%default\"] ", metavar = "logical"),
  
  make_option(c("--IC_threshold"), type = "numeric", default = 0.25, 
              help = "IC threshold to trim the motfis. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("--spike_IC_threshold"), type = "numeric", default = 0.25, 
              help = "IC threshold for the trimming spikes. [Default \"%default\"] ", metavar = "number"),
  
  make_option(c("--trim_values_output"), type = "character", default = "trim_values.txt", 
              help = "A path to a text file to output trim values [Default \"%default\"] ", metavar = "character")
  
  );

message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Output file prefix
tf.matrix.file.in    <- opt$matrix_file
tf.matrix.file.out   <- opt$output_file
logos.flag           <- opt$logos
format.in            <- opt$from
format.out           <- opt$to
rc.flag              <- opt$rc
trim.flag            <- opt$trim
trim.ic              <- as.numeric(opt$IC_threshold)
trim.spike.ic        <- as.numeric(opt$spike_IC_threshold)
trim.values.out.file <- opt$trim_values_output


## Mandatory input
if (!file.exists(tf.matrix.file.in)) {
  stop("Mandatory input file not found: --matrix_file")
}

if (is.null(tf.matrix.file.out)) {
  stop("Missing mandatory parameter: --output_file")
}

if (is.null(format.in)) {
  stop("Missing mandatory parameter: --from")
}

## Skip the trimming step when the IC threshold is not correct
## For example, if it is a character or numeric values out of the IC range (0-2)
if (trim.ic < 0 | is.character(trim.ic) | is.na(trim.ic) | trim.ic > 2) {
  warning("; Incorrect value in --IC_threshold argument. The value must be a number > 0 and < 2. The trimming step will not be applied.")
  trim.flag <- FALSE
}

if (trim.spike.ic < 0 | is.character(trim.spike.ic) | is.na(trim.spike.ic) | trim.spike.ic > 2) {
  warning("; Incorrect value in --spike_IC_threshold argument. The value must be a number > 0 and < 2. The trimming step will not be applied.")
  trim.flag <- FALSE
}


###########
## Debug ##
###########
## Example:
## Rscript convert-matrix.R -i data/OCT4_datasets/RSAT_OCT4_motifs.tf --from tf --to meme -o RSAT_OCT4_motifs.meme
# tf.matrix.file.in  <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/data/OCT4_datasets/RSAT_OCT4_motifs.tf"
# tf.matrix.file.in  <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/data/example/one_example_motif.tf"


# tf.matrix.file.in  <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/Debug/Ieva/pbm_example.jaspar"
# tf.matrix.file.out <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/Debug/Ieva/trimmed_motif.meme"
# format.in            <- "jaspar"
# format.out           <- "meme"

# rc.flag              <- TRUE
# trim.flag            <- TRUE
# trim.ic              <- 0.3
# trim.spike.ic        <- 0.6
# trim.values.out.file <- "trim_values.txt"
# logos.flag           <- TRUE


# tf.matrix.file.in  <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/data/example/tf2meme.tf"
# tf.matrix.file.out <- "/home/jaime/Documents/Personal/matrix-clustering_stand-alone/Debug/Ieva/tf2meme.meme"
# format.in          <- "tf"
# format.out         <- "meme"

## Format in/out test
fit <- check.supported.formats(motif.format = format.in)
fot <- check.supported.formats(motif.format = format.out)


#####################
## Read motif file ##
#####################
suppressWarnings(
  motifs.um <- read.motif.file(motif.file   = tf.matrix.file.in,
                               motif.format = format.in)
)
one.motif.input.flag <- ifelse(length(motifs.um) == 1, TRUE, FALSE)


#################
## Trim motifs ##
#################

if (trim.flag) {

  trimming_res <- trim.motifs.window(um                 = motifs.um,
                                     ic.threshold       = trim.ic,
                                     ic.spike.threshold = trim.spike.ic,
                                     window.k           = 1)
  
  motifs.um   <- trimming_res$trimmed_matrix
  trim_values <- trimming_res$trim_values
  
  dir.create(dirname(trim.values.out.file), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.frame(trim_values),
         file      = trim.values.out.file,
         col.names = TRUE,
         sep       = "\t")
}


##################
## RC conversion #
##################

# Rename one motif UniversalMotif object as the provided dirname
# Only when the input is 1 motif
if (one.motif.input.flag) {
  
  um.name    <- basename(fs:::path_ext_remove(tf.matrix.file.out))

  motifs.um[[1]]@name    <- um.name
  motifs.um[[1]]@altname <- um.name
 
}
motifs.um.rc <- universalmotif::motif_rc(motifs.um)

if (rc.flag) {
  message("; Generating reverse-complement")

  ## Add the suffix '_rc' before the file format
  ## Example: 'RSAT_OCT4_motifs.tf' becomes 'RSAT_OCT4_motifs_rc.tf' 
  tf.matrix.file.out.rc <- gsub(tf.matrix.file.out, pattern = "^(.+)(\\.\\D+)$", replacement = "\\1_rc\\2")
}


#######################
## Export motif file ##
#######################
write.motif.file(um.object    = motifs.um,
                 motif.format = format.out,
                 outfile.name = tf.matrix.file.out)

if (rc.flag) {
  message("; Exporting reverse-complement motifs in ", format.out, " format: ", tf.matrix.file.out.rc)
  write.motif.file(um.object    = motifs.um.rc,
                   motif.format = format.out,
                   outfile.name = tf.matrix.file.out.rc)
}


if (logos.flag) {
  
  message("; Exporting logos ")
  
  logos.dir <- file.path(dirname(tf.matrix.file.out), "logos")
  dir.create(logos.dir, recursive = TRUE, showWarnings = FALSE)

  export.logos(um        = motifs.um,
               outfolder = logos.dir,
               rev_tag   = FALSE)
  
  export.logos(um        = motifs.um.rc,
               outfolder = logos.dir,
               rev_tag   = TRUE)
  
}

message("; End of program")

# Rscript convert-matrix.R  -i data/example/one_example_motif.tf   --from tf --to jaspar --output_file results/convert-matrix_examples/DEBUG_IEVA/IEVA_MOTIF.jaspar --logos TRUE