#############################
## Load required libraries ##
#############################

## List of packages to install from CRAN
required.packages = c("furrr",          ## Run functions in parallel
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
  
  make_option(c("--from"), type = "character", default = NULL, 
              help = "Format of matrices in --matrix_file (Mandatory).", metavar = "character"),
  
  make_option(c("--to"), type = "character", default = "transfac", 
              help = "Output format of motifs. [Default: \"%default\" . Options: cluster-buster, cisbp, homer, jaspar, meme, tf, transfac, uniprobe]", metavar = "character"),
  
  make_option(c("--rc"), type = "logical", default = FALSE, 
              help = "Return the reverse-complement of the input motifs. [Default \"%default\"] ", metavar = "logical")
  
  );

message("; Reading arguments from command-line")
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


## Output file prefix
tf.matrix.file.in  <- opt$matrix_file
tf.matrix.file.out <- opt$output_file
format.in          <- opt$from
format.out         <- opt$to
rc.flag            <- opt$rc

## Mandatory input
if (!file.exists(tf.matrix.file.in)) {
  stop("Mandatory input file not found: ", tf.matrix.file.in)
}

if (is.null(tf.matrix.file.out)) {
  stop("Missing mandatory parameter: ", tf.matrix.file.out)
}

if (is.null(format.in)) {
  stop("Missing mandatory parameter: ", format.in)
}


###########
## Debug ##
###########
## Example:
## Rscript matrix-clustering.R -i data/OCT4_datasets/OCT4_motif_table.txt -o results/OCT4_motifs_example/OCT4_motif_analysis --number_of_workers 8 -q compare-matrices-quick/compare-matrices-quick --minimal_output FALSE -r ./R
tf.matrix.file.in <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/RSAT/matrix-clustering_stand-alone/data/OCT4_datasets/RSAT_OCT4_motifs.tf"
format.in         <- "tf"
format.out        <- "meme"
format.out        <- "homer"
format.out        <- "jaspar"
rc.flag           <- TRUE

## Format in/out test
fit <- check.supported.formats(motif.format = format.in)
fot <- check.supported.formats(motif.format = format.out)


#####################
## Read motif file ##
#####################
motifs.uo <- read.motif.file(motif.file   = tf.matrix.file.in,
                             motif.format = format.in)


##################
## RC conversion #
##################
if (rc.flag) {
  message("; Generating reverse-complement")
  motifs.uo.rc <- universalmotif::motif_rc(motifs.uo)
  
  ## Add the suffix '_rc' before the file format
  ## Example: 'RSAT_OCT4_motifs.tf' becomes 'RSAT_OCT4_motifs_rc.tf' 
  tf.matrix.file.out.rc <- gsub(tf.matrix.file.in, pattern = "^(.+)(\\.\\D+)$", replacement = "\\1_rc\\2")
}


#######################
## Export motif file ##
#######################
message("; Exporting motifs in ", format.out, " format: ", tf.matrix.file.out)

if (rc.flag) {
  message("; Exporting reverse-complement motifs in ", format.out, " format: ", tf.matrix.file.out.rc)
  write.motif.file(um.object    = motifs.uo.rc,
                   motif.format = format.out,
                   outfile.name = tf.matrix.file.out.rc)
}

message("; End of program")