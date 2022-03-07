# matrix-clustering_stand-alone

This is a stand-alone version of *RSAT matrix-clustering*. This version is faster and simplified compared to the original one but the graphical output is still under development.

*RSAT matrix-clustering* is a software for Transcription Factor binding motifs clustering and alignment. Here is a brief description of the method:

  - **Motif comparison**: The motifs are compared to each other using two comparison metrics (pearson correlation coeficient (*cor*) and a alignment-width correction (normalized pearson correlation (*Ncor*)).
  - **Hierarchical clustering**: The motifs are hierarchically clustered based in the values of a comparison metric (default = *Ncor*) .
  - **Tree partition**: the hierarchical tree is partitioned by calculating the average *cor* and *Ncor* values at each node, each time a node does not satisfy the thresholds (one value for *cor* and another for *Ncor*) the node is split in two clusters.
  - **Motif alignment**: for each cluster, the motifs are progressively aligned following the linkage order of the hierarchical tree, this ensures that each motif is aligned in relation to its most similar motif in the cluster. 


Originally, matrix-clustering was planned to be part of the RSAT suite (http://www.rsat.eu/) for motif analysis, we decided to create a portable stand-alone version that can be ran without installing the whole RSAT environment.


## Before starting

If you want to run the original version with all the graphical output, you can do it thorugh RSAT website: http://rsat-tagc.univ-mrs.fr/rsat/matrix-clustering_form.cgi

WARNING: This repository is under active development, so we expect many changes as long as you see this line.

The graphical output (interactive trees and heatmap will be added soon).


## Install required software


### Download this repository

```
git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
cd matrix-clustering_stand-alone
```

### R libraries

```R
##############################################################
## List of R packages requires in the 'minimal_output' mode ##
##############################################################
required.packages = c("dplyr",          ## Data manipulation
                      "data.table",     ## Read long matrices in a quick way
                      "furrr",          ## Run functions in parallel
                      "optparse",       ## Read command-line arguments
                      "purrr",          ## Iterations
                      "rcartocolor",    ## Cluster colors
                      "this.path",      ## Create relative paths
                      "tidyr",          ## Data manipulation
                      "dendsort",       ## To draw heatmap
                      "RColorBrewer",   ## Heatmap cell colors
                      "ape",            ## Export hclust tree in newick format
                      "RJSONIO",        ## Export hclust tree in JSON format


for (lib in required.packages) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
}


############################################
## Install required Bioconductor packages ##
############################################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("universalmotif")
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")
```

### Compile C dependencies

The motif comparison step is ran by `compare-matrices-quick`, a fast version of `RSAT compare-matrices` implemented in C (with less options but very fast).

This repository contains the script written in `C` but it needs to be compiled to generate the executable script that will be called inside `matrix-clustering`.

Assuming you are in the main directory, after cloning the repository:

```bash
cd compare-matrices-quick
make
```

The makefile contains the command to compile the `compare-matrices-quick.c` script, after running the makefile, be sure it generated the executable script

```
./compare-matrices-quick
``` 

This should print the help to run `compare-matrices-quick` and the exaplanation of the parameters, but don't worry you don't have to read it, this script will be called within the `R` scripts.


## How to run

Assuming you are in the root of the repository folder you can run the following example format.

```
Rscript matrix-clustering.R                          \
  -i data/OCT4_datasets/OCT4_motif_table.txt         \
  -o results/OCT4_motifs_example/OCT4_motif_analysis \
  --number_of_workers 8                              
```

In this example three motif collections are clustered. A ChIP-seq dataset was analyzed with three different motif discovery tools (RSAT peak-motifs, MEME-ChIP, and HOMER).


```

```

## Expected output

## Options


## How to cite this software?

If you use this software, please cite its publication: https://academic.oup.com/nar/article/45/13/e119/3862068

```
Castro-Mondragon JA, Jaeger S, Thieffry D, Thomas-Chollier M, and van Helden J. 2017 Nucleic Acids Research
```

## Acknowledgements

We thank the JASPAR curation team for their input to improve RSAT matrix-clustering; the RSAT developer teams for their constant support across many years of collaboration; and the users for their advices, suggestions and reporting bugs.

Special thanks to Ieva Rauluseviciute and Vipin Kumar (from Anthony Mathelier's lab) for testing this software, the discussions, ideas, suggestions of R libraries that make this script faster.
