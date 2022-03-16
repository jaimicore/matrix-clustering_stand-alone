# matrix-clustering_stand-alone

This is a stand-alone version of *RSAT matrix-clustering*. This version is faster and simplified compared to the original one but the graphical output is still under development.

*RSAT matrix-clustering* is a software for Transcription Factor binding motifs clustering and alignment. Here is a brief description of the method:

  - **Motif comparison**: The motifs are compared to each other using two comparison metrics (pearson correlation coeficient (*cor*) and a alignment-width correction (normalized pearson correlation (*Ncor*)).
  - **Hierarchical clustering**: The motifs are hierarchically clustered based in the values of a comparison metric (default = *Ncor*) .
  - **Tree partition**: the hierarchical tree is partitioned by calculating the average *cor* and *Ncor* values at each node, each time a node does not satisfy the thresholds (one value for *cor* and another for *Ncor*) the node is split in two clusters.
  - **Motif alignment**: for each cluster, the motifs are progressively aligned following the linkage order of the hierarchical tree, this ensures that each motif is aligned in relation to its most similar motif in the cluster. 

As in the original version of *RSAT matrix-clustering*, there is no limit in the input motif files (so far we have tried up to 900 input files). When users have two or more input files, some intersection statistics are calculated (e.g., overlap among input collections) visualized as heatmaps.

Originally, *RSAT matrix-clustering* was planned to be part of the [*RSAT* suite](http://www.rsat.eu/) for motif analysis, we decided to create a portable stand-alone version that can be ran without installing the whole *RSAT* environment and that can be easily integrated within pipelines.



## Before starting

If you want to run the original version with all the graphical output, you can do it thorugh [*RSAT* website](http://rsat-tagc.univ-mrs.fr/rsat/matrix-clustering_form.cgi) or alternatively, installing *RSAT* locally and run the command line version of *matrix-clustering*.

:warning: This repository is under active development, so we expect many changes as long as you see this line.

The graphical output (interactive trees and heatmaps will be added soon).



## :computer: Install required software


### Download this repository

```
git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
cd matrix-clustering_stand-alone
```

### R libraries

The following R/Bioconductor packages are required to run *RSAT matrix-clustering*, you can install them within `R` using the following commands

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
                      "flexclust)       ## Calculate adjusted rand index


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

The makefile contains the commands to compile the `compare-matrices-quick.c` script, after running the makefile, be sure it generated the executable script

```
./compare-matrices-quick
``` 

This should print the help to run `compare-matrices-quick` and the exaplanation of the parameters, but don't worry you don't have to read it, this script will be called within the `R` scripts.



## :arrow_forward: Qucikstart

Assuming you are in the root of the repository folder you can run the following examples. The input files are part of the repository, they are found in the folder `data`.

### Example 1

Clustering of 66 motifs separated in three motif collections (files). An [Oct4](https://doi.org/10.1016/j.cell.2008.04.043) ChIP-seq dataset was analyzed with three different motif discovery tools (*RSAT peak-motifs*, *MEME-ChIP*, and *HOMER*), the resulting motifs are used as input and we detected a cluster of Oct4 motifs, including the canonical motif, and other binding variants including homodimers and heterodimers, see [Fig 2 of the *RSAT matrix-clustering* paper](https://doi.org/10.1093/nar/gkx314) for a detailed explanation.

:hourglass_flowing_sand: Running time: ~1 minute

```bash
Rscript matrix-clustering.R                          \
  -i data/OCT4_datasets/OCT4_motif_table.txt         \
  -o results/OCT4_motifs_example/OCT4_motif_analysis \
  -w 8                              
```

### Example 2

We cluster the [*JASPAR 2022 plants*](https://jaspar.genereg.net/matrix-clusters/plants/) motif collection (656 motifs), we compare the resulting clusters detected by *RSAT matrix-clustering* against a user-provided reference annotation (in this case the Transcription Factor classes). We calculated the *Adjusted Rand Index* (ARI), a single-value metric (ranging from -1 to +1) indicating the proportion of consistent pairs between two classifications, in this example the ARI measures the proportion of motif pairs that are consistently classified between *RSAT matrix-clustering* results and the reference TF classes. We consider that a motif pair is consistently classified when the two motifs either belong to the same class and are co-clustered, or belong to different families and are not co-clustered.

The calculatio of ARI is a new functionality of *RSAT matrix-clustering*, in this example, using default paramters we obtained `ARI = 0.38`, changing parameters may increase/decrease the resulting ARI.

:hourglass_flowing_sand: Running time: ~5 minutes

```bash
Rscript matrix-clustering.R                         \
  -i data/JASPAR_2022/Jaspar_plants_motifs_tab.txt  \
  -o results/Jaspar_plants/Jaspar_plants            \
  -w 8                                              \
  -r data/JASPAR_2022/Jaspar_2022_plants_TF_fam.tab 
```



## :wrench: Arguments and Options


### Mandatory arguments

 - `-i` or `--matrix_file_table` : A text-delimited file where each line contain the following fields/columns. It does not expect a header, but it expects these columns in the indicated order.
 
    1. Motif file path
    2. Motif collection name
    3. Motif format. 

- `-o` or `--output_folder` : Folder to save the results.


### Comparison + Clustering arguments

- `-m` or `--comparison_metric` : Comparison metric used to build the hierarchical tree. Default: `Ncor`. [Options: `cor`, `Ncor`].
- `-l` or `--linkage_method` : Linkage/agglomeration method to build the hierarchical tree. Default: `average`. [Options: `average`, `complete`, `single`].
- `-c` or `--cor_th` : Pearson correlation lower threshold. Default: `0.75`. [Options: any value among `-1` and `+1`].
- `-n` or `--Ncor_th` : Normalized Pearson correlation lower threshold. Default: `0.55`. [Options: any value among `-1` and `+1`].


### Output files

- `-m` or `--minimal_output` : Only returns the aligned motifs and the alignment, clusters and motif description tables. Comparison results, plots and trees are not exported. Default : `FALSE`. [Options: `TRUE`, `FALSE`].
- `--export_newick` : Export hierarchical tree in Newick format. Default : `FALSE`. [Options: `TRUE`, `FALSE`].
- `--export_heatmap` : Export heatmap with clusters in PDF. Default : `FALSE`. [Options: `TRUE`, `FALSE`].


### Clusters vs Reference annotation

- `-r` or `--reference_cluster_annotation`: User defined cluster annotation tab. Calculates the Adjusted Rand Index (ARI) of the resulting clusters against the reference groups. If the input motifs are separated in many collections, concatenate all the annotations in a single file. A tab-delimited file with two columns, it does not expect a header, but it expects these columns in the indicated order:

    1. Motif ID
    2. Reference group

### Others

- `-w` or `--number_of_workers` : Number of cores to run in parallel. Default: `2`. [Options: depends in your machine].
- `--heatmap_color_palette` : Cell colors in clusters heatmap. Default: `RdGy`. [Options: any colorBrewer palette, see colorbrewer2.org for details ].
- `--color_palette_classes` : Number of classes to create color palette in clusters heatmap. Default: `11`. [Options: depends on the selected colorBrewer palette, see colorbrewer2.org for details ].
  
  
  
## :collision: Contact + Report issues 

This repository is maintained by [Jaime A Castro-Mondragon](https://jaimicore.github.io/)

:e-mail: j.a.c.mondragon@ncmm.uio.no

Twitter: [@jaimicore](https://twitter.com/jaimicore)

Use this space to report [issues](https://github.com/jaimicore/matrix-clustering_stand-alone/issues) related to this repository.



## :pushpin: To Do/Wishlist

- When calculating the ARI, implement an option to find an optimal threshold thorugh a grid search approach.
- Generate the interactive `html` output motif trees.
- Implement the option to annotate clusters.
- Export root motifs.
- Detect the central motif within each cluster.
- Export motif collection intersection stats.



## :tada: Acknowledgements

We thank the [*JASPAR*](https://jaspar.genereg.net/) curation team for their input to improve *RSAT matrix-clustering*; the [*RSAT developer team*](http://rsat-tagc.univ-mrs.fr/rsat/people.php) for their constant support across many years of collaboration; and the users for their advices, suggestions and reporting bugs :beetle:.

Special thanks to my colleagues Ieva Rauluseviciute (and her *gently reminders* :unamused: that pushed me to write this stand-alone version) and Vipin Kumar (both from [Anthony Mathelier's lab](https://mathelierlab.com/)) for testing this software, the discussions, ideas, and their suggestions of `R` libraries that make this script faster than the original version.



## :page_with_curl: How to cite this software?

If you use this software, please cite [its own publication](https://doi.org/10.1093/nar/gkx314) and/or the [latest *RSAT* publication](https://doi.org/10.1093/nar/gky317).

```
Castro-Mondragon JA, Jaeger S, Thieffry D, Thomas-Chollier M, and van Helden J. 2017. RSAT matrix-clustering: dynamic exploration and redundancy reduction of transcription factor binding motif collections. Nucleic Acids Research

Nguyen NTT, Contreras-Moreira B, et al. 2018. RSAT 2018: regulatory sequence analysis tools 20th anniversary. Nucleic Acids Research
```
