# matrix-clustering_stand-alone

This is a stand-alone version of *RSAT matrix-clustering*. This version is faster and with more features compared to the original one but the graphical output is still under development.

*RSAT matrix-clustering* is a software to cluster and align Transcription Factor Binding Motifs (*TFBM*s). Here is a brief description of the method:

  - **Motif comparison**: The motifs are compared to each other using two comparison metrics (pearson correlation coeficient (*cor*) and an alignment-width correction (normalized pearson correlation (*Ncor*)).
  - **Hierarchical clustering**: The motifs are hierarchically clustered based in the values of a comparison metric (default = *Ncor*) .
  - **Tree partition**: the hierarchical tree is partitioned by calculating the average *cor* and *Ncor* values at each node, each time a node does not satisfy any of these thresholds (one value for *cor* and another one for *Ncor*) the node is split in two clusters.
  - **Motif alignment**: for each cluster, the motifs are progressively aligned following the linkage order of the hierarchical tree, this ensures that each motif is aligned in relation to its closest motif in the cluster.
  - **Radial alignment (optional)**: all the motifs are forced to be aligned, the aligned logos are displayed in a radial (circular) tree. This option is useful to visualize entire motif collections. See example in the [*JASPAR* website](https://jaspar.genereg.net/matrix-clusters/nematodes/).

There is no limit in the input motif files (so far we have tried with up to 900 input files), each file is treated as a *collection*. When an analysis contains two or more *collections*, some intersection statistics are calculated (e.g., overlap among input collections) visualized as heatmaps.

*RSAT matrix-clustering* is part of the [*RSAT* suite](http://www.rsat.eu/) for motif analysis, we decided to create a portable stand-alone version that can be ran without installing the whole *RSAT* environment and that can be easily integrated within pipelines.

&nbsp;
&nbsp;


## Before starting

If you want to run the original version with all the graphical output, you can do it through the [*RSAT* website](https://rsat.france-bioinformatique.fr/fungi/matrix-clustering.cgi) or alternatively, installing *RSAT* locally and run the command line version of *matrix-clustering*.

:warning: This repository is under active development, so you can expect many changes as long as you see this line.

&nbsp;

## :wrench: Changes relative to the original version

- :placard: Functionality to create motif alignments in a radial (circular) representation (see **Example 2**). This representation allows to visualize entire motif collections (a motif browser) and highlight categories (TF classes, TF families, Motif collection, etc).

- :placard: Functionality to calculate how well the resulting clusters are similar to a user provided annotation (see **Example 3** for more details). This functionality could be used to select the parameters (thresholds in `cor` and `Ncor`) that maximizes the similarity to a user-provided annotation.

- :placard: Default thresholds are slightly different: `cor = 0.75` and `Ncor = 0.55`. To decide if a node in the hierarchical tree will be merged or split, we compute the average `cor` and `Ncor` of all the pairwise comparisons for all the motifs in a particualr node. We realized that the original version didn't considered all the pairwise comparisons, we corrected this problem, but now the original default thresholds are too permissive, so we updated them to obtain good results.

- :placard: We implemented a motif trimming algorithm that is robust to IC spikes, see [this reference](https://academic.oup.com/nar/article/52/D1/D174/7420101) for the details. Motifs can be trimmed before clustering, more of this in the **Extra** section. 

&nbsp;
&nbsp;

## :computer: Install required software


### Download this repository

```
git clone https://github.com/jaimicore/matrix-clustering_stand-alone.git
cd matrix-clustering_stand-alone
```

&nbsp;

### Docker container

The tools in this repository are available as a container. You can build your own container using requirement and Docker file in `docker` dorectory. Or you can use an [existing image](https://hub.docker.com/repository/docker/cbgr/matrix_clustering/general): `docker pull cbgr/matrix_clustering:1.0.0`.

&nbsp;

### R libraries

The following R/Bioconductor packages are required to run *RSAT matrix-clustering*, you can install them within `R` using the following commands

```R
# --------------------------- #
# List of required R packages #
# --------------------------- #
required.packages = c("dplyr",          # Data manipulation
                      "data.table",     # Read long matrices in a quick way
                      "furrr",          # Run functions in parallel
                      "optparse",       # Read command-line arguments
                      "purrr",          # Iterations
                      "rcartocolor",    # Cluster colors
                      "reshape2",       # Dataframe manipulation
                      "this.path",      # Create relative paths
                      "tidyr",          # Data manipulation
                      "dendsort",       # To draw heatmap
                      "ggplot2",
                      "ggseqlogo",      # Draw logos
                      "RColorBrewer",   # Heatmap cell colors
                      "ape",            # Export hclust tree in newick format
                      "RJSONIO",        # Export hclust tree in JSON format
                      "rjson",
                      "circlize",       # Required to draw heatmaps
                      "flexclust",      # Calculate adjusted rand index
                      "htmlwidgets",    # Save plotly output as html
                      "plotly",         # Interactive plots
                      "svglite",        # Easy export of ggplot content as svg
                      "jsonlite")       # To create the JSON file from the hclust outputs


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
BiocManager::install("ComplexHeatmap")
```

&nbsp;

### Compile C dependencies

The motif comparison step is ran by `compare-matrices-quick`, a fast version of `RSAT compare-matrices` implemented in `C` (with limited options but much faster).

This repository contains the script written in `C` but it needs to be compiled to generate the executable script that will be called inside `matrix-clustering.R`.

Assuming you are in the main directory, after cloning this repository:

```bash
cd compare-matrices-quick
make
```

The `makefile` contains the commands to compile the `compare-matrices-quick.c` script, after running the `makefile`, be sure the next executable script was created by running the following command:

```
./compare-matrices-quick
``` 

This should print the help to run `compare-matrices-quick` and the exaplanation of the parameters, but don't worry you don't have to read it, this script will be called within the `R` scripts.

&nbsp;
&nbsp;

## :arrow_forward: Quickstart

Assuming you are in the root of the repository folder you can run the following examples. The input files are part of the repository, they are found in the folder `data`.

### Example 0

The most common case to use this tool is when you have one motif collection as input. In this case we cluster a collection of motifs discovered by [*RSAT peak-motifs*](https://doi.org/10.1093/nar/gkr1104) in an [*Oct4*](https://doi.org/10.1016/j.cell.2008.04.043) ChIP-seq dataset.

Here we use the arguments `--M` (motif file path), `--C` (collection name), and `--F` (motif format) because the input is one motif file. In cases when you have two or more collections (input files), you can use the option `-i`, as shown in the other Examples.

:hourglass_flowing_sand: Running time: ~1 minute

```bash
Rscript matrix-clustering.R                  \
  --M data/OCT4_datasets/RSAT_OCT4_motifs.tf \
  --C RSAT_OCT4                              \
  --F tf                                     \
  -o results/OCT4_motifs/OCT4_motif_analysis \
  -w 8                              
```

### Example 1

Clustering of 66 motifs separated in three motif collections (files). An [*Oct4*](https://doi.org/10.1016/j.cell.2008.04.043) ChIP-seq dataset was analyzed with three different motif discovery tools ( [*RSAT peak-motifs*](https://doi.org/10.1093/nar/gkr1104), [*MEME-ChIP*](https://doi.org/10.1093/bioinformatics/btr189), and [*HOMER*](https://doi.org/10.1016/j.molcel.2010.05.004)), the resulting motifs are used as input and we detected a cluster of *Oct4* motifs, including the canonical motif, and other binding variants including homodimers and heterodimers, see [Fig 2 of the *RSAT matrix-clustering* paper](https://doi.org/10.1093/nar/gkx314) for a detailed explanation.

:hourglass_flowing_sand: Running time: ~1 minute

```bash
Rscript matrix-clustering.R                           \
  -i data/OCT4_datasets/OCT4_motif_table.txt          \
  -o results/OCT4_motifs_clusters/OCT4_motif_analysis \
  -w 8                              
```

&nbsp;


### Example 2

In this example we are reproducing the clustering of [JASPAR nematodes](https://jaspar.uio.no/matrix-clusters/nematodes/) (43 motifs corresponding to 12 TF classes). We use the option `--radial_tree = TRUE` to force the alignment of all motifs (as if they were in a single cluster), this alignment is displayed in a radial (circular) visualization.

Users can provide motif metadata with the option `-a`.

:hourglass_flowing_sand: Running time: ~2 minutes

```bash
Rscript matrix-clustering.R                           \
  -i data/JASPAR_2022/Jaspar_nematodes_motifs_tab.txt \
  -o results/JASPAR_nematodes_radial/JASPAR_nematodes \
  -a data/JASPAR_2022/JASPAR_nematodes_metadata.txt   \
  --radial_tree TRUE                                  \
  --title JASPAR_CORE_nematodes                       \
  -w 8

```

&nbsp;
&nbsp;


### Example 3

We cluster the [*JASPAR 2022 plants*](https://jaspar.uio.no/matrix-clusters/plants/) motif collection (656 motifs), we compare the resulting clusters detected by *RSAT matrix-clustering* against a user-provided reference annotation (in this case the Transcription Factor classes). We calculated the *Adjusted Rand Index* (*ARI*), a single-value metric (ranging from -1 to +1) indicating the proportion of consistent pairs between two classifications, in this example the *ARI* measures the proportion of motif pairs that are consistently classified between *RSAT matrix-clustering* results and the reference TF classes. We consider that a motif pair is consistently classified when the two motifs either belong to the same class and are grouped in the same cluster, or belong to different families and grouped in different clusters.

The calculation of *ARI* is a new functionality of *RSAT matrix-clustering*, in this example, using default paramters we obtained `ARI = 0.38`, changing parameters may increase/decrease the resulting *ARI*.

:hourglass_flowing_sand: Running time: ~5 minutes

```bash
Rscript matrix-clustering.R                         \
  -i data/JASPAR_2022/Jaspar_plants_motifs_tab.txt  \
  -o results/Jaspar_plants/Jaspar_plants            \
  -w 8                                              \
  --ARI TRUE                                        \
  -a data/JASPAR_2022/Jaspar_2022_plants_TF_fam.tab 
```

&nbsp;
&nbsp;

## :scroll: Input files

### Motifs (Mandatory)

This version of *RSAT matrix-clustering* relies on the R/Bioconductor package `universalmotif` for the motif manipulation steps.
This is the list of supported TF motif formats:

- cluster-buster
- cisbp
- homer
- jaspar
- meme
- transfac
- uniprobe

In case your motif format is not in this list, please contact me to add it.

&nbsp;

### Matrix file table (Mandatory, Option 1)

To avoid long commands when the input are many motif collections, we opted for a simple file format. The input file (`-i`) must be a tab-delimited file providing the following information (in the following order; no header):

1. Motif collection file path
2. Collection name (i.e. an alias given to the motif file)
3. Motif format: see above for the list of supported formats.

Each line in this table should correspond to a different motif file.

If a file path is duplicated it will be considered only once.

If a collection name is duplicated, the program will stop, collection names are needed to create unique motif IDs.

Input motifs may be in different formats.

Example:

```bash
data/OCT4_datasets/HOMER_OCT4_motifs.homer    HOMER_motifs	  homer
data/OCT4_datasets/MEME_OCT4_motifs.meme      MEME_motifs	  meme
data/OCT4_datasets/RSAT_OCT4_motifs.tf        RSAT_motifs	  tf
```

&nbsp;
&nbsp;

### Matrix file table (Mandatory, Option 2)

Alternatively, if you have one input file, instead of creating a matrix file table, you can pass the file path, collection name and motif format directly through the command-line (see the **Example 0**).

The three arguments (`--M`, `--C`, and `--F`) must be provided.

Example:

```bash
--M data/OCT4_datasets/RSAT_OCT4_motifs.tf  --C RSAT_motifs  --F tf
```

&nbsp;
&nbsp;

### Motif annotation table (Optional)

Users can provide a reference table that may be used for two purposes:

1. Compare the resulting clusters against a user-defined annotation (i.e., how good the resulting clusters resemble the reference classes in the annotation file). When this file is provided the *adjusted rand index* (*ARI*) will be calculated.
2. Add some color features (class numbers and background) in a radial tree (see this [example](https://jaspar.uio.no/matrix-clusters/plants/)). Note that this table is only used with the parameter `--radial_tree = TRUE`. 

The reference table (`-a`) must be a tab-delimited file providing at least the following columns (extra columns are ignored):

1. motif_id	
2. class
3. collection
4. url (this may be an empty column but the column name is expected)

The motif IDs and collectio names in this reference table must be the same IDs as in the motif files, if this is not the case the program will stop.

Example: 
```bash
motif_id  class       collection     url
MA1404.1	BBR/BPC     JASPAR_plants  https://jaspar.uio.no/matrix/MA1404.1
MA1403.1	BBR/BPC     JASPAR_plants  https://jaspar.uio.no/matrix/MA1403.1
MA1402.1	BBR/BPC     JASPAR_plants  https://jaspar.uio.no/matrix/MA1402.1
MA1197.1	CAMTA       JASPAR_plants  https://jaspar.uio.no/matrix/MA1197.1
MA0969.1	CAMTA       JASPAR_plants  https://jaspar.uio.no/matrix/MA0969.1
MA0970.1	CAMTA       JASPAR_plants  https://jaspar.uio.no/matrix/MA0970.1
MA0975.1	AP2/EREBP   JASPAR_plants  https://jaspar.uio.no/matrix/MA0975.1
MA0976.2	AP2/EREBP   JASPAR_plants  https://jaspar.uio.no/matrix/MA0976.2
MA1376.1	AP2/EREBP   JASPAR_plants  https://jaspar.uio.no/matrix/MA1376.1
```

&nbsp;
&nbsp;


## :crystal_ball: Example output:

This is the folder structure after running this software:

```bash
results
├── *_motifs
│   ├── individual_motifs_with_gaps
│   │   └── Two files per motif (direct and reverse orientation) in transfac format, these motifs are already aligned (may contain gaps).
│   │
│   ├── motifs_sep_by_cluster   (each folder contains the motifs belonging to a cluster)
│   │   └── Cluster_01
│   │   └── Cluster_02
│   │   └── ...
│   │   └── Cluster_N
│   │
│   └── root_motifs
│       └── Root_motifs.tf  (also referred as archetype motifs)
│
├── *_plots
│   ├── Clusters_vs_reference_contingency_table.pdf
│   └── Heatmap_clusters.pdf
│   
├── *_tables
│   ├── alignment_table.tab
│   ├── clusters.tab
│   ├── distance_table.tab
│   ├── pairwise_motif_comparison.tab
│   └── summary_table.tab
│
└── *_trees
    ├── tree.json
    ├── tree.newick
    └── tree.RData
```

&nbsp;

### Example 1

The main output of this analysis is an html file containing the logo forest (multiple hierarchical trees, each representing a cluster with aligned motifs).

<img src="data/images/Logo_forest.png" width="550px" align="center">

&nbsp;

The analysis produces the file named `alignment_table.tab` which contains one line per motif with its corresponding cluster name, orientation in the alignment, the number of upstream/downstream gaps, the aligned consensus, and the alignment width.

```bash
cluster	    id	                        name                consensus         rc_consensus                  strand  offset_up offset_down aligned_consensus    alignment_width

cluster_01  RSAT_positions_7nt_m1_n9    positions_7nt_m1    NNATTTGCATATGCAAATNN    NNATTTGCATATGCAAATNN    R       4         0   ----NNATTTGCATATGCAAATNN	24
cluster_01  MEME_MEME_ChIP_1_n1         MEME_ChIP_1         ATGYWAA                 TTWRCAT                 R       7         10  -------TTWRCAT----------	24
cluster_01  MEME_MEME_ChIP_15_n15       MEME_ChIP_15        TATGCAAAT               ATTTGCATA               R       6         9   ------ATTTGCATA---------	24
cluster_01  RSAT_local_words_7nt_m3_n7  local_words_7nt_m3  NNATATGCAAATNN          NNATTTGCATATNN          R       4         6   ----NNATTTGCATATNN------	24
cluster_01  RSAT_oligos_7nt_mkv5_m1_n1  oligos_7nt_mkv5_m1  NNATGCAAATNN            NNATTTGCATNN            R       4         8   ----NNATTTGCATNN--------	24
cluster_01  RSAT_local_words_7nt_m2_n6  local_words_7nt_m2  NHATTTGCATAACAAWNN      NNWTTGTTATGCAAATDN      D       4         2   ----NHATTTGCATAACAAWNN--	24
cluster_01  HOMER_homer_1_n1            homer_1             YWTTNWNATGCAAA          TTTGCATNWNAAWR          R       7         3   -------TTTGCATNWNAAWR---	24
cluster_01  RSAT_local_words_7nt_m4_n8  local_words_7nt_m4  NNATTGTTATGCATAACAATNN  NNATTGTTATGCATAACAATNN  D       0         2   NNATTGTTATGCATAACAATNN--	24
```


The file `clusters.tab` contains one line per cluster with the motifs IDs and names. 

```bash
cluster     id                                                              name
cluster_01  MEME_MEME_ChIP_3_n3,RSAT_oligos_7nt_mkv5_m4_n4,HOMER_homer_9_n9 MEME_ChIP_3,oligos_7nt_mkv5_m4,homer_9
cluster_02  HOMER_homer_6_n6,HOMER_homer_16_n16,HOMER_homer_18_n18          homer_6,homer_16,homer_18
cluster_03  MEME_MEME_ChIP_4_n4,MEME_MEME_ChIP_7_n7                         MEME_ChIP_4,MEME_ChIP_7
cluster_04  RSAT_oligos_7nt_mkv5_m3_n3                                      oligos_7nt_mkv5_m3
```

If the option `--export_heatmap TRUE` is indicated the file `Heatmap_clusters.pdf` will be generated. This is a heatmap of `N x N` where N is the number of motifs, each cell represents the motif similarity. The color annotation bar corresponds to the clusters. 

<img src="data/images/Heatmap_clusters.jpeg" width="750px" align="center">

&nbsp;


### Example 2

When the users activate the option `--radial_tree TRUE` all the motifs are forced to be aligned in a single cluster.
The output is an `html` document containing the code in [`D3`](https://d3js.org/) (a javascript library). Open this document in an internet browser to visualize the results.


You can zoom in/out using the mouse and change the motif orientation by clicking on the red buttons on the page top. Click on the `Hide/Show legend` button to ease the readability.

<img src="data/images/Radial_tree_JASPAR_nematodes.png" width="500px" align="center">

&nbsp;
&nbsp;

<img src="data/images/Radial_tree_JASPAR_nematodes_zoom.png" width="500px" align="center">

&nbsp;
&nbsp;

<img src="data/images/Radial_tree_JASPAR_nematodes_zoom_nolegend.png" width="500px" align="center">

&nbsp;
&nbsp;

:warning: It is possible that this `html` is not properly displayed by all browsers, we recommend to use [Firefox](https://www.mozilla.org/en-US/).

:warning: It is possible that this `html` have to be opened from a webserver and may need you have [Apache](https://httpd.apache.org/) ready to use. More information in the **Extra** section below.

&nbsp;
&nbsp;


### Example 3

When the users provide a reference annotation table (argument `-r` or `--reference_cluster_annotation`) the script will produce a contingency table comparing the resulting clusters and the reference groups, this table is visualized as a heatmap in the file `Clusters_vs_reference_contingency_table.pdf`.

<img src="data/images/Clusters_vs_reference_contingency_table.jpeg" width="700px" align="center">

&nbsp;
&nbsp;


## :hammer_and_wrench: Arguments and Options


### Mandatory arguments

&nbsp;

#### Option 1 : provide a matrix_file_table'. Use this when you have 2 or more input motif collections.
 - `-i` or `--matrix_file_table` : A text-delimited file where each line contain the following fields/columns. It does not expect a header, but it expects these columns in the indicated order.
 
    1. Motif file path
    2. Motif collection name
    3. Motif format.

&nbsp;

#### Option 2 : use these arguments when you have one input motif collections.

An alternative to the `--matrix_file_table` argument. If you have once inpt file you can pass the following three arguments instead of creating the motif file table.
The program expects the three `--M --C --F` arguments, if one of them is not provided it will stop.

  - `--M` : Path to motif file.
  - `--C` : Motif collection name.
  - `--F` : Motif format.

&nbsp;

**NOTE**: the `--matrix_file_table`  and  `--M --C --F` options are mutually exclusive. If both options are provided the program will stop.

&nbsp;

- `-o` or `--output_folder` : Folder to save the results.

&nbsp;
&nbsp;

### Comparison + Clustering arguments

- `-m` or `--comparison_metric` : Comparison metric used to build the hierarchical tree. Default: `Ncor`. [Options: `cor`, `Ncor`].
- `-l` or `--linkage_method` : Linkage/agglomeration method to build the hierarchical tree. Default: `average`. [Options: `average`, `complete`, `single`].
- `-c` or `--cor_th` : Pearson correlation lower threshold. Default: `0.75`. [Options: any value among `-1` and `+1`].
- `-n` or `--Ncor_th` : Normalized Pearson correlation lower threshold. Default: `0.55`. [Options: any value among `-1` and `+1`].
- `-W` or `--w_th` : Minimum number of aligned positions between a pair of motifs. Default: `5`.

&nbsp;

### Output files

- `--minimal_output` : Only returns the aligned motifs and the alignment, clusters and motif description tables. Comparison results, plots and trees are not exported. Default : `FALSE`. [Options: `TRUE`, `FALSE`].
- `--export_newick` : Export hierarchical tree in Newick format. Default : `FALSE`. [Options: `TRUE`, `FALSE`].
- `--export_heatmap` : Export heatmap with clusters in PDF. Default : `FALSE`. [Options: `TRUE`, `FALSE`].

&nbsp;

### Annotation table + ARI

- `-a` or `--annotation_table`: motif annotation tab. One line per motif, the proved color will be used as background in the radial tree, the class name text will be shown as an annotation layer (ring) in the radial tree. A tab-delimited file with the following columns (additional columns are ignored):

  1. motif_id	
  2. class
  3. collection
  4. url (may be an empty column but the header is expected)

- `--ARI` : when this option is `TRUE` and an annotation table is provided `--annotation_table` the program will calculate the ARI (partition similarity) among the resulting clusters and the provided annotation classes.

&nbsp;

### Radial trees

- `--radial_tree`: When this option is activated all the motifs are forced to be aligned in a single cluster. Default : `FALSE`. [Options: `TRUE`, `FALSE`].

&nbsp;

### Other

- `-w` or `--number_of_workers` : Number of cores to run in parallel. Default: `2`. [Options: depends in your machine].
- `--heatmap_color_palette` : Cell colors in clusters heatmap. Default: `RdGy`. [Options: any colorBrewer palette, see colorbrewer2.org for details ].
- `--color_palette_classes` : Number of classes to create color palette in clusters heatmap. Default: `11`. [Options: depends on the selected colorBrewer palette, see colorbrewer2.org for details ].
- `--title`: Set an analysis title. Default : `matrix-clustering`.

&nbsp;

### Optional output

- `--export_newick` : Export hierarchical tree in Newick format. Default : `FALSE`. [Options: `TRUE`, `FALSE`]
- `--export_heatmap` : Export heatmap with clusters in PDF. Default : `FALSE`. [Options: `TRUE`, `FALSE`]
- `-M` or `--minimal_output` : When TRUE only returns the alignment, clusters and motif description tables. Comparison results, plots and trees are not exported. Default : `FALSE`. [Options: `TRUE`, `FALSE`]


 &nbsp;
 &nbsp;
  
## :collision: Contact + Contributors + Report issues 

Contributors
- [Jaime A Castro-Mondragon](https://jaimicore.github.io/) 
- [Rafel Riudavets-Puig](https://github.com/rriupu)
- [Walter Santana-Garcia](https://github.com/santanaw)
- [Ieva Rauluseviciute](https://github.com/ievarau)

 &nbsp;
 &nbsp;

This repository is maintained mainly by [Jaime A Castro-Mondragon](https://jaimicore.github.io/) with important contributions from [Ieva Rauluseviciute](https://github.com/ievarau).

:e-mail: j.a.castro.mondragon@gmail.com

Bluesky: [@jaimicore](https://bsky.app/profile/jaimicore.bsky.social)

Use this space to report [issues](https://github.com/jaimicore/matrix-clustering_stand-alone/issues) related to this repository.


&nbsp;
&nbsp;

## :pushpin: To Do/Wishlist

- When calculating the ARI, implement an option to find an optimal threshold thorugh a grid search approach.
- Generate the interactive `html` output motif trees.
- Implement the option to annotate clusters.
- Trim root motifs
- Detect the central motif within each cluster.
- Export motif collection intersection stats.
- Add option to trim motifs during clustering (at the moment this is done with the *convert-matrix.R* script).

&nbsp;
&nbsp;

## :bangbang: Extra

This repository also contains the script `convert-matrix` which is a simplified version of the `RSAT convert-matrix` tool, for motif manipulation and format conversion.

For the moment this scripts has three main functions:

1. Motif format conversion, see above for the supported formats.
2. Export reverse-complement of the input motifs
3. Trim motifs (remove columns with low information content).
4. Export motif logos.


Simple motif conversion from `transfac` to `meme` format.
```bash
Rscript convert-matrix.R                    \
  -i data/OCT4_datasets/RSAT_OCT4_motifs.tf \
  --from tf --to jaspar                     \
  --output_file results/convert-matrix_examples/RSAT_OCT4_motifs.jaspar
```


Simple motif conversion from `transfac` to `meme` format with reverse-complement motifs.
The file with the reverse-complement motifs has the suffix `_rc` in its name.
In this example:
  - Output: `results/convert-matrix_examples/RSAT_OCT4_motifs.jaspar`
  - Output (RC) : `results/convert-matrix_examples/RSAT_OCT4_motifs_rc.jaspar`
```
Rscript convert-matrix.R                    \
  -i data/OCT4_datasets/RSAT_OCT4_motifs.tf \
  --from tf --to jaspar                     \
  --rc TRUE                                 \
  --output_file results/convert-matrix_examples/RSAT_OCT4_motifs.jaspar
```


Simple motif conversion from `transfac` to `meme` format after motif trimming.
```
Rscript convert-matrix.R                    \
  -i data/OCT4_datasets/RSAT_OCT4_motifs.tf \
  --from tf --to jaspar                     \
  --trim TRUE                               \
  --IC_threshold 0.25                       \
  --spike_IC_threshold 0.25                 \
  --trim_values_output results/convert-matrix_examples/RSAT_OCT4_motifs_trim_values.txt \
  --output_file results/convert-matrix_examples/RSAT_OCT4_motifs_trimmed.jaspar
```

In this figure we show the advantages of using a window-based approach to trim the motifs instead of using a single value, we use as example the [IRF7 motif from JASPAR](https://jaspar.genereg.net/matrix/MA0772.1/).

&nbsp;

<img src="data/images/Trim_motifs_example.png" width="500px" align="center">


&nbsp;
&nbsp;


### Install `apache`

Unfortunately, to visualize the content of the `html` file containing the radial tree it is required to have installed [apache2](https://httpd.apache.org/)
and open this html file as a localhost. If you don't do this, you will not see the html content.

To install `apache` you can follow this [instructions](https://ubuntu.com/tutorials/install-and-configure-apache#1-overview).

Once `apache` is installed in your computer:

1. Remove the folder `sudo rm -rf /var/www/html` 
2. Copy the result folder (including all directories) to `/var/www/` . You will need **sudo** permissions. For this example we use the results of the *Example 2*.

```
sudo rm -rf /var/www/html
sudo mkdir -p /var/www/html
sudo cp -fr results/JASPAR_nematodes_radial /var/www/html/JASPAR_nematodes_radial
```
   
3. Open your browser and type `localhost`. Now you can browse the files in  `/var/www/`
4. Search and open the file `*.html`


## :tada: Acknowledgements

We thank the [*JASPAR*](https://jaspar.genereg.net/) curation team for their input to improve *RSAT matrix-clustering*; the [*RSAT developer team*](https://rsat.eead.csic.es/plants/people.php) for their constant support across many years of collaboration; and the users for their advices, suggestions, merge requests and reporting bugs :beetle:.

Special thanks to my colleagues (now PhD) Ieva Rauluseviciute (and her *gently reminders* :unamused: that pushed me to write this stand-alone version), Vipin Kumar and Katalin Ferenc (from [Anthony Mathelier's lab](https://mathelierlab.com/)) for testing this software, the discussions, ideas, and their suggestions of `R` libraries that make this script faster than the original version.

&nbsp;
&nbsp;

## :page_with_curl: How to cite this software?

If you use this software, please cite the [*RSAT matrix-clustering* publication](https://doi.org/10.1093/nar/gkx314) and/or the [latest *RSAT* publication](https://doi.org/10.1093/nar/gkac312).

```
# RSAT matrix-clustering

@article{10.1093/nar/gkx314,
    author = {Castro-Mondragon, Jaime Abraham and Jaeger, Sébastien and Thieffry, Denis and Thomas-Chollier, Morgane and van Helden, Jacques},
    title = {RSAT matrix-clustering: dynamic exploration and redundancy reduction of transcription factor binding motif collections},
    journal = {Nucleic Acids Research},
    volume = {45},
    number = {13},
    pages = {e119-e119},
    year = {2017},
    month = {06},
    issn = {0305-1048},
    doi = {10.1093/nar/gkx314},
    url = {https://doi.org/10.1093/nar/gkx314},
    eprint = {https://academic.oup.com/nar/article-pdf/45/13/e119/25367647/gkx314.pdf},
}

# RSAT 2022

@article{10.1093/nar/gkac312,
    author = {Santana-Garcia, Walter and Castro-Mondragon, Jaime A and Padilla-Gálvez, Mónica and Nguyen, Nga Thi Thuy and Elizondo-Salas, Ana and Ksouri, Najla and Gerbes, François and Thieffry, Denis and Vincens, Pierre and Contreras-Moreira, Bruno and van Helden, Jacques and Thomas-Chollier, Morgane and Medina-Rivera, Alejandra},
    title = {RSAT 2022: regulatory sequence analysis tools},
    journal = {Nucleic Acids Research},
    volume = {50},
    number = {W1},
    pages = {W670-W676},
    year = {2022},
    month = {05},
    issn = {0305-1048},
    doi = {10.1093/nar/gkac312},
    url = {https://doi.org/10.1093/nar/gkac312},
    eprint = {https://academic.oup.com/nar/article-pdf/50/W1/W670/44378306/gkac312.pdf},
}
```
