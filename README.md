# <i>Salmonella enterica</i> pangenome
The repository with working scripts, figures, data, and supplementary materials for pangenome analysis of 1598 <i>Salmonella enterica</i> genomes.


## Contents

This repository contains scripts used for statistical analysis of pangenome, phylogeny, and phenotypical associations. Please consult the methods section in the paper for extra details:

Merkushova, A.V.; Shikov, A.E.; Nizhnikov, A.A.; Antonets, K.S. For Someone, You Are the Whole World: Host-Specificity of <i>Salmonella enterica</i>. Int. J. Mol. Sci. 2023, 24, 13670. https://doi.org/10.3390/ijms241813670

## Figures
Figures are available in the `pics/` directory. For the description, please consult the results section of the article.

## Supplemenatry
All supplementary material is located in the `supplementary/` directory.
`Supplementary_figures.docx` is a Microsoft Word document with a description of all supplementary figures.
`Supplementary_tables.xlsx` is an Excel table with all supplementary tables and description.

## Data
Analyzed data are included in the `data/` directory.

* `data/pangenome/` folder contains gene presence/absence table obtained using Panaroo.
* `data/trees/` folder contains phylogenetic inferences in Newick format based on core gene alignments (aligned with MAFFT).
* `data/pyseer/` folder includes Pyseer results attributed to particular hosts.

## Scripts
The `scripts/analysis` directory includes all code used for pangenome analysis:
* `IO_lib.py` is an ancillary script with functions for processing CSV files;
* `metadata_downloader.py` was applied to download metadata of studied assemblies;
* `split_fasta_to_chunks.py` - script for re-naming Panaroo-attributed gene codes to real protein accession numbers;
* `make_table_for_pysser.py` - script for generating a table with phenotypic traits used for Pyseer analysis;
* `create_matrix_for_enrichments.py` - script generates two distance matrices based on sequence alignments: one using the Jaccard metric and the other using the Shimkevich-Simpson metric
* `filter_virDB_hits.py` - script for filtering hits matching the VFDB database using identity and coverage;
* `topGO_enrichment.R` - script for GO enrichment analysis using topGO;

The `scripts/figures` directory includes all code used for building images:
* `Salmonella_trees_plots.R` - R script for generating phylogenetics tree figure;
* `Salmonella_statistics_plot.py` - script to build the picture that illustrates GC content, genome length and number of hypothetical proteins;
* `Salmonella_pangenome_plot.py` - script to plot the main characteristics of the reconstructed pangenome;
* `Salmonella_enrichment_plot.py` - was used to build results of testing GO terms using topGO.
