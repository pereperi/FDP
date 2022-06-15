# Final Degree Project

This repository includes the scripts for performing the analyses and preprocessing included in the Final Degree Project: 

Single-cell analysis of Plexiform Neurofibromas: Heterogeneity and tumor cell-of-origin


### ANALYSIS:

scRNA_PNF_corrected.R ------> Contains the full analysis of single-cell RNA-seq data from 3 pNFs.

scATAC_RNA.R ------> Contains the full analysis of the Multiome single-cell ATAC + Gene Expression data from the same 3 pNFs.

scRNA_merge.R ------> Contains the full analysis of single-cell RNA-seq data from iPSC-based in vitro 3D neurofibroma models, combined with the 3 pNFs.


### PREPROCESSING:

run_in_cluster_sc(...).sh ------> The scripts used to run the count pipeline from CellRanger for all of the samples.

run_in_cluster_(...)-arc.sh ------> The scripts used to run the count pipeline from CellRanger-ARC for all of the samples.

run_in_cluster_aggr.sh ------> The script used to run the aggr pipeline from CellRanger-ARC for Multiome data.

kb_processing.sh ------> The script used to generate spliced and unspliced count matrices using kallisto | bustools, to perform RNA velocity.

