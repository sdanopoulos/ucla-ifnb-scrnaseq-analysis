# IFN-beta Prenatal Lung Single-Cell RNA-seq Analysis

This repository contains the R analysis code associated with:

Excess type I interferon disrupts developmental crosstalk and promotes premature distalization in trisomy 21 lungs

## Analysis overview

The script performs:

- Loading and merging of eight 10x Genomics Cell Ranger datasets
- Quality-control filtering
- Removal of doublets using previously generated scDblFinder classifications
- SCTransform normalization
- Reciprocal PCA (RPCA) integration
- Clustering and UMAP visualization
- Cell-type annotation using Azimuth lung and fetal references
- Epithelial and mesenchymal lineage subclustering
- Cluster-marker identification
- IFN-beta versus control differential-expression analysis using MAST
- Cell-cycle scoring
- CellChat analysis of epithelial-mesenchymal communication
- Comparison of WNT, BMP, and FGF signaling between IFN and control conditions

## Samples

Four prenatal lung specimens were analyzed under paired control and IFN-beta
conditions:

- H29912_CTL and H29912_IFN
- H29921_CTL and H29921_IFN
- H29925_CTL and H29925_IFN
- H29934_CTL and H29934_IFN

## Repository contents

- `ucla.ifnb.explants.scrnaseq.2026.cc.R`: complete analysis script

## Required input

The analysis uses 10x Genomics Cell Ranger gene-expression matrices for eight
prenatal lung samples.

The file paths marked `path to data` and `path to output` in the script must be
replaced with local paths before running the analysis.

## Data availability

The sequencing data associated with this study are available through the NCBI
Gene Expression Omnibus under accession
[GSE325181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325181).

## Software requirements

The analysis was developed using R 4.4.3 and Seurat 5.3.0.
Required R packages include:

- Seurat
- CellChat
- tidyverse
- here
- patchwork
- dplyr
- ggplot2
- Azimuth
- MAST

Azimuth lung and fetal reference datasets are also required.

## Reproducibility

A random seed of 1234 is set at the beginning of the analysis. The script
writes the complete R session information to `sessionInfo.txt`.

## License

This code is distributed under the terms stated in the `LICENSE` file.

## Citation

This analysis code is permanently archived in Zenodo:

Danopoulos S, Chu C. *UCLA IFN-beta Prenatal Lung Single-Cell RNA-seq
Analysis*. Version 1.0.0. Zenodo. 2026.
DOI: 10.5281/zenodo.23108475

Additional citation metadata is provided in `CITATION.cff`.
