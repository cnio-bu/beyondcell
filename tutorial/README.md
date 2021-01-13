<img src="../.img/beyondcell.png" width="500">

## Installing beyondcell
We recommend installing **beyondcell** via gitlab using devtools:

```r
library("devtools")
devtools::install("bu_cnio/beyondcell")
```

See the DESCRIPTION file for a complete list of R dependencies. If the R dependencies are already installed, installation should finish promptly.

## Using beyondcell
For a correct analysis with **beyondcell**, users should follow these steps: 

 * Compute BCS
 * Compute Therapeutic clusters
 * Check clustering and look for unwanted sources of variation
 * Regress out unwanted sources of variation
 * Recompute UMAP
 * Obtain signature's statistics
 * Visualize the results


## Compute BCS
In order to correctly compute the scores, the transcriptomic data needs to be pre-processed. This means that proper cell-based quality control filters, as well as normalization, scaling and clustering of the data, should be applied prior to the analysis with **beyondcell**. The `bcCompute` function allows you to input either a pre-processed seurat object or a single cell matrix. Have in mind, that when a seurat object is used as an input, the `DefaultAssay` must be specified, both `SCT` and `RNA` assays are accepted.

```r
# Generate a gene set object
gs <- GenerateGenesets(SSc)
# Read single cell experiment
sc = readRDS(path_to_sc)
# Set Assay
DefaultAssay(sc) <- "RNA"
# Compute score for the SSc
bc <- bcScore(sc, gs, expr.thres = 0.1)
```

> TIP: we recommend to input cells with at least 1000-1500 genes detected.

## Compute Therapeutic clusters
The ouput of the `bcScore` computation is a `bc object`. The object contains the normalized and scaled **beyondcell** scores and switch point, as well as information concerning the parameters used for the analysis. The bc object can be used as an input for a dimensionality reduction and clustering analysis. With this analysis, cells can be classified into distinct **therapeutic clusters**, that represent sets of cells sharing a common response to a particular drug exposition. The Uniform Manifold Approximation and Projection (UMAP) will allow the visualization of the identified clusters. 

```r
# Generate a gene set object
bc <- bcUMAP(bc, pc = 5, res = 0.2, k.neighbors = 20)
```

## Check clustering
It is important to check whether any unwanted source of variation is guiding the clustering analysis. The `bcClusters` function allows us to colour the UMAP based on the metadata variables that migth be influencing the clustering. We recommend checking these sources of variation among others:

 * Number of detected genes per cell
 * Number of detected counts
 * Cell cycle phase
 * Batch

```r
# Visualize whether cells are clustered based on the number of genes detecter per each cell
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
# Visualize whether cells are clustered based on their cell cycle status
bcClusters(bc, UMAP = "beyondcell", idents = "Phase", factor.col = TRUE)
```
> TIP: the cell cycle information must be present in bc@meta.data and can be obtained using Seurat's function `CellCycleScoring`

## Regress out unwanted sources of variation
The `bcRegressOut` function will allow us to correct existing sources of variation. Have in mind that the number of detected genes per cell will *always* have an inpact in the final score.

```r
bc <- bcRegressOut(bc, vars.to.regress = c("nFeature_RNA"))
```
> TIP: is the regression step taking too long? Check the amount of NAs per cell of your bc@normalized matrix. You migth need to refine the filtering of your single cell experiment based on the amount of detected features.

## Recompute Therapeutic clusters
Once corrected, you will need to recompute the dimensionality reduction and clustering, in order to find the *true* **therapeutic clusters** present in your sample. 

```r
# Recompute UMAP
bc <- bcUMAP(bc, pc = 5, res = 0.2, add.DSS = FALSE, k.neighbors = 20) 
```

## Obtain signature's statistics
A summary table can be obtained using the `bcRanks`function. This table includes summary metrics such as: the switch point, mean, median, sd, variance, min, max, proportion of NaN and residuals' mean of each signature. This table aims to help you in the prioritization of drug candidates. 

```r
# Obtain condition-based statistics
bc <- bcRanks(bc, idents = "condition")
# Obtain therapeutic cluster-based statistics
bc <- bcRanks(bc, idents = "bc_clusters_res.0.2")
```
> TIP: We recommend prioritizing drugs taking into account both the switch point and residuals.

## Visualization of the results
**beyondcell** provides several visualization functions to help you better understand the results.

|**Function** | **Description** |
|--------------|--------------|
|**`bcClusters()`** |Returns a ggplot object with the UMAP reduction (either beyondcell's or Seurat's) colored by the specified metadata column.|
|**`bcSignatures()`** |Returns a list of patchwork objects containing ggplot2s with the desired UMAP reduction (either beyondcell's or Seurat's) colored by bcscores or gene expression values.|
|**`bcHistogram()`** |Drawns a histogram plot of bcscores for each signature of interest. The plot can be a single histogram (if idents = NULL) or a histogram for each level found in idents.|
|**`bcCellCycle()`** |Drawns, for each signature of interest, a violindot plot of the bcscores grouped by the cell cycle phase (G1, G2M or S). Note that this information must be present in `bc@meta.data` and can be obtained using Seurat's function `CellCycleScoring`.|
