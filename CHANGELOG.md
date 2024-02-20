# [2.2.1] - 2024-02-20
* **Fixed** dependency errors. Beyondcell requires `Seurat` >=4.0.0,<5.0.0 and `Matrix` ==1.6_1.1. Related to issue [#151](https://github.com/cnio-bu/beyondcell/issues/151).
* **Added** Spatial Transcriptomics tutorial example.

# [2.2.0] - 2023-07-06
* **Fixed** an issue ([#132](https://github.com/cnio-bu/beyondcell/issues/132)) where `bcScore`'s check for normalized data failed.
* **Fixed** an issue ([#135](https://github.com/cnio-bu/beyondcell/issues/135)) where `bcRanks` was not computing the residual's mean.
* **Updated** tutorials.
* **Removed** `bcRanks` "general" statistics computation.
* **Removed** `GetStatistics` function.

# [2.1.0] - 2023-03-31
* **Fixed** an error in `bcClusters` and `bcSignatures` when `spatial = TRUE` and `images = NULL`.
* **Fixed** title error in `bc4Squares`. The current title matches the corresponding `idents` level.
* **Fixed** error in `bc4Squares` when drug names were passed as `topnames`.
* **Updated** pre-loaded collections, which have been trimmed to 250 up and 250 down genes.
* **Updated** SSc collection, which only contains protein-coding genes.

# [2.0.0] - 2023-01-25
* **Added** functions to plot spatial transcriptomics data on top of selected images. 
* **Added** a function to load our collection of drug signatures: `GetCollection`.
* **Added** a set of functional pathways related to cancer disease: these are available using `include.pathways` in `GetCollection`.
* **Added** two new optional arguments to  `bc4Squares`: `x.cutoff` and `y.cutoff` to set custom quadrant thresholds. 
* **Added** two new optional arguments to `bcUMAP`: `npcs` and `seed` to allow users to both set the number of components to compute in the PCA step and reproduce a given UMAP by specifying the seed.
* **Changed** the KNN imputation methodology in `bcRegressOut`. Beyondcell now relies on **DMwR** instead of **bnstruct**.
* **Fixed** multiple issues in the package when only a single signature was being provided by the user. 
* **Updated** the SSc collection. A previous source of drug signatures (CCLE) has been replaced by those from the PRISM study. CTRP and GDSC derived signatures have been updated to the latest releases of the aforementioned projects' data. 
* **Updated** the drug annotation to harmonize the mechanism of action for the compounds in the _SSc_ collection.
* **Updated** all of our pre-loaded collections, which are stored as geneset objects. 
* **Updated** `GenerateGenesets`  to work on user-provided genesets stored in GMT format.
* **Updated** beyondcell's dependencies.
* **Updated** binned scales in `bcSignatures` to improve visualizations.
* **Removed** the `method` argument from `bcUMAP`. Only the UWOT library implementation of UMAP is available.

# [1.3.3] - 2022-01-12
* **Fixed** an error in `bcSubset` when filtering by the proportion of `NaN` values. Previously, the beyondcell object wasn't a subset when using `nan.sigs` or `nan.cells`. Now those signatures/cells below or equal to the `NaN` threshold are kept. Related to issue [#21](https://github.com/cnio-bu/beyondcell/issues/21).

# [1.3.2] - 2021-12-20
* **Fixed** an issue ([#17](https://github.com/cnio-bu/beyondcell/issues/17)) where if all the cells have a number of expressed genes < `expr.thres`, 
all the values for that signature will be `NaN` and so will be the switch point. 

# [1.3.1] - 2021-09-9
* **Added** group annotation to rank tables: high/low/differential sensitivity is now defined for each tested drug.

# [1.2.1] - 2021-06-29
* **Fixed** an issue ([#14](https://github.com/cnio-bu/beyondcell/issues/14)) where `bcSignatures` would fail when used with pathways

# [1.2.0] - 2021-06-8
* **Fixed** an issue where `bcMerge` would fail if the background matrix was empty. 

# [1.1.1] - 2021-05-20
* **Updated** pathways with the latest sources from MSigDB 7.3.
* **Added** a directory called `data-raw` with both the code and the data used to generate the pathway signatures. 

# [1.1.0] - 2021-04-14
* **Added** two new arguments in `bcUMAP`: `method` and `return.model`. `method` can be `"umap-learn"` or `"uwot"` (default). `return.model` is a logical argument that determines whether the function will return the uwot model (only valid when `method = "uwot"`).
* **Changed tutorials**: Typos have been corrected, chunks of code now specify the UMAP method (`"umap-learn"`) and figures have been redrawn
(they look the same as the originals).

# [1.0.0] - 2021-04-8
* Initial public release: subsequent Beyondcell releases will follow [Semantic Versioning guidelines](https://semver.org/).
