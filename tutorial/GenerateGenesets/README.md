# GenerateGenesets function



By default, `GenerateGenesets` returns a `geneset` with the `250` most upregulated and downregulated genes in each drug signature. You can change this behaviour by providing new values to `n.genes` and `mode`. Moreover, a small collection of functional pathways will be included in your `geneset` object. These pathways are related to the regulation of the epithelial-mesenchymal transition (EMT), cell cycle, proliferation, senescence and apoptosis. Note that `n.genes` and `mode` arguments do not affect to functional pathways.

```r
# Generate geneset object with one of the ready to use signature collections.
gset <- GenerateGenesets(PSc)
# Retrieve only the top 100 most upregulated genes in drug signatures (functional pathways remain unchanged)
up100 <- GenerateGenesets(PSc, n.genes = 100, mode = "up")
# You can deactivate the functional pathways option if you are not interested in evaluating them
nopath <- GenerateGenesets(PSc, include.pathways = FALSE)
```

Additionaly, you can computed a `geneset` from a pre-loaded PSc subset called DSS.

```r
# Generate geneset object with one of the ready to use signature collections
dss <- GenerateGenesets(DSS, include.pathways = FALSE)
```

Also, you can filter PSc, SSc and DDS objects by several fields (cap insensitive):

 * `drugs`: Drug name of interest (i.e sirolimus).
 * `IDs`: `sig_id` of the signature(s) of interest.
 * `MoA`: Desired mechanism of action of interest (i.e. MTOR INHIBITOR).
 * `targets`: Target gene of interest (i.e. MTOR).
 * `source`: `"LINCS"` (for PSc) or `"GDSC"`, `"CCLE"` and/or `"CTRP"` (for SSc)

```r
# Return a `geneset` with all sirolimus signatures, as well as signatures of sirolimus synonyms such as 
# rapamycin or BRD-K84937637
sirolimus <- GenerateGenesets(SSc, include.pathways = FALSE, filters = list(drugs = "sirolimus"))
# Return just a subset of sirolimus signatures
my_sigs <- GenerateGenesets(SSc, include.pathways = FALSE, filters = list(IDs = c("sig_2349", "sig_7409"))
# Return all MTOR INHIBITORS
MTORi <- GenerateGenesets(SSc, include.pathways = FALSE, filters = list(MoA = "MTOR INHIBITOR")
# Return all drugs targetting MTOR
mtor_targets <- GenerateGenesets(SSc, include.pathways = FALSE, filters = list(targets = "MTOR")
# Return only signatures derived from GDSC and CCLE
my_sources <- GenerateGenesets(SSc, include.pathways = FALSE, filters = list(source = c("GDSC", "CCLE"))
```

By calling `ListFilters` function, you can retrieve all the available values for a given field. The signatures that pass **ANY** of these filters are included in the final `geneset`.

```r
# Values for targets
ListFilters(entry = "targets")
# Geneset with all drugs taht target MTOR and sirolimus signatures
filter_combination <- GenerateGenesets(SSc, include.pathways = FALSE, 
                                       filters = list(drugs = "sirolimus", targets = "MTOR"))
```
You can check information about the pre-loaded signatures calling the object `drugInfo`. Also, each `geneset` object obtained using pre-loaded matrices contains a subset of `drugInfo` for the selected drugs.

```r
# drugInfo of the signatures of interest
gset@info
```

Finally, Beyondcell allows the user to input a GMT file containing the functional pathways/signatures of interest or a numeric matrix (containing a ranking criteria such as the t-statistic or logFoldChange).

 * **In case your input is a GMT file:** You must supply the path to the file. Take into account that the names of each gene set must end in `"_UP"` or `"_DOWN"` to specify its mode. In this case, `n.genes` and `mode` are deprecated.
 * **In case your input is a numeric matrix:** Make sure that rows correspond to genes and columns to signatures.
 
In both cases, `filters` argument is deprecated but you must indicate if the `comparison` that yielded your input was `"treated_vs_control"` or `"sensitive_vs_resistant"`.

```r
# Mock numeric matrix
m <- matrix(rnorm(500 * 25), ncol = 25, dimnames = list(rownames(PSc[[1]])[1:500], colnames(PSc[[1]])[1:25]))
num_matrix <- GenerateGenesets(m, n.genes = 100, mode = c("up", "down"), 
                               comparison = "treated_vs_control", include.pathways = TRUE)
```
