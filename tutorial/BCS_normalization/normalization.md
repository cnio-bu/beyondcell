---
title: "Beyondcell Score Normalization"
author: "Maria Jose Jimenez-Santos"
date: "`r Sys.Date()`"
output:
  word_document: default
  pdf_document: default
  html_document:
    df_print: paged
vignette: |
  %\VignetteIndexEntry{Vignette Title} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Beyondcell scores are normalised in order to **penalize** cells with a great 
number of **zeros** and/or cells with **outliers** (genes whose expression is 
much higher than the rest of genes within the same cell).

In this vignette we try to demonstrate the effect of this normalization using a 
mock single-cell experiment.

## Recapitulation on Beyondcell score calculation
As it was stated in Materials and Methods, to compute the Beyondcell score (BCS) 
for a mode $M = \{UP, DN\}$ we perform two steps:

1. Compute raw BCS
2. Normalise these scores

### Compute raw BCS
Suppose we have a single-cell expression matrix, which we call $X$, that has $n$ 
rows and $m$ columns. In this matrix, rows correspond to genes and columns to 
cells. Thus, we have a set of $n$ genes $I = \{i_1, i_2, ..., i_n\}$ and a set 
of $m$ cells $J = \{j_1, j_2, ..., j_m\}$. Also, we have a geneset 
$G_M = \{S_{M,1}, S_{M,2}, ..., S_{M,p}\}$, which is a set of set of genes 
(called signatures and denoted as $S_{M,p}$) per mode $M = \{UP, DN\}$.

Suppose we are interested in computing the raw BCS of a signature $S_{M}$. The
first step would be to calculate the common genes between $I$ and $S_{M}$.

$$ G_{M,S_M} = I \cap S_M = \{g_1, g_2, ..., g_q\}$$

Then, we subset $X$ to keep only the genes present in $G_{M,S_M}$. We call this 
subsetted matrix $Y$, which has $q$ rows and $m$ cells. The raw score for each
cell $j$ is the mean expression of the genes in $G_{M,S_M}$.

$$ raw_{M,S_M,j} = \frac{1}{|G_{M,S_M}|} \cdot \sum_{k=1}^q{y_{kj}} = \bar{y_j}$$

### Normalise raw BCS
Finally, we normalise the raw scores using the sum, mean and standard deviation 
of the total gene expression per cell.

$$ norm_{M,S_M,j} = raw_{M,S_M,j} \cdot f$$

The normalization factor $f$ can be decomposed as follows:

$$ f = \frac{\sum_{k=1}^q{y_{kj}} - \sqrt{\frac{\sum_{k=1}^q{(y_{kj} - \bar{y_j})^2}}{q-1}}}{\bar{y_j} + \sqrt{\frac{\sum_{k=1}^q{(y_{kj} - \bar{y_j})^2}}{q-1}}}$$

$$ f = \frac{sumexpr - sd}{mean + sd}$$

Thus, **the higher the standard deviation the lower $f$**. This results in a 
higher penalization of the raw score for cells with genes that are much more 
expressed than the rest of genes within the same cell. Moreover, **the lower **
**the mean and sumexpr, the lower $f$**, which further penalises cells with a 
great number of zeros.

## Example
```{r, message = FALSE}
library(scater)
library(ggplot2)
library(Seurat)
library(stringr)
library(proxyC)
library(dplyr)
# Random seed for reproducibility
set.seed(123)
```

We create a mock single-cell expression matrix $X$ with $n$ genes and $m$ cells.

```{r}
n <- 1000
m <- 50
X <- scater::mockSCE(ngenes = n, ncells = m)
head(assay(X)[, 1:5])
```

We draw the histogram of the counts to see its distribution. We can observe 
that there is a high proportion of zeros, as one would expect from a 
single-cell experiment.

```{r, fig.width = 7, fig.height = 4, message = FALSE}
df <- data.frame(counts = as.vector(assay(X)))
ggplot(df, aes(counts)) + geom_histogram()
```

We also print a summary of the raw counts.

```{r}
summary(df$counts)
```

We can observe than the max 5344 is and the mean and median are 200 and 38 
respectively.

Next, we are going to make up three cells in order to illustrate how the 
normalization process works. `Cell_0001` will have a high number of zeros, 
whereas `Cell_0002` will present few genes with an expression of 20-40 times the
maximum expression of our mock experiment. `Cell_0003` will combine these two 
conditions.

```{r}
# Cell_0001
assay(X)[sample(1:50, size = 30), 1] <- 0
# Cell_0002
assay(X)[sample(1:50, size = 5), 2] <- sample(100000:200000, size = 5, 
                                              replace = TRUE)
# Cell_0003
zero_idx <- sample(1:50, size = 30)
high_idx <- sample(which(!1:50 %in% zero_idx), size = 5)
assay(X)[zero_idx, 3] <- 0
assay(X)[high_idx, 3] <- sample(100000:200000, size = 5, replace = TRUE)
# Resulting X
head(assay(X)[, 1:5])
```

Then, we create a Seurat object

```{r}
X <- CreateSeuratObject(counts = assay(X))
```

Note that these counts are raw, and beyondcell requires normalized counts. So 
the next step is to normalize $X$.

```{r}
DefaultAssay(X) <- "RNA"
X <- NormalizeData(X, normalization.method = "LogNormalize", 
                   scale.factor = 10000)
# We extract the normalized counts and store them as a sparse matrix
X <- GetAssayData(X, slot = "data", assay = DefaultAssay(X))
head(X[, 1:5])
```

Now we are ready to compute the BCS. We are doing so for a signature of $q$ 
genes $S_M = \{g_1, g_2, ..., g_q\}$.

```{r}
q <- 50
SM <- paste0("Gene-", stringr::str_pad(1:q, 4, pad = "0"))
```

We subset $X$ to keep only the genes present in $S_M$. We call this subset $Y$.

```{r, fig.width = 7, fig.height = 4}
Y <- X[SM, ]
stats_df <- data.frame(Cells = colnames(Y), nzeros = colZeros(Y))
# Number of zeros per cell
ggplot(stats_df, aes(x = Cells, y = nzeros)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90)) + ylab("Number of zeros")
```

We can observe than `Cell_0001` and `Cell_0003` have much more zeros than the 
rest of the cells.

Then, we compute the raw score (mean of expression) of each cell.

```{r, fig.width = 7, fig.height = 4}
stats_df <- stats_df %>% mutate(raw = colMeans(Y))
# Raw score per cell
ggplot(stats_df, aes(x = Cells, y = raw)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Raw BCS")
```

In this case `Cell_0001` has the lowest raw BCS, followed by `Cell_0003`. In 
contrast, `Cell_0002` score is similar to the score of the rest of cells.

The next step is to compute the components of the normalization factor.

```{r, fig.width = 7, fig.height = 4}
stats_df <- stats_df %>% mutate(sumexpr = colSums(Y), sd = apply(Y, 2, sd))
# Sum of gene expression per cell
ggplot(stats_df, aes(x = Cells, y = sumexpr)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Sum of Gene Expression")
# Standard deviation per cell
ggplot(stats_df, aes(x = Cells, y = sd)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("SD")
```

On one hand, we can see that the sum of gene expression is proportional to the 
mean. On the other hand, in the second bar plot we observe that `Cell_0001` has 
the minimum and `Cell_0002` and `Cell_0003` have the maximum standard deviation.

Finally, we compute $f$ and compare the raw and the normalised BCS.

```{r, fig.width = 7, fig.height = 4}
stats_df <- stats_df %>% mutate(f = (sumexpr - sd)/(raw + sd))
stats_df <- stats_df %>% mutate(norm = raw * f)
# Raw score per cell
ggplot(stats_df, aes(x = Cells, y = raw)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Raw BCS")
# Normalised score per cell
ggplot(stats_df, aes(x = Cells, y = norm)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Normalised BCS")
```

When comparing these two last bar plots, we can observe than beyondcell's 
normalization step penalises cells with a great number of zeros (`Cell_0001`), 
cells with outliers (`Cell_0002`) or cells that satisfy these two conditions 
simultaneously (`Cell_0003`).

Numerically, we can compute the ratio of each score with respect to their mean 
and compare these two measures by calculating the log fold change.

```{r, fig.width = 7, fig.height = 4}
stats_df <- stats_df %>% mutate(log2FC = log2((norm/mean(norm))/(raw/mean(raw))))
# Log2Fold Change of normalized vs raw
ggplot(stats_df, aes(x = Cells, y = log2FC)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("log2FC of ratios") + ggtitle("Normalised vs Raw BCS")
```

As we can see, the BCS ratio is reduced drastically after normalization for 
cells 1 to 3. This means that these three cells are the ones most penalised by 
the normalization method.

```{r}
sessionInfo()
```
