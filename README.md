<img src="./.img/beyondcell.png" width="500">

[Package status](https://gitlab.com/bu_cnio/Beyondcell/commits/master)

## Introduction
**Beyondcell** is a methodology for the identification of drug vulnerabilities 
in single-cell RNA-seq (scRNA-seq) data. To this end, Beyondcell focuses on the 
analysis of drug-related commonalities between cells by classifying them into 
distinct Therapeutic Clusters (TCs).

## Workflow overview

**Beyondcell workflow.** Given two inputs, the scRNA-seq expression matrix and a 
collection of drug signatures, the methodology calculates a Beyondcell Score 
(BCS) for each drug-cell pair. The BCS ranges from 0 to 1 and measures the 
susceptibility of each cell to a given drug. The resulting BCS matrix can be 
used to determine the sample’s TCs. Furthermore, drugs are prioritized in a 
table and each individual drug score can be visualized in a UMAP.

![Beyondcell workflow](./.img/workflow_tutorial.png)

Depending on the evaluated signatures, the BCS represents the cell perturbation 
susceptibility (PSc) or the sensitivity to the drug effect (SSc). BCS can also 
be estimated from functional signatures  to evaluate each cell functional 
status.

![drug signatures](./.img/drug_signatures.png)

## Beyondcell's key applications
 * Analyse the intratumoural heterogeneity (ITH) of your experiment 
 * Classify your cells into TCs
 * Prioritize cancer treatments
 * If time points are available, identify the changes in drug tolerance
 * Identify mechanisms of resistance

## Installing Beyondcell
The Beyondcell algorithm is implemented in R (v. 4.0.0 or greater). We recommend 
running the installation via conda: 

```r
# Create a conda environment.
conda create -n beyondcell 
# Activate the environment.
conda activate beyondcell
# Install Beyondcell package and dependencies.
conda install -c bu_cnio r-beyondcell
```

## Results
We have validated Beyondcell in a population of MCF7-AA cells exposed to 500nM 
of bortezomib and collected at different time points: t0 (before treatment), 
t12, t48 and t96 (72h treatment followed by drug wash and 24h of recovery) 
obtained from *Ben-David U, et al., Nature, 2018*. We integrated all four 
conditions using the Seurat pipeline (left). After calculating the BCS for each 
cell using PSc, a clustering analysis was applied. Beyondcell was able to 
cluster the cells based on their treatment time point, to separate untreated 
cells from treated cells (center) and to recapitulate the changes arisen by the 
treatment with bortezomib (right). 

![results_golub](./.img/integrated_bendavid.png)


## How to run
For general instructions on running Beyondcell, check out the [analysis workflow](https://gitlab.com/bu_cnio/Beyondcell/-/tree/master/tutorial/analysis_workflow) and [visualization](https://gitlab.com/bu_cnio/Beyondcell/-/tree/master/tutorial/visualization) tutorials.


## Authors

 * Coral Fustero-Torre
 * María José Jiménez-Santos
 * Santiago García-Martín
 * Carlos Carretero-Puche
 * Luis García-Jimeno
 * Tomás Di Domenico
 * Gonzalo Gómez-López
 * Fátima Al-Shahrour


## Citation
 **Beyondcell: targeting cancer therapeutic heterogeneity in single-cell RNA-seq**\
Coral Fustero-Torre, María José Jiménez-Santos, Santiago García-Martín, Carlos Carretero-Puche, Luis García-Jimeno, Tomás Di Domenico, Gonzalo Gómez-López, Fátima Al-Shahrour.\
_bioRxiv_ 2021.04.08.438954; doi: https://doi.org/10.1101/2021.04.08.438954 

## Support
If you have any question regarding the use of Beyondcell, feel free to submit an [issue](https://gitlab.com/bu_cnio/Beyondcell/issues).
