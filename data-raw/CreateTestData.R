rm(list = ls())
library(Matrix)
set.seed(1)

# --- Functions ---
# Proportion of expressed genes per cell and signature.
proportion.expressed <- function(mtx, genelist) {
  sapply(genelist, FUN = function(x) {
    apply(mtx[x, ], 2, function(y) sum(y > 0)/length(x))
  })
}

# Creates a gmt from a list.
output.gmt <- function(l, filename) {
  invisible(lapply(1:length(l), function(i) {
    cat(c(names(l)[i], "na", l[[i]], "\n"), sep = "\t", file = filename,
        append = TRUE)
  }))
}

# --- Data ---
pbmc <- Seurat::Read10X("PBMC3K_subset/")

# --- Code ---
# Counts matrix.
counts <- as.matrix(pbmc)

# Genesets.
ssc <- beyondcell::GetCollection(SSc, n.genes = 100, include.pathways = TRUE)
dss <- beyondcell::GetCollection(DSS, n.genes = 100, include.pathways = FALSE)

# Genelists with up and down genes all together.
ssc.genelist <- lapply(ssc@genelist[c(1:100, 585:613)], FUN = function(x) unlist(x))
dss.genelist <- lapply(dss@genelist, FUN = function(x) unlist(x))
all.genelists <- c(ssc.genelist, dss.genelist)

# All unique genes.
all.genes <- unique(unlist(all.genelists))

# Subset and rename rows.
counts <- counts[1:length(all.genes), ]
rownames(counts) <- all.genes

# Number of appearances of each gene in each signature.
genematrix <- sapply(all.genelists, FUN = function(x) all.genes %in% x)
rownames(genematrix) <- all.genes

# In how many signatures each gene appears?
table(rowSums(genematrix))

# Keep genes that appear only in one signature.
one.gene.sig <- genematrix[rowSums(genematrix) == 1, ]

# Get the 2 signatures with most unique genes
unique.genes <- colSums(one.gene.sig)
table(unique.genes)
max.unique.genes <- as.numeric(tail(names(table(unique.genes)), n = 2))
sig.max.unique.genes <- unique.genes[unique.genes %in% max.unique.genes]
sig.max.unique.genes <- names(sort(sig.max.unique.genes, decreasing = TRUE))

# What is the proportion of expressed genes per cell and signature?
cellsigmatrix <- proportion.expressed(counts, genelist = all.genelists)
hist(cellsigmatrix)

# Get a distribution with the values !=0 in the original matrix
hist(counts)
table.counts <- table(counts)
values <- unlist(sapply(2:length(table.counts), FUN = function(i) {
  rep(as.numeric(names(table.counts)[i]), each = table.counts[i])
}))

# Update dynamically the counts matrix until all percentages > 0.08.
while (any(cellsigmatrix <= 0.1)) {
  zeroes <- which(counts == 0)
  counts[sample(zeroes, size = length(values), replace = FALSE)] <- values
  cellsigmatrix <- proportion.expressed(counts, genelist = all.genelists)
  hist(cellsigmatrix)
}

hist(counts)

# Make the expression values of all the genes in the signature with most unique 
# genes = 0, for all cells.
counts[all.genelists[[sig.max.unique.genes[1]]], ] <- 0
cellsigmatrix <- proportion.expressed(counts, genelist = all.genelists)
hist(cellsigmatrix)

# Make the expression values of 91% of the genes in the 2nd signature with most 
# unique genes = 0, for 50% of cells.
sig.2.unique <- all.genelists[[sig.max.unique.genes[2]]]
sig.2.unique <- sample(sig.2.unique, size = (length(sig.2.unique) * 0.9) + 1,
                       replace = FALSE)
cells.2 <- sample(1:ncol(counts), size = as.integer(ncol(counts) / 2), 
                  replace = FALSE)
counts[sig.2.unique, cells.2] <- 0
cellsigmatrix <- proportion.expressed(counts, genelist = all.genelists)
hist(cellsigmatrix)

# Check results.
names(which(colSums(cellsigmatrix == 0) > 0)) # Sigs with 0 expressed genes.
table(colSums(cellsigmatrix == 0)) # Number of cells with 0 expressed genes.

which(colSums(cellsigmatrix < 0.1) > 0)
table(colSums(cellsigmatrix < 0.1))

which(colSums(cellsigmatrix < 0.2) > 0)
table(colSums(cellsigmatrix < 0.2))

# Create gmts.
genesets100 <- head(unique(c(sig.max.unique.genes, names(ssc@genelist)[1:100])), 
                    n = 100)

gmt10 <- unlist(ssc@genelist[genesets100[2:11]], recursive = FALSE)
names(gmt10) <- gsub(pattern = "\\.", replacement = "_", names(gmt10))
names(gmt10)[6] <- gsub(pattern = "_down", replacement = "_dn", names(gmt10)[6])
names(gmt10)[c(7, 8, 9, 10)] <- toupper(names(gmt10)[c(7, 8, 9, 10)])
names(gmt10)

gmt10warning <- unlist(ssc@genelist[genesets100[1:10]], recursive = FALSE)
names(gmt10warning) <- gsub(pattern = "\\.", replacement = "_", 
                            names(gmt10warning))
gmt10duplicated <- c(gmt10[1:18], gmt10[1:2])
gmt10incorrect <- gmt10
names(gmt10incorrect)[19] <- gsub(pattern = "_up", replacement = "_bad", 
                                  names(gmt10incorrect)[19])

gmt20 <- unlist(ssc@genelist[genesets100[2:21]], recursive = FALSE)
names(gmt20) <- gsub(pattern = "\\.", replacement = "_", names(gmt20))
gmt100 <- unlist(ssc@genelist[genesets100[2:101]], recursive = FALSE)
names(gmt100) <- gsub(pattern = "\\.", replacement = "_", names(gmt100))

# Save.
Matrix::writeMM(obj = Matrix(counts), file = "../tests/testdata/matrix.mtx")
write.table(colnames(counts), file = "../tests/testdata/barcodes.tsv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(counts), file = "../tests/testdata/genes.tsv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

output.gmt(gmt10, filename = "../tests/testdata/correct10.gmt")
output.gmt(gmt10warning, filename = "../tests/testdata/score_warning10.gmt")
output.gmt(gmt10duplicated, filename = "../tests/testdata/duplicated10.gmt")
output.gmt(gmt10incorrect, filename = "../tests/testdata/incorrect_mode10.gmt")
output.gmt(gmt20, filename = "../tests/testdata/correct20.gmt")
output.gmt(gmt100, filename = "../tests/testdata/correct100.gmt")