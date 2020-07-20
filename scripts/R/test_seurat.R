library(SingleCellExperiment)
library(mclust)
library(Seurat)

# datasets <- c("biase", "yan1", "yan2", "goolam", "deng1", "deng2", "pollen", "kolodziejczyk", "patel", "fan", "treutlein")
#datasets <- c("yan1", "patel")
datasets <- c("treutlein")
n_trials <- 100
ari <- matrix(nrow = n_trials, ncol = length(datasets))
colnames(ari) <- datasets

for (dset in datasets)
{


  # DIFFERENT SETTINGS FOR DIFFERENT DATASETS
  ### biase
  if (dset == "biase")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/biase.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    counts <- counts[, 1:49]
    cells <- colData(sce)$cell_type1
    cells <-cells[1:49]
    nc <- 3
    res <- 1.2
  }

  else if (dset == "yan1")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/yan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
    res <- 2.8
  }

  else if (dset == "yan2")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/yan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type2
    nc <- 7
    res <- 2.9
  }

  else if (dset == "goolam")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/goolam.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- counts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
    res <- 1.0
  }

  else if (dset == "fan")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/fan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
    res <- 2.1
  }

  else if (dset == "treutlein")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/treutlein.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
    res <- 1.3
  }

  else if (dset == "patel")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/patel.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- logcounts(sce)
    # counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
    res <- 0.7
  }
  else if (dset == "deng1")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/deng-rpkms.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
    res <- 0.67
  }
  else if (dset == "deng2")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/deng-rpkms.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type2
    nc <- 10
    res <- 2.6
  }
  else if (dset == "pollen")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/pollen.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 11
    res <- 3.6
  }
  else if (dset == "kolodziejczyk")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/kolodziejczyk.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- counts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 3
    res <- 0.01
  }

  genes <- rownames(sce)

  # GENE FILTER
  if (dset == "patel")
  {

  }

  else
  {
    pct_dropout_max <- 90
    pct_dropout_min <- 10

    dropouts <- (rowSums(counts == 0) / dim(counts)[2]) * 100
    gene_filter <- (dropouts < pct_dropout_max) & (dropouts > pct_dropout_min)
    counts <- counts[gene_filter, ]
  }


  # NORMALIZE
  #for (i in 1:dim(counts)[2])
  #{
  #  nrm <- norm(counts[i,], "2")
  #  counts[i,] <- counts[i,]/nrm
  #}

  so <- CreateSeuratObject(counts = counts)
  so <- ScaleData(so)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
  so <- RunPCA(object = so, npcs = 30)
  so <- FindNeighbors(so, dims = 1:10)
  so <- FindClusters(so, resolution = res)
  
  # Finally Compute Clusters and Save the ARI
  for (i in 1:n_trials)
  {
    so <- FindClusters(so, resolution = res)
    ari[i, dset] <- adjustedRandIndex(cells, Idents(so))
  }
  
}

write.csv(ari, file = "seurat_results2.csv")

# Finally Compute the ARI
# print(adjustedRandIndex(cells, so@ident))