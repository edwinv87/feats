library(SingleCellExperiment)
library(mclust)
library(pcaReduce)


datasets <- c("biase", "yan1", "yan2", "goolam", "deng1", "deng2", "pollen", "kolodziejczyk", "patel", "fan", "treutlein")
# datasets <- c("yan1", "patel")
# datasets <- c("deng1")
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
    sce <- sce[, 1:49]
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 3
  }

  else if (dset == "yan1")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/yan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
  }

  else if (dset == "yan2")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/yan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type2
    nc <- 7
  }

  else if (dset == "goolam")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/goolam.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- counts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
  }

  else if (dset == "fan")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/fan.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
  }

  else if (dset == "treutlein")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/treutlein.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
  }

  else if (dset == "deng2")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/deng-rpkms.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type2
    nc <- 10
  }
  
  else if (dset == "deng1")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/deng-rpkms.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
  }

  else if (dset == "pollen")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/pollen.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 11
  }

  else if (dset == "kolodziejczyk")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/kolodziejczyk.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- counts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 3
  }

  else if (dset == "patel")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/patel.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- logcounts(sce)
    # counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 5
  }


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
  # for (i in 1:dim(counts)[2])
  # {
  #   nrm <- norm(counts[i,], "2")
  #   counts[i,] <- counts[i,]/nrm
  # }

  # Compute Clusters and Save ARI 
  Output_S <- PCAreduce(t(counts), nbt=100, q=30, method='S')

  for (i in 1:n_trials)
  {
    trial <- Output_S[[i]]
    ari[i, dset] <- adjustedRandIndex(cells, trial[, 30-(nc-2)])
  }

}


# Write the ARI values to file 
write.csv(ari, file = "pcaReduce_results_all.csv")