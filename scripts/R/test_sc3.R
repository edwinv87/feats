# R Script for testing SC3 clustering approach on various datasets
# Name: Edwin Vans
# Date: 19/09/2018

# Retested on 21/05/2020 with latest version of SC3

library(SingleCellExperiment)
library(SC3)
library(scater)
library(mclust)

datasets <- c("biase", "yan1", "yan2", "goolam", "deng1", "deng2", "pollen", "kolodziejczyk", "patel", "fan", "treutlein")
# datasets <- c("yan1", "patel")
# datasets <- c("yan1")
# datasets <- c("deng2", "pollen", "kolodziejczyk")
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
    sce <- sce[,1:49]
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

  else if (dset == "deng1")
  {
    path <- paste("D:/OneDrive/PhD Research/R-Datasets/deng-rpkms.rds", sep="", collapse = NULL)
    sce <- readRDS(path)
    counts <- normcounts(sce)
    counts <- log2(counts + 1)
    cells <- colData(sce)$cell_type1
    nc <- 6
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
    sce <- sce[gene_filter, ]
    logcounts(sce) <- as.matrix(counts)
  }



  # Finally Compute Clusters and Save the ARI
  for (i in 1:n_trials)
  {
    sce <- sc3(sce, ks = nc, biology = FALSE, gene_filter = FALSE)
    ari[i, dset] <- adjustedRandIndex(cells, colData(sce)[,paste("sc3_",nc,"_clusters", sep="", collapse = NULL)])
  }

} 

# Save ARI results to file
write.csv(ari, file = "sc3_results_all.csv")

