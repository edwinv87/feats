# FEATS

FEATS is a Python tool for performing the following downstream analysis on single-cell RNA-seq datasets:

1. Clustering
2. Estimating the number of clusters
3. Outlier detection
4. Batch correction and integration of data from multiple experiments

## Prerequisites

FEATS depends on the following packages

1. numpy
2. pandas
3. scikit-learn
4. scipy
5. [singlecelldata](https://edwinv87.github.io/singlecelldata/)

## Installation

The latest version of FEATS can be installed from PyPI:

`pip install feats`

## Documentation

The functional reference manual for FEATS is available [here](https://feats.readthedocs.io/en/latest/index.html).

### Examples

To use FEATS, please refer to the following example code presented in the Jupyter notebook.

1. [Clustering using FEATS](https://edwinv87.github.io/feats/docs/FEATS-Clustering.html)
2. [Performing outlier analysis](https://edwinv87.github.io/feats/docs/FEATS-Outlier-Detection.html)
3. [Performing batch correction](https://edwinv87.github.io/feats/docs/FEATS-Batch-Correction.html)

### Data

The data for the examples in this section is available [here](https://edwinvans.com/datasets/). The data is contained in subfolders in the datasets folder. The subfolders are named according to the dataset name. To load the data for the examples above, provide the path to the datasets folder on your local machine.

## Paper

The FEATS paper is published in the journal Briefings in Bioinformatics. It is available [here](https://academic.oup.com/bib/article-abstract/22/4/bbaa306/6025503)

## Contact

Contact the author on <info@edwinvans.com> to give feedback/suggestions for further improvements and to report issues.
