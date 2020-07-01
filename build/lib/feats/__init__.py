"""
FEATS
=====
FEATS is a new Python tool for performing the following downstream analysis on single-cell RNA-seq datasets:
1. Clustering
2. Estimating the number of clusters
3. Outlier detection
4. Batch correction and integration of data from multiple experiments

See https://github.com/edwinv87/feats for more information and documentation.
 
"""


from .feats_clustering import Cluster
from .feats_detectoutliers import DetectOutliers
from .feats_batchcorrection import IntegrateBatches, CorrectBatches, MergeBatches
from .feats_gapstatistics import GapStatistics
from .feats_filtering import LogFilter, HVGFilter, GeneFilter, CellFilter
from .feats_transformations import PCA
from .feats_utils import FeatureNormalize

__all__ = ['Cluster', 'DetectOutliers', 'IntegrateBatches', 'CorrectBatches', 'MergeBatches', 'GapStatistics', 'LogFilter', 'HVGFilter',
'GeneFilter', 'CellFilter', 'PCA', 'FeatureNormalize']