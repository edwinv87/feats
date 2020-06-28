import numpy as np

from sklearn.covariance import MinCovDet
from scipy.stats import chi2
from sklearn.decomposition import PCA


"""
============================
Method Name: DetectOutliers
============================

Method Description: This method implements the outlier detection algorithm. 


Arguments:
========== 

sc                  -   A single cell object which contains the data and metadata of genes and cells
n_clusters          -   The number of group in the clustered data

Returns:
========

sc                  -   The single cell object containing the outlier scores in the CellData assay. A column is added to the cell data 
                        assay containing the cluster labels. The column name is the method name 'FeatClust'. 

"""

def DetectOutliers(sc, cluster_label, red_dim = 2, outlier_prob_thres = 10**-4):

    # Store outlier probability in sc object 
    sc.addCellData(col_data = -np.log10(np.ones(sc.dim[1]) * outlier_prob_thres), col_name = 'Outlier_Thres')

    # First check if clustering has been performed
    if (sc.checkCellData(cluster_label) == False):
        raise ValueError("Clustering has not been done. Perform clustering first! ")

    else:
        print("Computing outliers . . .")
        # Get cluster labels
        labels = sc.getCellData(cluster_label)
        n_clusters = np.unique(labels)
        X = sc.getCounts()
        _, n_samples = X.shape

        # Sort according to F scores
        scores = sc.getGeneData('FEATS_F_Score')
        idx = np.argsort(scores, kind='mergesort')
        idx = idx[::-1] # Sort descending
        # X = X[idx[0:100], :]

        # PCA
        pc = PCA(n_components=red_dim)
        X_red = pc.fit_transform(X.T)
        X_red = X_red.T

        mcd = MinCovDet(assume_centered=True)
        #mcd = []
        #for i in range(n_clusters):
        #    mcd.append(MinCovDet(assume_centered=True))   # mcd object, to compute min cov determinant 

        oos = np.zeros(n_samples, dtype=bool)   # Out of sample estimates (bool), True if sample is not included 
                                                # in MCD computation

        squared_md = np.zeros(n_samples)        # Squared Mahalanobis Distance

        # For each cluster reduce the data and estimate the robust covariance
        for i in n_clusters:
            mask = (labels == i)

            # If number of samples is less than number of features in reduced data squared. 
            if (np.sum(mask) < red_dim**2):

                print("Number of samples is less than number of features squared.")
                print("Not performing outlier detection on cluster ", i)
                oos[mask] = False               # Set the samples as not an outlier
                squared_md[mask] = 0.0      # Set the mahalanobis distance as zero.

            else:
 
                cluster = X_red[:, mask]
                mcd.fit(cluster.T)          # Fit a minimum covariance determinant estimator
                # cluster_mu = mcd.location_
                # cluster_cov = mcd.covariance_
                squared_md[mask] = mcd.mahalanobis(cluster.T)
                oos[mask] = (mcd.support_ == False)

        outlier_score = -np.log10(chi2.sf(squared_md, red_dim))
        outliers = outlier_score > -np.log10(outlier_prob_thres)

        print ("Number of outliers = ", np.sum(outliers))
        print ("Number of points in out of sample = ", np.sum(oos))

        print("Saving outlier information in Single Cell object . . .")

        sc.addCellData(col_data = outliers, col_name = "FEATS_Outliers")
        sc.addCellData(col_data = squared_md, col_name = "FEATS_Msd")
        sc.addCellData(col_data = outlier_score, col_name = "Outlier_Score")
        sc.addCellData(col_data = oos, col_name = "FEATS_Oos")

    return sc