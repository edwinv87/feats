import numpy as np

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_samples
from sklearn.feature_selection import f_classif, SelectKBest, mutual_info_classif

from .feats_gapstatistics import GapStatistics
from .feats_utils import FeatureNormalize


# Compute Clustering Score
def ClusteringScore(X, labels):

    s_score = silhouette_samples(X, labels)
    s_average = np.mean(s_score)

    return s_average


def SelectTopNScores(clustering_score, n):

    idx = np.argsort(clustering_score, kind='mergesort')
    idx = np.flip(idx)
    #print("s_score array size: ", idx.size)

    if (idx.size < n):
        return idx
    else:
        return idx[0:n]



"""
============================
Method Name: AnovaHierarchical
============================

Method Description: This method implements the proposed FeatClust approach. 


Arguments:
========== 

sc_obj              -   A single cell object which contains the data and metadata of genes and cells
n_clusters          -   The number of group to cluster the data into. This is in the range 2 <= n_clusters < n. Default (10)
gene_filt           -   A boolean variable. Gene filtering is done if True. Default (True)
apply_normalization -   A boolean variable. Data is normalized if True. Default (True)
normalization       -   A string variable. The type of normalization to apply. Options are "l2", "mean" and "norm6". Default ("l2")
k                   -   The number of groups to divide the features of the dataset into. Valid range n_clusters <= k < d. DefaultDefault (10)


Returns:
========

sc_obj              -   The single cell object containing the cluster labels in the CellData assay. A column is added to the cell data 
                        assay containing the cluster labels. The column name is the method name 'FeatClust'. 

"""

# ANOVA Hierarchical Clustering
def AnovaHierarchical(
    sc,                                
    n_clusters, 
    classification_type,
    normalization,
    q
    ):

    method_name = "FEATS_"

    # Get expression matrix
    X = sc.getCounts()

    if (normalization == 'mean'):
        sc = FeatureNormalize(sc, 'mean')
        X_nrm = sc.getCounts()
        X_nrm = X_nrm.T
    elif (normalization == 'l2'):
        sc = FeatureNormalize(sc, 'l2')
        X_nrm = sc.getCounts()
        X_nrm = X_nrm.T
    elif (normalization == 'cosine'):
        sc = FeatureNormalize(sc, 'cosine')
        X_nrm = sc.getCounts()
        X_nrm = X_nrm.T
    else:
        X_nrm = X.T
    
    # FEATURE SELECTION STEP
    # Using Hierarchical Clustering, compute and store temporary clusters in sc object
    print("Computing Temporary Clusters . . .")
    model_hc = AgglomerativeClustering(n_clusters = n_clusters, linkage='ward')
    temp_labels = model_hc.fit_predict(X_nrm) + 1
    sc.addCellData(col_data = temp_labels, col_name = method_name + str(n_clusters) + "_Clusters_Temp")
    
    # Perform ANOVA analysis to find significant genes
    print("Performing Feature Selection . . .")
    if (classification_type == "f_classif"):
        feat_sel = SelectKBest(f_classif, k="all")
        
    elif (classification_type == "mutual_info_classif"):
        feat_sel = SelectKBest(mutual_info_classif, k="all")
    
    feat_sel.fit(X_nrm, temp_labels)
    feature_scores = feat_sel.scores_
    idx = np.argsort(feature_scores, kind='mergesort')
    idx = idx[::-1] # Sort descending

    # Save F-scores in the SC object
    sc.addGeneData(col_data = feature_scores, col_name = method_name + "F_Score")
    

    # CLUSTERING STEP
    clust_score = np.zeros(q.shape[0]) # for storing clustering score
    pred_labels = np.zeros([X_nrm.shape[0], q.shape[0]]) # for storing cluster labels

    # compute clusters and iterate
    print("Computing " + str(n_clusters) + " clusters using best features . . .")
    k = 0
    for j in q:
        X_red = X_nrm[:, idx[0:j]]
        pred_labels[:,k] = model_hc.fit_predict(X_red) + 1
        clust_score[k] = ClusteringScore(X.T, pred_labels[:,k])
        k = k + 1

    # FINAL CLUSTERING
    # Choose clustering with the max. clustering score
    mask = SelectTopNScores(clust_score, 1)
    final_labels = pred_labels[:, mask]
    final_labels = final_labels.astype(int)

    # SAVE FINAL CLUSTERING RESULT
    print("Saving final cluster labels in Single Cell object . . .")
    sc.addCellData(col_data = final_labels, col_name = method_name + str(n_clusters) + "_Clusters")
    sc.setCounts(X)

    return sc



def Cluster(    
    sc,                                
    n_clusters = "gap",
    n_clusters_max = 10, 
    classification_type = "f_classif",
    normalization = 'mean',
    q = None 
    ):

    """
    Clusters gene expression data stored in SingleCell object. Stores the cluster information
    (labels) in the celldata assay of the SingleCell object. The column name under which it
    stores the cluster information is 'FEATS_k_Clusters', where k is the number of clusters,
    which the method can estimate or is defined by the user as an int or an int list.

    Parameters
    ----------

    sc : SingleCell
        The SingleCell dataset containing gene expression data to be clustered.  
    
    n_clusters : int, list or str, optional
        This is an optional input for the number of clusters in the dataset. It is either an int,
        a Python list of ints or a str. If it is an int, the method will cluster the data into 
        n_clusters groups. If it is a Python list, the method will cluster the data into the number 
        of groups specified in the list and store the cluster information for all the values in the
        list. If it is the string 'gap', the method will use gap statistic to estimate the number 
        of clusters and cluster the data into the number of groups estimated. 
    
    n_clusters_max : int, optional
        The upper limit of the number of clusters when using gap statistic to estimate the 
        number of clusters. Ignored if n_clusters is not 'gap'. Default 10.
    
    normalization : str
        The normalization to use when normalizing the data. 

    Returns
    -------

    SingleCell
        The SingleCell object with the cluster information stored. 

    int
        The number of clusters in the dataset.  

    Raises
    ------
    
    ValueError
        If no common genes are found between the SingleCell datasets in the batches list. 

    """










    if (type(q) == type(None)):
        q = np.arange(1, int((5/100) * sc.dim[0]) + 1) # Default value of q

    if (type(n_clusters) is str):
        if (n_clusters == "gap"):
            if ((n_clusters_max != None) and (n_clusters_max >= 2)):
                n_clusters, _, _, _, _ = GapStatistics( 
                    sc,
                    n_clusters_max = n_clusters_max,
                    B = 500
                    )
                
                print ("Number of clusters estimated by gap statistics: ", n_clusters)

                sc = AnovaHierarchical(   
                    sc,
                    n_clusters = n_clusters,
                    classification_type = classification_type,
                    normalization = normalization,
                    q = q
                    )

            else:
                raise ValueError("If n_clusters is 'gap', then n_clusters_max should be an integer greater than or equal to 2.")


        else:
            raise ValueError("Unidentified string value for n_clusters. Supports only the string 'gap', to determine n_clusters using gap statistics")

    elif (type(n_clusters) is list):
        if (type(n_clusters[0]) is int):

            for n in n_clusters:

                sc = AnovaHierarchical(   
                    sc,
                    n_clusters = n,
                    classification_type = classification_type,
                    normalization = normalization,
                    q = q
                    )

        else:
            raise ValueError("Unknown type of values in list n_clusters. Only supports int types.")
    
    elif (type(n_clusters) is int):

        sc = AnovaHierarchical(   
            sc,
            n_clusters = n_clusters,
            classification_type = classification_type,
            normalization = normalization,
            q = q
            )

    else:

        raise ValueError("Unknown value for n_clusters. It should be an int, list of ints or the string 'gap'.")


    return sc, n_clusters