import numpy as np

from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering

def ComputeWk(X, labels, classes):
    """
    X       - (d x n) data matrix, where n is samples and d is dimentionality
    lables  - n dimentional vector which are class labels

    """
    Wk = 0

    for i in range(classes):
        mask = (labels == i)
        Wk = Wk + np.sum(np.sum((X[:, mask] - np.mean(X[:, mask], axis=1, keepdims=True))**2, axis=0))

    return Wk 



def GapStatistics(  sc_obj,
                    k_max,
                    B):

    """
    Computes the gap statistic and estimates the number of clusters in the
    gene expression dataset contained in SingleCell object. 

    Parameters
    ----------

    sc_obj : SingleCell
        The SingleCell object containing gene expression data and the metadata.  
    
    k_max : int
        The upper limit of the number of clusters. 
    
    B : int
        The number of reference datasets to generate to compute the gap quantities.

    Returns
    -------

    est_clusters : int
        The estimate of the number of clusters in the dataset. 

    Gap_k : list
        The gap statistic quantity for gap. The list contains gap values for each k,
        where k = 1, 2, ..., k_max.
    
    s_k : list 
        The gap statistic quantity for standard deviation. The list contains the  
        standard deviation for each k, where k = 1, 2, ..., k_max.
    
    W_k : list
        A gap statistic quantity for each k, where k = 1, 2, ..., k_max.

    w_bar : list
        A gap statistic quantity for each k, where k = 1, 2, ..., k_max.

    """

    k_clusters = np.arange(1, k_max + 1)
    Gap_k = np.zeros(k_clusters.shape)
    s_k = np.zeros(k_clusters.shape)
    W_k = np.zeros(k_clusters.shape)
    w_bar = np.zeros(k_clusters.shape)

    X = sc_obj.getCounts() # d x n

    pc = PCA(n_components=0.99)
    X_red = pc.fit_transform(X.T)

    Xmin = np.amin(X_red, axis = 0)
    Xmax = np.amax(X_red, axis = 0)

    # First cluster the sc data
    for idx in range(k_clusters.shape[0]):

        k = k_clusters[idx]

        hc = AgglomerativeClustering(n_clusters=k)

        if (k == 1):
            labels = np.zeros(X_red.shape[0])

        else:

            labels = hc.fit_predict(X_red)

        W_k[idx] = ComputeWk(X_red.T, labels, k)

        BW_k = np.zeros(B)
        for j in range(B):
            # Generate uniform random data with the same shape as X_red
            X_b = np.random.uniform(Xmin, Xmax, size=X_red.shape)
            # Cluster this data 
            if (k == 1):
                labels_b = np.zeros(X_b.shape[0])
            else:
                labels_b = hc.fit_predict(X_b)
            # Find W_k for the clustered data
            BW_k[j] = ComputeWk(X_b.T, labels_b, k)

        w_bar[idx] = (1/B) * np.sum(np.log(BW_k)) 
        Gap_k[idx] = w_bar[idx] - np.log(W_k[idx])
        sd_k = np.sqrt((1/B) * np.sum((np.log(BW_k) - w_bar[idx])**2))
        s_k[idx] = np.sqrt(1 + (1/B)) * sd_k
        print("For k = ", k, " Gap = ", Gap_k[idx], " s_k = ", s_k[idx])

    # Then determine the number of clusters
    for idx in range(1, k_clusters.shape[0] - 1):
        if (Gap_k[idx] > (Gap_k[idx + 1] - s_k[idx + 1])):
            est_clusters = k_clusters[idx]
            break

        else:
            est_clusters = k_clusters[idx + 1]


    return est_clusters, Gap_k, s_k, W_k, w_bar