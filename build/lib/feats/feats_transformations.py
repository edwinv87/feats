import numpy as np
import pandas as pd
from singlecelldata import SingleCell
from sklearn.metrics import pairwise_distances
from scipy.stats import spearmanr
from .feats_utils import FeatureNormalize


### Pairwise Distances 

def spearmans(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D, _ = spearmanr(X, axis = 1)
    D = 1 - D
    
    return D


def pearsons(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D = pairwise_distances(X, metric = "correlation")
    
    return D
    
    
def euclidean(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns D which is a (n x n) matrix of distances between samples
    
    """
    D = pairwise_distances(X, metric = "euclidean", n_jobs = -1)
    
    return D



### Kernels

def linear_kernel(X):
    """
    X is a (n x d) matrix where rows are samples
    Returns K which is a (n x n) kernel matrix
    
    """    
    K = np.dot(X, X.T)
    
    return K




### Utility Functions
# Converts Distances to Kernels
def dist_to_kernel(D):
    
    gamma = 1.0 / np.amax(D)
    S = np.exp(-D * gamma)
    
    return S



# Computes Gram Matrix from Kernel
def GramMatrix(K):
    
    N_one = np.ones(K.shape) / K.shape[0] 
    K_bar = K - np.dot(N_one, K) - np.dot(K, N_one) + np.dot(np.dot(N_one, K), N_one)
    
    return K_bar



def PCA(sc, n_comp = 1, dist_or_kernel = 'linear'):
    """
    Computes and stores the Principal Components of the gene expression data stored in 
    the SingleCell object.  

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  
    
    n_comp : int, optional
        The number of Principal Components to compute. Default 1. 

    dist_or_kernel : str, optional
        The distance metric or the kernel to compute. If a distance metric is passed,
        it computes the pairwise distance between the cells and then converts the distance
        metrics to kernels. If a kernel is passed, it computes the kernel. Valid values 
        are 'linear' (default) for linear kernel, 'spearmans' for Spearmans distance, 
        'euclidean' for Euclidean distance and 'pearsons' for Pearsons distance. 

    Returns
    -------

    sc : SingleCell
        The SingleCell object containing the dimensionnality reduced gene expression data. The 
        reduced dimensionality is n_comp. The gene names are removed and the features in the 
        reduced space are named 'PC1', 'PC2' and so on. 

    """

    sc = FeatureNormalize(sc, 'mean')
    X = sc.getCounts()

    if (dist_or_kernel == 'linear'):
        K = GramMatrix(linear_kernel(X.T))

    elif (dist_or_kernel == 'spearmans'):
        K = GramMatrix(dist_to_kernel(spearmans(X.T)))

    elif (dist_or_kernel == 'euclidean'):
        K = GramMatrix(dist_to_kernel(euclidean(X.T)))

    elif (dist_or_kernel == 'pearsons'):
        K = GramMatrix(dist_to_kernel(pearsons(X.T)))

    E_vec, E_val, _ = np.linalg.svd(K)

    #print(E_val)
    #print(E_vec[:,0])
    # Sort Eigenvector 
    idx = np.argsort(E_val, kind = 'mergesort')
    idx = np.flip(idx)
    E_val = E_val[idx]
    E_vec = E_vec[:, idx]

    # Remove zero eigenvectors
    idx2 = E_val > 0
    E_val = E_val[idx2]
    E_vec = E_vec[:, idx2]
    # print("Maximum components possible = ", E_val.size)

    # Scale eigenvectors so that np.dot(D[:,0].T, D[:, 0]) * E[0] = 1

    E_val = np.reshape(E_val, [1, E_val.size])
    # E_vec = E_vec / np.linalg.norm(E_vec, axis = 0, keepdims = True)
    E_vec = E_vec / np.sqrt(E_val)
    X_red = np.dot(E_vec[:, 0:n_comp].T, K)

    data = pd.DataFrame(X_red)
    celldata = sc.celldata
    pcs = []
    for i in range(n_comp):
        pcs.append('PC' + str(i+1))

    genedata = pd.DataFrame(pcs, index = data.index, columns = ['principal_components'])
    sc_red = SingleCell(dataset = sc.dataset + "_reduced", data = data, celldata = celldata, genedata = genedata)
    # print(E_vec[:,0])
    # print(X_red)
    return sc_red



