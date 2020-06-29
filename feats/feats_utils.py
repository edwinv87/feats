
import numpy as np
from sklearn.preprocessing import normalize

# Function to normalize the features of X
def FeatureNormalize(sc, norm):
    """
    Computes and stores the normalized gene expression data stored in 
    SingleCell object.  

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  
    
    norm : str
        The normalization to perform. Accepted values are 'l2', 'mean', 'norm6' and 
        'cosine'. 

    Returns
    -------

    sc : SingleCell
        The SingleCell object containing the normalized gene expression data.


    """
    

    # print ("Applying " + norm + " Normalization . . .")
    X = sc.getCounts()
    if (norm == "l2"):
        X_nrm = normalize(X.T, axis=0)
        X_nrm = X_nrm.T
        
    elif (norm == "mean"):
        mu = np.mean(X, axis=1, keepdims=True)
        sd = np.std(X, axis=1, keepdims=True)
        X_nrm = (X - mu)/sd
        
    elif (norm == "norm6"):
        min_f = np.amin(X, axis=1, keepdims=True)
        y = np.log(X + np.absolute(min_f) + 1)
        max_all = np.amax(y)
        X_nrm = y/max_all

    elif (norm == "cosine"):
        X_nrm = normalize(X.T, axis=1)
        X_nrm = X_nrm.T
    
    sc.setCounts(X_nrm)
    return sc