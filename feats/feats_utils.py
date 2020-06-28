
import numpy as np
from sklearn.preprocessing import normalize

# Function to normalize the features of X
def FeatureNormalize(sc, norm):
# X shape is (features, samples)
# norm is either l2, mean
# Returns a matrix X_nrm containing normalized values (same shape as X)
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