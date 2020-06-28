import numpy as np
import pandas as pd

from sklearn.metrics import pairwise_distances
from singlecelldata import SingleCell


def ComputeMNNPairs(X_ref, X, k):
    # X_ref, X - (d x n) gene expression matrix
    # k - number of nearest neighbours

    n1 = X_ref.shape[1]
    n2 = X.shape[1]
    mnn_pair1 = np.zeros([n1, n2], dtype = bool)
    mnn_pair2 = np.zeros([n1, n2], dtype = bool)

    # Pairwise distance between X_ref and X
    D = pairwise_distances(X = X_ref.T, Y = X.T, metric = "euclidean")

    for i in range(n1):
        idx = np.argsort(D[i, :])
        idx = idx[0:k]
        mnn_pair1[i, idx] = True

    for i in range(n2):
        idx = np.argsort(D[:, i])
        idx = idx[0:k]
        mnn_pair2[idx, i] = True  

    return np.logical_and(mnn_pair1, mnn_pair2), D


"""
===============================
Method Name: RemoveBatchEffects
===============================

Method Description: This method implements the MNN algorithm for removal of
batch effects in multiple datasets. 


Arguments:
========== 
batches - a list containing batches in the form of single cell object. The first batch will be taken as the reference batch.



Returns:
========


"""




def ReducedMean(X, idx):

    
    idx_unique = np.unique(idx)
    X_redmean = np.empty([X.shape[0], idx_unique.shape[0]])

    for i in range(idx_unique.shape[0]):

        idx_bool = (idx == idx_unique[i])

        if (np.sum(idx_bool) == 1):
            X_redmean[:, i] = np.squeeze(X[:, idx_bool])
        elif (np.sum(idx_bool) > 1):
            X_redmean[:, i] = np.squeeze(np.mean(X[:, idx_bool], axis=1, keepdims = True))

    return X_redmean, idx_unique


def LogspaceAdd(log_array):

    l = log_array.shape[0]
    cum = np.zeros(l)
    total = 0
    for i in range(l):

        if (i == 0):
            total = log_array[i]

        else:
            total = np.logaddexp(total, log_array[i])
        
        cum[i] = total

    return total, cum


def ComputeCorrectionVectors(V, X, tar_idx, sigma):

    ncells = X.shape[1]
    d = V.shape[0]
    nmnns = V.shape[1]
    weights = np.empty([nmnns, ncells])

    for i in range(nmnns):
        weights[i, :] =  -np.sum((X - X[:, tar_idx[i]].reshape(d,1))**2, axis=0)/(2*(sigma**2))


    top_weights = np.amax(weights, axis = 0)
    weights = np.exp(weights - top_weights)
    #weights = np.exp(weights)
    weights = weights / np.sum(weights, axis = 0)
    C_vecs = np.dot(V, weights) # Correction vectors in columns

    return C_vecs



def ComputeSpan(X, idx, svd_dim):

    # Center the data
    mu = np.mean(X[:, idx], axis = 1, keepdims=True)
    X_cen = X[:, idx] - mu
    U, _, _ = np.linalg.svd(X_cen)
    U = U[:, 0:svd_dim]

    return U



# This is a vectorized fuction
def SqDistToLine(ref, grad, point):

    # ref - (d x 1), the reference on which the correction is applied
    # grad - (d x 1), l2 normalized correction vector
    # point - (d x n), other points in the target batch
    d = ref.shape[0]
    n = point.shape[1]

    grad = grad.reshape(d, 1)
    
    w = ref.reshape(d, 1) - point
    scale = np.dot(grad.T, w)
    w = w - (scale * np.tile(grad, (1, n)))

    nrm = np.sum(w**2, axis = 0)

    return nrm



def AdjustShiftVariance(X_ref, X, Corr, sigma):

    # X_ref - (d x n1) matrix of reference batch 
    # X - (d x n2) gene expression for the target batch
    #
    #

    ncells = Corr.shape[1]

    scaling_factor = np.zeros(ncells)
    # Computing the l2 norm of correction vector and normalizing to unit l2 norm
    l2_norm = np.sum(Corr**2, axis = 0)
    Corr_norm = Corr / l2_norm

    # Projection of reference and target batch cells onto the correction vectors. (_ref represents reference batch quantities)
    Proj = np.dot(Corr_norm.T, X)
    Proj_ref = np.dot(Corr_norm.T, X_ref)

    for i in range(ncells):
        Nrm_same = SqDistToLine(X[:, i], Corr_norm[:, i], X)
        log_prob = -Nrm_same/sigma
        
        mask = (Proj[i, i] > Proj[i, :])

        prob, _ = LogspaceAdd(log_prob[mask])
        totalprob, _ = LogspaceAdd(log_prob)

        prob = prob - totalprob

        Nrm_other = SqDistToLine(X[:, i], Corr_norm[:, i], X_ref)
        log_prob_ref = -Nrm_other/sigma
        
        
        proj_ref_i = Proj_ref[i, :]

        proj_ref_i_idx = np.argsort(proj_ref_i)
        proj_ref_i = proj_ref_i[proj_ref_i_idx]
        log_prob_ref = log_prob_ref[proj_ref_i_idx]
        totalprob_ref, cumsum_ref = LogspaceAdd(log_prob_ref)

        target = prob + totalprob_ref

        mask2 = (cumsum_ref >= target)

        if (np.sum(mask2) < 2):
            ref_quan = np.flip(proj_ref_i)

        else:
            ref_quan = proj_ref_i[mask2]
        

        scaling_factor[i] = (ref_quan[0] - Proj[i, i])/l2_norm[i]
    
    return scaling_factor



def RemoveBatchEffects(ref_batch, tar_batch, k , sigma, svd_dim, adj_var):

    # Get the expression of the reference batch and normalize
    X_ref = ref_batch.getCounts()
    X_tar = tar_batch.getCounts()


    # Compute the mutual nearest neighbours (MNN's)
    mnn, _ = ComputeMNNPairs(X_ref, X_tar, k)
    ref_idxs, tar_idxs = np.nonzero(mnn)

    # Compute batch effect as the vector difference of the reference and target batches
    V = X_ref[:, ref_idxs] - X_tar[:, tar_idxs]

    # Take the average of correction vectors for the target cells
    V, tar_idx_uni = ReducedMean(V, tar_idxs)

    # Compute the correction vectors for all the cells in the target batch using Gaussian kernel smoothing
    C = ComputeCorrectionVectors(V, X_tar, tar_idx_uni, sigma)

    # Compute biological span and remove the batch vector components parallel to biological basis vectors
    if (svd_dim > 0):
        ref_idx_uni = np.unique(ref_idxs)
        U_ref = ComputeSpan(X_ref, ref_idx_uni, svd_dim)
        U_tar = ComputeSpan(X_tar, tar_idx_uni, svd_dim)

        for U in (U_tar, U_ref):
            bio_comp = np.dot(C.T, U)
            bio_comp = np.dot(bio_comp, U.T)
            C = C - bio_comp.T

    # Adjust shift variance
    if (adj_var):
        scaling = AdjustShiftVariance(X_ref, X_tar, C, sigma)
        # print(scaling)
        # print (np.sum(scaling))
        scaling = np.maximum(scaling, 1)
        # print(np.sum(scaling))
        C = scaling * C

    X_tar = X_tar + C

    # Modify the expression values of the ref and tar batches
    tar_batch.setCounts(X_tar)
    ref_batch.setCounts(X_ref)

    # Create a list 
    batches = [ref_batch, tar_batch]

    # Merge the reference and target batches and return 
    sc = MergeBatches(batches)

    return sc


def IntegrateBatches(batches, name_by):

    print ("Integrating Batches . . .")

    # Incorporate batch info to sc object 
    for i in range(len(batches)): 
        batch_no = [batches[i].dataset] * batches[i].dim[1]
        batches[i].addCellData(col_data = batch_no, col_name = 'batch')
    
    # Find common genes in all the datasets
    keep_genes = set()
    for i in range(len(batches)):
        gene_list = batches[i].getGeneData(name_by[i]).tolist()
        if len(keep_genes) == 0:
            keep_genes = set(gene_list)
        else:
            keep_genes &= set(gene_list)
        #print ("Num genes: ", len(keep_genes))

    # Warn user and exit if no common genes are found
    if len(keep_genes) == 0:
        raise ValueError('No common genes found in all datasets! Check if datasets have common genes.')


    # Print the number of common genes in all datasets
    print("Number of common genes in all batches: ", len(keep_genes))

    # Filter the sc object 
    ret_genes = np.array(sorted(keep_genes))
    for i in range(len(batches)):
        # Remove duplicate genes.
        genes = batches[i].getGeneData(name_by[i]).tolist()
        uniq_genes, uniq_idx = np.unique(genes, return_index=True)
        Xi = batches[i].getCounts()
        Xi = Xi[uniq_idx, :]

        sort_idx = np.argsort(uniq_genes)

        gene_idx = []
        for idx in sort_idx:
            if (uniq_genes[idx] in ret_genes):
                gene_idx.append(idx)

        assert(np.array_equal(uniq_genes[gene_idx], ret_genes))

        # Combine the dataframes of the batches together
        Xi = Xi[gene_idx, :]
        data = pd.DataFrame(Xi)
        celldata = batches[i].celldata
        genedata = pd.DataFrame(ret_genes, index = data.index, columns = ['gene_names'])
        batches[i] = SingleCell(dataset = batches[i].dataset, data = data, celldata = celldata, genedata = genedata)

    batches = MergeBatches(batches)

    return batches


def MergeBatches(batches):

    print("Merging Batches . . .")
    data_frames = [batches[0].data]
    cell_data_frames = [batches[0].celldata]
    merged_name = batches[0].dataset

    for i in range(1, len(batches)):
        data_frames.append(batches[i].data)
        cell_data_frames.append(batches[i].celldata)
        merged_name = merged_name + '+' + batches[i].dataset


    # Return the combined object with a batches in the celldata array indicating the batch number/dataset name
    merged_data = pd.concat(data_frames, axis = 1, ignore_index = True, sort = False)
    merged_celldata = pd.concat(cell_data_frames, axis = 0, ignore_index = True, sort = False)
    
    # Assuming that the genes in all the batches is the same, i.e., batches have been integrated using IntegrateBatches
    merged_genedata = batches[0].genedata
    #print(merged_data.shape)
    #print(merged_celldata.shape)

    # Create a new singlecell object 
    sc = SingleCell(dataset = merged_name, data = merged_data, celldata = merged_celldata, genedata = merged_genedata)

    return sc



def CorrectBatches(batches, correct_order, name_by = None, k = 20, sigma = 10, svd_dim = 0, adj_var = True):

    # Check the following:
    # 1. at least two batches must be present
    # 2. if 'integrated == False', name_by must be a list with the same length as batches

    if (len(correct_order) < 2):
        print('<correct_order> - specifies less than two batches are to be merged. Specify two or more batches to be corrected.')

    batch_info = batches.getCellData('batch')
    
    correct_order_split = correct_order[0].split('+')
    if (len(correct_order_split) == 1):
        idx = (batch_info == correct_order[0])
    else:
        idx = np.zeros(batches.dim[1], dtype = bool)
        for i in range(len(correct_order_split)):
            idx = np.logical_or(idx, (batch_info == correct_order_split[i]))
  
    ref_batch = batches[:, idx.tolist()]
    print(ref_batch.dim)
    ref_batch.dataset = correct_order[0]

    for i in range(1, len(correct_order)):
        correct_order_split = correct_order[i].split('+')
        if (len(correct_order_split) == 1):
            idx = (batch_info == correct_order[i])
        else:
            idx = np.zeros(batches.dim[1], dtype = bool)
            for j in range(len(correct_order_split)):
                idx = np.logical_or(idx, (batch_info == correct_order_split[j]))

        tar_batch = batches[:, idx.tolist()]
        tar_batch.dataset = correct_order[i]
        print('Correcting batches <', ref_batch.dataset, ' and ', tar_batch.dataset, '>' )
        ref_batch = RemoveBatchEffects( ref_batch, 
                                        tar_batch, 
                                        k = k, 
                                        sigma = sigma, 
                                        svd_dim = svd_dim, 
                                        adj_var = adj_var)

    return ref_batch