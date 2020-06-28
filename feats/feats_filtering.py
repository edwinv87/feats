
import numpy as np

'''
Find Highly Variable Genes




'''

def LogFilter(sc):
    X = sc.getCounts()
    X = np.log2(X + 1)
    sc.setCounts(X)

    return sc


def HVGFilter(sc, num_genes, name_by = 'gene_names'):

    # compute dispersion and store it in sc object (genedata assay)
    X = sc.getCounts()

    var_bat = np.var(X, axis = 1)
    mu_bat = np.mean(X, axis = 1)
    dispersion = var_bat/mu_bat
    idx = np.argsort(dispersion)
    idx = idx[::-1]

    # Print if gene names are available
    if (sc.checkGeneData(name_by)):
        gene_names = sc.getGeneData(name_by)
        print("The top ", num_genes, " highly variable genes are: ")
        print(gene_names[idx[0:num_genes]])

    hvg = np.zeros(dispersion.shape, dtype = bool)

    for i in range(num_genes):
        hvg[idx[i]] = True

    # Save the data
    sc.addGeneData(dispersion, "FEATS_dispersion")
    sc.addGeneData(var_bat, "FEATS_variance")
    sc.addGeneData(mu_bat, "FEATS_mean")
    sc.addGeneData(hvg, "FEATS_hvg")

    sc = sc[hvg, :]

    return sc


def GeneFilter(sc, min_cells, max_cells, expr_thres = 0):

    print("Applying Gene Filter . . .")
    X = sc.getCounts()
    
    expr_count = np.sum(X > expr_thres, axis = 1, keepdims = False)
    gene_filter = (expr_count <= max_cells) & (expr_count >= min_cells)


    print("Number of features remaining after gene filtering: ", np.sum(gene_filter))

    # If gene filtering removes all features:
    if (np.sum(gene_filter) != 0):
        # Save the gene filter mask
        # sc_obj.addGeneDataColumn(col_data = gene_filter, col_name = "gene_filter")
        sc = sc[gene_filter, :]

    else:
        print ("Gene filtering removed all the features!")
        print ("Set different gene filter parameters.")
        print ("Skipping Gene Filtering step...")

    return sc



def CellFilter(sc, min_genes, max_genes, expr_thres = 0):
    
    print("Applying Cell Filter . . .")

    X = sc.getCounts()

    expr_count = np.sum(X , axis = 0, keepdims = False)
    cell_filter = (expr_count <= max_genes) & (expr_count >= min_genes)

    print("Number of samples remaining after cell filtering: ", np.sum(cell_filter))

    # If cell filtering removes all samples:
    if (np.sum(cell_filter) != 0):
        # Save the cell filter mask
        # sc_obj.addCellDataColumn(col_data = cell_filter, col_name = "cell_filter")
        sc = sc[:, cell_filter]
    else:
        print ("Cell filtering removed all the samples!")
        print ("Set different cell filter parameters.")
        print ("Skipping Cell Filtering step...")

    return sc