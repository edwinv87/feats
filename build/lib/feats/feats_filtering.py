
import numpy as np


def LogFilter(sc):
    """
    Computes and stores the log-transformed gene expression data stored in 
    SingleCell object. Here Log2 is performed after adding a pseudocount of 1.

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  
    
    Returns
    -------

    sc : SingleCell
        The SingleCell object containing the log-transfprmed gene expression data.

    """

    X = sc.getCounts()
    X = np.log2(X + 1)
    sc.setCounts(X)

    return sc


def HVGFilter(sc, num_genes, name_by = 'gene_names'):
    """
    Filters/Selects the Highly Variable Genes in the SingleCell object by computing 
    the dispersion of each gene. Returns a sliced SingleCell object containing the top 
    Highly Variable Genes. Stores additional information such as dispersion ('FEATS_dispersion), 
    mean ('FEATS_mean') and variance ('FEATS_variance') in the genedata assay of SingleCell 
    object.

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  

    num_genes : int
        The number of Highly Variable Genes to select. 
    
    name_by : str
        The name of the column in SingleCell object that stores the gene names. This 
        is used to print the top Highly Variable Genes.
    
    Returns
    -------

    sc : SingleCell
        A sliced SingleCell object containing the top Highly Variable Genes. 

    """

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
    """
    Filters the lowly and/or highly expressed genes stored in the SingleCell object 
    based on the expression counts. Returns the sliced SingleCell object with the
    genes selected through the filtering criteria. This functions skips filtering and
    warns the user if the filtering criteria filters out all the genes.

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  

    min_cells : int
        The minimum number of cells in which the gene must be expressed with the 
        expression threshold. The value should be >= 0 and <= n, where n is the 
        number of cells the the dataset.  
    
    max_cells : int
        The maximum number of cells in which the gene must be expressed with the 
        expression threshold. The value should be >= 0 and <= n, where n is the 
        number of cells the the dataset.  

    expr_thres : int, optional
        The expression threshold. Default 0. The gene is considered not expressed if
        the expression count is <= this threshold. 
    
    Returns
    -------

    sc : SingleCell
        A sliced SingleCell object containing the filtered genes. 

    """


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
    """
    Filters the cells in which the genes are lowly and/or highly expressed. Returns 
    the sliced SingleCell object with the cells selected through the filtering criteria. 
    This functions skips filtering and warns the user if the filtering criteria filters 
    out all the samples/cells.

    Parameters
    ----------

    sc : SingleCell
        The SingleCell object containing gene expression data.  

    min_genes : int
        The minimum number of genes in the cell in which the genes must be expressed 
        with the expression threshold. The value should be >= 0 and <= d, where d is the 
        number of genes the the dataset.  
    
    max_genes : int
        The maximum number of genes in the cell in which the gene must be expressed 
        with the expression threshold. The value should be >= 0 and <= d, where d is the 
        number of genes the the dataset.  

    expr_thres : int, optional
        The expression threshold. Default 0. The gene is considered not expressed if
        the expression count is <= this threshold. 
    
    Returns
    -------

    sc : SingleCell
        A sliced SingleCell object containing the filtered cells. 

    """

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