# Import core libraries
import json
import pandas as pd
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score

# Import our libraries
from singlecelldata import SingleCell
from feats import Cluster, GeneFilter, LogFilter

# Python dictionary storing dataset specific parameters.
dataset_params = json.load(open("dset_param.txt"))
# print (dataset_params)

n_trials = 100
datasets = ['biase', 'yan1', 'yan2', 'goolam', 'deng1', 'deng2', 'pollen', 'kolodziejczyk', 'patel', 'fan', 'treutlein']
ari = np.zeros([n_trials, len(datasets)])

col = 0
for dataset in datasets:
    # Load Data
    dset_name = dataset_params[dataset]['dset_name']
    data_path = "../../New_Datasets/" + dset_name + '/' + dset_name + "_data.csv"
    celldata_path = "../../New_Datasets/" + dset_name + '/' + dset_name + "_celldata.csv"
    genedata_path = "../../New_Datasets/" + dset_name + '/' + dset_name + "_genedata.csv"

    data = pd.read_csv(data_path, index_col=0)
    celldata = pd.read_csv(celldata_path, index_col=0)
    genedata = pd.read_csv(genedata_path, index_col = 0)

    # Create a single cell object   
    sc = SingleCell(dataset,data,celldata,genedata)

    # Get parameters from the dictionary
    label = dataset_params[dataset]['label']
    nc = dataset_params[dataset]['nc']

    # Select only the first 49 cells for analysis of biase dataset
    if (dataset == "biase"):
        sc = sc[0:49]

    # Log Filtering
    if (dataset != "patel"):
        sc = LogFilter(sc)
    
    # Gene Filtering
    sc = GeneFilter(sc, min_cells = int((10/100) * sc.dim[1]), max_cells = int((90/100) * sc.dim[1]), expr_thres = 0)

    row = 0

    for i in range(n_trials):
        # Perform clustering 
        sc, _ = Cluster(    sc,
                            normalization="mean", 
                            k=nc)

        # Compute and print the Adjusted Rand Score
        ari[row, col] = adjusted_rand_score(sc.getNumericCellLabels(label), sc.getCellData("FEATS_" + str(nc) + "_Clusters"))
        print ("ARI for dataset " + dataset + " = ")
        print (ari[row, col])
        row = row + 1

    col = col + 1

df = pd.DataFrame(data = ari, columns = datasets)
df.to_excel("results/results-feats_clustering.xlsx")