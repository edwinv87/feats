# Import core libraries
import pandas as pd
import numpy as np
from time import time

from singlecelldata import SingleCell

# Import our libraries
from feats import Cluster

n_samples = np.array([1000,2000,3000,4000,5000,6000,7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000,16000,17000,18000,19000,20000])
#n_samples = np.array([1000,2000,3000])
n_trials = n_samples.size

clust_time = np.zeros([n_trials, 1])

col = 0

# Load Data
dset_name = "simulation"
data_path = "../../New_Datasets/" + dset_name + '/' + dset_name + "_data.csv"
data = pd.read_csv(data_path, index_col=0)

# Create a single cell object   
sc = SingleCell(dset_name,data)
sc.print()

nc = 3
row = 0

for i in range(n_trials):
    
    # Select n samples
    idx = np.random.randint(sc.dim[1], size=n_samples[i])
    sc_temp = sc[idx.tolist()]

    # Perform clustering 
    t0 = time()
    sc, _ = Cluster(    sc_temp,
                        normalization='l2',
                        k=nc)
    t1 = time()

    clust_time[row, col] = t1 - t0
    print("N = ", n_samples[i])
    print("Clustering time: ", clust_time[row, col])
    row = row + 1

    


df = pd.DataFrame(data = clust_time, columns = [dset_name])
df.to_excel("results/results-feats_comp_time_20000.xlsx")