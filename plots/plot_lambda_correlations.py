# Plots the correlations between lambda (HDRF parameter) and other hgraph properties

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

hgraphs_folder = '../resources/'
hgraph_files = ["sat14_itox_vc1130.cnf.dual.hgr","2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","pdb1HYS.mtx.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","sat14_E02F22.cnf.hgr","webbase-1M.mtx.hgr"]
experiment_prefix = '../test_default_'
partitioning = 'rHDRF'
num_processes = 12

lambda_column = 1

storePlots = False
image_format = 'pdf'

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

hgraph_col = []
lambda_col = []
hedges_col = []
vert_col = []
vert_hedge_ratio_col = []
cardinality_col = []
nnz_col = []

for i,hgraph in enumerate(hgraph_files):
    
    # load lambda
    part_history_file = experiment_prefix + hgraph + "_" + str(num_processes) + "_0.hgr_partition_history__" + str(num_processes)
    part_history = get_data_from_csv(part_history_file)
    lambda_value = part_history[-2][lambda_column]
    lambda_col.append(lambda_value)

    # load hgraph parameters
    hgraph_file = hgraphs_folder + hgraph
    hgraph_col.append(hgraph)
    with open(hgraph_file) as f:
        print("Loading graph: " + hgraph_file)
        # hyperedges and vertices
        header = f.readline().rstrip().split(" ") #skip header info
        hyperedges = int(header[0])
        vertices = int(header[1])
        hedges_col.append(hyperedges)
        vert_col.append(vertices)
        vert_hedge_ratio_col.append(vertices * 1.0 / hyperedges)
        # nonzeros and cardinality
        hedges = [x.rstrip().split(" ") for x in f]
        hedge_sizes = [len(x) for x in hedges]
        nnz_col.append(np.sum(hedge_sizes))
        cardinality_col.append(np.mean(hedge_sizes))
        #hedges = [int(vid.replace("\n","")) for hedge in hedges for vid in hedge] #np.array([int(x.replace("\n","")) for x in hedges]).flatten()
        #unique, vertex_counts = np.unique(hedges,return_counts=True)

# convert to dataframe
data = pd.DataFrame({'Hypergraph': hgraph_col,
                        'Lambda' : lambda_col, 
                        'Hyperedges' : hedges_col,
                        'Vertices' : vert_col,
                        'VH ratio' : vert_hedge_ratio_col,
                        'Cardinality' : cardinality_col,
                        'NNZ' : nnz_col})

# plot
#sns.jointplot(x='Hyperedges',y='Lambda',data=data,kind='reg')
sns.pairplot(data,hue='Hypergraph')

if storePlots:
        plt.savefig("lambda." + image_format,format=image_format,dpi=1000)
plt.show()