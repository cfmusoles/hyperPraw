## Compare clustering done by partitioning algorithms with ground truth clustering (for generated hypergraphs) 
## Can use NMI (normalised mutual information) or https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score

data_folder = '../results/'
ground_truth_cluster = 'random_power_law.hgr_clustering'
num_clusters = 12
experiment_prefix = 'test_default'
hgraph = 'random_power_law.hgr'
partitioning_candidates = ['rHDRF','zoltanVertex']


# load ground truth
ground_truth_clustering = np.genfromtxt(data_folder + ground_truth_cluster,skip_header=0,delimiter=",")
ground_truth_clustering = [int(i)  for i in ground_truth_clustering]

for partition in partitioning_candidates:
    #load partitioning scheme
    partitioning = np.genfromtxt(data_folder + experiment_prefix + '_' + hgraph + '_' + partition + '_partitioning__' + str(num_clusters),skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    similarity_score = adjusted_rand_score(ground_truth_clustering, partitioning)
    ami = adjusted_mutual_info_score(ground_truth_clustering, partitioning)
    
    #output results
    print("{}:\n Adjusted_rand {:.3f}\n Adjusted mutual info {:.3f}".format(partition,similarity_score,ami))


    