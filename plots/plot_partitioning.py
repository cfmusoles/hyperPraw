import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import collections

folder = "../results/hypergraph_generator/cluster_density/"
hgraph_files = ["custom.hgr"]#["2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","webbase-1M.mtx.hgr"]# ["atmosmodj.mtx.hgr","kkt_power.mtx.hgr"]
experiment_prefix = 'test_default'
partition = 'zoltanVertex'
num_processes = 12
expected_ratios = [1.0/num_processes for _ in range(num_processes)]
#expected_ratios = [0.05,0.05,0.05,0.05,0.05,0.05,0.1,0.05,0.20,0.25,0.05,0.05]
storePlot = True
image_format = 'pdf'
image_names = [partition + "_" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    #load partitioning scheme
    partitioning = np.genfromtxt(folder + experiment_prefix + '_' + hgraph + '_' + partition + '_partitioning__' + str(num_processes),skip_header=0,delimiter=",")
    #partitioning = np.genfromtxt(folder + hgraph + '_clustering',skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    counts = collections.Counter(partitioning)
    df = pd.DataFrame(list(counts.values()),columns=['Elements'],index=counts.keys())
    df = df.sort_index()
    df['Expected ratio'] = expected_ratios

    #boxplot of all partitions
    fig, ax = plt.subplots()
    df.boxplot(column='Elements',ax=ax)
    
    
    if storePlot:
        plt.savefig(image_names[i] + "_box_plot." + image_format,format=image_format)
    plt.show()

    #bar plot of each individual partition
    fig, ax = plt.subplots()
    
    total = sum(df['Elements'])
    df['Elements'] /= total
    df.plot.bar(ax=ax)
    ax.set_xlabel('Partitions')
    ax.set_ylabel('Elements ratio')

    if storePlot:
        plt.savefig(image_names[i] + "_ratios." + image_format,format=image_format)
    plt.show()

    

    