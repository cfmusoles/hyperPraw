## Plot element degree (vertex degree for edge partitioning, hedge degree for vertex partitioning) vs whether or not the element has been replicated (i.e. partitioned)
## Very relevant to demonstrate that HDRF actually prefers to replicate elements with high degree

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


hgraphs_folder = '../resources/'
#hgraph_files = ["atmosmodj.mtx.hgr","kkt_power.mtx.hgr","sat14_velev-vliw-uns-2.0-uq5.cnf.dual.hgr"]
hgraph_files = ["random_power_law.hgr"]#["2cubes_sphere.mtx.hgr","ABACUS_shell_hd.mtx.hgr","sparsine.mtx.hgr","sat14_10pipe_q0_k.cnf.primal.hgr","webbase-1M.mtx.hgr"]
experiment_folder = "../results/"
experiment_prefix = 'test_overlap'
partition = 'rHDRF'
num_processes = 12

storePlot = False
image_format = 'pdf'
image_names = ["rHDRF_" + str(i) for i,v in enumerate(hgraph_files)]

for i,hgraph in enumerate(hgraph_files):
    hgraph_file = hgraphs_folder + hgraph

    #load partitioning scheme
    partitioning = np.genfromtxt(experiment_folder + experiment_prefix + '_' + hgraph + '_' + partition + '_partitioning__' + str(num_processes),skip_header=0,delimiter=",")
    partitioning = [int(i)  for i in partitioning]

    # load pin degrees from hgraph file
    with open(hgraph_file) as f:
        f.readline() #skip header info
        elements = [x.rstrip().split(" ") for x in f]
        elements = [list( map(int,i) ) for i in elements]
        pin_degrees = [len(x) for x in elements]

    # create pandas frame with two columns:
    #   pin degree
    #   whether it is replicated
    replicated = []
    for element in elements:
        parts = set()
        for pins in element:
            parts.add(partitioning[pins-1])
        replicated.append(len(parts) > 1)


    df = pd.DataFrame({'Pin degrees' : pin_degrees, 'Replicated' : replicated})

    g = sns.catplot( x="Replicated", y="Pin degrees", data=df, legend=False, kind='violin')
    #ax = g.fig.get_axes()[0].set_yscale('log')

    if storePlot:
        plt.savefig(image_names[i]+ image_format,format=image_format,dpi=1000)
    plt.show()